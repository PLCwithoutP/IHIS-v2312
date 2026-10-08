#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./thin_time_dirs.sh [CASE_DIR] [N_KEEP] [--apply]
#
# Examples:
#   ./thin_time_dirs.sh . 1000
#   ./thin_time_dirs.sh . 1000 --apply
#
# Default behavior is DRY RUN.
# Add --apply to actually delete folders.

CASE_DIR="${1:-.}"
N_KEEP="${2:-1000}"
MODE="${3:-}"

if [[ ! -d "$CASE_DIR" ]]; then
    echo "Error: case directory '$CASE_DIR' does not exist."
    exit 1
fi

if ! [[ "$N_KEEP" =~ ^[0-9]+$ ]] || [[ "$N_KEEP" -lt 2 ]]; then
    echo "Error: N_KEEP must be an integer >= 2."
    exit 1
fi

APPLY=0
if [[ "$MODE" == "--apply" ]]; then
    APPLY=1
fi

cd "$CASE_DIR"

# Collect numeric time directories only.
# Accepts forms like:
#   0
#   1e-09
#   2.5e-06
#   0.001
mapfile -t ALL_TIMES < <(
    find . -maxdepth 1 -mindepth 1 -type d -printf '%f\n' \
    | awk '
        /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/ ||
        /^[.][0-9]+([eE][+-]?[0-9]+)?$/ {
            print
        }
    ' \
    | sort -g
)

if [[ "${#ALL_TIMES[@]}" -eq 0 ]]; then
    echo "No numeric time directories found."
    exit 0
fi

# Split out positive times for log spacing
mapfile -t POS_TIMES < <(
    printf '%s\n' "${ALL_TIMES[@]}" | awk '$1+0 > 0 {print}'
)

if [[ "${#POS_TIMES[@]}" -eq 0 ]]; then
    echo "Only time directory '0' exists. Nothing to thin."
    exit 0
fi

FIRST_POS="${POS_TIMES[0]}"
LAST_TIME="${ALL_TIMES[-1]}"

# If there are already fewer than or equal to requested count (+ keep 0),
# nothing useful to do.
if [[ "${#ALL_TIMES[@]}" -le $((N_KEEP + 1)) ]]; then
    echo "Found ${#ALL_TIMES[@]} time folders, which is already <= N_KEEP+1."
    echo "Nothing to thin."
    exit 0
fi

TMP_ALL="$(mktemp)"
TMP_KEEP="$(mktemp)"
trap 'rm -f "$TMP_ALL" "$TMP_KEEP"' EXIT

printf '%s\n' "${ALL_TIMES[@]}" > "$TMP_ALL"

# Build a keep-set:
# - always keep 0 if present
# - always keep latest time
# - keep nearest folders to N_KEEP log-spaced targets
awk -v nkeep="$N_KEEP" -v first="$FIRST_POS" -v last="$LAST_TIME" '
BEGIN {
    # Always keep 0 and latest
    keep["0"] = 1
    keep[last] = 1

    # Read all times into arrays later
}
{
    t[NR] = $1
    v[NR] = $1 + 0.0
    n = NR
}
END {
    # Mark all positive indices
    np = 0
    for (i = 1; i <= n; ++i) {
        if (v[i] > 0.0) {
            ++np
            p_idx[np] = i
        }
    }

    if (np == 0) {
        for (k in keep) print k
        exit
    }

    # If there are fewer positive times than requested, keep all positives
    if (np <= nkeep) {
        for (i = 1; i <= np; ++i) {
            keep[t[p_idx[i]]] = 1
        }
        for (k in keep) print k
        exit
    }

    log_first = log(first + 0.0)
    log_last  = log(last + 0.0)

    # Generate log-spaced targets and keep nearest existing folder
    for (j = 0; j < nkeep; ++j) {
        frac = (nkeep == 1 ? 0.0 : j / (nkeep - 1.0))
        target = exp(log_first + frac * (log_last - log_first))

        best_i = p_idx[1]
        best_d = (v[best_i] > target ? v[best_i] - target : target - v[best_i])

        for (i = 2; i <= np; ++i) {
            idx = p_idx[i]
            d = (v[idx] > target ? v[idx] - target : target - v[idx])
            if (d < best_d) {
                best_d = d
                best_i = idx
            }
        }

        keep[t[best_i]] = 1
    }

    for (k in keep) print k
}
' "$TMP_ALL" | sort -g -u > "$TMP_KEEP"

mapfile -t KEEP_TIMES < "$TMP_KEEP"

# Build delete list
mapfile -t DELETE_TIMES < <(
    awk '
        NR==FNR { keep[$1]=1; next }
        !($1 in keep) { print }
    ' "$TMP_KEEP" "$TMP_ALL"
)

echo "Case directory      : $(pwd)"
echo "Total time folders  : ${#ALL_TIMES[@]}"
echo "Requested keep count: $N_KEEP (+ 0 and latest if applicable)"
echo "Actual kept folders : ${#KEEP_TIMES[@]}"
echo "Folders to delete   : ${#DELETE_TIMES[@]}"
echo

if [[ "${#DELETE_TIMES[@]}" -eq 0 ]]; then
    echo "Nothing to delete."
    exit 0
fi

if [[ "$APPLY" -eq 0 ]]; then
    echo "DRY RUN. These folders would be deleted:"
    printf '  %s\n' "${DELETE_TIMES[@]}"
    echo
    echo "Re-run with --apply to actually remove them."
    exit 0
fi

echo "Deleting folders..."
for d in "${DELETE_TIMES[@]}"; do
    rm -rf -- "$d"
    echo "  deleted $d"
done

echo
echo "Done."
