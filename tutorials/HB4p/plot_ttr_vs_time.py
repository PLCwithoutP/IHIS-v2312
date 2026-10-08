#!/usr/bin/env python3
"""
Extract heat-bath time series from an IHIS (or nexaFoam) case into ONE CSV,
optionally with a quick-look plot.

Data source, in order of preference:
  1. postProcessing/<fo>/<startTime>/volFieldValue*.dat  (every time step; IHIS
     cases written by makeHeatBathCases.sh have the 'hbState' function object)
  2. numeric time directories  (only the written times, e.g. every 100th step)

The t = 0 row is taken from 0/initialConditions with the right key per field:
  TTR -> TTRValue,  TVib and TVib_<m> -> TVibValue,  p -> pressureValue
(override or extend with --initial-key FIELD=KEY). Fields without a key get NaN at t = 0.

CSV output:
    Time (s),TTR,TVib[,TVib_N2,...]

Examples:
    python3 plot_ttr_vs_time.py --case . --fields TTR TVib --output hB_case1.csv
    python3 plot_ttr_vs_time.py --case . --fields TTR TVib TVib_N2 TVib_O2 \
        --output hB_case3.csv --plot hB_case3_quicklook.png
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

NUM = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
TIME_RE = re.compile(rf"^{NUM}$")

DEFAULT_INITIAL_KEYS = {"TTR": "TTRValue", "p": "pressureValue"}


def initial_key_for(field: str, overrides: dict[str, str]) -> str | None:
    if field in overrides:
        return overrides[field]
    if field in DEFAULT_INITIAL_KEYS:
        return DEFAULT_INITIAL_KEYS[field]
    if field == "TVib" or field.startswith("TVib_"):
        return "TVibValue"
    return None


# --------------------------------------------------------------------------- #
# 0/initialConditions
# --------------------------------------------------------------------------- #
def read_initial_value(init_path: Path, key: str | None) -> float:
    if key is None or not init_path.exists():
        return math.nan
    text = init_path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(rf"\b{re.escape(key)}\b\s+({NUM})\s*;", text)
    return float(m.group(1)) if m else math.nan


# --------------------------------------------------------------------------- #
# Source 1: function-object output (every time step)
# --------------------------------------------------------------------------- #
def find_fo_files(case: Path, fo_name: str | None) -> list[Path]:
    pp = case / "postProcessing"
    if not pp.is_dir():
        return []
    fo_dirs = [pp / fo_name] if fo_name else sorted(d for d in pp.iterdir() if d.is_dir())
    files: list[Path] = []
    for fo in fo_dirs:
        if not fo.is_dir():
            continue
        for start in sorted((d for d in fo.iterdir() if d.is_dir() and TIME_RE.match(d.name)),
                            key=lambda d: float(d.name)):
            files.extend(sorted(start.glob("volFieldValue*.dat")))
        if files and fo_name is None:
            break  # first function object with data
    return files


def read_fo_series(files: list[Path], fields: list[str]) -> dict[float, dict[str, float]]:
    """Merge every volFieldValue*.dat (restarts: later files overwrite earlier times)."""
    rows: dict[float, dict[str, float]] = {}
    for f in files:
        col_of: dict[str, int] = {}
        for line in f.read_text(encoding="utf-8", errors="ignore").splitlines():
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                cols = s.lstrip("#").split()
                if cols and cols[0] == "Time":
                    col_of = {}
                    for i, c in enumerate(cols):
                        m = re.match(r"\w+\((.+)\)$", c)   # volAverage(TTR) -> TTR
                        col_of[m.group(1) if m else c] = i
                continue
            parts = s.split()
            try:
                t = float(parts[0])
            except ValueError:
                continue
            row = rows.setdefault(t, {})
            for fld in fields:
                i = col_of.get(fld)
                if i is not None and i < len(parts):
                    row[fld] = float(parts[i])
    return rows


# --------------------------------------------------------------------------- #
# Source 2: time directories (fallback)
# --------------------------------------------------------------------------- #
def parse_internal_field(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(rf"internalField\s+uniform\s+({NUM})\s*;", text)
    if m:
        return float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;",
                  text, flags=re.DOTALL)
    if m:
        vals = [float(v) for v in re.findall(NUM, m.group(2))]
        if len(vals) != int(m.group(1)):
            raise ValueError(f"{path}: expected {m.group(1)} values, found {len(vals)}")
        return sum(vals)/len(vals)
    raise ValueError(f"{path}: cannot parse internalField (macro or unsupported format)")


def read_timedir_series(case: Path, fields: list[str]) -> dict[float, dict[str, float]]:
    rows: dict[float, dict[str, float]] = {}
    for d in case.iterdir():
        if not (d.is_dir() and TIME_RE.match(d.name)):
            continue
        t = float(d.name)
        if t == 0.0:
            continue  # t = 0 comes from initialConditions (0/ fields hold $macros)
        for fld in fields:
            p = d / fld
            if p.exists():
                try:
                    rows.setdefault(t, {})[fld] = parse_internal_field(p)
                except ValueError:
                    pass
    return rows


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=".", help="Case directory (default: .)")
    ap.add_argument("--fields", nargs="+", default=["TTR", "TVib"],
                    help="Fields to extract (default: TTR TVib)")
    ap.add_argument("--fo", default=None,
                    help="Function-object name under postProcessing/ (default: first found, e.g. hbState)")
    ap.add_argument("--source", choices=["auto", "fo", "timedirs"], default="auto",
                    help="auto: function object if present, else time directories")
    ap.add_argument("--initial-key", action="append", default=[], metavar="FIELD=KEY",
                    help="Extra/override initialConditions key, e.g. Y_N2=N2_mass_fraction")
    ap.add_argument("--output", default="heatBath.csv", help="Merged CSV (default: heatBath.csv)")
    ap.add_argument("--plot", default=None, help="Optional quick-look PNG (log time axis)")
    ap.add_argument("--title", default="Heat bath", help="Quick-look title")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    case = Path(args.case)
    if not case.is_dir():
        raise SystemExit(f"Case directory not found: {case}")

    overrides = {}
    for kv in args.initial_key:
        fld, _, key = kv.partition("=")
        if not key:
            raise SystemExit(f"--initial-key expects FIELD=KEY, got '{kv}'")
        overrides[fld] = key

    rows: dict[float, dict[str, float]] = {}
    source = "none"
    if args.source in ("auto", "fo"):
        files = find_fo_files(case, args.fo)
        if files:
            rows, source = read_fo_series(files, args.fields), f"function object ({len(files)} file(s))"
    if not rows and args.source in ("auto", "timedirs"):
        rows, source = read_timedir_series(case, args.fields), "time directories"
    if not rows:
        raise SystemExit("No data found (no function-object output and no readable time directories).")

    # t = 0 from initialConditions, field by field
    init = case / "0" / "initialConditions"
    rows[0.0] = {f: read_initial_value(init, initial_key_for(f, overrides)) for f in args.fields}

    times = sorted(rows)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Time (s)", *args.fields])
        for t in times:
            w.writerow([repr(t), *(repr(rows[t].get(f, math.nan)) for f in args.fields)])

    missing = [f for f in args.fields if all(math.isnan(rows[t].get(f, math.nan)) for t in times if t > 0)]
    print(f"Source: {source}; {len(times)} rows; t = {times[1] if len(times) > 1 else 0:g} ... {times[-1]:g} s")
    print(f"t = 0: " + ", ".join(f"{f} = {rows[0.0][f]:g}" for f in args.fields))
    if missing:
        print(f"WARNING: no data for {missing} (not in the function object / time directories)")
    print(f"Saved CSV to: {out}")

    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))
        for f in args.fields:
            pts = [(t, rows[t][f]) for t in times
                   if t > 0 and f in rows[t] and math.isfinite(rows[t][f])]
            if pts:
                ax.plot(*zip(*pts), linewidth=1.8, label=f)
        ax.set_xscale("log")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature (K)")
        ax.set_title(args.title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=300)
        print(f"Saved figure to: {args.plot}")


if __name__ == "__main__":
    main()
