#!/usr/bin/env python3
"""
Plot TTR vs time from OpenFOAM timestep folders and export CSV.

Supports:
- normal field files like 1e-6/TTR, 2e-6/TTR, ...
- special initial value from 0/initialConditions using TTRValue

CSV output format:
    Time (s),Temperature_tr

Example:
    python3 plot_ttr_vs_time.py \
        --case "nexaFoam" . \
        --field TTR \
        --title "N2-N V-T Relaxation" \
        --output n2n_vt_ttr.png \
        --export-csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Iterable


TIME_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")


def is_time_folder(path: Path) -> bool:
    return path.is_dir() and TIME_RE.match(path.name) is not None


def parse_foam_scalar_field(field_path: Path) -> float:
    text = field_path.read_text(encoding="utf-8", errors="ignore")

    m_uniform = re.search(
        r"internalField\s+uniform\s+([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*;",
        text,
    )
    if m_uniform:
        return float(m_uniform.group(1))

    m_nonuniform = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if m_nonuniform:
        n_expected = int(m_nonuniform.group(1))
        body = m_nonuniform.group(2)

        values = [
            float(v)
            for v in re.findall(
                r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?",
                body,
            )
        ]

        if not values:
            raise ValueError(f"No scalar values found in nonuniform field: {field_path}")

        if len(values) != n_expected:
            raise ValueError(
                f"Expected {n_expected} values in {field_path}, found {len(values)}"
            )

        return sum(values) / len(values)

    raise ValueError(
        f"Could not parse internalField from {field_path}. "
        "Only uniform and nonuniform List<scalar> are supported."
    )


def parse_initial_conditions_value(initial_conditions_path: Path, key: str) -> float:
    text = initial_conditions_path.read_text(encoding="utf-8", errors="ignore")

    m = re.search(
        rf"\b{re.escape(key)}\b\s+([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*;",
        text,
    )
    if not m:
        raise ValueError(
            f"Could not find key '{key}' in {initial_conditions_path}"
        )

    return float(m.group(1))


def read_case_series(
    case_dir: Path,
    field_name: str,
    initial_dict_name: str = "initialConditions",
    initial_key: str = "TTRValue",
) -> tuple[list[float], list[float]]:
    if not case_dir.is_dir():
        raise FileNotFoundError(f"Case directory not found: {case_dir}")

    time_dirs = sorted(
        (p for p in case_dir.iterdir() if is_time_folder(p)),
        key=lambda p: float(p.name),
    )

    times: list[float] = []
    values: list[float] = []

    for tdir in time_dirs:
        time_value = float(tdir.name)
        field_path = tdir / field_name

        if field_path.exists():
            try:
                field_value = parse_foam_scalar_field(field_path)
                times.append(time_value)
                values.append(field_value)
                continue
            except ValueError:
                pass

        if abs(time_value) < 1e-30:
            initial_path = tdir / initial_dict_name
            if initial_path.exists():
                field_value = parse_initial_conditions_value(initial_path, initial_key)
                times.append(time_value)
                values.append(field_value)
                continue

    if not times:
        raise FileNotFoundError(
            f"No readable '{field_name}' files found in numeric time folders under {case_dir}, "
            f"and no '{initial_dict_name}' with key '{initial_key}' was usable at time 0."
        )

    paired = sorted(zip(times, values), key=lambda x: x[0])
    times = [t for t, _ in paired]
    values = [v for _, v in paired]
    return times, values


def read_csv_series(csv_path: Path) -> tuple[list[float], list[float]]:
    times: list[float] = []
    values: list[float] = []

    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        raise ValueError(f"CSV is empty: {csv_path}")

    start_idx = 0
    try:
        float(rows[0][0])
        float(rows[0][1])
    except Exception:
        start_idx = 1

    for row in rows[start_idx:]:
        if len(row) < 2:
            continue
        times.append(float(row[0]))
        values.append(float(row[1]))

    if not times:
        raise ValueError(f"No numeric data found in CSV: {csv_path}")

    return times, values


def positive_only(times: Iterable[float], values: Iterable[float]) -> tuple[list[float], list[float]]:
    t_out: list[float] = []
    v_out: list[float] = []
    for t, v in zip(times, values):
        if t > 0.0 and math.isfinite(t) and math.isfinite(v):
            t_out.append(t)
            v_out.append(v)
    return t_out, v_out


def sanitize_label(label: str) -> str:
    s = re.sub(r"\s+", "_", label.strip())
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    return s or "case"


def write_case_csv(csv_path: Path, times: list[float], values: list[float]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time (s)", "Temperature_tr"])
        for t, v in zip(times, values):
            writer.writerow([t, v])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot TTR vs time from OpenFOAM timestep folders and export CSV."
    )

    parser.add_argument(
        "--case",
        nargs=2,
        action="append",
        metavar=("LABEL", "CASE_DIR"),
        required=True,
        help="OpenFOAM case to plot. Can be repeated.",
    )

    parser.add_argument(
        "--csv",
        nargs=2,
        action="append",
        metavar=("LABEL", "CSV_FILE"),
        default=[],
        help="Optional reference CSV with columns: time,value. Can be repeated.",
    )

    parser.add_argument(
        "--field",
        default="TTR",
        help="Field file name inside each time directory. Default: TTR",
    )

    parser.add_argument(
        "--initial-dict",
        default="initialConditions",
        help="Dictionary name in 0/ used for initial value. Default: initialConditions",
    )

    parser.add_argument(
        "--initial-key",
        default="TTRValue",
        help="Key inside initial dictionary for TTR initial value. Default: TTRValue",
    )

    parser.add_argument(
        "--title",
        default="Air-5 Reacting Flow",
        help="Plot title.",
    )

    parser.add_argument(
        "--ylabel",
        default="Temperature (K)",
        help="Y-axis label.",
    )

    parser.add_argument(
        "--xlabel",
        default="Time (s)",
        help="X-axis label.",
    )

    parser.add_argument(
        "--output",
        default="TTR_vs_time.png",
        help="Output figure path.",
    )

    parser.add_argument(
        "--export-csv",
        action="store_true",
        help="Export one CSV per case.",
    )

    parser.add_argument(
        "--csv-out-dir",
        default=".",
        help="Directory for exported case CSV files. Default: current directory",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plot window in addition to saving.",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))

    csv_out_dir = Path(args.csv_out_dir)

    for label, case_dir_str in args.case:
        case_dir = Path(case_dir_str)
        times, values = read_case_series(
            case_dir,
            args.field,
            initial_dict_name=args.initial_dict,
            initial_key=args.initial_key,
        )

        if args.export_csv:
            csv_name = f"{sanitize_label(label)}.csv"
            csv_path = csv_out_dir / csv_name
            write_case_csv(csv_path, times, values)
            print(f"Saved CSV to: {csv_path}")

        plot_times, plot_values = positive_only(times, values)
        plt.plot(plot_times, plot_values, label=label, linewidth=2)

    for label, csv_path_str in args.csv:
        csv_path = Path(csv_path_str)
        times, values = read_csv_series(csv_path)
        plot_times, plot_values = positive_only(times, values)
        plt.plot(plot_times, plot_values, label=label)

    plt.xscale("log")
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.title(args.title)
    plt.grid(True, which="both", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)

    print(f"Saved figure to: {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
