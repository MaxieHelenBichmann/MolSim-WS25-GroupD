#!/usr/bin/env python3
"""
Plot a diffusion diagram Var(t) from a CSV file with rows:
<iteration>,<var>

time is computed as:
t = iteration * delta_t + start_time

Use command (matplotlib required):
python ./plot_diffusion.py ./<FILE>.csv --delta-t <DELTA_T> --start-time <START_TIME> --out <FILE>.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt


def read_iteration_var_csv(path: Path) -> Tuple[List[float], List[float]]:
    """
    Reads the first two columns from a CSV.
    - Skips empty lines
    - Skips comment lines starting with '#'
    - If the first non-empty line is non-numeric, treats it as a header and skips it
    """
    iterations: List[float] = []
    vars_: List[float] = []

    def try_parse_row(row: List[str]) -> Optional[Tuple[float, float]]:
        if len(row) < 2:
            return None
        a, b = row[0].strip(), row[1].strip()
        if not a or not b:
            return None
        try:
            return float(a), float(b)
        except ValueError:
            return None

    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            joined = "".join(row).strip()
            if not joined:
                continue
            if joined.lstrip().startswith("#"):
                continue

            parsed = try_parse_row(row)
            if parsed is None:
                # likely a header; skip and continue
                continue

            it, v = parsed
            iterations.append(it)
            vars_.append(v)

    if not iterations:
        raise ValueError(
            f"No numeric data found in '{path}'. Expected rows like: <iteration>,<var>"
        )

    return iterations, vars_


def main() -> None:
    p = argparse.ArgumentParser(description="Plot diffusion diagram Var(t) from CSV.")
    p.add_argument("csv", type=Path, help="Input CSV file (<iteration>,<var>)")
    p.add_argument("--delta-t", type=float, required=True, help="Time step size Δt")
    p.add_argument("--start-time", type=float, default=0.0, help="Start time offset t0")
    p.add_argument(
        "--x",
        choices=["time", "iteration"],
        default="time",
        help="Use computed time or raw iteration on the x-axis (default: time)",
    )
    p.add_argument(
        "--title", default="Diffusion Var(t)", help="Plot title (default: Diffusion Var(t))"
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output image path (e.g., plot.png or plot.pdf). If omitted, shows a window.",
    )
    p.add_argument("--dpi", type=int, default=200, help="DPI for saved images (default: 200)")
    args = p.parse_args()

    iterations, var_values = read_iteration_var_csv(args.csv)

    if args.x == "time":
        x_vals = [it * args.delta_t + args.start_time for it in iterations]
        x_label = "Time"
    else:
        x_vals = iterations
        x_label = "Time Step (Iteration)"

    plt.figure(figsize=(9, 4.8))
    plt.plot(x_vals, var_values, linewidth=1.5)
    plt.ylim(bottom=0)
    plt.title(args.title)
    plt.xlabel(x_label)
    plt.ylabel("Diffusion Var(t)")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    plt.tight_layout()

    if args.out is None:
        plt.show()
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved plot to: {args.out}")


if __name__ == "__main__":
    main()