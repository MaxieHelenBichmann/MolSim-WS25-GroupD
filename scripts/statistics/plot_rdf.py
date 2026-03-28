#!/usr/bin/env python3
"""
Plot RDF/density curves from a CSV with rows:
<iteration>,<distance>,<density>

Groups by iteration and plots selected iterations as separate curves.
Legend labels can show:
t = iteration * delta_t + start_time

Use command (matplotlib required):
python plot_rdf.py <FILE>.csv --delta-t <DELTA_T> --start-time <START_TIME> --out <FILE>.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import matplotlib.pyplot as plt


Row = Tuple[float, float, float]  # (iteration, distance, density)


def _try_parse_row(row: List[str]) -> Optional[Row]:
    if len(row) < 3:
        return None
    a, b, c = row[0].strip(), row[1].strip(), row[2].strip()
    if not a or not b or not c:
        return None
    try:
        return float(a), float(b), float(c)
    except ValueError:
        return None


def read_iter_dist_density_csv(path: Path) -> Dict[float, List[Tuple[float, float]]]:
    """
    Returns dict: iteration -> list of (distance, density)
    Skips empty lines, comments starting with '#', and non-numeric header rows.
    """
    data: Dict[float, List[Tuple[float, float]]] = {}

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

            parsed = _try_parse_row(row)
            if parsed is None:
                # header or malformed; skip
                continue

            it, dist, dens = parsed
            data.setdefault(it, []).append((dist, dens))

    if not data:
        raise ValueError(
            f"No numeric data found in '{path}'. Expected rows like: <iteration>,<distance>,<density>"
        )

    return data


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def main() -> None:
    p = argparse.ArgumentParser(description="Plot RDF/density curves from CSV.")
    p.add_argument(
        "csv", type=Path, help="Input CSV file (<iteration>,<distance>,<density>)"
    )
    p.add_argument(
        "--delta-t", type=float, default=1.0, help="Δt for time labels (default: 1.0)"
    )
    p.add_argument(
        "--start-time",
        type=float,
        default=0.0,
        help="Start time offset t0",
    )
    p.add_argument(
        "--x",
        choices=["time", "iteration"],
        default="time",
        help="Legend labels show computed time or iteration (default: time)",
    )

    p.add_argument(
        "--iterations",
        type=str,
        default=None,
        help="Comma-separated iterations to plot (default: all)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output image path (png/pdf). If omitted, shows window.",
    )
    p.add_argument(
        "--dpi", type=int, default=200, help="DPI for saved images (default: 200)"
    )
    args = p.parse_args()

    data = read_iter_dist_density_csv(args.csv)
    all_iterations = sorted(data.keys())

    if args.iterations:
        selected = parse_float_list(args.iterations)
        missing = [it for it in selected if it not in data]
        if missing:
            raise ValueError(f"These iterations are not present in the CSV: {missing}")
        selected_iterations = selected
    else:
        selected_iterations = all_iterations

    plt.figure(figsize=(9, 4.8))

    for it in selected_iterations:
        pts = sorted(data[it], key=lambda x: x[0])
        dists = [d for d, _ in pts]
        dens = [y for _, y in pts]

        if args.x == "time":
            t = it * args.delta_t + args.start_time
            lbl = f"Time {t:g}"
        else:
            lbl = f"Iter {it:g}"

        plt.plot(dists, dens, linewidth=1.5, label=lbl)

    plt.title("Radial Distribution Function")
    plt.xlabel(
        f"Distance\n RDF with i ∈ [0,{max(dists):.2f}] and δr = {(dists[1]-dists[0]):.2f}"
    )
    plt.ylabel("Local Density")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    plt.legend()
    plt.ylim(bottom=0)

    plt.tight_layout()

    if args.out is None:
        plt.show()
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved plot to: {args.out}")


if __name__ == "__main__":
    main()
