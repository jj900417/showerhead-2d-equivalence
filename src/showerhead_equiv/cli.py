"""Batch command line: geometry CSV (+ reactions CSV) -> equivalent 2D table.

    showerhead-equiv examples/geometry_example.csv \
        --reactions examples/reactions_example.csv -o result.csv

Geometry columns: label, r_mm, N, L_mm, H_mm [, ratio_override]
Reaction columns: name, gamma
"""

from __future__ import annotations

import argparse
import csv
import sys

from .core import Ring, equivalent_2d


def read_rings(path: str) -> list[Ring]:
    with open(path, newline="", encoding="utf-8") as f:
        rings = []
        for i, row in enumerate(csv.DictReader(f), start=1):
            override = (row.get("ratio_override") or "").strip()
            rings.append(
                Ring(
                    label=(row.get("label") or "").strip() or f"row {i}",
                    r=float(row["r_mm"]),
                    N=int(row["N"]),
                    L=float(row["L_mm"]),
                    H=float(row["H_mm"]),
                    ratio_override=float(override) if override else None,
                )
            )
    return rings


def read_reactions(path: str) -> list[tuple[str, float]]:
    with open(path, newline="", encoding="utf-8") as f:
        return [
            ((row.get("name") or "").strip() or f"reaction {j}", float(row["gamma"]))
            for j, row in enumerate(csv.DictReader(f), start=1)
        ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="showerhead-equiv", description=__doc__.splitlines()[0])
    p.add_argument("geometry", help="geometry CSV, one row per hole ring")
    p.add_argument("--reactions", help="reactions CSV (name, gamma)")
    p.add_argument("-o", "--output", help="output CSV (default: stdout)")
    args = p.parse_args(argv)

    rings = read_rings(args.geometry)
    reactions = read_reactions(args.reactions) if args.reactions else []

    header = ["label", "R_prime_mm", "R_o_mm", "R_i_mm", "wall_area_ratio_3D_to_2D", "residual_rel"]
    header += [f"gamma_eff[{name}]" for name, _ in reactions]

    out = open(args.output, "w", newline="", encoding="utf-8") if args.output else sys.stdout
    try:
        w = csv.writer(out)
        w.writerow(header)
        for ring in rings:
            eq = equivalent_2d(ring)
            w.writerow(
                [ring.label, f"{eq.r_prime:.8g}", f"{eq.ro:.8g}", f"{eq.ri:.8g}",
                 f"{eq.area_ratio:.8g}", f"{eq.residual_rel:.3e}"]
                + [f"{eq.gamma_eff(g):.6g}" for _, g in reactions]
            )
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
