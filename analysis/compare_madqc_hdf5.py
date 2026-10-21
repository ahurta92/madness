"""Inspect native MADNESS HDF5 archives and compare structured pymra fields."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def classify_hdf5(path: Path) -> str:
    path = Path(path)
    if not path.exists():
        raise ValueError(f"HDF5 path does not exist: {path}")
    try:
        import h5py
    except ImportError as exc:
        raise ValueError(f"h5py unavailable: {exc}") from exc
    with h5py.File(path, "r") as handle:
        keys = set(handle.keys())
    if "madness_archive" in keys:
        return "native_madness_archive"
    if {"coeffs", "keys", "meta"}.issubset(keys):
        return "pymra_structured"
    return "unknown_hdf5"


def compare_fields(coarse_path: Path, fine_path: Path, *, ngrid: int = 64,
                   cell: tuple[float, float] = (-20.0, 20.0)) -> dict:
    if classify_hdf5(coarse_path) != "pymra_structured":
        return {"available": False, "reason": "coarse archive is not a pymra structured export"}
    if classify_hdf5(fine_path) != "pymra_structured":
        return {"available": False, "reason": "fine archive is not a pymra structured export"}
    try:
        from pymra import eval_points, read_function
    except ImportError as exc:
        return {"available": False, "reason": f"pymra unavailable: {exc}"}
    lo, hi = cell
    axis = np.linspace(lo + (hi-lo)/(2*ngrid), hi - (hi-lo)/(2*ngrid), ngrid)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    coarse = np.asarray(eval_points(read_function(str(coarse_path)), points), dtype=float)
    fine = np.asarray(eval_points(read_function(str(fine_path)), points), dtype=float)
    # The exported H2 orbital is a restricted alpha orbital; rho=2*psi^2 is a
    # ground-state density proxy, not a response-state density.
    coarse_rho, fine_rho = 2*coarse**2, 2*fine**2
    delta = coarse_rho - fine_rho
    dv = ((hi-lo)/ngrid)**3
    radius = np.sqrt(np.sum(points*points, axis=1))
    inner = radius < 3.0
    return {
        "available": True,
        "field_construction": "closed-shell ground-state proxy rho=2*psi_alpha^2",
        "grid": {"n_per_axis": ngrid, "points": int(ngrid**3), "cell": [lo, hi]},
        "coarse_density_l2": float(np.sqrt(np.sum(coarse_rho**2)*dv)),
        "fine_density_l2": float(np.sqrt(np.sum(fine_rho**2)*dv)),
        "density_difference_l2": float(np.sqrt(np.sum(delta**2)*dv)),
        "density_difference_max_abs": float(np.max(np.abs(delta))),
        "density_difference_pointwise_rms": float(np.sqrt(np.mean(delta**2))),
        "inner_r_lt_3_l2": float(np.sqrt(np.sum(delta[inner]**2)*dv)),
        "outer_r_ge_3_l2": float(np.sqrt(np.sum(delta[~inner]**2)*dv)),
        "limitations": [
            "sampled common-grid comparison, not an exact MRA inner product",
            "ground-state orbital-derived density proxy, not response density",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--grid", type=int, default=64)
    args = parser.parse_args()
    result = compare_fields(args.coarse, args.fine, ngrid=args.grid)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    return 0 if result.get("available", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
