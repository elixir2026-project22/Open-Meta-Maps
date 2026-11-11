"""Reproducibility across seeds and stability when a few reactions are removed.

These compare two position maps on shared nodes (no layout algorithm changes here).
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import procrustes

Positions = dict[str, tuple[float, float]]


def position_displacement(a: Positions, b: Positions) -> dict[str, float | None]:
    """After Procrustes alignment: mean / median / p90 residual distance on common nodes."""
    common = sorted(set(a) & set(b))
    if len(common) < 3:
        return {
            "common_nodes": float(len(common)),
            "mean_displacement": None,
            "median_displacement": None,
            "p90_displacement": None,
            "procrustes_disparity": None,
        }
    xa = np.array([a[n] for n in common], float)
    xb = np.array([b[n] for n in common], float)
    _, xb_fit, disparity = procrustes(xa, xb)
    # procrustes returns normalised matrices; report disparity and raw centroid-aligned residuals
    xa0 = xa - xa.mean(0)
    xb0 = xb - xb.mean(0)
    # best rotation (no scale) via Kabsch
    h = xb0.T @ xa0
    u, _, vt = np.linalg.svd(h)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    resid = np.linalg.norm(xa0 - xb0 @ r, axis=1)
    return {
        "common_nodes": float(len(common)),
        "mean_displacement": float(resid.mean()),
        "median_displacement": float(np.median(resid)),
        "p90_displacement": float(np.percentile(resid, 90)),
        "procrustes_disparity": float(disparity),
    }


def reproducibility(maps: list[Positions]) -> dict[str, float | None]:
    """Pairwise displacement summary across several seed layouts."""
    if len(maps) < 2:
        return {"pairs": 0.0, "mean_displacement_mean": None, "procrustes_mean": None}
    rows = []
    for i in range(len(maps)):
        for j in range(i + 1, len(maps)):
            rows.append(position_displacement(maps[i], maps[j]))
    means = [r["mean_displacement"] for r in rows if r["mean_displacement"] is not None]
    procs = [r["procrustes_disparity"] for r in rows if r["procrustes_disparity"] is not None]
    return {
        "pairs": float(len(rows)),
        "mean_displacement_mean": float(np.mean(means)) if means else None,
        "procrustes_mean": float(np.mean(procs)) if procs else None,
    }
