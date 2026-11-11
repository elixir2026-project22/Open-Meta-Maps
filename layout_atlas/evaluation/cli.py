"""CLI: score Open Meta Maps Escher output with the unified evaluation path.

  python3 project.py evaluate --gem MODEL.xml --map DIR/map.json [--out DIR/evaluation.json]
  python3 project.py evaluate --gem MODEL.xml --map DIR/map.json --reference REF.sbgn
  python3 project.py evaluate --batch MANIFEST.json --out results/evaluation/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

# layout_atlas is on sys.path when launched via project.py / cli.py
from evaluation.score import score_paths


def _dump(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str))


def _flat_row(report: dict) -> dict:
    b = report.get("benchmark") or {}
    g = (report.get("geometry") or {}).get("core") or {}
    cross = g.get("crossings") or {}
    coll = g.get("collisions") or {}
    return {
        "model": report.get("model"),
        "placement_share": (report.get("coverage") or {}).get("placement_share"),
        "runtime_seconds": report.get("runtime_seconds"),
        "crossings": b.get("crossings"),
        "crossing_density": b.get("crossing_density"),
        "overlapping_pairs": b.get("overlapping_pairs"),
        "overlap_area_share": b.get("overlap_area_share"),
        "participant_consistency": b.get("participant_consistency"),
        "species_in_other_compartment": b.get("species_in_other_compartment"),
        "long_edge_share": b.get("long_edge_share"),
        "edge_length_deviation": b.get("edge_length_deviation"),
        "aspect_ratio": b.get("aspect_ratio") if b.get("aspect_ratio") is not None else b.get("aspect_ratio_bbox"),
        "area_per_node": b.get("area_per_node") if b.get("area_per_node") is not None else b.get("area_per_node_bbox"),
        "stress": b.get("stress"),
        "procrustes": b.get("procrustes"),
        "distance_correlation": b.get("distance_correlation"),
        "knn_preservation": b.get("knn_preservation"),
        "geom_crossings_total": cross.get("total"),
        "geom_node_overlap_pairs": coll.get("node_overlap_pairs"),
        "geom_label_overlap_pairs": coll.get("label_overlap_pairs"),
        "gdmetrix": (report.get("dependencies") or {}).get("gdmetrix"),
    }


def run_one(args) -> int:
    if not args.gem or not args.map:
        print("evaluate needs --gem and --map", file=sys.stderr)
        return 2
    report = score_paths(
        args.gem,
        args.map,
        reference=args.reference,
        compare_map=args.compare_map,
        include_geometry=not args.no_geometry,
    )
    out = Path(args.out) if args.out else Path(args.map).with_name("evaluation.json")
    _dump(out, report)
    row = _flat_row(report)
    print(json.dumps(row, indent=2))
    print(f"wrote {out}")
    return 0


def run_batch(args) -> int:
    """MANIFEST.json: list of {name, gem, map, reference?} or a directory of compose results."""
    manifest_path = Path(args.batch)
    out_dir = Path(args.out or "results/evaluation")
    out_dir.mkdir(parents=True, exist_ok=True)
    if manifest_path.is_dir():
        jobs = []
        for map_json in sorted(manifest_path.glob("*/map.json")):
            # gem may be recorded in report.json
            gem = None
            report = map_json.with_name("report.json")
            if report.exists():
                rj = json.loads(report.read_text())
                gem = (rj.get("model") or {}).get("path") or rj.get("gem")
            upload = map_json.parent / "upload"
            if gem is None and upload.exists():
                files = list(upload.glob("model.*"))
                gem = str(files[0]) if files else None
            if gem:
                jobs.append(dict(name=map_json.parent.name, gem=gem, map=str(map_json)))
    else:
        jobs = json.loads(manifest_path.read_text())
    rows = []
    reports = []
    for job in jobs:
        name = job.get("name") or Path(job["gem"]).stem
        print(f"== {name}", flush=True)
        try:
            report = score_paths(
                job["gem"],
                job["map"],
                reference=job.get("reference"),
                compare_map=job.get("compare_map"),
                include_geometry=not args.no_geometry,
            )
        except Exception as ex:  # noqa: BLE001
            report = dict(model=name, error=f"{type(ex).__name__}: {ex}")
            print(f"  FAILED: {report['error']}", file=sys.stderr)
        _dump(out_dir / f"{name}.json", report)
        reports.append(report)
        if "error" not in report:
            rows.append(_flat_row(report))
    if rows:
        fields = list(rows[0])
        with (out_dir / "summary.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"summary: {out_dir / 'summary.csv'}")
    _dump(out_dir / "all.json", reports)
    return 0 if all("error" not in r for r in reports) else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Evaluate Open Meta Maps Escher layouts (Project 22 metrics).")
    p.add_argument("--gem", help="SBML / model file")
    p.add_argument("--map", help="Escher map.json")
    p.add_argument("--reference", help="Optional reference: SBGN, Escher JSON, or {id:[x,y]} JSON")
    p.add_argument("--compare-map", help="Optional second map.json for stability / reproducibility")
    p.add_argument("--batch", help="Manifest JSON or a folder of compose results")
    p.add_argument("--out", help="Output JSON path (single) or directory (batch)")
    p.add_argument("--no-geometry", action="store_true", help="Skip layout_quality final-geometry block")
    args = p.parse_args(argv)
    if args.batch:
        return run_batch(args)
    return run_one(args)


if __name__ == "__main__":
    raise SystemExit(main())
