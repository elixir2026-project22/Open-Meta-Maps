"""EC-off vs EC-on compose comparison for the three example models (iML1515 skipped).

  python3 layout_atlas/ec_extension_eval.py --run     # compose both settings
  python3 layout_atlas/ec_extension_eval.py --summarize DIR

Does not change the default; writes a JSON summary for the recommendation note.
"""
import argparse
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import *

MODELS = ['e_coli_core', 'iIT341', 'iMM904']


def _lq(report):
    fg = ((report.get('layout_quality') or {}).get('final_geometry') or {}).get('core') or {}
    crossings = (fg.get('crossings') or {}).get('total')
    col = fg.get('collisions') or {}
    length = fg.get('length_p90') or {}
    return dict(
        crossings=crossings,
        node_overlap_pairs=col.get('node_overlap_pairs'),
        label_overlap_pairs=col.get('label_overlap_pairs'),
        curves_through_nodes=col.get('curves_through_nodes'),
        length_p90_px=length.get('px'),
        length_p90_spacings=length.get('in_reaction_spacings'),
    )


def _evidence(report):
    ev = report.get('evidence') or {}
    return dict(
        layout_provenance=ev.get('layout_provenance') or {},
        identity_evidence=ev.get('identity_evidence') or {},
        transfer_mode=ev.get('transfer_mode') or {},
        chemical_compatibility=ev.get('chemical_compatibility') or {},
    )


def summarize_job(folder):
    folder = Path(folder)
    report = read(folder / 'report.json')
    prov = read(folder / 'reaction_provenance.json') if (folder / 'reaction_provenance.json').exists() else {}
    # selected structures by match_level / source
    mods = read(folder / 'modules.json') if (folder / 'modules.json').exists() else []
    levels = Counter(m.get('match_level') for m in mods)
    sources = Counter(m.get('source') for m in mods)
    # ec_topology reaction sample
    ec_rxns = []
    for rid, rec in (prov.items() if isinstance(prov, dict) else []):
        if not isinstance(rec, dict):
            continue
        iden = rec.get('identity_evidence') or {}
        st = rec.get('structure') or {}
        if st.get('match_level') == 'ec_topology' or iden.get('kind') == 'ec_only':
            ec_rxns.append(dict(
                reaction=rid,
                provenance=rec.get('provenance'),
                layout_provenance=rec.get('layout_provenance'),
                shared_ec=(iden.get('shared_ec') or []),
                identity_kind=iden.get('kind'),
                chemical=(rec.get('chemical_compatibility') or {}).get('status'),
                difference_types=(rec.get('chemical_compatibility') or {}).get('difference_types') or [],
                source=(st.get('source') or (rec.get('transfer') or {}).get('source_kind')),
                source_reaction=(rec.get('transfer') or {}).get('source_reaction'),
            ))
    return dict(
        model=report.get('model', {}).get('id') or folder.name,
        reactions=report.get('model', {}).get('reactions') or report.get('covered', {}).get('total'),
        seconds=report.get('seconds'),
        covered=report.get('covered'),
        provenance=report.get('provenance'),
        levels=report.get('levels'),
        modules=report.get('modules'),
        module_match_levels=dict(levels),
        module_sources=dict(sources),
        ec_screening=report.get('ec_screening') or {},
        evidence=_evidence(report),
        layout_quality=_lq(report),
        ec_topology_reactions=ec_rxns,
        n_ec_topology=len(ec_rxns),
    )


def compose(gem, out, ec):
    cmd = [sys.executable, str(ROOT / 'project.py'), 'compose', '--gem', str(gem), '--out', str(out),
           '--source-paths', '--layout-style', 'expert']
    if ec:
        cmd.append('--ec')
    print('RUN', ' '.join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode:
        raise SystemExit(r.returncode)


def run_all(out_root):
    out_root = Path(out_root)
    rows = {}
    for name in MODELS:
        gem = EXAMPLES / f'{name}.xml'
        for ec, tag in [(False, 'ec_off'), (True, 'ec_on')]:
            dest = out_root / f'{name}_{tag}'
            compose(gem, dest, ec=ec)
            rows[f'{name}_{tag}'] = summarize_job(dest)
    dump(out_root / 'ec_extension_summary.json', dict(models=MODELS, jobs=rows))
    return rows


def sample_suspicious(ec_rxns, limit=24):
    """Flag likely-wrong EC placements: transport/exchange ids, empty/promiscuous EC patterns, conflicts."""
    suspicious, ok = [], []
    for r in ec_rxns:
        rid = r['reaction']
        reasons = []
        if rid.startswith(('EX_', 'DM_', 'SK_', 'sink_')) or 't2' in rid.lower() or rid.endswith('tex') or 'abc' in rid.lower():
            reasons.append('transport_or_exchange_id')
        if r.get('chemical') == 'conflict':
            reasons.append('equation_conflict')
        ecs = r.get('shared_ec') or []
        if any(e.endswith('.-') or e.count('.') < 3 for e in ecs):
            reasons.append('incomplete_ec')
        # very broad dehydrogenase / transferase classes often promiscuous when alone
        if len(ecs) == 1 and ecs[0].startswith(('1.1.1.', '2.7.1.', '3.6.3.', '7.')):
            if ecs[0] in {'1.1.1.1', '2.7.1.40', '3.6.3.17', '7.1.2.2'}:
                reasons.append('broad_or_transporter_ec')
        if reasons:
            suspicious.append(dict(r, flags=reasons))
        else:
            ok.append(r)
    # prefer diverse sample
    suspicious.sort(key=lambda x: (x['reaction'], str(x.get('shared_ec'))))
    ok.sort(key=lambda x: x['reaction'])
    return dict(suspicious=suspicious[:limit], ok_sample=ok[:limit], n_suspicious=len(suspicious), n_ok=len(ok))


def comparison_table(jobs):
    rows = []
    for name in MODELS:
        off, on = jobs.get(f'{name}_ec_off'), jobs.get(f'{name}_ec_on')
        if not off or not on:
            continue
        ev_off, ev_on = off['evidence'], on['evidence']
        scr = on.get('ec_screening') or {}
        kept = sum(v.get('proposals_kept', 0) for v in scr.values() if isinstance(v, dict))
        rejected = sum(v.get('proposals_rejected_incompatible', 0) for v in scr.values() if isinstance(v, dict))
        rows.append(dict(
            model=name,
            off=dict(modules=off['modules'], covered=off['covered'], provenance=off['provenance'],
                     layout=off['layout_quality'], seconds=off['seconds'],
                     exact=(ev_off['layout_provenance'] or {}).get('expert_exact', 0),
                     adapted=(ev_off['layout_provenance'] or {}).get('expert_adapted', 0),
                     topological=(ev_off['layout_provenance'] or {}).get('expert_topological', 0)),
            on=dict(modules=on['modules'], covered=on['covered'], provenance=on['provenance'],
                    layout=on['layout_quality'], seconds=on['seconds'],
                    exact=(ev_on['layout_provenance'] or {}).get('expert_exact', 0),
                    adapted=(ev_on['layout_provenance'] or {}).get('expert_adapted', 0),
                    topological=(ev_on['layout_provenance'] or {}).get('expert_topological', 0),
                    functional=(ev_on['layout_provenance'] or {}).get('expert_functional', 0),
                    ec_topology_reactions=on['n_ec_topology'],
                    ec_proposals_kept=kept, ec_proposals_rejected=rejected,
                    ec_screening=scr),
            suspicious=sample_suspicious(on['ec_topology_reactions']),
        ))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', action='store_true')
    p.add_argument('--out', default=str(ROOT / 'results/ec_extension_eval'))
    p.add_argument('--summarize', help='existing output directory')
    a = p.parse_args()
    out = Path(a.summarize or a.out)
    if a.run:
        run_all(out)
    jobs = {}
    summary_path = out / 'ec_extension_summary.json'
    if summary_path.exists():
        jobs = read(summary_path).get('jobs', {})
    else:
        for name in MODELS:
            for tag in ('ec_off', 'ec_on'):
                folder = out / f'{name}_{tag}'
                if (folder / 'report.json').exists():
                    jobs[f'{name}_{tag}'] = summarize_job(folder)
        dump(summary_path, dict(models=MODELS, jobs=jobs))
    table = comparison_table(jobs)
    dump(out / 'ec_extension_comparison.json', dict(comparison=table))
    for row in table:
        print(f"\n=== {row['model']} ===")
        o, n = row['off'], row['on']
        print(f"  modules {o['modules']} -> {n['modules']}; seconds {o['seconds']} -> {n['seconds']}")
        print(f"  expert_module_rxns {o['covered'].get('expert_module_reactions')} -> {n['covered'].get('expert_module_reactions')}; "
              f"ec_topology {o['covered'].get('ec_topology',0)} -> {n['covered'].get('ec_topology',0)}")
        print(f"  layout_prov exact/adapted/topo/func: "
              f"{o['exact']}/{o['adapted']}/{o['topological']}/0 -> {n['exact']}/{n['adapted']}/{n['topological']}/{n['functional']}")
        print(f"  EC proposals kept/rejected: {n['ec_proposals_kept']}/{n['ec_proposals_rejected']}")
        print(f"  crossings {o['layout'].get('crossings')} -> {n['layout'].get('crossings')}; "
              f"label_ov {o['layout'].get('label_overlap_pairs')} -> {n['layout'].get('label_overlap_pairs')}; "
              f"node_ov {o['layout'].get('node_overlap_pairs')} -> {n['layout'].get('node_overlap_pairs')}; "
              f"p90 {o['layout'].get('length_p90_px')} -> {n['layout'].get('length_p90_px')}")
        s = row['suspicious']
        print(f"  ec_topology reactions {n['ec_topology_reactions']}; suspicious flags {s['n_suspicious']}/{s['n_suspicious']+s['n_ok']}")
        for x in s['suspicious'][:8]:
            print(f"    ! {x['reaction']} ec={x['shared_ec']} flags={x['flags']} chem={x['chemical']} src={x.get('source_reaction')}")


if __name__ == '__main__':
    main()
