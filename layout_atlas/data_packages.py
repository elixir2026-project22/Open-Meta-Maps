"""Data packages: export the curated Atlas as one archive; fetch, verify and install optional release snapshots.

The default Atlas is already in the repository: `data/open_core.sqlite` plus the in-tree extensions under
`data/extensions/` (BRENDA, Human-maps, WikiPathways geometry, Metabolism Regulation, PANTHER, iCH360). The only
optional package is `atlas-full` (a full local `atlas.sqlite` snapshot). This command downloads a package, checks its
SHA-256 against `data/packages.json` and extracts it into `data/`. Nothing is overwritten silently: an existing file is
moved to `<name>.bak` first, and symlinks are never touched.

  python3 project.py data list                      # packages, size, licence, installed or not
  python3 project.py data fetch atlas-full          # download + verify + install
  python3 project.py data fetch all                 # all packages
  python3 project.py data fetch all --from DIR      # install from a local directory (or a base URL / file:// URL) instead of the default release
  python3 project.py data export                    # whole Atlas, unique after dedup -> dist/openmetamaps-atlas-<date>.tar.gz
  python3 project.py data export --source panther --licence CC-BY-4.0 --all-windows --out atlas.tar.gz

`data export` writes modules/<source>/<module_id>.json (Escher JSON, same format as the Atlas Browser download),
manifest.json (per-module provenance, licence, canonical reaction ids and also_in for identical reaction sets in other
sources), README.md and THIRD_PARTY_NOTICES.md. Filters are repeatable. The archive is byte-identical across re-runs.
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import tempfile
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
DATA = BASE.parent / 'data'
MANIFEST = DATA / 'packages.json'
REGISTRY = DATA / 'atlas_registry.json'
PROVENANCE = ('lineage', 'quality', 'source_sha256', 'original_instance', 'extraction', 'atlas_schema_version')


def load():
    return json.loads(MANIFEST.read_text())


def installed(pkg):
    return (DATA / pkg['marker']).exists()


def cmd_list(_args):
    m = load()
    print(f"Open Meta Maps data packages ({m['version']}) - base URL {m['base_url']}\n")
    for p in m['packages']:
        print(f"  {p['name']:15s} {p['bytes'] / 1e6:7.1f} MB  {'installed' if installed(p) else 'not installed'}")
        print(f"      {p['title']}")
        print(f"      licence / terms: {p['licence']}\n")


def _open(url_or_path, base_hint):
    s = str(url_or_path)
    if '://' in s:
        return urllib.request.urlopen(s)
    return open(s, 'rb')


def _download(pkg, source, dest):
    if source is None:
        loc = load()['base_url'] + pkg['file']
    elif '://' in source:
        loc = source.rstrip('/') + '/' + pkg['file']
    else:
        loc = str(Path(source) / pkg['file'])
    print(f"  downloading {loc}")
    h = hashlib.sha256()
    with _open(loc, None) as r, open(dest, 'wb') as w:
        total = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            w.write(chunk)
            h.update(chunk)
            total += len(chunk)
    if h.hexdigest() != pkg['sha256']:
        raise SystemExit(f"checksum mismatch for {pkg['file']}: expected {pkg['sha256']}, got {h.hexdigest()}; nothing was installed")
    print(f"  checksum ok ({total / 1e6:.1f} MB)")


def _safe_members(tf):
    root = DATA.resolve()
    for m in tf.getmembers():
        target = (DATA / m.name).resolve()
        if not str(target).startswith(str(root) + os.sep) and target != root:
            raise SystemExit(f'unsafe path in package: {m.name}')
        if m.issym() or m.islnk() or m.isdev():
            raise SystemExit(f'unexpected special file in package: {m.name}')
        yield m


def _install(pkg, archive):
    with tarfile.open(archive, 'r:gz') as tf:
        members = list(_safe_members(tf))
        backed_up = 0
        for m in members:
            if not m.isfile():
                continue
            target = DATA / m.name
            if target.is_symlink():
                raise SystemExit(f'{target} is a symlink (a developer checkout that already links to this data); not touching it')
            if target.exists():
                bak = target.with_name(target.name + '.bak')
                if bak.exists():
                    raise SystemExit(f'{bak} already exists; move it away first')
                target.rename(bak)
                backed_up += 1
        tf.extractall(DATA, members=members)
    print(f"  installed {len([m for m in members if m.isfile()])} files into {DATA.relative_to(BASE.parent)}" + (f" ({backed_up} existing file(s) kept as .bak)" if backed_up else ''))
    print(f"  notice: {DATA / pkg['notice']}")


def cmd_fetch(args):
    m = load()
    names = [a for a in args if not a.startswith('--')]
    source = None
    if '--from' in args:
        i = args.index('--from')
        source = args[i + 1]
        names = [a for a in names if a != source]
    if not names:
        raise SystemExit('usage: python3 project.py data fetch <name|all> [--from DIR_OR_URL]')
    pk = {p['name']: p for p in m['packages']}
    todo = list(pk) if names == ['all'] else names
    for n in todo:
        if n not in pk:
            raise SystemExit(f"unknown package '{n}'; known: {', '.join(pk)}")
    for n in todo:
        p = pk[n]
        print(f"[{n}] {p['title']}\n  terms: {p['licence']}")
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / p['file']
            _download(p, source, archive)
            _install(p, archive)


def _atlas_rows(ab):
    """One row per licence-admitted Atlas module: manifest card, store, dedup key and provenance."""
    rows = []
    for store, m in ab.iter_modules():
        card = ab.module_card(m, store=store)
        if card is None:
            continue
        prov = {k: m[k] for k in PROVENANCE if m.get(k) is not None}
        if 'original_instance' in prov:  # repository-relative, whatever checkout built the Atlas
            prov['original_instance'] = re.sub(r'^/.*?/(?=data/)', '', prov['original_instance'])
        rows.append(dict(card, source=ab.source_key(m), key=ab.reaction_set(m), provenance=prov))
    return rows


def _select(ab, rows, sources, licences, all_windows):
    for flag, wanted, known in (('--source', sources, {r['source'] for r in rows}),
                                ('--licence', licences, {r['license'] for r in rows})):
        unknown = sorted(set(wanted) - known)
        if unknown:
            raise SystemExit(f"unknown {flag} {', '.join(unknown)}; known: {', '.join(sorted(known))}")
    picked = [r for r in rows if (not sources or r['source'] in sources) and (not licences or r['license'] in licences)]
    if not all_windows:
        best = {}
        for r in picked:
            rank = (-r['size'],) + ab._compactness(r)
            if r['key'] not in best or rank < best[r['key']][0]:
                best[r['key']] = (rank, r)
        picked = [r for _, r in best.values()]
    return sorted(picked, key=lambda r: (r['source'], r['id']))


def _record(ab, row, peers):
    rec = ab.manifest_row(row)
    rec['file'] = f"modules/{row['source']}/{ab.module_file(row['id'])}"
    rec['canonical_reaction_ids'] = sorted(row['key'])
    rec['provenance'] = row['provenance']
    rec['also_in'] = sorted((dict(module_id=p['id'], source=p['source'], map_name=p['name'])
                             for p in peers if p['source'] != row['source']),
                            key=lambda d: (d['source'], d['module_id']))
    return rec


def _archive_readme(ab, version, records, all_windows):
    by_source = Counter(r['source'] for r in records)
    by_licence = Counter(r['license'] for r in records)
    lines = [
        f'# Open Meta Maps curated-structure Atlas ({version})',
        '',
        f"{len(records):,} curated modules ({'every cut window' if all_windows else 'unique after dedup'}) "
        f"from {ab._count_label(len(by_source), 'source')}, cut from published pathway maps.",
        '',
        '- `modules/<source>/<module_id>.json`: one Escher map per module; curated geometry is uniformly scaled and translated, never rotated or mirrored.',
        '- `manifest.json`: per module: source map, licence, attribution text, reaction count, canonical reaction ids, '
        'extraction provenance and `also_in` (modules with the identical reaction set in other sources).',
        '- `THIRD_PARTY_NOTICES.md`: source credits and licences.',
        '',
        '| Source | Modules |',
        '|---|---:|',
        *[f'| {k} | {v:,} |' for k, v in sorted(by_source.items())],
        '',
        '| Licence | Modules |',
        '|---|---:|',
        *[f'| {k} | {v:,} |' for k, v in sorted(by_licence.items())],
        '',
        'CC BY modules require the `attribution_text` from the manifest when reused.',
        '',
    ]
    return '\n'.join(lines)


def export_atlas(out, sources=(), licences=(), all_windows=False):
    """Write the Atlas (or a filtered part) as a deterministic .tar.gz. Returns a summary dict."""
    import atlas_browser as ab
    registry = json.loads(REGISTRY.read_text())
    version = registry['generated']
    rows = _atlas_rows(ab)
    peers = defaultdict(list)
    for r in rows:
        peers[r['key']].append(r)
    records, files, missing = [], [], []
    for row in _select(ab, rows, sources, licences, all_windows):
        try:
            text = ab.module_json(ab.load_module(row['id'], store=row['store']))
        except ValueError:
            missing.append(row['id'])
            continue
        rec = _record(ab, row, peers[row['key']])
        records.append(rec)
        files.append((rec['file'], text.encode()))
    manifest = dict(
        format_version=1,
        generated=version,
        policy=ab.MANIFEST_POLICY,
        atlas_registry_sha256=ab.registry_fingerprint(),
        selection=dict(sources=sorted(sources), licences=sorted(licences),
                       dedup='every cut window' if all_windows else 'unique after dedup (identical canonical reaction sets)'),
        module_count=len(records),
        modules=records,
        missing=missing,
    )
    files += [
        ('README.md', _archive_readme(ab, version, records, all_windows).encode()),
        ('THIRD_PARTY_NOTICES.md', (BASE.parent / 'THIRD_PARTY_NOTICES.md').read_bytes()),
        ('manifest.json', (json.dumps(manifest, ensure_ascii=False, indent=2) + '\n').encode()),
    ]
    mtime = int(datetime.fromisoformat(version).replace(tzinfo=timezone.utc).timestamp())
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'wb') as f, gzip.GzipFile(filename='', mode='wb', fileobj=f, compresslevel=6, mtime=mtime) as gz, \
            tarfile.open(fileobj=gz, mode='w', format=tarfile.PAX_FORMAT) as tf:
        for name, data in sorted(files):
            info = tarfile.TarInfo(name)
            info.size, info.mtime, info.mode = len(data), mtime, 0o644
            tf.addfile(info, io.BytesIO(data))
    return dict(path=out, modules=len(records), sources=ab._count_label(len({r['source'] for r in records}), 'source'),
                licences=dict(Counter(r['license'] for r in records)), bytes=out.stat().st_size)


def cmd_export(args):
    p = argparse.ArgumentParser(prog='python3 project.py data export',
                                description='Export the curated Atlas as one .tar.gz (Escher JSON per module + manifest + notices).')
    p.add_argument('--out', type=Path, help='archive path (default: dist/openmetamaps-atlas-<registry date>.tar.gz)')
    p.add_argument('--source', action='append', default=[], help='only this source (repeatable), e.g. panther')
    p.add_argument('--licence', '--license', dest='licence', action='append', default=[],
                   help='only this licence (repeatable), e.g. CC0-1.0')
    p.add_argument('--all-windows', action='store_true', help='every cut window instead of unique after dedup')
    a = p.parse_args(args)
    out = a.out or BASE.parent / 'dist' / f"openmetamaps-atlas-{json.loads(REGISTRY.read_text())['generated']}.tar.gz"
    t0 = time.time()
    s = export_atlas(out, a.source, a.licence, a.all_windows)
    lic = ', '.join(f'{k} {v:,}' for k, v in sorted(s['licences'].items()))
    print(f"{s['modules']:,} modules from {s['sources']} ({lic}) -> {s['path']} "
          f"({s['bytes'] / 1e6:.1f} MB, {time.time() - t0:.0f} s)")


def main(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return
    if argv[0] == 'list':
        cmd_list(argv[1:])
    elif argv[0] == 'fetch':
        cmd_fetch(argv[1:])
    elif argv[0] == 'export':
        cmd_export(argv[1:])
    else:
        raise SystemExit(f'unknown subcommand {argv[0]!r}; use list, fetch or export')


if __name__ == '__main__':
    main(sys.argv[1:])
