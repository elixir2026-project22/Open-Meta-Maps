"""Atlas Browser: compact search index, filters, Escher export and licence manifest.

Index fields are derived generically from module payloads and registry/provenance
(no per-source hardcoding). Curated coordinates are never rotated or mirrored.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
import sqlite3
import time
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path

from core import DATA, ROOT, dump, read
from atlas_registry import ATTRIBUTION, DEFAULT_EXTENSIONS, OPEN_CORE, unpack

INDEX_PATH = DATA / 'atlas_browser_index.json.gz'
INDEX_SCHEMA = 1
# Admit CC0 / CC BY / MIT only (reject BY-SA, NC, ND, etc.).
ALLOWED_LICENSE = re.compile(
    r'(?i)^(CC0(?:-1\.0)?(?:-data)?|CC-?BY(?:-4\.0)?|MIT)$'
)
PAGE_DEFAULT = 24
SOURCE_COLORS = {
    'escher': '#2d8f6f',
    'wikipathways': '#3b77b8',
    'reactome': '#a8559a',
    'human_maps': '#c45c26',
    'brenda': '#8a6d3b',
    'metabolismregulation': '#3d7a5a',
    'panther': '#5b6b8a',
    'ich360': '#128781',
}

_index = None
_payload_cache = {}


def license_ok(lic):
    s = str(lic or '').strip()
    if not s or re.search(r'(?i)(BY-SA|BY-NC|BY-ND|NC-SA|GPL)', s):
        return False
    return bool(ALLOWED_LICENSE.match(s))


def _extension_paths():
    """Default extension module dumps, same membership as database.modules()."""
    out = []
    for name in DEFAULT_EXTENSIONS:
        path = DATA / 'extensions' / name / 'modules.json.gz'
        if path.is_file():
            out.append(path)
    return out


def iter_modules():
    """Yield (store_key, module) from open_core.sqlite + default extensions."""
    if OPEN_CORE.is_file():
        db = sqlite3.connect(OPEN_CORE)
        try:
            for (blob,) in db.execute('SELECT payload FROM modules ORDER BY id'):
                yield 'core', unpack(blob)
        finally:
            db.close()
    for path in _extension_paths():
        key = path.parent.name
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            for m in json.load(f):
                yield key, m


def _norm(s):
    return re.sub(r'\s+', ' ', str(s or '').strip().lower())


def _tokens(*parts):
    out = []
    for p in parts:
        if p is None:
            continue
        t = _norm(p)
        if t:
            out.append(t)
            out.append(re.sub(r'[^a-z0-9.:+_-]+', '', t))
    return [x for x in out if x]


def _reaction_ec_bridge():
    path = DATA / 'identifier_bridge.json'
    if not path.is_file():
        return {}
    try:
        return read(path).get('reaction_ec') or {}
    except Exception:
        return {}


def _collect_ids(module, bridge_ec):
    """Gather searchable reaction / EC / metabolite tokens from a module payload."""
    rx, ec, met = set(), set(), set()
    for r in module.get('reactions') or []:
        rid = r.get('id')
        if rid:
            rx.update(_tokens(rid))
        for ns, vals in (r.get('references') or {}).items():
            ns_l = _norm(ns)
            items = vals if isinstance(vals, (list, tuple)) else [vals]
            for v in items:
                if v is None:
                    continue
                tok = _tokens(v)
                rx.update(tok)
                rx.update(_tokens(f'{ns_l}:{v}'))
                if 'ec' in ns_l:
                    ec.update(_tokens(v))
            if rid and ns_l in ('bigg.reaction', 'bigg'):
                for code in (bridge_ec.get(str(rid)) or {}):
                    ec.update(_tokens(code))
        key = r.get('key')
        if key:
            rx.update(_tokens(str(key).split(':', 1)[-1]))
    for n in module.get('nodes') or []:
        if n.get('role') != 'metabolite':
            continue
        ent = n.get('entity')
        if ent:
            base = str(ent).split('@', 1)[0]
            met.update(_tokens(base, ent, n.get('label')))
        elif n.get('label'):
            met.update(_tokens(n.get('label')))
    return sorted(rx), sorted(ec), sorted(met)


def record_from_module(module, bridge_ec=None, attribution_fallback=None, store='core'):
    bridge_ec = bridge_ec if bridge_ec is not None else {}
    lic = module.get('license') or module.get('licence') or ''
    if not license_ok(lic):
        return None
    attr = module.get('attribution') or module.get('source_citation') or ''
    if not attr:
        attr = (attribution_fallback or {}).get(module.get('source')) or ATTRIBUTION.get(module.get('source'), '')
    rx, ec, met = _collect_ids(module, bridge_ec)
    name = module.get('name') or ''
    map_id = module.get('map_id') or ''
    lineage = module.get('lineage') or ''
    source = module.get('source') or ''
    map_tokens = sorted(set(_tokens(name, map_id, lineage, source, *name.split('.'))))
    return dict(
        id=module['id'],
        store=store,
        source=source,
        name=name,
        map_id=map_id,
        lineage=lineage,
        species=module.get('species') or '',
        license=lic,
        attribution=attr,
        size=int(module.get('size') or len(module.get('reactions') or []) or 0),
        quality=module.get('quality') or '',
        rx=rx,
        ec=ec,
        met=met,
        map=map_tokens,
    )


def build_index(path=None, modules=None):
    """Build and write the compact browser index. Returns summary dict."""
    path = Path(path or INDEX_PATH)
    bridge = _reaction_ec_bridge()
    records = []
    skipped = 0
    sources = Counter()
    licences = Counter()
    t0 = time.time()
    stream = modules if modules is not None else iter_modules()
    for item in stream:
        if isinstance(item, tuple) and len(item) == 2 and isinstance(item[1], dict):
            store, m = item
        else:
            store, m = 'core', item
        rec = record_from_module(m, bridge, store=store)
        if rec is None:
            skipped += 1
            continue
        records.append(rec)
        sources[rec['source']] += 1
        licences[rec['license']] += 1
    records.sort(key=lambda r: (-r['size'], r['source'], r['name'], r['id']))
    payload = dict(
        schema_version=INDEX_SCHEMA,
        generated=str(date.today()),
        module_count=len(records),
        skipped_license=skipped,
        sources=dict(sources),
        licences=dict(licences),
        modules=records,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(payload, ensure_ascii=False, separators=(',', ':')).encode()
    with gzip.open(path, 'wb', compresslevel=6) as f:
        f.write(raw)
    return dict(
        path=str(path.relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path),
        modules=len(records),
        skipped_license=skipped,
        bytes=path.stat().st_size,
        uncompressed=len(raw),
        seconds=round(time.time() - t0, 2),
        sources=dict(sources),
        licences=dict(licences),
    )


def load_index(path=None, force=False):
    global _index
    path = Path(path or INDEX_PATH)
    if _index is not None and not force and _index.get('_path') == str(path):
        return _index
    if not path.is_file():
        build_index(path)
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        data = json.load(f)
    data['_path'] = str(path)
    data['_by_id'] = {m['id']: m for m in data['modules']}
    _index = data
    return data


def ensure_index(path=None):
    path = Path(path or INDEX_PATH)
    if not path.is_file():
        return build_index(path)
    return dict(path=str(path), modules=None, bytes=path.stat().st_size, status='present')


def _match_list(values, query):
    q = _norm(query)
    if not q:
        return True
    for v in values:
        if q in v:
            return True
    return False


def search(filters=None, page=1, per_page=PAGE_DEFAULT, index=None):
    """Filter modules. All filter values are optional substring / equality matches."""
    filters = filters or {}
    idx = index or load_index()
    rows = idx['modules']
    src = _norm(filters.get('source') or filters.get('database') or '')
    lic = _norm(filters.get('license') or filters.get('licence') or '')
    ec = filters.get('ec') or filters.get('ec_number') or ''
    rx = filters.get('reaction') or filters.get('reaction_id') or ''
    met = filters.get('metabolite') or filters.get('metabolite_id') or ''
    mmap = filters.get('map') or filters.get('map_id') or filters.get('model') or ''
    q = filters.get('q') or filters.get('query') or ''

    hits = []
    for m in rows:
        if src and src != _norm(m['source']):
            continue
        if lic and lic != _norm(m['license']):
            continue
        if ec and not _match_list(m['ec'], ec):
            continue
        if rx and not _match_list(m['rx'], rx):
            continue
        if met and not _match_list(m['met'], met):
            continue
        if mmap and not _match_list(m['map'], mmap):
            continue
        if q:
            blob = ' '.join([m['name'], m['map_id'], m['lineage'], m['source'], m['species'],
                             ' '.join(m['rx'][:40]), ' '.join(m['ec'][:20]), ' '.join(m['met'][:40])])
            if _norm(q) not in _norm(blob):
                continue
        hits.append(m)

    per_page = max(1, min(int(per_page or PAGE_DEFAULT), 100))
    page = max(1, int(page or 1))
    total = len(hits)
    pages = max(1, (total + per_page - 1) // per_page)
    if page > pages:
        page = pages
    start = (page - 1) * per_page
    slice_ = hits[start:start + per_page]
    return dict(
        total=total,
        page=page,
        per_page=per_page,
        pages=pages,
        results=[public_card(m) for m in slice_],
        filters={k: v for k, v in filters.items() if v},
    )


def public_card(m):
    return dict(
        id=m['id'],
        source=m['source'],
        name=m['name'],
        map_id=m['map_id'],
        species=m.get('species') or '',
        license=m['license'],
        attribution=m.get('attribution') or '',
        size=m['size'],
        quality=m.get('quality') or '',
        thumbnail=f"/api/atlas/modules/{m['id']}/thumbnail.svg",
    )


def meta(index=None):
    idx = index or load_index()
    return dict(
        schema_version=idx.get('schema_version'),
        generated=idx.get('generated'),
        module_count=idx.get('module_count'),
        sources=sorted((idx.get('sources') or {}).items()),
        licences=sorted((idx.get('licences') or {}).items()),
        filters=['q', 'ec', 'reaction', 'metabolite', 'map', 'source', 'license'],
    )


def load_module(module_id, index=None):
    idx = index or load_index()
    card = (idx.get('_by_id') or {}).get(module_id)
    store = (card or {}).get('store') or 'core'
    if store == 'core':
        if not OPEN_CORE.is_file():
            return None
        db = sqlite3.connect(OPEN_CORE)
        try:
            row = db.execute('SELECT payload FROM modules WHERE id=?', (module_id,)).fetchone()
        finally:
            db.close()
        return unpack(row[0]) if row else None
    path = DATA / 'extensions' / store / 'modules.json.gz'
    if not path.is_file():
        return None
    key = str(path)
    if key not in _payload_cache:
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            _payload_cache[key] = {m['id']: m for m in json.load(f)}
    return _payload_cache[key].get(module_id)


def module_detail(module_id, index=None):
    idx = index or load_index()
    card = (idx.get('_by_id') or {}).get(module_id)
    mod = load_module(module_id)
    if mod is None:
        return None
    if not license_ok(mod.get('license') or mod.get('licence') or (card or {}).get('license')):
        return None
    reactions = []
    for r in mod.get('reactions') or []:
        reactions.append(dict(
            id=r.get('id'),
            key=r.get('key'),
            references=r.get('references') or {},
            reversible=bool(r.get('reversible')),
            coefficients=r.get('coefficients') or {},
        ))
    bridge = _reaction_ec_bridge()
    rec = card or record_from_module(mod, bridge) or {}
    return dict(
        id=mod['id'],
        source=mod.get('source'),
        name=mod.get('name'),
        map_id=mod.get('map_id'),
        lineage=mod.get('lineage'),
        species=mod.get('species'),
        license=mod.get('license') or mod.get('licence'),
        attribution=rec.get('attribution') or mod.get('attribution') or ATTRIBUTION.get(mod.get('source'), ''),
        size=rec.get('size') or mod.get('size'),
        quality=mod.get('quality'),
        source_sha256=mod.get('source_sha256'),
        original_instance=mod.get('original_instance'),
        layout_origin=mod.get('layout_origin'),
        extraction=mod.get('extraction'),
        core_family=mod.get('core_family'),
        variant=mod.get('variant'),
        atlas_schema_version=mod.get('atlas_schema_version'),
        reactions=reactions,
        ports=mod.get('ports') or [],
        node_count=len(mod.get('nodes') or []),
        edge_count=len(mod.get('edges') or []),
        thumbnail=f"/api/atlas/modules/{mod['id']}/thumbnail.svg",
    )


def thumbnail_svg(module, width=220, height=140):
    """Reuse the Map Generator thumbnail style; uniform scale, no rotate/mirror."""
    from tool_pipeline import thumbnail_svg as _thumb
    color = SOURCE_COLORS.get(module.get('source'), '#128781')
    return _thumb(module, color, width=width, height=height) or (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">'
        f'<text x="12" y="{height // 2}" fill="#5b7480" font-size="12">No geometry</text></svg>'
    )


def module_to_escher(module):
    """Export one curated module as Escher JSON. Coordinates are translated only (no flip/rotate)."""
    nodes_in = module.get('nodes') or []
    edges = module.get('edges') or []
    reactions_in = {str(r.get('key') or ('r:' + str(r.get('id')))): r for r in (module.get('reactions') or [])}
    # Prefer stored xy; drop nodes without coordinates.
    usable = [n for n in nodes_in if n.get('xy') is not None]
    if not usable:
        raise ValueError('module has no drawable coordinates')
    xs = [n['xy'][0] for n in usable]
    ys = [n['xy'][1] for n in usable]
    # Uniform scale to a readable canvas, preserving aspect (no mirror).
    min_x, min_y = min(xs), min(ys)
    span = max(max(xs) - min_x, max(ys) - min_y, 1e-9)
    scale = 400.0 / span
    pad = 40.0

    def xy(n):
        return (pad + (n['xy'][0] - min_x) * scale, pad + (n['xy'][1] - min_y) * scale)

    node_ids = {}
    nodes = {}
    for i, n in enumerate(usable):
        nid = str(i)
        node_ids[n['key']] = nid
        x, y = xy(n)
        if n.get('role') == 'reaction':
            nodes[nid] = dict(node_type='midmarker', x=x, y=y)
        else:
            ent = str(n.get('entity') or n.get('label') or n['key']).split('@', 1)[0]
            nodes[nid] = dict(
                node_type='metabolite', x=x, y=y, bigg_id=ent,
                name=n.get('label') or ent, node_is_primary=True,
                label_x=x + 8, label_y=y - 8,
            )

    reactions = {}
    seg_i = 0
    for ri, n in enumerate(usable):
        if n.get('role') != 'reaction':
            continue
        rk = n['key']
        rx = reactions_in.get(rk) or reactions_in.get(rk.split(':', 1)[-1]) or {}
        rid = str(rx.get('id') or rk.split(':', 1)[-1])
        x, y = xy(n)
        segs = {}
        mets = []
        coeffs = rx.get('coefficients') or {}
        # Fall back to edge flow when coefficients are absent.
        seen_met = set()
        for e in edges:
            a, b = e.get('a'), e.get('b')
            if a != rk and b != rk:
                continue
            mk = b if a == rk else a
            if mk not in node_ids or mk in seen_met:
                continue
            seen_met.add(mk)
            met_node = next((x for x in usable if x['key'] == mk), None)
            ent = str((met_node or {}).get('entity') or mk).split('@', 1)[0]
            flow = e.get('flow')
            if flow == 'substrate':
                coef = -abs(float(e.get('coefficient') or coeffs.get(ent) or coeffs.get((met_node or {}).get('entity')) or 1))
            elif flow == 'product':
                coef = abs(float(e.get('coefficient') or coeffs.get(ent) or coeffs.get((met_node or {}).get('entity')) or 1))
            else:
                coef = float(e.get('coefficient') or coeffs.get(ent) or 1)
                if coef == 0:
                    coef = 1.0
            mets.append(dict(bigg_id=ent, coefficient=coef))
            segs[str(seg_i)] = dict(from_node_id=node_ids[rk], to_node_id=node_ids[mk], b1=None, b2=None)
            seg_i += 1
        if not mets and coeffs:
            for ent, coef in coeffs.items():
                mets.append(dict(bigg_id=str(ent).split('@', 1)[0], coefficient=float(coef)))
        reactions[str(ri)] = dict(
            bigg_id=rid, name=rid, reversibility=bool(rx.get('reversible')),
            label_x=x + 6, label_y=y - 6, gene_reaction_rule='', genes=[],
            metabolites=mets, segments=segs,
        )

    width = pad * 2 + (max(xs) - min_x) * scale
    height = pad * 2 + (max(ys) - min_y) * scale
    header = dict(
        schema='https://escher.github.io/escher/jsonschema/1-0-0#',
        homepage='https://escher.github.io',
        map_id=module.get('id'),
        map_name=module.get('name') or module.get('id'),
        map_description='Open Meta Maps Atlas module export; curated geometry (no rotate/mirror)',
        metadata=dict(
            atlas_module_id=module.get('id'),
            source=module.get('source'),
            license=module.get('license') or module.get('licence'),
            attribution=module.get('attribution') or ATTRIBUTION.get(module.get('source'), ''),
            map_id=module.get('map_id'),
            source_sha256=module.get('source_sha256'),
        ),
    )
    body = dict(nodes=nodes, reactions=reactions, text_labels={}, canvas=dict(x=0, y=0, width=width, height=height))
    return [header, body]


def attribution_text(module_or_card):
    """CC BY (and other) attribution line for manifests."""
    lic = module_or_card.get('license') or ''
    attr = module_or_card.get('attribution') or ATTRIBUTION.get(module_or_card.get('source'), '')
    name = module_or_card.get('name') or module_or_card.get('id')
    src = module_or_card.get('source') or ''
    base = attr or f'{src} module {name}'
    if re.search(r'(?i)CC-?BY', lic):
        return f'{base}. Licensed under {lic}. Attribution required under Creative Commons Attribution.'
    if re.search(r'(?i)^MIT\b', lic):
        return f'{base}. Licensed under {lic}.'
    if re.search(r'(?i)^CC0', lic):
        return f'{base}. Dedicated to the public domain under {lic}.'
    return f'{base}. Licence: {lic}.'


def download_manifest(module_ids, index=None, fmt='json'):
    """Build provenance/licence manifest rows for selected modules."""
    idx = index or load_index()
    by = idx.get('_by_id') or {m['id']: m for m in idx['modules']}
    rows = []
    missing = []
    for mid in module_ids:
        card = by.get(mid)
        if not card:
            missing.append(mid)
            continue
        rows.append(dict(
            module_id=card['id'],
            source=card['source'],
            map_name=card['name'],
            map_id=card['map_id'],
            species=card.get('species') or '',
            reaction_count=card['size'],
            license=card['license'],
            attribution=card.get('attribution') or '',
            attribution_text=attribution_text(card),
        ))
    if fmt == 'csv':
        buf = io.StringIO()
        fields = ['module_id', 'source', 'map_name', 'map_id', 'species', 'reaction_count',
                  'license', 'attribution', 'attribution_text']
        w = csv.DictWriter(buf, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
        return buf.getvalue(), rows, missing
    payload = dict(
        format_version=1,
        generated=str(date.today()),
        policy='Per-module source licences preserved; CC BY requires attribution_text',
        modules=rows,
        missing=missing,
    )
    return json.dumps(payload, ensure_ascii=False, indent=2) + '\n', rows, missing


def build_download_zip(module_ids, index=None):
    """Zip of one Escher JSON per module plus JSON and CSV manifests."""
    idx = index or load_index()
    manifest_json, rows, missing = download_manifest(module_ids, index=idx, fmt='json')
    manifest_csv, _, _ = download_manifest(module_ids, index=idx, fmt='csv')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr('manifest.json', manifest_json)
        zf.writestr('manifest.csv', manifest_csv)
        for row in rows:
            mod = load_module(row['module_id'])
            if mod is None:
                missing.append(row['module_id'])
                continue
            try:
                data = module_to_escher(mod)
            except Exception:
                missing.append(row['module_id'])
                continue
            safe = re.sub(r'[^A-Za-z0-9._-]+', '_', row['module_id'])[:80]
            zf.writestr(f'modules/{safe}.json', json.dumps(data, ensure_ascii=False, indent=2) + '\n')
        if missing:
            zf.writestr('missing.txt', '\n'.join(missing) + '\n')
    return buf.getvalue(), rows, missing


def page_html():
    """Static Atlas Browser page (English), same visual language as Map Generator."""
    import ui
    css_extra = '''<style>
.ab-toolbar{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:8px 0 12px}
.ab-toolbar .sel{font-size:13px;color:var(--mute)}
.ab-card{position:relative;cursor:pointer}
.ab-card.sel{outline:2px solid var(--accent)}
.ab-card .chk{position:absolute;top:10px;right:10px;z-index:2}
.ab-card img.th,.ab-card object.th{display:block;width:100%;height:140px;object-fit:contain;background:var(--bg);border-radius:8px;margin-bottom:10px;border:1px solid var(--line)}
.ab-empty{color:var(--mute);padding:24px 0}
.ab-detail{display:none}.ab-detail.on{display:block}
.ab-list.on{display:none}
.ab-rx{max-height:280px;overflow:auto;font:12px/1.45 ui-monospace,Menlo,monospace}
button.primary{background:var(--accent);color:var(--accent-ink);border-color:var(--accent);font-weight:600}
.filters .grow{flex:1;min-width:min(280px,100%)}
</style>'''
    body = css_extra + '''
<div id="listView" class="ab-list">
<form class="filters" id="filters" autocomplete="off">
<label class="grow">Search<input type="search" name="q" id="q" placeholder="name, id, free text"></label>
<label>EC number<input name="ec" id="ec" placeholder="1.1.1.1"></label>
<label>Reaction ID<input name="reaction" id="reaction" placeholder="BiGG / Rhea / source"></label>
<label>Metabolite ID<input name="metabolite" id="metabolite" placeholder="g6p / chebi:…"></label>
<label>Map / model<input name="map" id="map" placeholder="map id or name"></label>
<label>Source database<select name="source" id="source"><option value="">All</option></select></label>
<label>Licence<select name="license" id="license"><option value="">All</option></select></label>
<button type="submit" class="primary">Search</button>
<button type="button" id="reset">Reset</button>
</form>
<div class="ab-toolbar">
<label class="sel"><input type="checkbox" id="selPage"> Select all on page</label>
<span class="count" id="count">Loading index…</span>
<button type="button" id="dlSel" disabled>Download selected</button>
</div>
<div class="grid" id="grid"></div>
<nav class="pager" id="pager" aria-label="Pages"></nav>
</div>
<div id="detailView" class="ab-detail">
<p class="toolbar"><button type="button" class="lnk" id="back">← Back to results</button>
 · <button type="button" id="dlOne">Download this module</button></p>
<div class="panel" id="detail"></div>
</div>
<script>
(function(){
var $=function(id){return document.getElementById(id)};
var state={page:1,per:24,selected:{},lastIds:[],detail:null};
function esc(s){return String(s==null?'':s).replace(/[&<>"']/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]})}
function params(){var f=$('filters'),o={};['q','ec','reaction','metabolite','map','source','license'].forEach(function(k){var v=f[k].value.trim();if(v)o[k]=v});o.page=state.page;o.per_page=state.per;return o}
function qs(o){return Object.keys(o).map(function(k){return encodeURIComponent(k)+'='+encodeURIComponent(o[k])}).join('&')}
function selCount(){return Object.keys(state.selected).length}
function syncDl(){$('dlSel').disabled=!selCount();$('dlSel').textContent=selCount()?'Download selected ('+selCount()+')':'Download selected'}
function loadMeta(){
 fetch('/api/atlas/meta').then(function(r){return r.json()}).then(function(m){
  var s=$('source'),l=$('license');
  (m.sources||[]).forEach(function(x){var o=document.createElement('option');o.value=x[0];o.textContent=x[0]+' ('+x[1]+')';s.appendChild(o)});
  (m.licences||[]).forEach(function(x){var o=document.createElement('option');o.value=x[0];o.textContent=x[0]+' ('+x[1]+')';l.appendChild(o)});
  search();
 }).catch(function(){$('count').textContent='Could not load Atlas index. Is the tool server running?'})}
function card(m){
 var checked=state.selected[m.id]?' checked':'';
 return '<article class="card ab-card'+(state.selected[m.id]?' sel':'')+'" data-id="'+esc(m.id)+'">'
  +'<label class="chk" onclick="event.stopPropagation()"><input type="checkbox" data-id="'+esc(m.id)+'"'+checked+'></label>'
  +'<img class="th" loading="lazy" alt="" src="'+esc(m.thumbnail)+'">'
  +'<h2 class="dn">'+esc(m.name||m.id)+'</h2>'
  +'<p class="meta">'+esc(m.source)+' · '+esc(m.name||m.map_id)+'</p>'
  +'<div class="tags"><span class="tag src-'+esc(m.source)+'">'+esc(m.source)+'</span>'
  +'<span class="tag">'+esc(m.size)+' reactions</span>'
  +'<span class="tag lic">'+esc(m.license)+'</span></div></article>'}
function render(data){
 state.lastIds=(data.results||[]).map(function(m){return m.id});
 $('grid').innerHTML=(data.results||[]).map(card).join('')||'<p class="ab-empty">No modules match these filters.</p>';
 $('count').textContent=data.total.toLocaleString()+' modules · page '+data.page+' / '+data.pages;
 var p=$('pager'),html='',i,start=Math.max(1,data.page-3),end=Math.min(data.pages,data.page+3);
 if(data.page>1)html+='<a href="#" data-p="'+(data.page-1)+'">‹</a>';
 for(i=start;i<=end;i++)html+=i===data.page?'<span>'+i+'</span>':'<a href="#" data-p="'+i+'">'+i+'</a>';
 if(data.page<data.pages)html+='<a href="#" data-p="'+(data.page+1)+'">›</a>';
 p.innerHTML=html;
 $('selPage').checked=state.lastIds.length>0&&state.lastIds.every(function(id){return state.selected[id]});
 syncDl();
}
function search(){
 $('count').textContent='Searching…';
 fetch('/api/atlas/search?'+qs(params())).then(function(r){return r.json()}).then(render)
  .catch(function(){$('count').textContent='Search failed.'})}
$('filters').onsubmit=function(e){e.preventDefault();state.page=1;search()};
$('reset').onclick=function(){['q','ec','reaction','metabolite','map'].forEach(function(k){$(k).value=''});$('source').value='';$('license').value='';state.page=1;search()};
$('pager').onclick=function(e){var a=e.target.closest('[data-p]');if(!a)return;e.preventDefault();state.page=+a.getAttribute('data-p');search()};
function toggleSel(cb){if(!cb)return;if(cb.checked)state.selected[cb.getAttribute('data-id')]=1;else delete state.selected[cb.getAttribute('data-id')];syncDl();var art=cb.closest('.ab-card');if(art)art.classList.toggle('sel',cb.checked)}
$('grid').onclick=function(e){
 var cb=e.target.closest('input[type=checkbox][data-id]');
 if(cb){toggleSel(cb);return}
 var card=e.target.closest('.ab-card');if(card)openDetail(card.getAttribute('data-id'))};
$('grid').onchange=function(e){var cb=e.target.closest('input[type=checkbox][data-id]');if(cb)toggleSel(cb)};
$('selPage').onchange=function(){
 state.lastIds.forEach(function(id){if($('selPage').checked)state.selected[id]=1;else delete state.selected[id]});
 Array.prototype.forEach.call(document.querySelectorAll('.ab-card'),function(art){
  var id=art.getAttribute('data-id');var on=!!$('selPage').checked;art.classList.toggle('sel',on);
  var cb=art.querySelector('input[type=checkbox]');if(cb)cb.checked=on});
 syncDl()};
function openDetail(id){
 fetch('/api/atlas/modules/'+encodeURIComponent(id)).then(function(r){return r.json()}).then(function(d){
  if(d.error){alert(d.error.message_en||'Not found');return}
  state.detail=d;
  $('listView').classList.add('on');$('detailView').classList.add('on');
  var rx=(d.reactions||[]).map(function(r){
   var refs=Object.keys(r.references||{}).map(function(ns){return ns+':'+(r.references[ns]||[]).join(',')}).join('; ');
   return '<tr><td>'+esc(r.id||'')+'</td><td>'+esc(refs)+'</td></tr>'}).join('');
  $('detail').innerHTML='<div class="detail-grid"><div><img class="th" style="height:240px;width:100%" src="'+esc(d.thumbnail)+'" alt="">'
   +'</div><div><dl class="kv">'
   +'<dt>Source</dt><dd>'+esc(d.source)+'</dd>'
   +'<dt>Map</dt><dd>'+esc(d.name)+' <span class="muted">'+esc(d.map_id)+'</span></dd>'
   +'<dt>Species</dt><dd>'+esc(d.species||'')+'</dd>'
   +'<dt>Licence</dt><dd>'+esc(d.license)+'</dd>'
   +'<dt>Attribution</dt><dd>'+esc(d.attribution||'')+'</dd>'
   +'<dt>Reactions</dt><dd>'+esc(d.size)+'</dd>'
   +'<dt>Quality</dt><dd>'+esc(d.quality||'')+'</dd>'
   +'<dt>SHA-256</dt><dd class="muted">'+esc(d.source_sha256||'')+'</dd>'
   +'</dl></div></div>'
   +'<h2>Member reactions</h2><div class="tablewrap"><table><thead><tr><th>ID</th><th>References</th></tr></thead><tbody class="ab-rx">'+rx+'</tbody></table></div>'
   +(d.extraction?'<details class="det"><summary>Extraction provenance</summary><pre>'+esc(JSON.stringify(d.extraction,null,2))+'</pre></details>':'');
 }).catch(function(){alert('Could not load module')})}
$('back').onclick=function(){$('listView').classList.remove('on');$('detailView').classList.remove('on');state.detail=null};
function download(ids){
 fetch('/api/atlas/download',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({ids:ids})})
  .then(function(r){if(!r.ok)return r.json().then(function(j){throw j});return r.blob()})
  .then(function(b){var a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='atlas_modules.zip';a.click();URL.revokeObjectURL(a.href)})
  .catch(function(e){alert((e&&e.error&&e.error.message_en)||'Download failed')})}
$('dlSel').onclick=function(){download(Object.keys(state.selected))};
$('dlOne').onclick=function(){if(state.detail)download([state.detail.id])};
loadMeta();
})();
</script>'''
    dest = ROOT / 'layout_atlas' / 'results' / 'atlas_browser'
    dest.mkdir(parents=True, exist_ok=True)
    lead = ui.bi('Search curated Atlas modules by EC, reaction, metabolite, source map or licence, then download Escher JSON with a provenance manifest.')
    ui.write(dest / 'index.html', 'Atlas Browser', body, section='atlas', lead=lead)
    return dest / 'index.html'
