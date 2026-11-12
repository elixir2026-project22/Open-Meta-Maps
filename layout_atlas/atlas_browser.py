"""Atlas Browser: compact search index, filters, Escher export and licence manifest.

Index fields are derived generically from module payloads and registry/provenance
(no per-source hardcoding). Curated coordinates are never rotated or mirrored.
The search index is built on first use under data/cache/ and rebuilt when
data/atlas_registry.json changes.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import sqlite3
import time
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape

from core import DATA, ROOT, dump, read, sha
from atlas_registry import ATTRIBUTION, DEFAULT_EXTENSIONS, OPEN_CORE, REGISTRY_PATH, unpack

CACHE_DIR = DATA / 'cache' / 'atlas_browser'
INDEX_PATH = CACHE_DIR / 'index.json.gz'
META_PATH = CACHE_DIR / 'meta.json'
INDEX_SCHEMA = 2
# Admit CC0 / CC BY / MIT only (reject BY-SA, NC, ND, etc.).
ALLOWED_LICENSE = re.compile(
    r'(?i)^(CC0(?:-1\.0)?(?:-data)?|CC-?BY(?:-4\.0)?|MIT)$'
)
PAGE_DEFAULT = 24
CANON_NS = ('bigg.reaction', 'metanetx.reaction', 'rhea', 'ec-code')
QUALITY_INFO = {
    'semantic_transfer': (
        'Curated reaction layout',
        'Geometry from a curated source map with reaction identities attached.',
    ),
    'geometry_only': (
        'Geometry only',
        'Drawn layout without reaction identifiers (shape cut from the map).',
    ),
    'annotated_topology': (
        'Annotated topology',
        'Reaction topology with annotations; the layout may be automatic.',
    ),
}
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
SOURCE_LABELS = {
    'escher': 'Escher',
    'wikipathways': 'WikiPathways',
    'reactome': 'Reactome',
    'human_maps': 'Human maps',
    'brenda': 'BRENDA',
    'metabolismregulation': 'Metabolism Regulation',
    'panther': 'PANTHER',
    'ich360': 'iCH360',
}
Q_KIND_LABELS = {
    'ec': 'Matched as EC number',
    'reaction': 'Matched as reaction ID',
    'metabolite': 'Matched as metabolite ID',
    'map': 'Matched as map / model',
    'text': 'Matched as free text',
}

_index = None
_payload_cache = {}


def license_ok(lic):
    s = str(lic or '').strip()
    if not s or re.search(r'(?i)(BY-SA|BY-NC|BY-ND|NC-SA|GPL)', s):
        return False
    return bool(ALLOWED_LICENSE.match(s))


def attribution_string(value, source=''):
    """Normalise attribution to plain text (some sources store a dict)."""
    if isinstance(value, dict):
        parts = [value.get('title') or value.get('attribution') or value.get('text') or '']
        if value.get('paper'):
            parts.append(str(value['paper']))
        if value.get('doi'):
            parts.append('doi:' + str(value['doi']))
        if value.get('url'):
            parts.append(str(value['url']))
        text = '; '.join(p for p in parts if p)
        return text or ATTRIBUTION.get(source, '')
    text = str(value or '').strip()
    return text or ATTRIBUTION.get(source, '')


def quality_label(code):
    label, tip = QUALITY_INFO.get(code or '', (None, None))
    if label:
        return dict(code=code or '', label=label, tip=tip)
    if not code:
        return dict(code='', label='Unknown', tip='No quality tag on this module.')
    return dict(code=code, label=str(code).replace('_', ' ').capitalize(), tip='Stored quality tag from the Atlas payload.')


def _extension_paths():
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


def _vals(refs, ns):
    vals = (refs or {}).get(ns) or []
    if isinstance(vals, (str, int, float)):
        vals = [vals]
    return [str(v) for v in vals if v is not None and str(v).strip()]


def canonical_reaction_ids(module):
    """Stable chemical identity tokens for one module (generic over namespaces)."""
    out = []
    for r in module.get('reactions') or []:
        refs = r.get('references') or {}
        picked = None
        for ns in CANON_NS:
            vals = _vals(refs, ns)
            if vals:
                picked = f'{ns}:{vals[0]}'
                break
        if picked is None:
            for ns in sorted(refs):
                if ns in CANON_NS:
                    continue
                vals = _vals(refs, ns)
                if vals:
                    picked = f'{ns}:{vals[0]}'
                    break
        if picked is None and r.get('id'):
            picked = f"id:{r['id']}"
        if picked:
            out.append(picked)
    return tuple(sorted(set(out)))


def reaction_set_signature(module):
    ids = canonical_reaction_ids(module)
    if not ids:
        # Geometry-only: do not collapse distinct drawings together.
        return ''
    return hashlib.sha1('\n'.join(ids).encode()).hexdigest()[:16]


def _collect_ids(module, bridge_ec):
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
    source = module.get('source') or ''
    attr = attribution_string(
        module.get('attribution') or module.get('source_citation') or '',
        source,
    )
    if not attr and attribution_fallback:
        attr = attribution_string(attribution_fallback.get(source), source)
    rx, ec, met = _collect_ids(module, bridge_ec)
    name = module.get('name') or ''
    map_id = module.get('map_id') or ''
    lineage = module.get('lineage') or ''
    map_tokens = sorted(set(_tokens(name, map_id, lineage, source, *name.split('.'))))
    rx_canon = list(canonical_reaction_ids(module))
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
        size=int(module.get('size') or len(module.get('reactions') or []) or len(rx_canon) or 0),
        quality=module.get('quality') or '',
        rx_sig=reaction_set_signature(module),
        rx_canon=rx_canon,
        rx=rx,
        ec=ec,
        met=met,
        map=map_tokens,
    )


def registry_fingerprint():
    path = Path(REGISTRY_PATH) if REGISTRY_PATH else DATA / 'atlas_registry.json'
    if not path.is_file():
        return ''
    return sha(path.read_bytes())


def build_index(path=None, modules=None):
    """Build and write the compact browser index under data/cache/. Returns summary dict."""
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
    # Larger / more complete modules first.
    records.sort(key=lambda r: (-r['size'], -len(r.get('rx_canon') or []), r['source'], r['name'], r['id']))
    unique_rx = len({r.get('rx_sig') for r in records if r.get('rx_sig')})
    # Stage B README unique-after-dedup (per-source nest dedupe); see atlas_overlap.json.
    overlap_unique = None
    overlap_path = ROOT / 'docs' / 'images' / 'atlas' / 'atlas_overlap.json'
    if overlap_path.is_file():
        try:
            overlap_unique = json.loads(overlap_path.read_text(encoding='utf-8')).get('unique_reaction_sets_global')
        except Exception:
            overlap_unique = None
    reg_sha = registry_fingerprint()
    payload = dict(
        schema_version=INDEX_SCHEMA,
        generated=str(date.today()),
        registry_sha256=reg_sha,
        module_count=len(records),
        unique_reaction_sets=unique_rx,
        readme_cut_modules=len(records) + skipped,
        readme_unique_reaction_sets=overlap_unique,
        skipped_license=skipped,
        sources=dict(sources),
        licences=dict(licences),
        modules=records,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(payload, ensure_ascii=False, separators=(',', ':')).encode()
    with gzip.open(path, 'wb', compresslevel=6) as f:
        f.write(raw)
    dump(META_PATH if path == INDEX_PATH else path.with_name('meta.json'), dict(
        schema_version=INDEX_SCHEMA,
        generated=payload['generated'],
        registry_sha256=reg_sha,
        modules=len(records),
        unique_reaction_sets=unique_rx,
        readme_cut_modules=payload['readme_cut_modules'],
        readme_unique_reaction_sets=overlap_unique,
        bytes=path.stat().st_size,
    ))
    global _index
    _index = None
    return dict(
        path=str(path.relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path),
        modules=len(records),
        unique_reaction_sets=unique_rx,
        readme_cut_modules=payload['readme_cut_modules'],
        readme_unique_reaction_sets=overlap_unique,
        skipped_license=skipped,
        bytes=path.stat().st_size,
        uncompressed=len(raw),
        seconds=round(time.time() - t0, 2),
        registry_sha256=reg_sha,
        sources=dict(sources),
        licences=dict(licences),
    )


def _index_is_fresh(path=None):
    path = Path(path or INDEX_PATH)
    if not path.is_file():
        return False
    meta_path = META_PATH if path == INDEX_PATH else path.with_name('meta.json')
    try:
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            head = json.load(f)
        if int(head.get('schema_version') or 0) != INDEX_SCHEMA:
            return False
        want = registry_fingerprint()
        got = head.get('registry_sha256') or ''
        if want and got and want != got:
            return False
        if meta_path.is_file():
            meta = read(meta_path)
            if want and meta.get('registry_sha256') and meta['registry_sha256'] != want:
                return False
        return True
    except Exception:
        return False


def load_index(path=None, force=False):
    """Load index from disk. `force` only bypasses the in-memory cache (does not rebuild)."""
    global _index
    path = Path(path or INDEX_PATH)
    if _index is not None and not force and _index.get('_path') == str(path):
        return _index
    if path == INDEX_PATH:
        ensure_index(path)
    elif not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        data = json.load(f)
    data['_path'] = str(path)
    data['_by_id'] = {m['id']: m for m in data['modules']}
    _index = data
    return data


def ensure_index(path=None):
    path = Path(path or INDEX_PATH)
    if _index_is_fresh(path):
        return dict(path=str(path), modules=None, bytes=path.stat().st_size, status='present')
    return build_index(path)


def _match_list(values, query):
    q = _norm(query)
    if not q:
        return True
    for v in values:
        if q in v:
            return True
    return False


def _truthy(v):
    if isinstance(v, bool):
        return v
    return str(v or '').strip().lower() in ('1', 'true', 'yes', 'on', 'all')


# EC: 1.1.1.1 / 1.1.1.- / 1.1.1.100
_EC_RE = re.compile(r'^\d{1,2}\.\d{1,2}\.\d{1,2}(?:\.(?:\d{1,3}|-))?$')
# BiGG-like metabolite: g6p_c, glc__D_e
_MET_RE = re.compile(r'(?i)^(?:chebi[:_]\d+|[a-z][a-z0-9]*(?:__[a-z0-9]+)?_[a-z]\d?)$')
# Model / map tokens (BiGG GEM ids, e_coli_core, RECON1, Human-GEM, dotted map names).
_MODEL_RE = re.compile(
    r'(?i)^(?:e_coli_core|recon\d*|human[_-]?gem|yeast[_-]?gem|'
    r'i[a-z]{1,3}\d{3,4}|[a-z][a-z0-9]*(?:_[a-z0-9]+)+)$'
)
_MAP_DOT_RE = re.compile(r'(?i)^[A-Za-z][A-Za-z0-9_-]*\.[A-Za-z0-9_.-]+$')
# Namespaced reaction ids (Rhea / Reactome / MetaNetX) and bare BiGG-like tokens.
_RX_NS_RE = re.compile(r'(?i)^(rhea|reactome|metanetx\.reaction):\S+$')
_RX_RE = re.compile(r'(?i)^[A-Za-z][A-Za-z0-9_.-]{1,40}$')


def classify_query(q):
    """Detect whether free-text looks like an EC, reaction, metabolite, or model/map id."""
    s = str(q or '').strip()
    if not s:
        return 'text'
    if _EC_RE.match(s):
        return 'ec'
    if s.lower().startswith('chebi:') or _MET_RE.match(s):
        return 'metabolite'
    if _MODEL_RE.match(s) or _MAP_DOT_RE.match(s):
        return 'map'
    if _RX_NS_RE.match(s):
        return 'reaction'
    # Prefer reaction for ALLCAPS / id-like tokens; avoid swallowing plain English.
    if _RX_RE.match(s) and (' ' not in s) and (s.isupper() or re.search(r'\d', s)):
        return 'reaction'
    return 'text'


_UNKNOWN_SPECIES = re.compile(
    r'(?i)^(unspecified|unknown|n/?a|none|pathway-level\b.*|.*organism assignment not.*|'
    r'.*multi-organism.*)$'
)


def species_display(value):
    """Return a user-facing species string, or None when unknown/placeholder."""
    sp = str(value or '').strip()
    if not sp or _UNKNOWN_SPECIES.match(sp):
        return None
    return sp


def _quality_rank(code):
    order = {'semantic_transfer': 0, 'annotated_topology': 1, 'geometry_only': 2}
    return order.get(code or '', 9)


def _group_key(m):
    sig = m.get('rx_sig') or ''
    return sig if sig else f"id:{m['id']}"


def source_label(code):
    return SOURCE_LABELS.get(code, code or '')


def _map_key(m):
    return (m.get('source') or '', m.get('name') or m.get('map_id') or m.get('id') or '')


def _query_matches(m, ec='', rx='', met='', mmap='', q=''):
    """Return (matched_id_list, hit_count) for ranking and card subtitles."""
    matched = []
    seen = set()

    def add(token):
        t = str(token or '').strip()
        if not t:
            return
        key = _norm(t)
        if key in seen:
            return
        # Prefer the bare id over a namespaced duplicate (e.g. HEX1 over bigg.reaction:hex1).
        bare = key.split(':', 1)[-1]
        if ':' in key and bare in seen:
            return
        if ':' not in key:
            # Replace any prior namespaced form of the same bare id.
            for i, prev in enumerate(list(matched)):
                if _norm(prev).endswith(':' + key) or _norm(prev) == key:
                    matched[i] = t
                    seen.discard(_norm(prev))
                    seen.add(key)
                    return
        seen.add(key)
        matched.append(t)

    def collect(values, query):
        qn = _norm(query)
        if not qn:
            return
        # Prefer the user's typed token for the subtitle (preserve casing).
        hit = False
        for v in values or []:
            if qn in _norm(v) or _norm(v).endswith(':' + qn):
                hit = True
                break
        if hit:
            add(query)
            return
        for v in values or []:
            if qn in _norm(v):
                add(v)

    collect(m.get('ec'), ec)
    collect(m.get('rx'), rx)
    collect(m.get('met'), met)
    collect(m.get('map'), mmap)
    qn = _norm(q)
    if qn:
        for bucket in (m.get('ec'), m.get('rx'), m.get('met'), m.get('map')):
            hit = False
            for v in bucket or []:
                if qn in _norm(v) or _norm(v).endswith(':' + qn):
                    hit = True
                    break
            if hit:
                add(q)
        # Free-text name hit still counts toward ranking even without an id token.
        blob = ' '.join([m.get('name') or '', m.get('map_id') or '', m.get('lineage') or ''])
        if qn in _norm(blob) and not matched:
            add(q)
    return matched, len(matched)


def _count_label(n, singular, plural=None):
    """English count label with simple pluralisation (1 reaction, 2 reactions)."""
    n = int(n or 0)
    word = singular if n == 1 else (plural or (singular + 's'))
    return f'{n} {word}'


def _compactness(m):
    """Lower is more compact when size and hit counts already tie."""
    return _quality_rank(m.get('quality')), m.get('id') or ''


def _rep_rank(m):
    """Largest query-matching window first; then more hits; then compactness."""
    return (-int(m.get('size') or 0), -int(m.get('_hits') or 0)) + _compactness(m)


def _group_by_source_map(hits, show_all_windows=False):
    """One representative per source map; other windows kept for the badge."""
    by_map = {}
    order = []
    for m in hits:
        key = _map_key(m)
        if key not in by_map:
            by_map[key] = []
            order.append(key)
        by_map[key].append(m)
    cards = []
    for key in order:
        members = by_map[key]
        members.sort(key=_rep_rank)
        # Windows list for the badge: largest first.
        windows_sorted = sorted(members, key=lambda r: (-int(r.get('size') or 0), r.get('id') or ''))
        if show_all_windows:
            for m in members:
                cards.append((m, [m], windows_sorted))
        else:
            cards.append((members[0], windows_sorted, windows_sorted))
    return cards


def _collapse_cross_source(map_cards):
    """Collapse identical reaction sets across sources into one card with also_in."""
    by_sig = {}
    for rep, windows, all_on_map in map_cards:
        sig = rep.get('rx_sig') or f"id:{rep['id']}"
        by_sig.setdefault(sig, []).append((rep, windows, all_on_map))
    out = []
    for sig, group in by_sig.items():
        # Prefer the largest / highest-hit representative as the visible card.
        group.sort(key=lambda t: _rep_rank(t[0]))
        rep, windows, all_on_map = group[0]
        peers = [t[0] for t in group]
        out.append((rep, windows, all_on_map, peers))
    out.sort(key=lambda t: (-int(t[0].get('size') or 0), -int(t[0].get('_hits') or 0),
                            -len({_map_key(p) for p in t[3]}), t[0].get('name') or '', t[0].get('id') or ''))
    return out


def search(filters=None, page=1, per_page=PAGE_DEFAULT, index=None):
    """Filter modules and group by source map (one card per map).

    Free-text `q` auto-routes EC / reaction / metabolite / model-map shaped
    tokens. Default listing collapses windows on the same source map and links
    identical reaction sets across sources via also_in.
    """
    filters = filters or {}
    idx = index or load_index()
    rows = idx['modules']
    src = _norm(filters.get('source') or filters.get('database') or '')
    lic = _norm(filters.get('license') or filters.get('licence') or '')
    ec = filters.get('ec') or filters.get('ec_number') or ''
    rx = filters.get('reaction') or filters.get('reaction_id') or ''
    met = filters.get('metabolite') or filters.get('metabolite_id') or ''
    mmap = filters.get('map') or filters.get('map_id') or filters.get('model') or ''
    q_raw = (filters.get('q') or filters.get('query') or '').strip()
    show_variants = _truthy(filters.get('variants') or filters.get('show_variants'))
    q_kind = classify_query(q_raw) if q_raw else None
    q = q_raw
    # Route free-text into the matching field when the dedicated field is empty.
    if q and q_kind == 'ec' and not ec:
        ec, q = q, ''
    elif q and q_kind == 'reaction' and not rx:
        rx, q = q, ''
    elif q and q_kind == 'metabolite' and not met:
        met, q = q, ''
    elif q and q_kind == 'map' and not mmap:
        mmap, q = q, ''

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
        matched, n_hits = _query_matches(m, ec=ec, rx=rx, met=met, mmap=mmap, q=q_raw)
        m = dict(m)
        m['_hits'] = n_hits
        m['_matched'] = matched
        hits.append(m)

    module_total = len(hits)
    map_cards = _group_by_source_map(hits, show_all_windows=False)
    source_map_total = len(map_cards)
    if show_variants:
        # Expand every window as its own card (still carry peer maps for also_in).
        expanded = []
        collapsed = _collapse_cross_source(map_cards)
        peer_by_sig = {}
        for rep, _w, _a, peers in collapsed:
            peer_by_sig[rep.get('rx_sig') or f"id:{rep['id']}"] = peers
        for rep, windows, all_on_map in map_cards:
            for w in windows:
                sig = w.get('rx_sig') or f"id:{w['id']}"
                expanded.append((w, [w], all_on_map, peer_by_sig.get(sig) or [w]))
        items = expanded
        items.sort(key=lambda t: (-int(t[0].get('_hits') or 0), -int(t[0].get('size') or 0),
                                  t[0].get('name') or '', t[0].get('id') or ''))
    else:
        items = _collapse_cross_source(map_cards)

    per_page = max(1, min(int(per_page or PAGE_DEFAULT), 100))
    page = max(1, int(page or 1))
    total = len(items)
    pages = max(1, (total + per_page - 1) // per_page)
    if page > pages:
        page = pages
    start = (page - 1) * per_page
    slice_ = items[start:start + per_page]
    results = [
        public_card(rep, windows=windows, peers=peers, matched=rep.get('_matched') or [])
        for rep, windows, _all_on_map, peers in slice_
    ]
    return dict(
        total=total,
        module_total=module_total,
        source_map_total=source_map_total,
        page=page,
        per_page=per_page,
        pages=pages,
        deduped=not show_variants,
        q_kind=q_kind,
        q_kind_label=Q_KIND_LABELS.get(q_kind) if q_kind else None,
        count_label=f"{source_map_total:,} source maps · {module_total:,} modules",
        results=results,
        filters={k: v for k, v in filters.items() if v not in (None, '')},
    )


def public_card(m, windows=None, peers=None, matched=None, variants=None):
    """Build a result card. `windows` = same source map; `peers` = cross-source same rx set."""
    windows = windows or [m]
    peers = peers or [m]
    matched = list(matched or m.get('_matched') or [])
    # Back-compat for older call sites that pass variants=.
    if variants and windows == [m]:
        windows = variants
        peers = variants
    peer_maps = sorted({_map_key(v) for v in peers})
    also = None
    if len(peer_maps) > 1:
        sources = sorted({v.get('source') for v in peers if v.get('source')})
        n_src, n_maps = len(sources), len(peer_maps)
        also = dict(
            modules=len(peers),
            sources=n_src,
            maps=n_maps,
            label=f"also in: {_count_label(n_src, 'source')} / {_count_label(n_maps, 'map')}",
        )
    q = quality_label(m.get('quality'))
    sp = species_display(m.get('species'))
    src = m.get('source') or ''
    match_bits = ', '.join(matched[:6])
    subtitle = source_label(src)
    if match_bits:
        subtitle = f"{subtitle} · matches {match_bits}"
    size = int(m.get('size') or 0)
    windows_sorted = sorted(
        windows,
        key=lambda w: (-int(w.get('size') or 0), w.get('id') or ''),
    )
    return dict(
        id=m['id'],
        source=src,
        source_label=source_label(src),
        name=m['name'],
        map_id=m['map_id'],
        subtitle=subtitle,
        matched=matched,
        species=sp,
        license=m['license'],
        attribution=attribution_string(m.get('attribution'), src),
        size=size,
        size_label=_count_label(size, 'reaction'),
        quality=m.get('quality') or '',
        quality_label=q['label'],
        quality_tip=q['tip'],
        rx_sig=m.get('rx_sig') or '',
        also_in=also,
        windows_count=len(windows_sorted),
        windows=[dict(
            id=w['id'], source=w.get('source') or '', name=w.get('name') or '',
            map_id=w.get('map_id') or '', size=int(w.get('size') or 0),
            size_label=_count_label(w.get('size'), 'reaction'),
            license=w.get('license') or '',
            thumbnail=f"/api/atlas/modules/{w['id']}/thumbnail.svg",
        ) for w in windows_sorted],
        variants=[dict(
            id=v['id'], source=v['source'], name=v['name'], map_id=v.get('map_id') or '',
            license=v.get('license') or '',
            attribution=attribution_string(v.get('attribution'), v.get('source') or ''),
            size=v.get('size') or 0,
        ) for v in peers],
        thumbnail=f"/api/atlas/modules/{m['id']}/thumbnail.svg",
    )


def meta(index=None):
    idx = index or load_index()
    return dict(
        schema_version=idx.get('schema_version'),
        generated=idx.get('generated'),
        module_count=idx.get('module_count'),
        registry_sha256=idx.get('registry_sha256'),
        sources=sorted((idx.get('sources') or {}).items()),
        licences=sorted((idx.get('licences') or {}).items()),
        filters=['q', 'source', 'license', 'variants'],
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
    mod = load_module(module_id, index=idx)
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
    q = quality_label(mod.get('quality') or rec.get('quality'))
    # Same reaction-set siblings for provenance listing.
    sig = rec.get('rx_sig') or ''
    siblings = []
    if sig:
        for m in idx['modules']:
            if m.get('rx_sig') == sig:
                siblings.append(dict(
                    id=m['id'], source=m['source'], name=m['name'], map_id=m.get('map_id') or '',
                    license=m.get('license') or '', attribution=m.get('attribution') or '', size=m.get('size') or 0,
                ))
    else:
        siblings = [dict(
            id=mod['id'], source=mod.get('source'), name=mod.get('name'), map_id=mod.get('map_id') or '',
            license=mod.get('license') or '', attribution=rec.get('attribution') or '', size=rec.get('size') or 0,
        )]
    layout = mod.get('layout_origin') or {}
    layout_label = None
    if isinstance(layout, dict) and layout.get('status'):
        st = str(layout['status'])
        layout_label = st.replace('_', ' ')
    return dict(
        id=mod['id'],
        source=mod.get('source'),
        name=mod.get('name'),
        map_id=mod.get('map_id'),
        lineage=mod.get('lineage'),
        species=species_display(mod.get('species')),
        license=mod.get('license') or mod.get('licence'),
        attribution=attribution_string(rec.get('attribution') or mod.get('attribution'), mod.get('source') or ''),
        size=rec.get('size') or mod.get('size'),
        quality=mod.get('quality'),
        quality_label=q['label'],
        quality_tip=q['tip'],
        source_sha256=mod.get('source_sha256'),
        original_instance=mod.get('original_instance'),
        layout_origin=layout,
        layout_origin_label=layout_label,
        extraction=mod.get('extraction'),
        core_family=mod.get('core_family'),
        variant=mod.get('variant'),
        atlas_schema_version=mod.get('atlas_schema_version'),
        reactions=reactions,
        ports=mod.get('ports') or [],
        node_count=len(mod.get('nodes') or []),
        edge_count=len(mod.get('edges') or []),
        also_in=(lambda n_src, n_maps: dict(
            modules=len(siblings),
            sources=n_src,
            maps=n_maps,
            label=f"{_count_label(n_src, 'source')} / {_count_label(n_maps, 'map')}",
        ))(len({s['source'] for s in siblings}), len({(s['source'], s['name']) for s in siblings}))
        if len(siblings) > 1 else None,
        variants=[dict(s, attribution=attribution_string(s.get('attribution'), s.get('source') or '')) for s in siblings],
        thumbnail=f"/api/atlas/modules/{mod['id']}/thumbnail.svg",
    )


def thumbnail_svg(module, width=240, height=160):
    """Curated drawing preview matching docs/images/atlas/atlas_examples style.

    Whole drawing visible (never cropped): viewBox = content bbox + 6% pad with
    preserveAspectRatio=xMidYMid meet. Grey letterboxing comes from the CSS
    frame (.ab-th), not from a mismatched SVG fill. Edges/nodes use screen-pixel
    sizes via vector-effect=non-scaling-stroke (~1.2px / ~2.5px).
    """
    usable = [n for n in (module.get('nodes') or []) if n.get('xy') is not None]
    if not usable:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'width="{width}" height="{height}" preserveAspectRatio="xMidYMid meet">'
            f'<text x="12" y="{height // 2}" fill="#5b7480" font-size="12">No geometry</text></svg>'
        )
    # Drop exact coordinate duplicates that pile into unreadable clumps.
    seen_xy = set()
    deduped = []
    for n in usable:
        key = (round(float(n['xy'][0]), 4), round(float(n['xy'][1]), 4), n.get('role'))
        if key in seen_xy:
            continue
        seen_xy.add(key)
        deduped.append(n)
    usable = deduped or usable
    xs = [float(n['xy'][0]) for n in usable]
    ys = [float(n['xy'][1]) for n in usable]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    bw = (max_x - min_x) or 1.0
    bh = (max_y - min_y) or 1.0
    # 6% content pad, then grow the viewBox a little more so ~2.5px nodes are
    # not clipped at the frame after meet.
    pad = 0.06 * max(bw, bh)
    vb_x, vb_y = min_x - pad, min_y - pad
    vb_w, vb_h = bw + 2 * pad, bh + 2 * pad
    meet = min(width / vb_w, height / vb_h) or 1.0
    r_rx = round(2.5 / meet, 6)
    r_met = round(2.5 / meet, 6)
    extra = max(r_rx, r_met) + (1.2 / meet)
    vb_x -= extra
    vb_y -= extra
    vb_w += 2 * extra
    vb_h += 2 * extra
    meet = min(width / vb_w, height / vb_h) or 1.0
    r_rx = round(2.5 / meet, 6)
    r_met = round(2.5 / meet, 6)

    def pt(n):
        return (round(float(n['xy'][0]), 4), round(float(n['xy'][1]), 4))

    by_key = {n['key']: n for n in usable}
    # Transparent SVG; .ab-th grey shows through meet letterboxing (no side bands).
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb_x:.4f} {vb_y:.4f} {vb_w:.4f} {vb_h:.4f}" '
        f'width="{width}" height="{height}" preserveAspectRatio="xMidYMid meet" '
        f'role="img" aria-label="{xml_escape(str(module.get("name") or module.get("id") or "module"))}">'
    ]
    ve = 'vector-effect="non-scaling-stroke"'
    for e in module.get('edges') or []:
        a, b = by_key.get(e.get('a')), by_key.get(e.get('b'))
        if not a or not b:
            continue
        (ax, ay), (bx, by) = pt(a), pt(b)
        parts.append(
            f'<line x1="{ax}" y1="{ay}" x2="{bx}" y2="{by}" stroke="#111" stroke-width="1.2" '
            f'stroke-linecap="round" {ve}/>'
        )
    for n in usable:
        x, y = pt(n)
        title = xml_escape(str(n.get('label') or n.get('entity') or ''))
        if n.get('role') == 'reaction':
            parts.append(
                f'<circle cx="{x}" cy="{y}" r="{r_rx}" fill="#111111"><title>{title}</title></circle>'
            )
        else:
            parts.append(
                f'<circle cx="{x}" cy="{y}" r="{r_met}" fill="#ffffff" stroke="#111111" '
                f'stroke-width="1.2" {ve}><title>{title}</title></circle>'
            )
    parts.append('</svg>')
    return ''.join(parts)


def module_to_escher(module):
    """Export one curated module as Escher JSON. Coordinates are translated only (no flip/rotate)."""
    nodes_in = module.get('nodes') or []
    edges = module.get('edges') or []
    reactions_in = {str(r.get('key') or ('r:' + str(r.get('id')))): r for r in (module.get('reactions') or [])}
    usable = [n for n in nodes_in if n.get('xy') is not None]
    if not usable:
        raise ValueError('module has no drawable coordinates')
    xs = [n['xy'][0] for n in usable]
    ys = [n['xy'][1] for n in usable]
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
    idx = index or load_index()
    manifest_json, rows, missing = download_manifest(module_ids, index=idx, fmt='json')
    manifest_csv, _, _ = download_manifest(module_ids, index=idx, fmt='csv')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr('manifest.json', manifest_json)
        zf.writestr('manifest.csv', manifest_csv)
        for row in rows:
            mod = load_module(row['module_id'], index=idx)
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
.ab-card{position:relative;cursor:pointer;display:flex;flex-direction:column;min-width:0;overflow:hidden}
.ab-card.sel{outline:2px solid var(--accent)}
.ab-card .chk{position:absolute;top:10px;right:10px;z-index:2}
.ab-th{position:relative;display:block;width:100%;aspect-ratio:3/2;height:auto;margin-bottom:10px;border-radius:8px;border:1px solid var(--line);background:#fafafa;overflow:hidden}
.ab-th::before{content:"Loading…";position:absolute;inset:0;display:flex;align-items:center;justify-content:center;font-size:12px;color:var(--mute);pointer-events:none;z-index:1}
.ab-th.is-ready::before,.ab-th.is-error::before{content:none}
.ab-th.is-error::after{content:"Preview unavailable";position:absolute;inset:0;display:flex;align-items:center;justify-content:center;font-size:12px;color:var(--mute);background:#fafafa}
.ab-th img.th{position:absolute;inset:0;display:block;width:100%!important;height:100%!important;max-width:none!important;max-height:none!important;object-fit:contain!important;object-position:center;opacity:0;transition:opacity .15s ease}
.ab-th.is-ready img.th{opacity:1}
.ab-empty{color:var(--mute);padding:24px 0}
.ab-detail{display:none}.ab-detail.on{display:block}
.ab-list.on{display:none}
.ab-rx{max-height:280px;overflow:auto;font:12px/1.45 ui-monospace,Menlo,monospace}
button.primary{background:var(--accent);color:var(--accent-ink);border-color:var(--accent);font-weight:600}
.filters .grow{flex:1;min-width:min(280px,100%)}
.ab-also{font-size:12px;color:var(--mute);margin:0 0 6px}
.ab-detect{font-size:12px;color:var(--mute);margin:4px 0 0;min-height:1.2em}
.ab-wins{font-size:12px;margin:0 0 6px;min-width:0}
.ab-wins button{font:inherit;font-size:12px;border:1px solid var(--line);background:var(--soft);color:var(--ink);border-radius:999px;padding:1px 10px;cursor:pointer}
.ab-wins ul.winlist{display:block;margin:6px 0 0;padding:4px 8px 4px 18px;color:var(--mute);max-height:4.8em;overflow-x:hidden;overflow-y:auto;overscroll-behavior:contain;list-style:disc;border:1px solid var(--line);border-radius:6px;background:var(--bg);box-sizing:border-box}
.ab-wins ul.winlist li{margin:0 0 2px;padding:0;line-height:1.35;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.tip{border-bottom:1px dotted var(--mute);cursor:help}
</style>'''
    body = css_extra + r'''
<div id="listView" class="ab-list">
<form class="filters" id="filters" autocomplete="off">
<label class="grow">Search<input type="search" name="q" id="q" placeholder="1.1.1.100, HEX1, g6p_c, iJO1366"><span class="ab-detect" id="qDetect"></span></label>
<label>Source database<select name="source" id="source"><option value="">All</option></select></label>
<label>Licence<select name="license" id="license"><option value="">All</option></select></label>
<button type="submit" class="primary">Search</button>
<button type="button" id="reset">Reset</button>
</form>
<div class="ab-toolbar">
<label class="sel"><input type="checkbox" id="selPage"> Select all on page</label>
<label class="sel"><input type="checkbox" id="showVariants"> Show all windows</label>
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
var state={page:1,per:24,selected:{},lastIds:[],detail:null,cards:{}};
function esc(s){return String(s==null?'':s).replace(/[&<>"']/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]})}
function params(){var f=$('filters'),o={};['q','source','license'].forEach(function(k){var v=f[k].value.trim();if(v)o[k]=v});o.page=state.page;o.per_page=state.per;if($('showVariants').checked)o.variants='1';return o}
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
var TH_PLACEHOLDER='data:image/svg+xml,'+encodeURIComponent('<svg xmlns="http://www.w3.org/2000/svg" width="240" height="160"><rect width="100%" height="100%" fill="#f3f4f6"/></svg>');
var thumbIO=null;
function ensureThumbIO(){
 if(thumbIO)return thumbIO;
 if(!('IntersectionObserver' in window))return null;
 thumbIO=new IntersectionObserver(function(entries){
  entries.forEach(function(en){
   if(!en.isIntersecting)return;
   loadThumb(en.target);
   thumbIO.unobserve(en.target);
  });
 },{root:null,rootMargin:'200px 0px',threshold:0.01});
 return thumbIO;
}
function loadThumb(wrap){
 if(!wrap||wrap.getAttribute('data-loaded'))return;
 var img=wrap.querySelector('img.th');if(!img)return;
 var src=wrap.getAttribute('data-src')||img.getAttribute('data-src')||'';
 if(!src){wrap.classList.add('is-error');wrap.setAttribute('data-loaded','1');return}
 wrap.setAttribute('data-loaded','1');
 img.onload=function(){wrap.classList.add('is-ready');wrap.classList.remove('is-error')};
 img.onerror=function(){img.removeAttribute('src');wrap.classList.remove('is-ready');wrap.classList.add('is-error')};
 img.src=src;
}
function wireThumbs(){
 var wraps=document.querySelectorAll('.ab-th[data-src]');
 var io=ensureThumbIO();
 Array.prototype.forEach.call(wraps,function(wrap,i){
  // Eagerly load first two rows (~12) so above-the-fold never stays empty.
  if(i<12||!io){loadThumb(wrap);return}
  io.observe(wrap);
 });
 // Fallback: any wrap already near the viewport (IO can miss mid-render).
 Array.prototype.forEach.call(wraps,function(wrap){
  if(wrap.getAttribute('data-loaded'))return;
  var r=wrap.getBoundingClientRect();
  if(r.top<window.innerHeight+240&&r.bottom>-40)loadThumb(wrap);
 });
}
if(!window.__abThumbScroll){
 window.__abThumbScroll=1;
 window.addEventListener('scroll',function(){
  Array.prototype.forEach.call(document.querySelectorAll('.ab-th[data-src]:not([data-loaded])'),function(wrap){
   var r=wrap.getBoundingClientRect();
   if(r.top<window.innerHeight+240&&r.bottom>-40)loadThumb(wrap);
  });
 },{passive:true});
}
function reactionLabel(n){n=+n||0;return n===1?'1 reaction':n+' reactions'}
function card(m){
 var checked=state.selected[m.id]?' checked':'';
 var also=m.also_in?'<p class="ab-also">'+esc(m.also_in.label)+'</p>':'';
 var wins='';
 if(m.windows_count&&m.windows_count>1){
  wins='<div class="ab-wins"><button type="button" class="winbtn" data-id="'+esc(m.id)+'">'+esc(m.windows_count)+' windows</button><ul class="winlist" id="wins-'+esc(m.id)+'" hidden></ul></div>';
 }
 var thumb=m.thumbnail||'';
 var sub=m.subtitle||((m.source_label||m.source||'')+(m.matched&&m.matched.length?' · matches '+m.matched.join(', '):''));
 var sizeLab=m.size_label||reactionLabel(m.size);
 return '<article class="card ab-card'+(state.selected[m.id]?' sel':'')+'" data-id="'+esc(m.id)+'">'
  +'<label class="chk" onclick="event.stopPropagation()"><input type="checkbox" data-id="'+esc(m.id)+'"'+checked+'></label>'
  +'<div class="ab-th" data-src="'+esc(thumb)+'"><img class="th" alt="" decoding="async" data-src="'+esc(thumb)+'" src="'+TH_PLACEHOLDER+'"></div>'
  +'<h2 class="dn">'+esc(m.name||m.id)+'</h2>'
  +'<p class="meta">'+esc(sub)+'</p>'+also+wins
  +'<div class="tags"><span class="tag src-'+esc(m.source)+'">'+esc(m.source_label||m.source)+'</span>'
  +'<span class="tag">'+esc(sizeLab)+'</span>'
  +'<span class="tag lic">'+esc(m.license)+'</span></div></article>'}
function render(data){
 state.lastIds=(data.results||[]).map(function(m){return m.id});
 state.cards={};(data.results||[]).forEach(function(m){state.cards[m.id]=m});
 $('grid').innerHTML=(data.results||[]).map(card).join('')||'<p class="ab-empty">No modules match these filters.</p>';
 wireThumbs();
 $('qDetect').textContent=data.q_kind_label||detectKind($('q').value)||'';
 var label=data.count_label||((data.source_map_total||data.total||0).toLocaleString()+' source maps · '+(data.module_total||0).toLocaleString()+' modules');
 $('count').textContent=label+' · page '+data.page+' / '+data.pages;
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
function detectKind(q){
 q=String(q||'').trim();if(!q)return '';
 if(/^\d+(\.\d+){1,3}(\.-)?$/i.test(q))return 'Matched as EC number';
 if(/^chebi:/i.test(q)||/^[a-z][a-z0-9]*_[a-z0-9]+$/i.test(q))return 'Matched as metabolite ID';
 if(/^[A-Za-z][A-Za-z0-9_-]*(\.[A-Za-z0-9_.-]+)?$/.test(q)&&(q.indexOf('.')>=0||/^i[A-Z0-9]/i.test(q)||/_(core|GEM)$/i.test(q)||/^[ei][A-Z]{2}\d/i.test(q)))return 'Matched as map / model';
 if(/^(rhea|reactome|metanetx\.reaction):/i.test(q))return 'Matched as reaction ID';
 if(/^[A-Za-z][A-Za-z0-9_.-]{1,40}$/.test(q)&&q.indexOf(' ')<0&&(q===q.toUpperCase()||/\d/.test(q)))return 'Matched as reaction ID';
 return 'Matched as free text'}
$('q').addEventListener('input',function(){$('qDetect').textContent=detectKind($('q').value)});
$('filters').onsubmit=function(e){e.preventDefault();state.page=1;search()};
$('reset').onclick=function(){$('q').value='';$('source').value='';$('license').value='';$('showVariants').checked=false;$('qDetect').textContent='';state.page=1;search()};
$('showVariants').onchange=function(){state.page=1;search()};
$('pager').onclick=function(e){var a=e.target.closest('[data-p]');if(!a)return;e.preventDefault();state.page=+a.getAttribute('data-p');search()};
function toggleSel(cb){if(!cb)return;if(cb.checked)state.selected[cb.getAttribute('data-id')]=1;else delete state.selected[cb.getAttribute('data-id')];syncDl();var art=cb.closest('.ab-card');if(art)art.classList.toggle('sel',cb.checked)}
$('grid').onclick=function(e){
 var win=e.target.closest('.winbtn');
 if(win){e.preventDefault();e.stopPropagation();
  var id=win.getAttribute('data-id'),list=$('wins-'+id),card=state.cards[id];
  if(!list||!card)return;
  if(list.hidden){
   var rows=(card.windows||[]).slice().sort(function(a,b){return (+b.size||0)-(+a.size||0)||String(a.id||'').localeCompare(String(b.id||''))});
   var seen={},uniq=[];
   rows.forEach(function(w){var s=+w.size||0;if(seen[s])return;seen[s]=1;uniq.push(w)});
   list.innerHTML=uniq.map(function(w){
    var lab=w.size_label||reactionLabel(w.size);
    return '<li><a href="#" class="tomod" data-id="'+esc(w.id)+'" title="'+esc(w.id)+'">'+esc(lab)+'</a></li>';
   }).join('');
   list.hidden=false;
  }else list.hidden=true;
  return}
 var tomod=e.target.closest('a.tomod');
 if(tomod){e.preventDefault();e.stopPropagation();openDetail(tomod.getAttribute('data-id'));return}
 var cb=e.target.closest('input[type=checkbox][data-id]');
 if(cb){toggleSel(cb);return}
 var card=e.target.closest('.ab-card');if(card)openDetail(card.getAttribute('data-id'))};
$('grid').onchange=function(e){var cb=e.target.closest('input[type=checkbox][data-id]');if(cb)toggleSel(cb)};
$('selPage').onchange=function(){
 state.lastIds.forEach(function(id){if($('selPage').checked)state.selected[id]=1;else delete state.selected[id]});
 Array.prototype.forEach.call(document.querySelectorAll('.ab-card'),function(art){
  var on=!!$('selPage').checked;art.classList.toggle('sel',on);
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
  var qTip=d.quality_tip?(' title="'+esc(d.quality_tip)+'"'):'';
  var qLab=d.quality_label||d.quality||'';
  var also='';
  if(d.also_in&&d.variants&&d.variants.length>1){
   also='<h2>Also drawn in</h2><p class="meta">'+esc(d.also_in.label)+'</p><div class="tablewrap"><table><thead><tr><th>Source</th><th>Map</th><th>Licence</th><th>Attribution</th></tr></thead><tbody>'
    +d.variants.map(function(v){return '<tr><td>'+esc(v.source)+'</td><td>'+esc(v.name)+'</td><td>'+esc(v.license)+'</td><td>'+esc(v.attribution||'')+'</td></tr>'}).join('')
    +'</tbody></table></div>'}
  var spRow=d.species?'<dt>Species</dt><dd>'+esc(d.species)+'</dd>':'';
  $('detail').innerHTML='<div class="detail-grid"><div><img class="th" style="height:240px;width:100%;object-fit:contain;background:var(--bg);border:1px solid var(--line);border-radius:8px" src="'+esc(d.thumbnail)+'" alt="">'
   +'</div><div><dl class="kv">'
   +'<dt>Source</dt><dd>'+esc(d.source)+'</dd>'
   +'<dt>Map</dt><dd>'+esc(d.name)+' <span class="muted">'+esc(d.map_id)+'</span></dd>'
   +spRow
   +'<dt>Licence</dt><dd>'+esc(d.license)+'</dd>'
   +'<dt>Attribution</dt><dd>'+esc(d.attribution||'')+'</dd>'
   +'<dt>Reactions</dt><dd>'+esc(d.size)+'</dd>'
   +'<dt>Layout quality</dt><dd><span class="tip"'+qTip+'>'+esc(qLab)+'</span></dd>'
   +(d.layout_origin_label?'<dt>Layout origin</dt><dd>'+esc(d.layout_origin_label)+'</dd>':'')
   +'<dt>Source file hash</dt><dd class="muted">'+esc(d.source_sha256||'')+'</dd>'
   +'</dl></div></div>'
   +'<h2>Member reactions</h2><div class="tablewrap"><table><thead><tr><th>ID</th><th>References</th></tr></thead><tbody class="ab-rx">'+rx+'</tbody></table></div>'
   +also
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
    lead = ui.bi('Search curated Atlas modules by EC, reaction, metabolite, source map or licence, then download Escher JSON with a provenance manifest. Identical reaction sets are collapsed by default.')
    ui.write(dest / 'index.html', 'Atlas Browser', body, section='atlas', lead=lead)
    return dest / 'index.html'
