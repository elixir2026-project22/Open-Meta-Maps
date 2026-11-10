"""Cross-check the decoded BRENDA SVG against BRENDA's public SPARQL pathway graph (read-only, one query)."""
import json, re, urllib.parse, urllib.request
from collections import Counter, defaultdict
from core import *
import brenda_library as bl

ENDPOINT = 'https://sparql.dsmz.de/api/brenda'
QUERY = '''SELECT ?pw ?label ?rx WHERE { ?pw a <https://purl.dsmz.de/schema/Pathway> . ?pw <http://www.w3.org/2000/01/rdf-schema#label> ?label .
OPTIONAL { ?pw <https://purl.dsmz.de/schema/hasReaction> ?rx } }'''
SNAP = bl.EXT / 'sparql_pathways.json'


def fetch():
    url = ENDPOINT + '?' + urllib.parse.urlencode({'query': QUERY})
    req = urllib.request.Request(url, headers={'Accept': 'application/sparql-results+json', 'User-Agent': 'open-pathway-maps-validation'})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = json.load(r)
    rows = [{k: v['value'] for k, v in b.items()} for b in data['results']['bindings']]
    dump(SNAP, dict(endpoint=ENDPOINT, query=QUERY, retrieved='2026-11-09', rows=rows))
    return rows


def norm(s):
    return re.sub(r'[^a-z0-9]+', '-', s.lower()).strip('-')


def run(refresh=False):
    rows = fetch() if refresh or not SNAP.exists() else read(SNAP)['rows']
    sp = defaultdict(lambda: dict(label=None, reactions=set()))
    for r in rows:
        e = sp[r['pw']]
        e['label'] = r['label']
        if r.get('rx'):
            e['reactions'].add(r['rx'])
    d = bl.decode()
    svg = {}
    for slug, p in d['pathways'].items():
        ens = [n for n in d['nodes'].values() if n['pathway'] == slug and n['kind'] == 'nodeEN']
        svg[slug] = dict(title=p.get('title') or slug, reactions=len(ens), tokens={n['reaction_token'] for n in ens if n.get('reaction_token')})
    by_label = defaultdict(list)
    for pw, e in sp.items():
        by_label[norm(e['label'])].append(pw)
    exact = 0; matched = []; unmatched_svg = []
    for slug, s in svg.items():
        keys = [norm(s['title']), slug]
        hit = next((k for k in keys if k in by_label), None)
        if hit and len(by_label[hit]) == 1:
            matched.append((slug, by_label[hit][0]));
        else:
            unmatched_svg.append(slug)
    diffs = [(slug, svg[slug]['reactions'], len(sp[pw]['reactions'])) for slug, pw in matched]
    same = sum(a == b for _, a, b in diffs)
    # numeric id overlap between svg reaction tokens and SPARQL reaction ids (e.g. BR37041 vs reaction/I/36918)
    sp_numbers = {pw: {m[1] for x in e['reactions'] if (m := re.search(r'/([IS])/(\d+)$', x))} for pw, e in sp.items()}
    overlap = [(len({re.sub(r'\D', '', t) for t in svg[slug]['tokens']} & sp_numbers[pw]), len(svg[slug]['tokens'])) for slug, pw in matched]
    result = dict(sparql_pathways=len(sp), sparql_reaction_links=sum(len(e['reactions']) for e in sp.values()), svg_pathways=len(svg),
                  svg_enzyme_nodes=sum(s['reactions'] for s in svg.values()), matched_by_name=len(matched), unmatched_svg_pathways=unmatched_svg,
                  unmatched_sparql_pathways=[e['label'] for pw, e in sp.items() if pw not in {b for _, b in matched}],
                  matched_with_equal_reaction_count=same,
                  total_abs_count_difference=sum(abs(a - b) for _, a, b in diffs),
                  largest_differences=sorted(((s, a, b) for s, a, b in diffs), key=lambda x: -abs(x[1] - x[2]))[:8],
                  svg_tokens_total=sum(b for _, b in overlap), svg_tokens_numeric_found_in_sparql_ids=sum(a for a, _ in overlap),
                  note='SVG reaction tokens (rIDBR/BS) and SPARQL reaction URIs (reaction/I|S/n) use related but not identical identifiers; numeric overlap is reported, not assumed.')
    dump(bl.RESULTS / 'sparql_validation.json', result)
    return result


if __name__ == '__main__':
    print(json.dumps(run(refresh=True), indent=1))
