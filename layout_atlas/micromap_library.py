"""Independent, provenance-preserving NC-SA MicroMap layout extension."""
import gzip
import itertools
import sqlite3
import html
from core import *
import database
import celldesigner

EXT=DATA/'extensions/micromap'
RESULTS=OUT/'micromap_extension'
CONFIG=dict(larger_window_roots=384,window_sizes=[4,8,12,24],
            max_module_nodes=100,variant_threshold=.25,variant_representatives=12)


def extract(diagram):
    rx={r['key']:r for r in diagram['events']};g=nx.Graph()
    def core(r):
        return [k for k in r['visible'] if diagram['glyphs'][k]['entity'] and chem(diagram['glyphs'][k]['entity']) not in CURRENCY]
    for key,r in rx.items():
        g.add_node(key,role='reaction')
        for a in core(r):g.add_node(a,role='metabolite');g.add_edge(key,a)
    subsets={(key,) for key in rx};pairs=set()
    for k in sorted(g):
        if k in rx:continue
        pairs.update(tuple(sorted(p)) for p in itertools.combinations(sorted(g[k]),2))
    subsets.update(pairs)
    roots=sorted(rx);step=max(1,math.ceil(len(roots)/CONFIG['larger_window_roots']))
    for root in roots[::step]:
        order=[k for k in nx.bfs_tree(g,root) if k in rx][:max(CONFIG['window_sizes'])]
        for n in CONFIG['window_sizes']:
            if len(order)>=n:subsets.add(tuple(sorted(order[:n])))
    all_uses=defaultdict(set)
    for record in diagram['raw_reactions']:
        for link in record['glyph_links']:
            if link['valid']:all_uses[link['alias']].add('r:'+record['id'])
    modules=[];excluded=[]
    for rs in sorted(subsets):
        keys=set(rs)
        for r in rs:keys.update(core(rx[r]))
        if len(keys)>CONFIG['max_module_nodes']:excluded.append(dict(reactions=list(rs),reason='node_budget'));continue
        sub=g.subgraph(keys)
        if len(keys)<=len(rs) or not nx.is_connected(sub):excluded.append(dict(reactions=list(rs),reason='empty_or_disconnected_core'));continue
        nodes=[]
        for k in sorted(keys):
            if k in rx:nodes.append(dict(key=k,role='reaction',xy=rx[k]['xy'],entity=rx[k]['id'],observed=False,position_origin=rx[k]['position_origin']))
            else:
                n=diagram['glyphs'][k];nodes.append(dict(key=k,role='metabolite',xy=n['xy'],entity=n['entity'],observed=True,position_origin='observed',source_alias=k,label=n['label']))
        original={n['key']:n['xy'] for n in nodes};normalized=normalize(original)
        for n in nodes:n['xy']=normalized[n['key']]
        ports=[dict(key=k,outside_reactions=sorted(all_uses[k]-set(rs)),unparsed_outside_reactions=sorted(n for n in all_uses[k]-set(rs) if n not in rx),scope='all_drawn_source_alias_links') for k in sorted(keys-set(rs)) if all_uses[k]-set(rs)]
        edges=[dict(a=r,b=a,flow='substrate' if rx[r]['visible'][a]<0 else 'product',coefficient=rx[r]['visible'][a]) for r in rs for a in core(rx[r])]
        # Raw absolute coordinates and paths are not averaged or overwritten.
        modules.append(dict(id=uid(['micromap-instance',rs]),source='micromap',species=diagram['species'],
            lineage=diagram['lineage'],map_id=diagram['id'],name=diagram['name'],
            license=diagram['license'],source_citation=diagram['citation'],partition=diagram['partition'],admitted_to_open_core=False,
            quality='semantic_transfer',size=len(rs),nodes=nodes,edges=edges,
            reactions=[rx[k] for k in rs],attachments=[l for k in rs for l in rx[k]['all_glyph_links'] if l['kind'] in ('reactantLink','productLink')],
            ports=ports,source_sha256=diagram['source_sha256'],original_instance=diagram['path'],
            original_coordinates=original,review_status='not_required',
            compartment_evidence='labels only; raw SBML compartment is default',
            biochemical_visibility='full SBML equation and observed primary/side glyph identities; only primary glyphs in core'))
    return modules,excluded


def build():
    EXT.mkdir(parents=True,exist_ok=True);RESULTS.mkdir(parents=True,exist_ok=True)
    path=DATA/'references/micromap/MicroMap.xml'
    assert sha(path.read_bytes())==read(str(path)+'.source.json')['sha256']
    d=celldesigner.decode(path)
    metadata=read(DATA/'references/micromap/dataset-v4.0.json')['data']
    fields={r['typeName']:r['value'] for r in metadata['metadataBlocks']['citation']['fields']}
    d['citation']=dict(doi='10.7910/DVN/FZKMJ8',version='4.0',title=fields['title'],authors=[r['authorName']['value'] for r in fields['author']],license=metadata['license'],source_map=d['id'],source_sha256=d['source_sha256'],changes='Source modules extracted and normalized; target layouts adapted and algorithmically completed')
    dump(EXT/'SOURCE.json',d['citation'])
    with gzip.open(EXT/'observations.json.gz','wt') as f:json.dump(d,f,ensure_ascii=False)
    mods,excluded=extract(d);families=defaultdict(list);familyrows=[];variant_timeout=0
    db=sqlite3.connect(EXT/'atlas.sqlite')
    db.executescript('''CREATE TABLE IF NOT EXISTS sources(id TEXT PRIMARY KEY,payload BLOB);
    CREATE TABLE IF NOT EXISTS modules(id TEXT PRIMARY KEY,size INTEGER,core_family TEXT,variant TEXT,payload BLOB);
    CREATE TABLE IF NOT EXISTS reactions(module_id TEXT,reaction_key TEXT,source_label TEXT,PRIMARY KEY(module_id,reaction_key));
    CREATE INDEX IF NOT EXISTS family_idx ON modules(core_family);
    CREATE INDEX IF NOT EXISTS reaction_idx ON reactions(source_label);
    DELETE FROM sources;DELETE FROM modules;DELETE FROM reactions;''')
    db.execute('INSERT INTO sources VALUES(?,?)',(d['id'],database.pack({k:v for k,v in d.items() if k not in ('events','glyphs','raw_reactions','raw_species')})))
    for i,m in enumerate(mods):
        gm=database.graph(m);h=nx.weisfeiler_lehman_graph_hash(gm,node_attr='role');bucket=(len(gm),gm.number_of_edges(),h)
        family=None
        for f in families[bucket]:
            if database.mappings(gm,f['graph'],1):family=f;break
        if family is None:
            family=dict(id=uid(['micromap-core',bucket,m['id']]),graph=gm,variants=[],instances=[]);families[bucket].append(family);familyrows.append(family)
        variant=None
        for v in family['variants'][:CONFIG['variant_representatives']]:
            mappings=database.mappings(gm,v['graph'],4)
            if mappings and min(database.geometry_error(m,v['representative'],mp) for mp in mappings)<CONFIG['variant_threshold']:variant=v;break
        if variant is None:
            variant=dict(id=uid([family['id'],m['id']]),graph=gm,representative=m,instances=[]);family['variants'].append(variant)
        family['instances'].append(m['id']);variant['instances'].append(m['id'])
        m['core_family']=family['id'];m['variant']=variant['id']
        arr=np.array([n['xy'] for n in m['nodes']]);m['geometry_fingerprint']=uid(np.round(np.sort(np.linalg.norm(arr[:,None]-arr[None,:],axis=2).ravel()),3).tolist())
        db.execute('INSERT INTO modules VALUES(?,?,?,?,?)',(m['id'],m['size'],m['core_family'],m['variant'],database.pack(m)))
        for r in m['reactions']:db.execute('INSERT INTO reactions VALUES(?,?,?)',(m['id'],r['key'],r['id']))
        if i%2000==0:print('MicroMap index',i,'/',len(mods),flush=True)
    db.commit();db.close()
    with gzip.open(EXT/'modules.json.gz','wt') as f:json.dump(mods,f,ensure_ascii=False)
    dump(EXT/'families.json',dict(families=[dict(id=f['id'],instances=f['instances'],variants=[dict(id=v['id'],template=v['representative']['id'],instances=v['instances']) for v in f['variants']]) for f in familyrows]))
    summary=dict(source_map=d['id'],source_sha256=d['source_sha256'],dataset_version='4.0',dataset_doi='10.7910/DVN/FZKMJ8',
        license=d['license'],partition=d['partition'],admitted_to_open_core=False,
        independent_designs=1,source_reactions=len(d['raw_reactions']),decoded_reactions=len(d['events']),
        glyphs=len(d['glyphs']),parser_audit=d['audit'],modules=len(mods),sizes=dict(Counter(m['size'] for m in mods)),
        families=len(familyrows),variants=sum(len(f['variants']) for f in familyrows),
        singleton_families=sum(len(f['instances'])==1 for f in familyrows),
        excluded_module_instances=len(excluded),config=CONFIG,
        code_sha256={p:sha((BASE/p).read_bytes()) for p in ['celldesigner.py','micromap_library.py']},
        limits=['One master design, not thousands of independent layouts or species.',
                'Only nonbranch direct geometry decoded; source controls are center-line coordinates.',
                'Process center recovered from source segment, never observed expert coordinate supervision.',
                'Subsystem spreadsheet contains names only; no fabricated reaction membership.',
                'Source paths retained; target exports use algorithmic routing.'])
    dump(RESULTS/'library_summary.json',summary);dump(RESULTS/'excluded.json',dict(reactions=d['excluded'],modules=excluded))
    catalog(mods,summary)
    return mods,summary


def load():
    with gzip.open(EXT/'modules.json.gz','rt') as f:return json.load(f)


def catalog(mods,summary):
    folder=RESULTS/'catalog';folder.mkdir(parents=True,exist_ok=True)
    # Representatives selected deterministically across sizes; not an expert endorsement.
    chosen=[]
    for size in sorted({m['size'] for m in mods}):chosen.extend(sorted([m for m in mods if m['size']==size],key=lambda m:m['id'])[:6])
    page='<!doctype html><meta charset="utf-8"><style>body{font:16px system-ui;margin:35px;color:#17384d}article{display:inline-block;vertical-align:top;width:420px;margin:18px}svg{width:400px;height:260px}small{display:block;word-break:break-all}</style><h1>MicroMap research layout extension</h1><p>CC BY-NC-SA4.0. One independent master design. Not in the open core. Deterministic first6 IDs per size; observed glyphs and decoded center-line controls, algorithm-derived process centers.</p>'
    for m in chosen:
        coords=m['original_coordinates'];arr=np.array(list(coords.values()));mn=arr.min(0)-30;extent=np.ptp(arr,axis=0)+60;nodes={n['key']:n for n in m['nodes']}
        svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{mn[0]} {mn[1]} {extent[0]} {extent[1]}">']
        for r in m['reactions']:
            pts=r['path_centerline'];svg.append('<polyline points="'+' '.join(f'{p[0]},{p[1]}' for p in pts)+'" fill="none" stroke="#93aab5" stroke-width="2"/>')
        for k,v in coords.items():svg.append(f'<circle cx="{v[0]}" cy="{v[1]}" r="5" fill="'+('#ce7c34' if nodes[k]['role']=='reaction' else '#128781')+'"/>')
        svg.append('</svg>');(folder/(m['id']+'.svg')).write_text(''.join(svg))
        page+=f'<article><h2>{m["size"]} reactions</h2>'+''.join(svg)+f'<small>Instance {m["id"]}<br>Family {m["core_family"]}<br>Variant {m["variant"]}<br>Source MicroMap master v4.0 / CC BY-NC-SA4.0</small></article>'
    (folder/'index.html').write_text(page)

if __name__=='__main__':print(build()[1])
