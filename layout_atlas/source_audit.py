"""Reference-source inspection without promoting unresolved licenses or lineage."""
import xml.etree.ElementTree as ET
from core import *

def run():
 refs=read(DATA/'acquisition.json')['references'];rows=[]
 for r in refs:
  if not r.get('path','').endswith('.json'):continue
  path=BASE/r['path'];a=read(path)
  if not isinstance(a,list) or len(a)!=2:continue
  nodes=a[1].get('nodes',{});rx=a[1].get('reactions',{});rows.append(dict(source=r['source'],file=r['path'],sha256=sha(path.read_bytes()),map_name=a[0].get('map_name'),reaction_ids=sorted(v.get('bigg_id','') for v in rx.values()),reaction_count=len(rx),glyph_count=len(nodes),source_annotations=a[0],element_lineage_classification={k:'unknown_original_vs_moved_vs_added_vs_algorithmic' for k in nodes},license_status='file-level conflict unresolved; reference partition only',training_allowed=False))
 historic=DATA/'references/yeast_original/Yeast8.xml';summary={}
 if historic.exists():
  root=ET.parse(historic).getroot();count=Counter(local(e) for e in root.iter());summary=dict(path=str(historic.relative_to(BASE)),sha256=sha(historic.read_bytes()),reaction_elements=count['reaction'],species_elements=count['species'],species_aliases=count['speciesAlias'],complex_aliases=count['complexSpeciesAlias'],note='CellDesigner source recovered. No verified glyph correspondence to script-modified Escher exports yet; elemental lineage remains unknown.')
 result=dict(escher_reference_files=len(rows),by_source=dict(Counter(r['source'] for r in rows)),historic_yeast=summary,human_reviews=0,rules_history=read(DATA/'references/rules_history.json'),files=rows);dump(OUT/'reference_source_audit.json',result);return {k:v for k,v in result.items() if k!='files'}
if __name__=='__main__':print(run())
