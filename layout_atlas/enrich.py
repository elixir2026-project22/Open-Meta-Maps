"""Retain candidate cross-references, not unqualified chemical equivalence classes."""
import csv
from core import *
import database

def run():
 links=defaultdict(set)
 for row in csv.DictReader((DATA/'rhea/rhea2xrefs.tsv').open(),delimiter='\t'):
  if row['DB']=='REACTOME':links[row['ID'].rsplit('.',1)[0]].add(row['MASTER_ID'])
 db=database.connect();count=0
 for key,payload in db.execute("SELECT id,payload FROM modules WHERE source='reactome'").fetchall():
  m=database.unpack(payload)
  for r in m['reactions']:
   ids=set()
   for event in r.get('references',{}).get('reactome',[]):ids.update(links.get(event.rsplit('.',1)[0] if '.' in event else event,set()))
   if ids:
    r['references']['rhea']=sorted(ids);r['reference_candidate_evidence']=dict(source='Rhea 143 rhea2xrefs.tsv',method='public Reactome cross-reference; not forced equivalence',candidate_count=len(ids));count+=1
    for value in ids:db.execute('INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)',(key,r['key'],'rhea',value))
  db.execute('UPDATE modules SET payload=? WHERE id=?',(database.pack(m),key))
 db.commit();db.close();dump(OUT/'reaction_reference_enrichment.json',dict(module_reaction_records_with_rhea_candidates=count,reactome_reference_ids=len(links),chemical_equivalence_policy='Candidates require visible participant and full equation checks during retrieval'))
if __name__=='__main__':run()
