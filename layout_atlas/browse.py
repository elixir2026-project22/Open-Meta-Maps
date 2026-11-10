"""Rebuild every browsable HTML page from saved results only (no experiments, no data acquisition)."""
import re
from core import *
import ui,pages,navigation,gallery,micromap_report,toolpage

def regenerate():
 done=[]
 done.append(('catalog',gallery.make()))
 (OUT/'composition/index.html').write_text(pages.composition_index(read(OUT/'composition/summary.json')))
 micromap_report.generate()
 g=read(OUT/'geometry_p0/e_coli_core/summary.json');d=OUT/'geometry_p0/e_coli_core';(d/'index.html').write_text(pages.geometry_index(g['results'],d))
 for folder in sorted((OUT/'models').glob('*/*/')):
  if not (folder/'summary.json').exists() or not (folder/'views.json').exists():continue
  s=read(folder/'summary.json');rows=read(folder/'views.json')
  if s.get('status')!='ok':continue
  (folder/'index.html').write_text(pages.view_index(folder,s['model'],s['condition'],s['reactions'],s['strict_equation_reactions'],s['reference_compatible_reactions'],s['module_reactions'],s['bridge_reactions'],rows))
 ui.rewrap(OUT/'cases/index.html','Deterministically selected cases',section='cases')
 ui.rewrap(OUT/'micromap_extension/catalog/index.html','MicroMap research layout extension',section='micromap');ui.localize_micromap(OUT/'micromap_extension/catalog/index.html')
 done.append(toolpage.make());done.append(navigation.generate());return done
if __name__=='__main__':print(regenerate())
