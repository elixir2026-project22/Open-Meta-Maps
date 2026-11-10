"""Core stages and reproducible, explicitly separated research experiments."""
import argparse,json
from core import *

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','acquire','archive','routing','build','catalog','sensitivity','ranking','tasks','learning','models','geometry','micromap','brenda','ec','tool','compose','deploy','verify','report','all']);p.add_argument('--model');p.add_argument('--port',type=int,default=8765);p.add_argument('--gem');p.add_argument('--out');p.add_argument('--ec',action='store_true');p.add_argument('--no-library',action='store_true');p.add_argument('--library',choices=['windows','units_sub'],default='windows');p.add_argument('--reaction-style',choices=['line','marker'],default='line');p.add_argument('--completion',choices=['spaced','classic'],default='spaced');p.add_argument('--source-paths',action='store_true');p.add_argument('--keep-running',action='store_true');p.add_argument('--layout',choices=['global','tiled'],default='global');p.add_argument('--orientation',choices=['off','sides','full'],default='full');p.add_argument('--block-layout',choices=['auto','single'],default='single');p.add_argument('--instances',choices=['auto','off'],default='auto');p.add_argument('--context',choices=['auto','off'],default='auto');p.add_argument('--layout-style',choices=['legacy','expert'],default='legacy');p.add_argument('--completion-context',choices=['auto','off'],default='auto');p.add_argument('--edge-routing',choices=['libavoid','hybrid','none'],default='libavoid');p.add_argument('--condition',default='all_geometry');p.add_argument('--ranker-source',choices=['community','reactome','wikipathways']);args=p.parse_args()
 def stage(name):
  if name=='compose':
   if args.gem:
    import compose_cli;return compose_cli.run(args.gem,args.out or 'results/local_mvp',ec=args.ec,no_library=args.no_library,library=args.library,reaction_style=args.reaction_style,completion=args.completion,source_paths=args.source_paths,layout=args.layout,orientation=args.orientation,block_layout=args.block_layout,instances=args.instances,context=args.context,edge_routing=args.edge_routing,completion_context=args.completion_context,layout_style=args.layout_style)
   import composition_experiment;return composition_experiment.run(args.model)
  if name=='micromap':
   import micromap_experiment,micromap_index;result=micromap_experiment.run(args.model);micromap_index.publish(result);return result
  if name=='tool':
   import tool_server;return tool_server.main(args.port,keep_running=args.keep_running)
  if name=='brenda':
   import subprocess;return subprocess.run([sys.executable,str(BASE/'brenda_library.py')],check=True)
  if name=='ec':
   import ec_experiment;return ec_experiment.run([args.model] if args.model else None)
  if name=='geometry':
   import geometry_experiment;return geometry_experiment.run(args.model or 'e_coli_core')
  if name=='freeze':return freeze()
  if name=='acquire':
   import acquire;return acquire.acquire()
  if name=='archive':
   import archive;archive.inspect_archive();return archive.extract_metabolism()
  if name=='routing':
   import routing;return routing.build()
  if name=='build':
   import database,enrich,diagnostics;database.build();enrich.run();diagnostics.run()
   import source_audit;return source_audit.run()
  if name=='sensitivity':
   import sensitivity;return sensitivity.run()
  if name=='catalog':
   import gallery;return gallery.make()
  if name in ('ranking','tasks','models'):
   import database;mods=database.modules()
   if name=='ranking':
    import experiments;return experiments.ranking(mods)
   if name=='tasks':
    import completion;return completion.prepare(mods)
   import experiments
   if not args.model:return experiments.model_experiments(mods)
   model=next(m for m in experiments.MODELS if m[0]==args.model);return experiments.run_model(*model,mods,condition=args.condition,ranker_source=args.ranker_source)
  if name=='learning':
   import neural;return neural.evaluate()
  if name=='deploy':
   import end_to_end;return end_to_end.run(name=args.model or 'e_coli_core',condition=args.condition)
  if name=='verify':
   import verify;return verify.verify()
  if name=='report':
   import report;return report.generate()
 if args.stage=='all':
  for name in ['freeze','acquire','archive','routing','build','catalog','sensitivity','ranking','tasks','learning','models','geometry','micromap','compose','deploy','verify','report']:
   print('STAGE',name,flush=True);stage(name)
 else:stage(args.stage)
if __name__=='__main__':main()
