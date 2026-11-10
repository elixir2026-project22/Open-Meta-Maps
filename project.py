#!/usr/bin/env python3
"""Stable root entry point for the active atlas."""
from pathlib import Path
import sys,subprocess
ROOT=Path(__file__).resolve().parent

def main():
 args=sys.argv[1:]
 if not args or args in (['-h'],['--help']):
  print('Open Meta Maps — current atlas v0.3\n\nRead: README.md; docs/\n\nCommands:\n  python3 project.py setup        # one-time: build the edge routers (needs clang++ and network)\n  python3 project.py verify       # integrity and tests\n  python3 project.py report       # rebuild report from saved results\n  python3 project.py models --model e_coli_core\n  python3 project.py compose      # context composition search\n  python3 project.py compose --gem examples/e_coli_core.xml --out results/local_mvp/   # GEM in -> map + report out (add --ec for extended EC mode, --no-library for a from-scratch layout)\n  python3 project.py data list    # optional data packages (MicroMap, BRENDA, reference maps, full Atlas); `data fetch all` installs them\n  python3 project.py tool         # local tool page: upload a GEM, map it to expert structures, build skeleton + map (127.0.0.1:8765)\n  python3 project.py brenda       # BRENDA pathway-map extension (trial, layout origin unverified)\n  python3 project.py ec           # EC-number channel: added coverage vs complete-equation matching\n  python3 project.py geometry     # fixed-candidate geometry ablation\n  python3 project.py deploy       # EGNN deployment example\n  python3 project.py all          # full pipeline, includes data acquisition\n\nFull stage options: python3 layout_atlas/cli.py --help')
  return
 if args and args[0]=='data':
  sys.path.insert(0,str(ROOT/'layout_atlas'));import data_packages;data_packages.main(args[1:]);return
 if args==['setup']:
  # one-time build of the libavoid edge routers (needs clang++ and network access for the pinned source)
  sys.path.insert(0,str(ROOT/'layout_atlas'));import routing;routing.build();routing.build_labels();print('routers built');return
 result=subprocess.run([sys.executable,str(ROOT/'layout_atlas/cli.py'),*args],cwd=ROOT);raise SystemExit(result.returncode)
if __name__=='__main__':main()
