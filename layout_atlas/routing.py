"""Pinned libavoid build and a narrow numeric stdin bridge; no shell interpolation."""
import subprocess,zipfile,io,urllib.request
from core import *
from acquire import fetch

def build():
 folder=BASE/'build';binary=folder/'avoid_router'
 if binary.exists():return binary
 meta=json.loads(fetch('https://api.github.com/repos/mjwybrow/adaptagrams/commits/master','routing/commit.json'));commit=meta['sha'];raw=fetch('https://codeload.github.com/mjwybrow/adaptagrams/zip/'+commit,'routing/source.zip')
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  for n in z.namelist():
   parts=n.split('/')[1:]
   if len(parts)>2 and parts[0]=='cola' and parts[1] in ('libavoid','libvpsc') and not n.endswith('/'):
    target=folder/'source'/Path(*parts);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(n))
 cpp=folder/'route.cpp';cpp.write_text('''#include <iostream>
#include <vector>
#include "libavoid/libavoid.h"
int main(){int n,m;std::cin>>n>>m;Avoid::Router router(Avoid::OrthogonalRouting);std::vector<Avoid::ShapeRef*> shapes;for(int i=0;i<n;++i){double x,y;std::cin>>x>>y;Avoid::Rectangle rectangle(Avoid::Point(x-6,y-6),Avoid::Point(x+6,y+6));auto s=new Avoid::ShapeRef(&router,rectangle,i+1);auto pin=new Avoid::ShapeConnectionPin(s,1,.5,.5,true,0,Avoid::ConnDirAll);pin->setExclusive(false);shapes.push_back(s);}std::vector<Avoid::ConnRef*> connections;for(int i=0;i<m;++i){int a,b;std::cin>>a>>b;auto c=new Avoid::ConnRef(&router,Avoid::ConnEnd(shapes[a],1),Avoid::ConnEnd(shapes[b],1),n+i+1);c->setRoutingType(Avoid::ConnType_Orthogonal);connections.push_back(c);}router.processTransaction();for(int i=0;i<m;++i){const auto &p=connections[i]->displayRoute();std::cout<<i<<" "<<p.size();for(unsigned j=0;j<p.size();++j)std::cout<<" "<<p.ps[j].x<<" "<<p.ps[j].y;std::cout<<"\\n";}return 0;}
''')
 sources=sorted((folder/'source/cola/libavoid').glob('*.cpp'));cmd=['clang++','-std=c++11','-O1','-DNDEBUG','-I'+str(folder/'source/cola'),str(cpp),*map(str,sources),'-o',str(binary)];result=subprocess.run(cmd,text=True,capture_output=True,timeout=180);dump(OUT/'routing_build.json',dict(commit=commit,command=cmd,exit_code=result.returncode,stderr=result.stderr[-12000:],source_license='See source/cola/libavoid/COPYING; no relicensing',status='built' if result.returncode==0 else 'failed'))
 if result.returncode:raise RuntimeError(result.stderr[-1500:])
 return binary

def avoid(positions,edges,obstacles=None):
 if obstacles is not None:return avoid_labels(positions,edges,obstacles)
 binary=BASE/'build/avoid_router'
 if not binary.exists():raise RuntimeError('libavoid not built; run routing.py')
 keys=sorted(positions);idx={k:i for i,k in enumerate(keys)};lines=[f'{len(keys)} {len(edges)}']+[f'{positions[k][0]} {positions[k][1]}' for k in keys]+[f'{idx[a]} {idx[b]}' for a,b in edges];r=subprocess.run([str(binary)],input='\n'.join(lines)+'\n',text=True,capture_output=True,timeout=45,check=True);paths={}
 for line in r.stdout.splitlines():
  parts=line.split();i=int(parts[0]);count=int(parts[1]);assert len(parts)==2+2*count;paths[i]=[[float(parts[j]),float(parts[j+1])] for j in range(2,len(parts),2)]
 assert len(paths)==len(edges) and all(len(p)>=2 for p in paths.values());return paths
if __name__=='__main__':
 try:print('Built',build(),flush=True)
 except Exception as e:print('Build failed:',e,flush=True);raise


def build_labels():
 """Same pinned libavoid sources, with explicit label rectangles."""
 folder=BASE/'build';binary=folder/'avoid_labels_router';cpp=folder/'route_labels.cpp'
 source=(folder/'route.cpp').read_text()
 source=source.replace('int n,m;std::cin>>n>>m;', 'int n,m,k;std::cin>>n>>m>>k;')
 source=source.replace('std::vector<Avoid::ConnRef*> connections;', 'for(int j=0;j<k;++j){double x0,y0,x1,y1;std::cin>>x0>>y0>>x1>>y1;Avoid::Rectangle rectangle(Avoid::Point(x0,y0),Avoid::Point(x1,y1));new Avoid::ShapeRef(&router,rectangle,n+m+j+1);}std::vector<Avoid::ConnRef*> connections;')
 if binary.exists() and cpp.exists() and cpp.read_text()==source:return binary
 cpp.write_text(source)
 sources=sorted((folder/'source/cola/libavoid').glob('*.cpp'))
 cmd=['clang++','-std=c++11','-O1','-DNDEBUG','-I'+str(folder/'source/cola'),str(cpp),*map(str,sources),'-o',str(binary)]
 result=subprocess.run(cmd,text=True,capture_output=True,timeout=180)
 dump(OUT/'routing_labels_build.json',dict(source_commit=read(OUT/'routing_build.json')['commit'],command=cmd,exit_code=result.returncode,stderr=result.stderr[-12000:],source_sha256=sha(source.encode()),binary_sha256=sha(binary.read_bytes()) if binary.exists() and result.returncode==0 else None))
 if result.returncode:raise RuntimeError(result.stderr[-1500:])
 return binary


def avoid_labels(positions,edges,obstacles):
 binary=build_labels();keys=sorted(positions);idx={k:i for i,k in enumerate(keys)}
 assert all(len(b)==4 and b[0]<b[2] and b[1]<b[3] for b in obstacles)
 lines=[f'{len(keys)} {len(edges)} {len(obstacles)}']+[f'{positions[k][0]} {positions[k][1]}' for k in keys]+[' '.join(map(str,b)) for b in obstacles]+[f'{idx[a]} {idx[b]}' for a,b in edges]
 result=subprocess.run([str(binary)],input='\n'.join(lines)+'\n',text=True,capture_output=True,timeout=60,check=True);paths={}
 for line in result.stdout.splitlines():
  parts=line.split();i=int(parts[0]);count=int(parts[1]);assert len(parts)==2+2*count;paths[i]=[[float(parts[j]),float(parts[j+1])] for j in range(2,len(parts),2)]
 assert len(paths)==len(edges) and all(len(p)>=2 for p in paths.values())
 return paths
