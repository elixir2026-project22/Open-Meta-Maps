"""Bounded HTTP-range access to the official fixed Reactome v97 ZIP archive."""
import io,urllib.request,zipfile,json
from core import *
from acquire import fetch
class Remote(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.transferred=0
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos)
  if n<=0:return b''
  if n>64*1024*1024 or self.transferred+n>128*1024*1024:raise ValueError('Archive range extraction budget 128MiB')
  request=urllib.request.Request(self.url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}','Accept-Encoding':'identity'})
  with urllib.request.urlopen(request,timeout=40) as r:
   if r.status!=206:raise RuntimeError('Server did not honor Range; refusing multi-gigabyte download')
   raw=r.read(n+1)
  if len(raw)!=n:raise ValueError('Unexpected range length')
  self.pos+=n;self.transferred+=n;return raw

def inspect_archive():
 meta=json.loads(fetch('https://zenodo.org/api/records/21383214','reactome/zenodo-v97.json'));file=meta['files'][0];remote=Remote(file['links']['self'],file['size']);result=dict(record=meta['id'],version=meta['metadata']['version'],checksum=file.get('checksum'),bytes=file['size'])
 try:
  with zipfile.ZipFile(remote) as z:
   members=z.namelist();result['members']=members
   wanted=[n for n in members if n.endswith(('ReactomePathways.txt','ReactomePathwaysRelation.txt'))]
   for name in wanted:
    raw=z.read(name);path=DATA/'reactome/archive-v97'/Path(name).name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
   result.update(status='range_indexed',extracted=wanted,transferred=remote.transferred)
 except Exception as e:result.update(status='unavailable_for_bounded_extraction',error=str(e),transferred=remote.transferred)
 dump(OUT/'reactome_archive.json',result);print(result['status'],result.get('error',''),flush=True)
if __name__=='__main__':inspect_archive()

def extract_metabolism():
 import struct,zlib,concurrent.futures,time
 meta=read(DATA/'reactome/zenodo-v97.json');file=meta['files'][0];remote=Remote(file['links']['self'],file['size'])
 relation=(DATA/'reactome/archive-v97/ReactomePathwaysRelation.txt').read_text();edges=[line.split() for line in relation.splitlines() if len(line.split())==2];g=nx.DiGraph(edges);wanted=nx.descendants(g,'R-HSA-1430728')|{'R-HSA-1430728'}
 with zipfile.ZipFile(remote) as z:
  infos={i.filename:i for i in z.infolist()};ids=sorted(r for r in wanted if '97/diagram/'+r+'.json' in infos and '97/diagram/'+r+'.graph.json' in infos)
 entries=[];errors=[]
 def get_member(name):
  dest=DATA/'reactome/archive-v97/diagram'/Path(name).name;info=infos[name]
  if dest.exists():raw=dest.read_bytes()
  else:
   time.sleep(1.55)
   r=Remote(file['links']['self'],file['size']);r.seek(info.header_offset);chunk=r.read(min(info.compress_size+30+len(info.filename.encode())+1024,file['size']-info.header_offset));header=chunk[:30];fields=struct.unpack('<4s5H3I2H',header);fn,extra=fields[-2:];compressed=chunk[30+fn+extra:30+fn+extra+info.compress_size];raw=zlib.decompress(compressed,-15) if info.compress_type==zipfile.ZIP_DEFLATED else compressed if info.compress_type==zipfile.ZIP_STORED else None
   if raw is None:raise ValueError('Unsupported ZIP compression')
   dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
  if len(raw)!=info.file_size or zlib.crc32(raw)&0xffffffff!=info.CRC:raise ValueError('Archive member CRC/length mismatch')
  json.loads(raw);return dict(path=str(dest.relative_to(BASE)),sha256=sha(raw),zip_crc32=info.CRC,entry=name,record=meta['id'],release='v97')
 def get(rid):
  try:
   a=get_member('97/diagram/'+rid+'.json');b=get_member('97/diagram/'+rid+'.graph.json');return dict(id=rid,layout_path=a['path'],graph_path=b['path'],layout_sha256=a['sha256'],graph_sha256=b['sha256'],archive_record=meta['id'],release='v97',members=[a,b]),None
  except Exception as e:return None,dict(id=rid,error=str(e))
 with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
  for ok,err in pool.map(get,ids):
   if ok:entries.append(ok)
   if err:errors.append(err)
 acq=read(DATA/'acquisition.json');acq['live_discovery_diagrams']=len(acq['reactome_diagrams']);acq['reactome_diagrams']=entries;acq['reactome_release_archive_status']='Official Zenodo record 21383214 v97; selected members CRC32 and SHA256 verified; whole 4.7GB archive checksum not verified';acq['failures'].extend(dict(source='reactome_archive_member',**e) for e in errors);dump(DATA/'acquisition.json',acq);dump(OUT/'reactome_archive_extraction.json',dict(record=21383214,version='v97',metabolic_descendants=len(wanted),eligible_diagrams=len(ids),downloaded=len(entries),errors=errors));print('Fixed Reactome v97 diagrams',len(entries),'of',len(ids),flush=True)
