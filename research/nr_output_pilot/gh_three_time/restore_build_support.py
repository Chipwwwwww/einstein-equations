"""Restore selected build dependencies from previously pinned official layers."""
import hashlib,json,pathlib,tarfile,urllib.request,shutil,time,os
base=pathlib.Path('/workspace/scratch/spectre-pilot');out=base/'gh-executable-gate/build-support';out.mkdir(exist_ok=True)
manifest=json.load(open(base/'image-amd64-manifest.json'));records=[];size=0
indices=[0,2,3,4,5,6,11,12,13,14,15,16,17,18,19,20,21,26,27,28]
def keep(n):
 if '..' in pathlib.PurePosixPath(n).parts:return False
 if n.startswith(('usr/include/','usr/local/include/','usr/lib/gcc/','usr/share/cmake-','sxscollaboration/charm/')):return True
 if n.startswith('usr/lib/llvm-21/'):
  return '/lib/clang/' in n or pathlib.PurePosixPath(n).name in {'clang','clang++','clang-21','lld','ld.lld','llvm-ar','llvm-ranlib'}
 if n.startswith('usr/local/lib/') and n.endswith('.a'):return True
 if n.startswith('sxscollaboration/spectre/build/'):
  return '/CMakeFiles/' in n and not n.endswith(('.o','.a')) or n.endswith('.a')
 return False
for i in indices:
 layer=manifest['layers'][i];t0=time.monotonic()
 token=json.load(urllib.request.urlopen('https://auth.docker.io/token?service=registry.docker.io&scope=repository:sxscollaboration/spectre:pull',timeout=30))['token']
 req=urllib.request.Request('https://registry-1.docker.io/v2/sxscollaboration/spectre/blobs/'+layer['digest'],headers={'Authorization':'Bearer '+token})
 class Reader:
  def __init__(self,s):self.s=s;self.h=hashlib.sha256();self.bytes=0
  def read(self,n=-1):
   b=self.s.read(n);self.h.update(b);self.bytes+=len(b);return b
 r=Reader(urllib.request.urlopen(req,timeout=120));count=0;kept=0
 with tarfile.open(fileobj=r,mode='r|gz') as tar:
  for m in tar:
   n=m.name.removeprefix('./').lstrip('/')
   if not keep(n) or not (m.isfile() or m.isdir() or m.issym() or m.islnk()):continue
   if '.wh.' in pathlib.PurePosixPath(n).name:raise RuntimeError('Whiteout needs handling '+n)
   if m.isfile():
    size+=m.size;kept+=m.size
    if size>2*1024**3 or shutil.disk_usage(base).free<6*1024**3:raise RuntimeError('Build support extraction resource guard')
   m.name=n
   if m.islnk() and not (out/m.linkname.lstrip('/')).exists():continue
   tar.extract(m,out,filter='fully_trusted');count+=1
 while r.read(1024**2):pass
 if 'sha256:'+r.h.hexdigest()!=layer['digest']:raise RuntimeError('Layer hash mismatch')
 rec={'layer':i,'digest':layer['digest'],'compressed_bytes':r.bytes,'selected_bytes':kept,'entries':count,'seconds':round(time.monotonic()-t0,2)};records.append(rec);print(json.dumps(rec),flush=True)
 (base/'gh-executable-gate/restored-layers.json').write_text(json.dumps(records,indent=2)+'\n')
print('DONE selected bytes',size,flush=True)
