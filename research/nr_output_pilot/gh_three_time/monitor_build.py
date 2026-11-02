import json,os,pathlib,shutil,subprocess,time
base=pathlib.Path('/workspace/scratch/spectre-pilot/gh-executable-gate');env=os.environ.copy()
for k in ['DOCKER_HOST','DOCKER_CONTEXT','DOCKER_TLS','DOCKER_TLS_VERIFY','DOCKER_CERT_PATH']:env.pop(k,None)
cmd=['docker','--host=unix:///var/run/docker.sock']
with open(base/'resource-samples.jsonl','a') as log:
 while True:
  sample={'unix':time.time(),'free_disk_bytes':shutil.disk_usage(base).free}
  r=subprocess.run(cmd+['exec','spectre-output-pilot','/bin/bash','-c','cat /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory.events'],env=env,capture_output=True,text=True);sample['cgroup']=r.stdout.strip();sample['docker_exit']=r.returncode;log.write(json.dumps(sample)+'\n');log.flush()
  if sample['free_disk_bytes']<4*1024**3:
   (base/'resource-stop.txt').write_text('Disk guard: less than 4 GiB free; stopped the build container.\n');subprocess.run(cmd+['stop','--time','1','spectre-output-pilot'],env=env);break
  if 'exit_code=' in (base/'build-timing.txt').read_text():break
  time.sleep(5)
