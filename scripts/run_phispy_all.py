#!/usr/bin/env python3
import csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
BASE="/home/bioinfo/cdprophage"
PHISPY=os.path.expanduser("~/.local/bin/phispy")
GBFF=os.path.join(BASE,"genomes","gbff")
OUT=os.path.join(BASE,"phispy")
LOG=os.path.join(BASE,"scaleup_241","logs")
os.makedirs(LOG,exist_ok=True)
env=os.environ.copy()
tools=os.path.join(BASE,"tools","root")
env["PATH"]=os.path.join(tools,"usr","bin")+os.pathsep+os.path.expanduser("~/.local/bin")+os.pathsep+env["PATH"]
env["LD_LIBRARY_PATH"]=os.path.join(tools,"usr","lib","ncbi-blast+")+os.pathsep+os.path.join(tools,"usr","lib","x86_64-linux-gnu")
accs=[r["assembly_accession"] for r in csv.DictReader(open(os.path.join(BASE,"genome_metadata.csv")))]
def run(acc):
    outdir=os.path.join(OUT,acc); coord=os.path.join(outdir,"prophage_coordinates.tsv")
    logfile=os.path.join(LOG,"phispy_"+acc+".log")
    if os.path.exists(coord): return acc,"skip",""
    gbff=os.path.join(GBFF,acc+".gbff")
    if not os.path.exists(gbff): return acc,"missing",""
    cmd=[PHISPY,"-o",outdir,"--threads","2","--quiet",gbff]
    try:
        with open(logfile,"w") as log:
            proc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=2400)
        if os.path.exists(coord): return acc,"ok",""
        return acc,"failed",f"exit={proc.returncode}"
    except subprocess.TimeoutExpired: return acc,"timeout",""
    except Exception as exc: return acc,"error",str(exc)
res=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    futs={pool.submit(run,a):a for a in accs}
    for i,f in enumerate(as_completed(futs),1):
        acc,status,note=f.result(); res.append((acc,status,note))
        print(f"[{i}/{len(accs)}] {acc}\t{status}\t{note}",flush=True)
summary={}
for _,s,_ in res: summary[s]=summary.get(s,0)+1
print("SUMMARY",summary,flush=True)
