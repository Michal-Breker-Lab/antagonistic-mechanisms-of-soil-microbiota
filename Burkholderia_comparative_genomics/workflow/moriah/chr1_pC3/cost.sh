#!/bin/bash
# Total billed core-hours across every job this rebuild has submitted.
for j in 45902096 45902097 45902685 45902696 45902725 45903743; do
  sacct -j $j --format=JobID,State,Elapsed,AllocCPUS -n -P 2>/dev/null | grep '\.batch'
done > /tmp/allacct.$$
python3 - /tmp/allacct.$$ <<'PY'
import sys
ch=0.0; n=0
for ln in open(sys.argv[1]):
    f=ln.strip().split("|")
    if len(f)<4: continue
    el=f[2]; cpus=int(f[3] or 0)
    d=0
    if "-" in el: d,el = el.split("-"); d=int(d)
    p=el.split(":")
    try: secs=d*86400+int(p[0])*3600+int(p[1])*60+float(p[2])
    except Exception: secs=0
    ch += secs/3600*cpus; n+=1
print(f"  billed task-steps : {n}")
print(f"  core-hours        : {ch:,.1f}")
print(f"  cost @ $0.008     : ${ch*0.008:,.2f}")
PY
rm -f /tmp/allacct.$$
