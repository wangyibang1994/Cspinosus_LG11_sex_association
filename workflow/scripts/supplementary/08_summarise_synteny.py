#!/usr/bin/env python3
"""Select the orthologous target chromosome and interval from minimap2 PAF."""
from __future__ import annotations
import csv, sys
from collections import defaultdict

if len(sys.argv) != 4:
    raise SystemExit("Usage: 08_summarise_synteny.py input.paf species output_prefix")
paf, species, prefix = sys.argv[1:]
rows=[]
with open(paf) as handle:
    for line in handle:
        f=line.rstrip().split("\t")
        if len(f)<12: continue
        qn,ql,qs,qe,strand,tn,tl,ts,te,nm,aln,mq=f[:12]
        rec=dict(species=species,query=qn,query_length=int(ql),query_start=int(qs),
                 query_end=int(qe),strand=strand,target=tn,target_length=int(tl),
                 target_start=int(ts),target_end=int(te),n_matches=int(nm),
                 alignment_length=int(aln),mapq=int(mq))
        identity=rec["n_matches"]/max(rec["alignment_length"],1)
        if rec["alignment_length"]>=2000 and rec["mapq"]>=20 and identity>=0.70: rows.append(rec)
if not rows: raise SystemExit(f"No informative alignments in {paf}")
score=defaultdict(int)
for r in rows: score[r["target"]]+=r["n_matches"]
best=max(score,key=score.get)
z=[r for r in rows if r["target"]==best]
# Remove isolated repeat-derived anchors using a weighted-median target midpoint.
ordered=sorted(z,key=lambda r:(r["target_start"]+r["target_end"])/2)
half=sum(r["n_matches"] for r in ordered)/2; acc=0; centre=0
for r in ordered:
    acc+=r["n_matches"]
    if acc>=half:
        centre=(r["target_start"]+r["target_end"])//2; break
z=[r for r in z if abs((r["target_start"]+r["target_end"])//2-centre)<=2_500_000]
start=max(0,min(r["target_start"] for r in z)-100_000)
end=min(z[0]["target_length"],max(r["target_end"] for r in z)+100_000)
with open(prefix+".blocks.tsv","w",newline="") as h:
    w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(prefix+".orthologous_interval.tsv","w",newline="") as h:
    w=csv.writer(h,delimiter="\t"); w.writerow(["species","target","start0","end","length","matched_bases"])
    w.writerow([species,best,start,end,end-start,sum(r['n_matches'] for r in z)])
with open(prefix+".region.txt","w") as h: h.write(f"{best}:{start+1}-{end}\n")
print(species,best,start+1,end,score[best],sep="\t")
