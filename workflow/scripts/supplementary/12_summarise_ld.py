#!/usr/bin/env python3
"""Summarise PLINK pairwise LD by physical-distance bins."""
from __future__ import annotations
import csv, statistics, sys

if len(sys.argv)<3 or (len(sys.argv)-2)%2:
    raise SystemExit('Usage: 12_summarise_ld.py output.tsv label1 file1.ld [label2 file2.ld ...]')
out=sys.argv[1]; pairs=list(zip(sys.argv[2::2],sys.argv[3::2]))
edges=[0,10_000,25_000,50_000,100_000,250_000,500_000,700_001]
rows=[]
for label,path in pairs:
    bins=[[] for _ in range(len(edges)-1)]
    with open(path) as h:
        header=h.readline().split(); ix={k:i for i,k in enumerate(header)}
        for line in h:
            f=line.split()
            try: dist=abs(int(f[ix['BP_B']])-int(f[ix['BP_A']])); r2=float(f[ix['R2']])
            except (ValueError,KeyError,IndexError): continue
            for i in range(len(edges)-1):
                if edges[i]<=dist<edges[i+1]: bins[i].append(r2); break
    for i,z in enumerate(bins):
        rows.append({'group':label,'distance_start_bp':edges[i],'distance_end_bp':edges[i+1]-1,
                     'n_pairs':len(z),'mean_r2':sum(z)/len(z) if z else 'NA',
                     'median_r2':statistics.median(z) if z else 'NA',
                     'proportion_r2_ge_0.8':sum(x>=0.8 for x in z)/len(z) if z else 'NA'})
with open(out,'w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)
