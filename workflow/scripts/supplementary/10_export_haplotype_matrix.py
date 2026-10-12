#!/usr/bin/env python3
"""Export discovery-oriented high-delta phased haplotypes for tree analysis."""
from __future__ import annotations
import csv, gzip, sys

if len(sys.argv)!=5:
    raise SystemExit("Usage: 10_export_haplotype_matrix.py phased.vcf.gz marker_stats.tsv haplotypes.tsv output.tsv")
vcf, marker_file, haplotype_file, output = sys.argv[1:]
selected=set()
with open(marker_file) as h:
    for r in csv.DictReader(h,delimiter='\t'):
        d=float(r['delta_AF'])
        if abs(d)>=0.30:
            selected.add((r['CHROM'],int(r['POS']),r['REF'],r['ALT']))
meta={}
with open(haplotype_file) as h:
    for r in csv.DictReader(h,delimiter='\t'): meta[(r['sample_id'],int(r['haplotype']))]=r

samples=[]; markers=[]; hap=[]
with gzip.open(vcf,'rt') as h:
    for line in h:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip().split('\t')[9:]
            hap=[[[],[]] for _ in samples]
            continue
        f=line.rstrip().split('\t'); key=(f[0],int(f[1]),f[3],f[4])
        if key not in selected: continue
        gi=f[8].split(':').index('GT'); markers.append(f"{f[0]}:{f[1]}:{f[3]}:{f[4]}")
        for i,cell in enumerate(f[9:]):
            gt=cell.split(':')[gi]
            a,b=gt.split('|') if '|' in gt else ('.','.')
            hap[i][0].append(a); hap[i][1].append(b)
with open(output,'w',newline='') as h:
    fields=['tip_id','sample_id','haplotype','cohort','gonadal_sex','haplotype_class']+markers
    w=csv.writer(h,delimiter='\t'); w.writerow(fields)
    for i,s in enumerate(samples):
        for k in (0,1):
            mm=meta[(s,k+1)]
            cls=mm['haplotype_class']
            w.writerow([f"{s}_H{k+1}",s,k+1,mm['cohort'],mm['gonadal_sex'],cls]+hap[i][k])
print(f"Wrote {2*len(samples)} haplotypes across {len(markers)} markers")
