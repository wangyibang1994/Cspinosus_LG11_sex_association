#!/usr/bin/env python3
"""Summarise homologous-interval size, TE coverage, k-mer conservation and duplications."""
from __future__ import annotations
import csv, glob, os, sys
from collections import defaultdict

if len(sys.argv)!=3: raise SystemExit('Usage: 15_summarise_comparative.py analysis_dir output_prefix')
root,prefix=sys.argv[1:]

def fasta_len(path):
    return sum(len(x.strip()) for x in open(path) if not x.startswith('>'))
def merged_bp(iv):
    if not iv: return 0
    iv=sorted(iv); total=0; a,b=iv[0]
    for x,y in iv[1:]:
        if x<=b: b=max(b,y)
        else: total+=b-a+1; a,b=x,y
    return total+b-a+1

species=('Csp','Cluc','Ccuc'); rows=[]
for sp in species:
    fa=f'{root}/intervals/{sp}.orthologous.fa'; length=fasta_len(fa)
    rm=glob.glob(f'{root}/repeats/{sp}/*.out')
    iv=defaultdict(list); n=defaultdict(int)
    if rm:
        for line in open(rm[0],errors='replace'):
            f=line.split()
            if len(f)<11 or not f[0].isdigit(): continue
            try: a,b=int(f[5]),int(f[6])
            except ValueError: continue
            cls=f[10].split('/')[0]
            iv['All_repeats'].append((a,b)); iv[cls].append((a,b)); n['All_repeats']+=1; n[cls]+=1
    for cls in sorted(iv):
        bp=merged_bp(iv[cls]); rows.append([sp,length,cls,n[cls],bp,bp/length])
with open(prefix+'.repeat_summary.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['species','interval_length','repeat_class','n_elements','covered_bp','covered_fraction']); w.writerows(rows)

krows=[]
for sp in species:
    path=f'{root}/kmers/{sp}.core83.exact.tsv'; hits=defaultdict(int)
    if os.path.exists(path):
        for line in open(path):
            f=line.split('\t');
            if f: hits[f[0]]+=1
    krows.append([sp,len(hits),sum(hits.values())])
with open(prefix+'.kmer_conservation.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['species','core83_with_exact_hit','total_exact_alignments']); w.writerows(krows)

drows=[]
for sp in species:
    path=f'{root}/self_align/{sp}.self_blast.tsv'; seen=set(); vals=[]
    if os.path.exists(path):
        for line in open(path):
            f=line.split('\t');
            if len(f)<8: continue
            pident=float(f[2]); ln=int(f[3]); q=sorted((int(f[4]),int(f[5]))); s=sorted((int(f[6]),int(f[7])))
            if ln<1000 or pident<90 or q==s: continue
            key=tuple(sorted((tuple(q),tuple(s))))
            if key in seen: continue
            seen.add(key); vals.append((pident,ln,*q,*s))
    drows.append([sp,len(vals),sum(v[1] for v in vals),max((v[1] for v in vals),default=0)])
with open(prefix+'.segmental_duplication_summary.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['species','n_nonredundant_alignments','sum_alignment_bp','max_alignment_bp']); w.writerows(drows)

grows=[]
for sp in species:
    path=f'{root}/proteins/{sp}.priority6.miniprot.gff'; hits=set()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if len(f)<9 or f[2] not in {'mRNA','transcript'}: continue
            attrs=f[8]
            hits.add(attrs)
    grows.append([sp,len(hits)])
with open(prefix+'.priority_gene_mapping_summary.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['species','n_priority_transcript_models']); w.writerows(grows)
