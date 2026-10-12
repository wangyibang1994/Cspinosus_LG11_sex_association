#!/usr/bin/env python3
"""Region-based GFF3 annotation for pre-ranked high-delta-AF SNPs.

This reports genomic context, gene overlap and nearest-gene distance.  It does
not infer synonymous/non-synonymous consequences; that requires a transcript
effect predictor with the reference sequence and CDS phases.
"""
import csv, sys
from bisect import bisect_left

snps, gff, out, genes_out = sys.argv[1:5]

def attrs(s):
    return {x.split("=",1)[0]: x.split("=",1)[1] for x in s.split(";") if "=" in x}

genes, exons, cds = [], [], []
with open(gff) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"): continue
        c, source, typ, start, end, score, strand, phase, at = line.rstrip("\n").split("\t")
        a=attrs(at); rec=(int(start),int(end),strand,a.get("ID",a.get("Parent","")),a.get("Parent",""))
        if typ=="gene": genes.append(rec)
        elif typ=="exon": exons.append(rec)
        elif typ=="CDS": cds.append(rec)
genes.sort(); starts=[g[0] for g in genes]

def overlaps(items,p): return [x for x in items if x[0]<=p<=x[1]]
def nearest(p):
    i=bisect_left(starts,p); cand=genes[max(0,i-1):min(len(genes),i+1)]
    best=min(cand,key=lambda x: 0 if x[0]<=p<=x[1] else min(abs(p-x[0]),abs(p-x[1])))
    d=0 if best[0]<=p<=best[1] else min(abs(p-best[0]),abs(p-best[1]))
    return best,d

all_rows=[]
with open(snps) as fh:
    for r in csv.DictReader(fh):
        p=int(r["POS"]); gg,d=nearest(p); gene_hit=overlaps(genes,p); exon_hit=overlaps(exons,p); cds_hit=overlaps(cds,p)
        if cds_hit: context="CDS"
        elif exon_hit: context="exon_nonCDS_or_UTR"
        elif gene_hit: context="intron"
        elif d<=10000:
            # direction relative to the nearest gene is reported only descriptively.
            context="near_gene_10kb"
        else: context="intergenic"
        r.update({"Functional_context":context,"Overlapping_gene":";".join(x[3] for x in gene_hit) or "", "Nearest_gene":gg[3],"Nearest_gene_distance_bp":d,"Nearest_gene_strand":gg[2],"Nearest_gene_start":gg[0],"Nearest_gene_end":gg[1]})
        all_rows.append(r)
fields=list(all_rows[0].keys())
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(all_rows)

summary={}
for r in all_rows: summary[r["Functional_context"]]=summary.get(r["Functional_context"],0)+1
with open(genes_out,"w",newline="") as fh:
    w=csv.writer(fh);w.writerow(["Functional_context","N_SNP"])
    for k,v in sorted(summary.items()):w.writerow([k,v])
