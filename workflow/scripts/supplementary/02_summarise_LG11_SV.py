#!/usr/bin/env python3
"""Summarise sex-stratified DELLY genotypes and read evidence.

Input is the wide table emitted by 01_delly_LG11_joint.lsf.sh. Each record has
six fixed fields followed by seven fields per sample: SAMPLE, GT, FT, DV, DR, RV, RR.
"""
from __future__ import annotations
import csv, math, sys
from collections import defaultdict

if len(sys.argv) != 4:
    raise SystemExit("Usage: 02_summarise_LG11_SV.py wide.tsv sample_sex.tsv output.tsv")
wide, sex_file, output = sys.argv[1:]
sex = {}
with open(sex_file) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        sex[row["Sample"]] = row["Group"]

def present(gt: str) -> bool:
    if gt in {"./.", ".|.", "."}: return False
    return any(a not in {"0", "."} for a in gt.replace("|", "/").split("/"))

def called(gt: str) -> bool:
    return gt not in {"./.", ".|.", "."} and "." not in gt.replace("|", "/").split("/")

def safe_int(x: str) -> int:
    try: return int(x)
    except (TypeError, ValueError): return 0

def logchoose(n: int, k: int) -> float:
    if k < 0 or k > n: return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    r1, r2, col1, n = a+b, c+d, a+c, a+b+c+d
    lo, hi = max(0, col1-r2), min(r1, col1)
    def prob(x):
        return math.exp(logchoose(col1,x)+logchoose(n-col1,r1-x)-logchoose(n,r1))
    observed = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo,hi+1) if prob(x) <= observed + 1e-12))

rows=[]
with open(wide) as handle:
    for line in handle:
        f=line.rstrip("\n").split("\t")
        if len(f) < 13 or (len(f)-6) % 7: continue
        chrom,pos,vid,svtype,end,qual=f[:6]
        calls=[]
        for i in range(6,len(f),7):
            sample,gt,ft,dv,dr,rv,rr=f[i:i+7]
            calls.append((sample,gt,ft,safe_int(dv),safe_int(dr),safe_int(rv),safe_int(rr)))
        counts=defaultdict(lambda: {"carrier":0,"n":0,"pe":0,"sr":0})
        for sample,gt,ft,dv,dr,rv,rr in calls:
            group=sex.get(sample)
            if group not in {"F","M"}: continue
            if ft != "PASS": continue
            if not called(gt): continue
            counts[group]["n"] += 1
            counts[group]["carrier"] += present(gt)
            counts[group]["pe"] += dv
            counts[group]["sr"] += rv
        a=counts["M"]["carrier"]; b=counts["M"]["n"]-a
        c=counts["F"]["carrier"]; d=counts["F"]["n"]-c
        start=int(pos); stop=int(end) if end not in {".",""} else start
        overlap_candidate = chrom=="LG11" and stop>=4_800_000 and start<=5_500_000
        overlap_kmer_cluster = chrom=="LG11" and stop>=5_168_000 and start<=5_336_000
        rows.append({
            "CHROM":chrom,"POS":start,"END":stop,"ID":vid,"SVTYPE":svtype,"QUAL":qual,
            "male_carriers":a,"male_total":counts["M"]["n"],
            "female_carriers":c,"female_total":counts["F"]["n"],
            "delta_carrier_frequency":a/max(counts["M"]["n"],1)-c/max(counts["F"]["n"],1),
            "fisher_p":fisher_two_sided(a,b,c,d),
            "male_variant_PE":counts["M"]["pe"],"female_variant_PE":counts["F"]["pe"],
            "male_variant_SR":counts["M"]["sr"],"female_variant_SR":counts["F"]["sr"],
            "overlaps_4.80_5.50Mb":overlap_candidate,
            "overlaps_5.168_5.336Mb_kmer_cluster":overlap_kmer_cluster,
        })
# Benjamini-Hochberg adjustment across all tested SVs.
order=sorted(range(len(rows)), key=lambda i: rows[i]["fisher_p"])
qvals=[1.0]*len(rows); running=1.0; m=len(rows)
for rank0 in range(m-1,-1,-1):
    i=order[rank0]; rank=rank0+1
    running=min(running, rows[i]["fisher_p"]*m/rank)
    qvals[i]=min(1.0,running)
for i,q in enumerate(qvals): rows[i]["fisher_q_BH"]=q
rows.sort(key=lambda r:(r["fisher_p"],-abs(r["delta_carrier_frequency"]),r["POS"]))
fields=list(rows[0]) if rows else ["CHROM","POS","END","ID","SVTYPE"]
with open(output,"w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t")
    writer.writeheader(); writer.writerows(rows)
print(f"Wrote {len(rows)} SVs to {output}")
