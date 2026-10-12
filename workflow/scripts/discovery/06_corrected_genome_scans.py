#!/usr/bin/env python3
import argparse, gzip, math, os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("vcf_dir", help="Directory containing LG1.filtered.snps.vcf.gz ... LG24.filtered.snps.vcf.gz")
parser.add_argument("metadata", help="TSV with columns record/sample/group or equivalent three-column header")
parser.add_argument("out_dir")
parser.add_argument("--window", type=int, default=200000)
parser.add_argument("--step", type=int, default=50000)
args = parser.parse_args()
vcf_dir = args.vcf_dir
meta_file = args.metadata
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

sex = {}
with open(meta_file) as fh:
    next(fh)
    for line in fh:
        _, sample, group = line.rstrip().split("\t")
        sex[sample] = group

site_out = open(f"{out_dir}/corrected_site_metrics.tsv", "w")
site_out.write("CHROM\tPOS\tAF_F\tAF_M\tDELTA_AF\tED\tFST_NUM\tFST_DEN\tFST_HUDSON\n")
sites = defaultdict(list)

def gt_dosage(gt):
    gt = gt.replace("|", "/")
    if gt == "0/0": return 0
    if gt in ("0/1", "1/0"): return 1
    if gt == "1/1": return 2
    return None

for chrom_i in range(1,25):
    chrom = f"LG{chrom_i}"
    path = f"{vcf_dir}/{chrom}.filtered.snps.vcf.gz"
    if not os.path.exists(path):
        raise SystemExit(f"Missing {path}")
    samples = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                fields=line.rstrip().split("\t"); samples=fields[9:]
                idxF=[i for i,s in enumerate(samples) if sex[s]=="F"]
                idxM=[i for i,s in enumerate(samples) if sex[s]=="M"]
                continue
            f=line.rstrip().split("\t"); pos=int(f[1]); fmt=f[8].split(":"); gi=fmt.index("GT")
            dos=[]
            for q in f[9:]:
                a=q.split(":"); dos.append(gt_dosage(a[gi] if gi<len(a) else "."))
            fvals=[dos[i] for i in idxF if dos[i] is not None]; mvals=[dos[i] for i in idxM if dos[i] is not None]
            if len(fvals)<2 or len(mvals)<2: continue
            pF=sum(fvals)/(2*len(fvals)); pM=sum(mvals)/(2*len(mvals)); delta=abs(pF-pM); ed=math.sqrt(2)*delta
            nF=2*len(fvals); nM=2*len(mvals)
            num=(pF-pM)**2 - pF*(1-pF)/(nF-1) - pM*(1-pM)/(nM-1)
            den=pF*(1-pM)+pM*(1-pF)
            fst=num/den if den>0 else float("nan")
            sites[chrom].append((pos,num,den,delta,ed))
            site_out.write(f"{chrom}\t{pos}\t{pF:.8g}\t{pM:.8g}\t{delta:.8g}\t{ed:.8g}\t{num:.8g}\t{den:.8g}\t{fst:.8g}\n")
site_out.close()

win=args.window; step=args.step
with open(f"{out_dir}/corrected_window_metrics.tsv","w") as out:
    out.write("CHROM\tSTART\tEND\tN_SNP\tWEIGHTED_FST\tMEAN_ED\tMAX_DELTA_AF\tN_DELTA_AF_GE_0.3\n")
    for ci in range(1,25):
        chrom=f"LG{ci}"; a=sites[chrom]
        if not a: continue
        maxpos=a[-1][0]; left=0; right=0
        for start in range(0,maxpos+1,step):
            end=start+win
            while left<len(a) and a[left][0]<start: left+=1
            if right<left:right=left
            while right<len(a) and a[right][0]<end:right+=1
            q=a[left:right]
            if not q: continue
            sn=sum(x[1] for x in q); sd=sum(x[2] for x in q); fst=sn/sd if sd>0 else float("nan")
            meaned=sum(x[4] for x in q)/len(q); maxd=max(x[3] for x in q); nd=sum(x[3]>=.3 for x in q)
            out.write(f"{chrom}\t{start}\t{end}\t{len(q)}\t{fst:.8g}\t{meaned:.8g}\t{maxd:.8g}\t{nd}\n")

with open(f"{out_dir}/README.txt","w") as out:
    out.write("Corrected-label sensitivity scans. Hudson-type per-site FST components were summed within 200-kb windows stepped by 50 kb. ED=sqrt((pF_ref-pM_ref)^2+(pF_alt-pM_alt)^2)=sqrt(2)*|deltaAF|.\n")
