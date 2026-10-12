#!/usr/bin/env python3
"""Quantify a discovery-defined male-associated haplotype and replication.

No sex labels enter phasing. Discovery sex is used only to orient alleles and
define the high-score haplotype class; external gonadal sex is used only for
the final evaluation.
"""
from __future__ import annotations
import csv, gzip, math, sys

if len(sys.argv) != 5:
    raise SystemExit("Usage: 04_haplotype_sex_association.py phased.vcf.gz discovery_sex.tsv external_sex.tsv out_prefix")
vcf, discovery_file, external_file, prefix = sys.argv[1:]

def read_sex(path, sample_col, sex_col):
    out={}
    with open(path) as h:
        d=csv.DictReader(h, delimiter='\t')
        for r in d: out[r[sample_col]]=r[sex_col][0].upper()
    return out

disc=read_sex(discovery_file,'Sample','Group')
ext=read_sex(external_file,'sample_id','sex')

samples=[]; variants=[]; phased=[]
with gzip.open(vcf,'rt') as h:
    for line in h:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip().split('\t')[9:]; continue
        f=line.rstrip().split('\t'); fmt=f[8].split(':'); gi=fmt.index('GT')
        g=[]
        for cell in f[9:]:
            gt=cell.split(':')[gi]
            if '|' not in gt: g.append((None,None))
            else:
                a,b=gt.split('|'); g.append((None if a=='.' else int(a),None if b=='.' else int(b)))
        variants.append((f[0],int(f[1]),f[3],f[4])); phased.append(g)

# Discovery allele-frequency contrast; retain strongly differentiated markers.
selected=[]; directions=[]; marker_stats=[]
for j,var in enumerate(variants):
    af={s:[] for s in ('F','M')}
    for i,sample in enumerate(samples):
        if sample not in disc: continue
        af[disc[sample]].extend(a for a in phased[j][i] if a is not None)
    if not af['F'] or not af['M']: continue
    pf=sum(af['F'])/len(af['F']); pm=sum(af['M'])/len(af['M']); delta=pm-pf
    marker_stats.append((*var,pf,pm,delta))
    if abs(delta)>=0.30:
        selected.append(j); directions.append(1 if delta>0 else 0)
if len(selected)<5: raise SystemExit(f"Only {len(selected)} markers have |delta AF| >= 0.30")

hap_rows=[]; scores=[]
for i,sample in enumerate(samples):
    for hap in (0,1):
        match=[]
        for j,direction in zip(selected,directions):
            allele=phased[j][i][hap]
            if allele is not None: match.append(allele==direction)
        score=sum(match)/len(match) if match else float('nan')
        scores.append(score)
        hap_rows.append({'sample_id':sample,'cohort':'discovery' if sample in disc else 'external',
                         'gonadal_sex':disc.get(sample,ext.get(sample,'')),'haplotype':hap+1,
                         'male_allele_score':score,'n_informative_markers':len(match)})

# Deterministic one-dimensional k-means (K=2) on haplotype scores.
vals=[x for x in scores if not math.isnan(x)]; centres=[min(vals),max(vals)]
for _ in range(100):
    groups=[[x for x in vals if abs(x-centres[k]) <= abs(x-centres[1-k])] for k in (0,1)]
    new=[sum(g)/len(g) if g else centres[k] for k,g in enumerate(groups)]
    if max(abs(new[k]-centres[k]) for k in (0,1))<1e-12: break
    centres=new
male_associated_cluster=max(range(2),key=lambda k: centres[k])
for r in hap_rows:
    x=r['male_allele_score']; k=min(range(2),key=lambda z:abs(x-centres[z]))
    # Statistical phasing does not establish chromosome identity or copy number.
    # Use evidence-bounded labels rather than implying a resolved X/Y system.
    r['haplotype_class']='male-associated' if k==male_associated_cluster else 'other'

with open(prefix+'.haplotypes.tsv','w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=hap_rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(hap_rows)
with open(prefix+'.marker_stats.tsv','w',newline='') as h:
    fields=['CHROM','POS','REF','ALT','AF_female_discovery','AF_male_discovery','delta_AF']
    w=csv.writer(h,delimiter='\t'); w.writerow(fields); w.writerows(marker_stats)

per_sample=[]
for sample in samples:
    rr=[r for r in hap_rows if r['sample_id']==sample]
    per_sample.append({'sample_id':sample,'cohort':rr[0]['cohort'],'gonadal_sex':rr[0]['gonadal_sex'],
                       'max_male_allele_score':max(r['male_allele_score'] for r in rr),
                       'n_male_associated_haplotypes':sum(r['haplotype_class']=='male-associated' for r in rr),
                       'male_associated_carrier':sum(r['haplotype_class']=='male-associated' for r in rr)>0})
with open(prefix+'.individuals.tsv','w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=per_sample[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(per_sample)

with open(prefix+'.summary.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['cohort','sex','n','male_associated_carriers','carrier_frequency'])
    for cohort in ('discovery','external'):
        for sx in ('F','M'):
            z=[r for r in per_sample if r['cohort']==cohort and r['gonadal_sex']==sx]
            ncar=sum(r['male_associated_carrier'] for r in z)
            w.writerow([cohort,sx,len(z),ncar,ncar/len(z) if z else 'NA'])
print(f"Informative markers: {len(selected)}; haplotype centres: {centres}")
