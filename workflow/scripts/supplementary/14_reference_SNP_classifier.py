#!/usr/bin/env python3
"""Project the unsexed reference juvenile onto discovery-sex SNP centroids."""
from __future__ import annotations
import csv, gzip, math, re, sys

if len(sys.argv)!=7:
    raise SystemExit('Usage: 14_reference_SNP_classifier.py discovery.vcf.gz sex.tsv panel.csv ref.mpileup per_marker.tsv summary.tsv')
vcf, sexfile, panelfile, pilefile, marker_out, summary_out=sys.argv[1:]
sex={}
with open(sexfile) as h:
    for r in csv.DictReader(h,delimiter='\t'): sex[r['Sample']]=r['Group']
panel={}
with open(panelfile) as h:
    for r in csv.DictReader(h): panel[(r['CHROM'],int(r['POS']))]=(r['REF'].upper(),r['ALT'].upper())

samples=[]; means={}; sample_dos={}
with gzip.open(vcf,'rt') as h:
    for line in h:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip().split('\t')[9:]; continue
        f=line.rstrip().split('\t'); key=(f[0],int(f[1]))
        if key not in panel or (f[3].upper(),f[4].upper())!=panel[key]: continue
        gi=f[8].split(':').index('GT'); vals={'F':[],'M':[]}; bysample={}
        for s,cell in zip(samples,f[9:]):
            gt=cell.split(':')[gi].replace('|','/')
            if '.' in gt: continue
            ds=sum(int(a) for a in gt.split('/')); vals[sex[s]].append(ds); bysample[s]=ds
        means[key]={'F':sum(vals['F'])/len(vals['F']),'M':sum(vals['M'])/len(vals['M'])}
        sample_dos[key]=bysample

def base_counts(ref,bases):
    c={b:0 for b in 'ACGT'}; i=0
    while i<len(bases):
        x=bases[i]
        if x=='^': i+=2; continue
        if x=='$': i+=1; continue
        if x in '+-':
            j=i+1
            while j<len(bases) and bases[j].isdigit(): j+=1
            n=int(bases[i+1:j]); i=j+n; continue
        if x in '.,': c[ref]+=1
        elif x.upper() in c: c[x.upper()]+=1
        i+=1
    return c

refdos={}; depth={}
with open(pilefile) as h:
    for line in h:
        f=line.rstrip().split('\t'); key=(f[0],int(f[1]))
        if key not in panel: continue
        ref,alt=panel[key]; counts=base_counts(ref,f[4]); dp=sum(counts.values())
        if dp:
            refdos[key]=2*counts.get(alt,0)/dp; depth[key]=dp

rows=[]; sq={'F':0.0,'M':0.0}; used=0
for key,(ref,alt) in panel.items():
    if key not in means or key not in refdos: continue
    r={'CHROM':key[0],'POS':key[1],'REF':ref,'ALT':alt,'depth':depth[key],
       'reference_alt_dosage':refdos[key],'female_mean_dosage':means[key]['F'],
       'male_mean_dosage':means[key]['M']}
    r['closer_sex']='M' if abs(r['reference_alt_dosage']-r['male_mean_dosage']) < abs(r['reference_alt_dosage']-r['female_mean_dosage']) else 'F'
    for s in ('F','M'): sq[s]+=(r['reference_alt_dosage']-means[key][s])**2
    rows.append(r); used+=1

# Leave-one-out discovery benchmark using the frozen panel and nearest sex centroid.
correct=total=0
for sample in samples:
    if sample not in sex: continue
    dist={'F':0.0,'M':0.0}; n=0
    for key in panel:
        if key not in sample_dos or sample not in sample_dos[key]: continue
        x=sample_dos[key][sample]
        for group in ('F','M'):
            z=[v for s,v in sample_dos[key].items() if sex[s]==group and s!=sample]
            if z: dist[group]+=(x-sum(z)/len(z))**2
        n+=1
    if n:
        pred=min(dist,key=dist.get); correct+=pred==sex[sample]; total+=1

with open(marker_out,'w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)
pred=min(sq,key=sq.get) if used else 'NA'
with open(summary_out,'w',newline='') as h:
    w=csv.writer(h,delimiter='\t'); w.writerow(['metric','value'])
    w.writerow(['n_panel_markers',len(panel)]); w.writerow(['n_reference_markers_scored',used])
    w.writerow(['female_centroid_squared_distance',sq['F']]); w.writerow(['male_centroid_squared_distance',sq['M']])
    w.writerow(['reference_predicted_profile',pred])
    w.writerow(['discovery_leave_one_out_correct',correct]); w.writerow(['discovery_leave_one_out_total',total])
    w.writerow(['discovery_leave_one_out_accuracy',correct/total if total else 'NA'])
