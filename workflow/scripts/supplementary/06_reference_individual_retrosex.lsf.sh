#!/usr/bin/env bash
#BSUB -J ref_retrosex
#BSUB -q normal
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=3000]"
#BSUB -M 48000
#BSUB -W 72:00
#BSUB -o 10_logs/ref_retrosex.%J.out
#BSUB -e 10_logs/ref_retrosex.%J.err
set -euo pipefail

project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/03_reference_individual"
raw=${REFERENCE_JUVENILE_FASTQ_DIR:?Set REFERENCE_JUVENILE_FASTQ_DIR to the reference-juvenile FASTQ directory}
kmcenv=${KMC_BIN_DIR:?Set KMC_BIN_DIR to a directory containing kmc and kmc_tools}
ref="$project/02_reference/genome.fa"
samtools=${SAMTOOLS:-samtools}
bcftools=${BCFTOOLS:-bcftools}
bwa=${BWA:-bwa}
mkdir -p "$out/tmp/kmc" "$project/10_logs"

# The original A/B/C Illumina WGS libraries from the 2023 reference juvenile
# already reside on the server; use them directly rather than downloading GSA copies.
find "$raw" -maxdepth 1 -type f -name 'C_spinosus_[ABC]_illumina_R[12].fq.gz' | sort \
  > "$out/reference_WGS_fastq.list"
[[ $(wc -l < "$out/reference_WGS_fastq.list") -eq 6 ]]

# Alignment-free test against the panel frozen in the discovery cohort.
awk -F, 'NR>1{gsub(/"/,"",$1); gsub(/"/,"",$2); print ">"$1"\n"$2}' \
  "$out/frozen_core83_kmers.csv" > "$out/core83.fa"
"$kmcenv/kmc" -k31 -t8 -m24 -ci1 -cs1000000 -fa "$out/core83.fa" \
  "$out/core83" "$out/tmp/kmc"
"$kmcenv/kmc" -k31 -t16 -m40 -ci1 -cs1000000 -fq "@$out/reference_WGS_fastq.list" \
  "$out/reference_WGS_k31" "$out/tmp/kmc"
"$kmcenv/kmc_tools" -t8 simple "$out/core83" "$out/reference_WGS_k31" \
  intersect "$out/reference_core83_hits" -ocright
"$kmcenv/kmc_tools" transform "$out/reference_core83_hits" dump \
  "$out/reference_core83_hits.tsv"
awk 'BEGIN{OFS="\t"; print "metric","value"}
  $2>=3{n++} $1=="AACATACACGTATATGCTGATGCAATGCACG"{k202=$2}
  END{print "core83_markers_count_ge3",n+0; print "K202_read_count",k202+0}' \
  "$out/reference_core83_hits.tsv" > "$out/reference_kmer_profile.tsv"

# Targeted LG11 alignment. Reads are aligned to the full genome, but only LG11
# alignments are retained, avoiding a large whole-genome BAM.
mapfile -t r1 < <(find "$raw" -maxdepth 1 -name 'C_spinosus_[ABC]_illumina_R1.fq.gz' | sort)
mapfile -t r2 < <(find "$raw" -maxdepth 1 -name 'C_spinosus_[ABC]_illumina_R2.fq.gz' | sort)
"$bwa" mem -t 16 -R '@RG\tID:reference_juvenile\tSM:reference_juvenile\tPL:ILLUMINA' \
  "$ref" <(gzip -dc "${r1[@]}") <(gzip -dc "${r2[@]}") | \
  "$samtools" view -h - | awk 'BEGIN{OFS="\t"} /^@/ || $3=="LG11"' | \
  "$samtools" view -b - | "$samtools" sort -@ 8 -o "$out/reference_juvenile.LG11.bam"
"$samtools" index "$out/reference_juvenile.LG11.bam"

awk -F, 'BEGIN{OFS="\t"} NR>1{gsub(/"/,"",$1); gsub(/"/,"",$2); print $1,$2-1,$2}' \
  "$out/frozen_top50_discovery_snps.csv" > "$out/frozen_top50.bed"
"$samtools" mpileup -q 20 -Q 20 -l "$out/frozen_top50.bed" -f "$ref" \
  "$out/reference_juvenile.LG11.bam" > "$out/reference_juvenile.top50.mpileup"
"$bcftools" mpileup -f "$ref" -R "$out/frozen_top50.bed" -q 20 -Q 20 \
  -a FORMAT/AD,FORMAT/DP -Ou "$out/reference_juvenile.LG11.bam" | \
  "$bcftools" call -m -Oz -o "$out/reference_juvenile.top50.vcf.gz"
"$bcftools" index -f "$out/reference_juvenile.top50.vcf.gz"
"$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n' \
  "$out/reference_juvenile.top50.vcf.gz" > "$out/reference_juvenile.top50_genotypes.tsv"
python3 "$project/09_scripts/14_reference_SNP_classifier.py" \
  "$project/06_filtering/LG11_4.8_5.5Mb.filtered.snps.vcf.gz" \
  "$project/00_metadata/sample_sex.corrected.tsv" \
  "$out/frozen_top50_discovery_snps.csv" "$out/reference_juvenile.top50.mpileup" \
  "$out/reference_juvenile.top50_centroid_markers.tsv" \
  "$out/reference_juvenile.genetic_profile_summary.tsv"

{
  printf 'source_project\tCRA009849\n'
  printf 'WGS_runs\tCRR660781,CRR660782,CRR660783\n'
  printf 'biosamples\tSAMC1092750,SAMC1092751,SAMC1092752\n'
  "$kmcenv/kmc" --version 2>&1 | head -1
  "$bwa" 2>&1 | grep Version | head -1
  "$samtools" --version | head -1
  "$bcftools" --version | head -1
} > "$out/software_and_data_sources.txt"
