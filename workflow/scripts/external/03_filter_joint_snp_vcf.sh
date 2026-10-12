#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 3 ]]; then echo "Usage: $0 REF RAW_COHORT_VCF OUT_PREFIX" >&2; exit 64; fi
ref=$1; raw=$2; p=$3
gatk VariantFiltration -R "$ref" -V "$raw" -O "${p}.hardfiltered.vcf.gz" \
  --filter-name QD2 --filter-expression 'QD < 2.0' \
  --filter-name FS60 --filter-expression 'FS > 60.0' \
  --filter-name SOR3 --filter-expression 'SOR > 3.0' \
  --filter-name MQ40 --filter-expression 'MQ < 40.0' \
  --filter-name MQRankSum12_5 --filter-expression 'MQRankSum < -12.5' \
  --filter-name ReadPosRankSum8 --filter-expression 'ReadPosRankSum < -8.0'
gatk SelectVariants -R "$ref" -V "${p}.hardfiltered.vcf.gz" --exclude-filtered true \
  --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O "${p}.pass.biallelic.vcf.gz"
gatk VariantFiltration -R "$ref" -V "${p}.pass.biallelic.vcf.gz" -O "${p}.gtfiltered.vcf.gz" \
  --genotype-filter-name lowDP --genotype-filter-expression 'DP < 4' \
  --genotype-filter-name lowGQ --genotype-filter-expression 'GQ < 20'
gatk SelectVariants -R "$ref" -V "${p}.gtfiltered.vcf.gz" --set-filtered-gt-to-nocall \
  --max-nocall-fraction 0.1 --exclude-non-variants --select 'AF >= 0.01 && AF <= 0.99' \
  -O "${p}.final.vcf.gz"
