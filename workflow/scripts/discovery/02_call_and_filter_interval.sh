#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 5 ]]; then echo "Usage: $0 REF REGION OUT_PREFIX BAM..." >&2; exit 64; fi
ref=$1; region=$2; prefix=$3; shift 3
raw="${prefix}.raw.vcf.gz"; out="${prefix}.filtered.snps.vcf.gz"
bcftools mpileup --threads "${BCF_THREADS:-8}" -f "$ref" -r "$region" -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD -Ou "$@" |
  bcftools call --threads "${BCF_THREADS:-8}" -m -v -Ou |
  bcftools norm --threads "${BCF_THREADS:-8}" -f "$ref" -m -any -Oz -o "$raw"
bcftools index -f -t "$raw"
bcftools view --threads "${BCF_THREADS:-8}" -m2 -M2 -v snps -Ou "$raw" |
  bcftools filter --threads "${BCF_THREADS:-8}" -e 'QUAL<30 || INFO/DP<300 || INFO/DP>3000' -Ou |
  bcftools +fill-tags --threads "${BCF_THREADS:-8}" -Ou -- -t MAF,F_MISSING |
  bcftools view --threads "${BCF_THREADS:-8}" -i 'F_MISSING<=0.10 && MAF>=0.05' -Oz -o "$out"
bcftools index -f "$out"
