#!/usr/bin/env bash
#BSUB -J lg11_phase
#BSUB -q normal
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=3000]"
#BSUB -M 24000
#BSUB -W 12:00
#BSUB -o 10_logs/lg11_phase.%J.out
#BSUB -e 10_logs/lg11_phase.%J.err
set -euo pipefail

sexproj=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
popproj=${CSPINOSUS_EXTERNAL_PROJECT:?Set CSPINOSUS_EXTERNAL_PROJECT to the external-project directory}
env=${HAPLOTYPE_BIN_DIR:?Set HAPLOTYPE_BIN_DIR to a directory containing bcftools and beagle}
out="$sexproj/12_supplementary_analyses/02_LG11_haplotypes"
disc="$sexproj/06_filtering/LG11_4.8_5.5Mb.filtered.snps.vcf.gz"
ext="$popproj/20_sex_association_validation/LG11_equivalent_4.8_5.5Mb.vcf.gz"
mkdir -p "$out/isec_shared" "$sexproj/10_logs"

printf 'CM056831.1\tLG11\n' > "$out/external.rename_chrs.tsv"

# Normalise, retain biallelic SNPs, assign allele-specific stable IDs and remove
# duplicate records before phasing.
"$env/bcftools" norm -f "$sexproj/02_reference/genome.fa" -m-any "$disc" -Ou | \
  "$env/bcftools" view -m2 -M2 -v snps -Ou | \
  "$env/bcftools" annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou | \
  "$env/bcftools" norm -d exact -Oz -o "$out/discovery.normalised.vcf.gz"
"$env/bcftools" index -f "$out/discovery.normalised.vcf.gz"

"$env/bcftools" annotate --rename-chrs "$out/external.rename_chrs.tsv" "$ext" -Ou | \
  "$env/bcftools" norm -f "$sexproj/02_reference/genome.fa" -m-any -Ou | \
  "$env/bcftools" view -m2 -M2 -v snps -Ou | \
  "$env/bcftools" annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou | \
  "$env/bcftools" norm -d exact -Oz -o "$out/external.normalised.vcf.gz"
"$env/bcftools" index -f "$out/external.normalised.vcf.gz"

# Phase each cohort separately for cohort-specific architecture.
"$env/beagle" gt="$out/discovery.normalised.vcf.gz" out="$out/discovery.phased" \
  chrom=LG11:4800000-5500000 nthreads=8 seed=20260922
"$env/beagle" gt="$out/external.normalised.vcf.gz" out="$out/external.phased" \
  chrom=LG11:4800000-5500000 nthreads=8 seed=20260922

# Build an exact shared-allele set. The two cohorts are then jointly phased
# without sex labels solely to put replicated haplotypes in a common space.
"$env/bcftools" isec -p "$out/isec_shared" \
  "$out/discovery.normalised.vcf.gz" "$out/external.normalised.vcf.gz"
"$env/bcftools" view "$out/isec_shared/0002.vcf" -Oz -o "$out/discovery.shared.vcf.gz"
"$env/bcftools" view "$out/isec_shared/0003.vcf" -Oz -o "$out/external.shared.vcf.gz"
"$env/bcftools" index -f "$out/discovery.shared.vcf.gz"
"$env/bcftools" index -f "$out/external.shared.vcf.gz"
"$env/bcftools" merge -m none -Oz -o "$out/combined120.shared.unphased.vcf.gz" \
  "$out/discovery.shared.vcf.gz" "$out/external.shared.vcf.gz"
"$env/bcftools" index -f "$out/combined120.shared.unphased.vcf.gz"
"$env/beagle" gt="$out/combined120.shared.unphased.vcf.gz" \
  out="$out/combined120.shared.phased" chrom=LG11:4800000-5500000 \
  nthreads=8 seed=20260922

{
  "$env/beagle" 2>&1 | head -2
  "$env/bcftools" --version | head -2
  printf 'discovery_variants\t'; "$env/bcftools" view -H "$out/discovery.normalised.vcf.gz" | wc -l
  printf 'external_variants\t'; "$env/bcftools" view -H "$out/external.normalised.vcf.gz" | wc -l
  printf 'exact_shared_variants\t'; "$env/bcftools" view -H "$out/discovery.shared.vcf.gz" | wc -l
  printf 'beagle_seed\t20260922\n'
} > "$out/software_and_variant_counts.txt"
