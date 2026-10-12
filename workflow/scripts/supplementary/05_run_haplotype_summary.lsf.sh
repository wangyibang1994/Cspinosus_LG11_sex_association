#!/usr/bin/env bash
#BSUB -J lg11_hapsum
#BSUB -q normal
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000
#BSUB -W 02:00
#BSUB -o 10_logs/lg11_hapsum.%J.out
#BSUB -e 10_logs/lg11_hapsum.%J.err
set -euo pipefail
project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/02_LG11_haplotypes"
python3 "$project/09_scripts/04_haplotype_sex_association.py" \
  "$out/combined120.shared.phased.vcf.gz" \
  "$project/00_metadata/sample_sex.corrected.tsv" \
  "$out/external60_final_sex.tsv" \
  "$out/combined120_shared"
