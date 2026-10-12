#!/usr/bin/env bash
#BSUB -J lg11_ldtree
#BSUB -q normal
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=4000]"
#BSUB -M 16000
#BSUB -W 08:00
#BSUB -o 10_logs/lg11_ldtree.%J.out
#BSUB -e 10_logs/lg11_ldtree.%J.err
set -euo pipefail
project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/02_LG11_haplotypes"
env=${HAPLOTYPE_BIN_DIR:?Set HAPLOTYPE_BIN_DIR to a directory containing Rscript and plink}
vcf="$out/combined120.shared.phased.vcf.gz"
mkdir -p "$out/ld"

python3 "${EXPORT_HAPLOTYPE_MATRIX_PY:?Set EXPORT_HAPLOTYPE_MATRIX_PY}" "$vcf" \
  "$out/combined120_shared.marker_stats.tsv" "$out/combined120_shared.haplotypes.tsv" \
  "$out/combined120_shared.haplotype_matrix.tsv"
"$env/Rscript" "${HAPLOTYPE_TREE_R:?Set HAPLOTYPE_TREE_R}" \
  "$out/combined120_shared.haplotype_matrix.tsv" "$out/combined120_shared.haplotypes"

"$env/plink" --vcf "$vcf" --double-id --allow-extra-chr --make-bed \
  --out "$out/ld/combined120" >/dev/null
awk 'BEGIN{OFS="\t"} NR>1{print $2,$2}' "$project/00_metadata/sample_sex.corrected.tsv" > "$out/ld/discovery_all.keep"
awk 'BEGIN{OFS="\t"} NR>1 && $3=="F"{print $2,$2}' "$project/00_metadata/sample_sex.corrected.tsv" > "$out/ld/discovery_F.keep"
awk 'BEGIN{OFS="\t"} NR>1 && $3=="M"{print $2,$2}' "$project/00_metadata/sample_sex.corrected.tsv" > "$out/ld/discovery_M.keep"
awk 'BEGIN{OFS="\t"} NR>1 && substr($2,1,1)=="F"{print $1,$1}' "$out/external60_final_sex.tsv" > "$out/ld/external_F.keep"
awk 'BEGIN{OFS="\t"} NR>1 && substr($2,1,1)=="M"{print $1,$1}' "$out/external60_final_sex.tsv" > "$out/ld/external_M.keep"
for group in discovery_F discovery_M external_F external_M; do
  "$env/plink" --bfile "$out/ld/combined120" --allow-extra-chr \
    --keep "$out/ld/${group}.keep" --r2 --ld-window 999999 \
    --ld-window-kb 700 --ld-window-r2 0 --out "$out/ld/$group" >/dev/null
done
python3 "${SUMMARISE_LD_PY:?Set SUMMARISE_LD_PY}" "$out/ld/sex_stratified_LD_summary.tsv" \
  discovery_F "$out/ld/discovery_F.ld" discovery_M "$out/ld/discovery_M.ld" \
  external_F "$out/ld/external_F.ld" external_M "$out/ld/external_M.ld"
