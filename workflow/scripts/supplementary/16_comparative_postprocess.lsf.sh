#!/usr/bin/env bash
#BSUB -J lg11_compost
#BSUB -q normal
#BSUB -n 6
#BSUB -R "span[hosts=1] rusage[mem=2000]"
#BSUB -M 12000
#BSUB -W 04:00
#BSUB -o 10_logs/lg11_compost.%J.out
#BSUB -e 10_logs/lg11_compost.%J.err
set -euo pipefail
project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/04_comparative_genomics"
blastn=${BLASTN:-blastn}
mkdir -p "$out/self_align" "$out/summary"
for spec in Csp Cluc Ccuc; do
  "$blastn" -task megablast -query "$out/intervals/${spec}.orthologous.fa" \
    -subject "$out/intervals/${spec}.orthologous.fa" -dust yes -word_size 28 \
    -perc_identity 90 -evalue 1e-20 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
    > "$out/self_align/${spec}.self_blast.tsv"
done
python3 "${SUMMARISE_COMPARATIVE_PY:?Set SUMMARISE_COMPARATIVE_PY}" "$out" "$out/summary/LG11_comparative"
