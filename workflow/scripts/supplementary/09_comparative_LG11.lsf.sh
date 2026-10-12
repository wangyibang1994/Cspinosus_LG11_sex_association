#!/usr/bin/env bash
#BSUB -J lg11_compare
#BSUB -q normal
#BSUB -n 12
#BSUB -R "span[hosts=1] rusage[mem=3000]"
#BSUB -M 36000
#BSUB -W 48:00
#BSUB -o 10_logs/lg11_compare.%J.out
#BSUB -e 10_logs/lg11_compare.%J.err
set -euo pipefail

project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/04_comparative_genomics"
samtools=${SAMTOOLS:-samtools}
minimap2=${MINIMAP2:-minimap2}
miniprot=${MINIPROT:-miniprot}
blastn=${BLASTN:-blastn}
rmask=${REPEATMASKER:-RepeatMasker}
mkdir -p "$out"/{queries,synteny,intervals,proteins,repeats,kmers,self_align}

luc=$(find "$out/genomes/GCA_965644535.1" -name '*_genomic.fna' -type f | head -1)
cuc=$(find "$out/genomes/GCA_971798655.1" -name '*_genomic.fna' -type f | head -1)
[[ -s "$luc" && -s "$cuc" ]]
"$samtools" faidx "$luc"; "$samtools" faidx "$cuc"

# C. spinosus query intervals: broad structural interval and complete LG11.
"$samtools" faidx "$project/02_reference/genome.fa" LG11:4800000-5500000 > "$out/queries/Csp_LG11_4.8_5.5Mb.fa"
"$samtools" faidx "$project/02_reference/genome.fa" LG11 > "$out/queries/Csp_LG11.fa"

for spec in Cluc Ccuc; do
  [[ "$spec" == Cluc ]] && target="$luc" || target="$cuc"
  "$minimap2" -x asm20 -c --secondary=no -t 12 "$target" "$out/queries/Csp_LG11.fa" \
    > "$out/synteny/${spec}.whole_LG11.paf"
  "$minimap2" -x asm20 -c --secondary=no -t 12 "$target" "$out/queries/Csp_LG11_4.8_5.5Mb.fa" \
    > "$out/synteny/${spec}.candidate.paf"
python3 "${SUMMARISE_SYNTENY_PY:?Set SUMMARISE_SYNTENY_PY to 08_summarise_synteny.py}" \
    "$out/synteny/${spec}.candidate.paf" "$spec" "$out/synteny/${spec}"
  region=$(cat "$out/synteny/${spec}.region.txt")
  "$samtools" faidx "$target" "$region" > "$out/intervals/${spec}.orthologous.fa"
done
cp "$out/queries/Csp_LG11_4.8_5.5Mb.fa" "$out/intervals/Csp.orthologous.fa"

# Map priority proteins to each orthologous interval to assess local gene order.
awk '/^>/{keep=($0 ~ /gene=Csp11G001(190|280|320|370|380|400)([^0-9]|$)/)} keep' \
  "$out/proteins/Cspinosus.final_gene.pep" > "$out/proteins/Csp.priority6.pep"
for spec in Csp Cluc Ccuc; do
  "$miniprot" --gff -t 12 "$out/intervals/${spec}.orthologous.fa" \
    "$out/proteins/Csp.priority6.pep" > "$out/proteins/${spec}.priority6.miniprot.gff"
done

# Use the same RepeatMasker library and parameters for all three intervals.
for spec in Csp Cluc Ccuc; do
  mkdir -p "$out/repeats/$spec"
  "$rmask" -species teleostei -pa 12 -gff -dir "$out/repeats/$spec" \
    "$out/intervals/${spec}.orthologous.fa"
done

# Test whether frozen male-associated 31-mers are conserved exactly.
awk -F, 'NR>1{gsub(/"/,"",$1); gsub(/"/,"",$2); print ">"$1"\n"$2}' \
  "$project/12_supplementary_analyses/03_reference_individual/frozen_core83_kmers.csv" \
  > "$out/kmers/core83.fa"
for spec in Csp Cluc Ccuc; do
  "$blastn" -task blastn-short -query "$out/kmers/core83.fa" \
    -subject "$out/intervals/${spec}.orthologous.fa" -dust no -word_size 11 \
    -perc_identity 100 -qcov_hsp_perc 100 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
    > "$out/kmers/${spec}.core83.exact.tsv"
done

# Sequence-level duplication screening is performed by the post-processing
# script with reciprocal BLAST alignments collapsed to non-redundant pairs.

{
  "$minimap2" --version
  "$miniprot" --version
  "$blastn" -version | head -1
  "$rmask" -version | head -2
  printf 'C_lucerna_accession\tGCA_965644535.1\n'
  printf 'C_cuculus_accession\tGCA_971798655.1\n'
} > "$out/software_and_accessions.txt"
