#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 6 ]]; then echo "Usage: $0 REF SAMPLE R1 R2 OUT_BAM METRICS" >&2; exit 64; fi
ref=$1; sample=$2; r1=$3; r2=$4; out=$5; metrics=$6
mkdir -p "$(dirname "$out")" "$(dirname "$metrics")"
sorted="${out%.markdup.bam}.sorted.bam"
bwa mem -M -t "${BWA_THREADS:-12}" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$ref" "$r1" "$r2" |
  samtools sort -@ "${SORT_THREADS:-4}" -m 2G -o "$sorted" -
samtools index "$sorted"
gatk MarkDuplicates -I "$sorted" -O "$out" -M "$metrics" --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
samtools quickcheck -v "$out"
