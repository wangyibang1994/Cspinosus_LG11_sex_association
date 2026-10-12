#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 5 ]]; then echo "Usage: $0 REF SAMPLE R1 R2 OUT_BAM" >&2; exit 64; fi
ref=$1; sample=$2; r1=$3; r2=$4; out=$5
mkdir -p "$(dirname "$out")"
tmp="${out%.bam}.sorttmp"
bwa mem -t "${BWA_THREADS:-16}" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$ref" "$r1" "$r2" |
  samtools sort -@ "${SORT_THREADS:-4}" -T "$tmp" -o "$out" -
samtools index "$out"
samtools quickcheck -v "$out"
samtools flagstat "$out" > "${out%.bam}.flagstat.txt"
