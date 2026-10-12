#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 4 ]]; then echo "Usage: $0 REF SAMPLE_BAM OUT_GVCF TMPDIR" >&2; exit 64; fi
ref=$1; bam=$2; out=$3; tmp=$4
mkdir -p "$tmp" "$(dirname "$out")"
gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmp" HaplotypeCaller \
  -R "$ref" -I "$bam" -O "$out" -ERC GVCF --native-pair-hmm-threads 8
