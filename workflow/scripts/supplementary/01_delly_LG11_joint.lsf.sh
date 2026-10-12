#!/usr/bin/env bash
#BSUB -J lg11_delly
#BSUB -q normal
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=4000]"
#BSUB -M 32000
#BSUB -W 24:00
#BSUB -o 10_logs/lg11_delly.%J.out
#BSUB -e 10_logs/lg11_delly.%J.err
set -euo pipefail

project=${CSPINOSUS_DISCOVERY_PROJECT:?Set CSPINOSUS_DISCOVERY_PROJECT to the discovery-project directory}
out="$project/12_supplementary_analyses/01_LG11_SV"
delly=${DELLY:-delly}
bcftools=${BCFTOOLS:-bcftools}
ref="$project/02_reference/genome.fa"
mkdir -p "$out" "$project/10_logs"

mapfile -t bams < <(find "$project/04_mapping/bam" -maxdepth 1 -name '*.sorted.bam' -type f | sort)
[[ ${#bams[@]} -eq 60 ]] || { echo "Expected 60 BAMs, found ${#bams[@]}" >&2; exit 1; }

samtools=${SAMTOOLS:-samtools}
subset="$out/region_bams"
mkdir -p "$subset"
export subset samtools
printf '%s\n' "${bams[@]}" | xargs -n1 -P8 -I{} bash -c '
  bam="$1"; base=$(basename "$bam" .sorted.bam)
  "$samtools" view -bh -o "$subset/${base}.LG11_4.0_6.3Mb.bam" "$bam" LG11:4000000-6300000
  "$samtools" index "$subset/${base}.LG11_4.0_6.3Mb.bam"
' _ {}
mapfile -t region_bams < <(find "$subset" -maxdepth 1 -name '*.bam' -type f | sort)
[[ ${#region_bams[@]} -eq 60 ]] || { echo "Expected 60 regional BAMs, found ${#region_bams[@]}" >&2; exit 1; }

"$delly" sr -t ALL -h 8 -g "$ref" \
  -o "$out/LG11_4.2_6.1Mb.joint.raw.bcf" "${region_bams[@]}"
"$delly" filter -f germline -o "$out/LG11_4.2_6.1Mb.joint.germline.bcf" \
  "$out/LG11_4.2_6.1Mb.joint.raw.bcf"
"$bcftools" index -f "$out/LG11_4.2_6.1Mb.joint.germline.bcf"
"$bcftools" view -r LG11:4200000-6100000 -Oz -o "$out/LG11_4.2_6.1Mb.joint.germline.vcf.gz" \
  "$out/LG11_4.2_6.1Mb.joint.germline.bcf"
"$bcftools" index -f "$out/LG11_4.2_6.1Mb.joint.germline.vcf.gz"

"$bcftools" query -f '%CHROM\t%POS\t%ID\t%SVTYPE\t%END\t%QUAL[\t%SAMPLE\t%GT\t%FT\t%DV\t%DR\t%RV\t%RR]\n' \
  "$out/LG11_4.2_6.1Mb.joint.germline.vcf.gz" > "$out/LG11_SV_genotypes.wide.tsv"

{
  "$delly" --version
  "$bcftools" --version | head -2
  printf 'region\tLG11:4200000-6100000\n'
  printf 'n_bams\t%s\n' "${#bams[@]}"
} > "$out/software_and_parameters.txt"
