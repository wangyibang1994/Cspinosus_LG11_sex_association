#!/usr/bin/env python3
"""One label permutation for genome-wide calibration of the sex-association scan.

The calculation exactly matches the discovery scan: corrected labels, all 60
individuals, Hudson-style weighted FST, mean ED, 200-kb windows and 50-kb
steps.  It writes only one maximum per metric, avoiding unnecessary large
intermediate site tables.
"""

from __future__ import annotations

import gzip
import math
import random
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 5:
    raise SystemExit("Usage: genome_wide_permutation_calibration.py <vcf_dir> <metadata.tsv> <seed> <output.tsv>")

vcf_dir, metadata, seed_text, output_path = sys.argv[1:]
vcf_dir, metadata, output_path = map(Path, (vcf_dir, metadata, output_path))
seed = int(seed_text)

sex = {}
with metadata.open() as handle:
    next(handle)
    for line in handle:
        _, sample, group = line.rstrip().split("\t")
        sex[sample] = group
female_samples = {sample for sample, group in sex.items() if group == "F"}
male_samples = {sample for sample, group in sex.items() if group == "M"}
if len(female_samples) != 30 or len(male_samples) != 30:
    raise SystemExit("Expected 30 female-labelled and 30 male-labelled discovery samples")

rng = random.Random(seed)
all_samples = sorted(sex)
permuted_female_samples = set(rng.sample(all_samples, len(female_samples)))


def dosage(gt: str):
    gt = gt.replace("|", "/")
    return {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}.get(gt)


sites = defaultdict(list)
for chrom_i in range(1, 25):
    chrom = f"LG{chrom_i}"
    vcf_path = vcf_dir / f"{chrom}.filtered.snps.vcf.gz"
    if not vcf_path.exists():
        raise SystemExit(f"Missing {vcf_path}")
    with gzip.open(vcf_path, "rt") as handle:
        sample_names = None
        female_indices = None
        male_indices = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip().split("\t")[9:]
                female_indices = [i for i, sample in enumerate(sample_names) if sample in permuted_female_samples]
                male_indices = [i for i, sample in enumerate(sample_names) if sample not in permuted_female_samples]
                if len(female_indices) != 30 or len(male_indices) != 30:
                    raise SystemExit(f"Unexpected sample count in {vcf_path}")
                continue
            fields = line.rstrip().split("\t")
            fmt = fields[8].split(":")
            gt_i = fmt.index("GT")
            values = []
            for sample_field in fields[9:]:
                entries = sample_field.split(":")
                values.append(dosage(entries[gt_i]) if gt_i < len(entries) else None)
            fvals = [values[i] for i in female_indices if values[i] is not None]
            mvals = [values[i] for i in male_indices if values[i] is not None]
            if len(fvals) < 2 or len(mvals) < 2:
                continue
            p_f = sum(fvals) / (2 * len(fvals))
            p_m = sum(mvals) / (2 * len(mvals))
            delta = abs(p_f - p_m)
            ed = math.sqrt(2) * delta
            n_f, n_m = 2 * len(fvals), 2 * len(mvals)
            numerator = (p_f - p_m) ** 2 - p_f * (1 - p_f) / (n_f - 1) - p_m * (1 - p_m) / (n_m - 1)
            denominator = p_f * (1 - p_m) + p_m * (1 - p_f)
            sites[chrom].append((int(fields[1]), numerator, denominator, delta, ed))

win, step = 200000, 50000
maxima = {
    "WEIGHTED_FST": (-float("inf"), None),
    "MEAN_ED": (-float("inf"), None),
    "N_DELTA_AF_GE_0.3": (-float("inf"), None),
}
for chrom_i in range(1, 25):
    chrom = f"LG{chrom_i}"
    records = sites[chrom]
    left = right = 0
    for start in range(0, records[-1][0] + 1, step):
        end = start + win
        while left < len(records) and records[left][0] < start:
            left += 1
        right = max(right, left)
        while right < len(records) and records[right][0] < end:
            right += 1
        window = records[left:right]
        if not window:
            continue
        fst_den = sum(item[2] for item in window)
        values = {
            "WEIGHTED_FST": sum(item[1] for item in window) / fst_den if fst_den > 0 else float("nan"),
            "MEAN_ED": sum(item[4] for item in window) / len(window),
            "N_DELTA_AF_GE_0.3": sum(item[3] >= 0.3 for item in window),
        }
        for metric, value in values.items():
            if math.isfinite(value) and value > maxima[metric][0]:
                maxima[metric] = (value, (chrom, start, end))

output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w") as handle:
    handle.write("seed\tmetric\tmax_value\tchrom\tstart\tend\n")
    for metric, (value, location) in maxima.items():
        chrom, start, end = location
        handle.write(f"{seed}\t{metric}\t{value:.10g}\t{chrom}\t{start}\t{end}\n")
