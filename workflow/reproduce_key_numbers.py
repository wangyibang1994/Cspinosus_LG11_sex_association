#!/usr/bin/env python3
"""Recompute central manuscript counts from deposited source data.

This script is deliberately dependency-light and serves as the release smoke
test. It does not replace raw-read reanalysis.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        yield from csv.DictReader(handle)


def wilson(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [centre - half, centre + half]


root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parents[1]).resolve()
sd = root / "source_data"
out = root / "results" / "key_numbers.json"

windows = list(rows(sd / "Figures3_4_window_metrics.csv"))
maxima = {}
for metric in ("WEIGHTED_FST", "MEAN_ED", "N_DELTA_AF_GE_0.3"):
    winner = max(windows, key=lambda x: float(x[metric]))
    maxima[metric] = {
        "value": float(winner[metric]),
        "chrom": winner["CHROM"],
        "start": int(float(winner["START"])),
        "end": int(float(winner["END"])),
    }

pca = list(rows(sd / "Figure5_LG11_PCA.csv"))
heatmap = list(rows(sd / "Figure6_top120_genotypes_long.csv"))
kmers = list(rows(sd / "FigureS2_strict_per_kmer_consistency.csv"))
validation = list(rows(sd / "Figure8_external60_final_sex_validation.csv"))

sex_counts = Counter(r["true_sex"] for r in validation)
snp_correct = sum(r["snp_correct"] == "Yes" for r in validation)
conclusive = [r for r in validation if r["kmer_conclusive"] == "Yes"]
kmer_correct = sum(r["kmer_strict_correct"] == "Yes" for r in conclusive)
inconclusive = [r["sample"] for r in validation if r["kmer_conclusive"] != "Yes"]

result = {
    "genome_scan_maxima": maxima,
    "local_pca_samples": len(pca),
    "heatmap_samples": len({r["sample"] for r in heatmap}),
    "heatmap_marker_ranks": len({r["rank"] for r in heatmap}),
    "core_kmers": len(kmers),
    "core_kmer_rule_check": {
        "all_males_ge_29": all(int(r["n_male_individuals"]) >= 29 for r in kmers),
        "all_females_eq_0": all(int(r["n_female_individuals"]) == 0 for r in kmers),
    },
    "external_validation": {
        "n": len(validation),
        "sex_counts": dict(sex_counts),
        "snp_correct": snp_correct,
        "snp_accuracy": snp_correct / len(validation),
        "snp_accuracy_wilson95": wilson(snp_correct, len(validation)),
        "kmer_conclusive_n": len(conclusive),
        "kmer_conclusive_correct": kmer_correct,
        "kmer_conclusive_accuracy": kmer_correct / len(conclusive),
        "kmer_conclusive_accuracy_wilson95": wilson(kmer_correct, len(conclusive)),
        "kmer_coverage": len(conclusive) / len(validation),
        "kmer_inconclusive_samples": inconclusive,
    },
}

assert len(pca) == 60
assert len({r["sample"] for r in heatmap}) == 60
assert len({r["rank"] for r in heatmap}) == 120
assert len(kmers) == 83
assert len(validation) == 60
assert sex_counts == {"Female": 37, "Male": 23}
assert snp_correct == 54
assert len(conclusive) == 58 and kmer_correct == 56
assert set(inconclusive) == {"BY1", "BY5"}

out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, ensure_ascii=False))
print(f"WROTE {out}")
