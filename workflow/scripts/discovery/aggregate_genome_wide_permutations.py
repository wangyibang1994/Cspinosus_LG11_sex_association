#!/usr/bin/env python3
"""Aggregate permutation maxima and calculate family-wise empirical P values."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: aggregate_genome_wide_permutations.py <primary_window_metrics.csv> <replicate_dir> <summary.csv> <null_long.csv>"
    )

primary_path, replicate_dir, summary_path, null_path = map(Path, sys.argv[1:])
metrics = ("WEIGHTED_FST", "MEAN_ED", "N_DELTA_AF_GE_0.3")

observed = {}
locations = {}
with primary_path.open(newline="") as handle:
    for row in csv.DictReader(handle):
        for metric in metrics:
            value = float(row[metric])
            if metric not in observed or value > observed[metric]:
                observed[metric] = value
                locations[metric] = (row["CHROM"], int(row["START"]), int(row["END"]))

by_metric: dict[str, list[dict[str, str]]] = defaultdict(list)
for path in sorted(replicate_dir.glob("perm_*.tsv")):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 3 or {row["metric"] for row in rows} != set(metrics):
        raise SystemExit(f"Malformed replicate: {path}")
    for row in rows:
        by_metric[row["metric"]].append(row)

if any(len(by_metric[metric]) != 500 for metric in metrics):
    raise SystemExit("Expected exactly 500 valid replicates per metric")

null_path.parent.mkdir(parents=True, exist_ok=True)
with null_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["seed", "metric", "max_value", "chrom", "start", "end"])
    writer.writeheader()
    for metric in metrics:
        writer.writerows(sorted(by_metric[metric], key=lambda row: int(row["seed"])))

summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w", newline="") as handle:
    columns = [
        "metric", "observed_genome_max", "observed_chrom", "observed_start", "observed_end",
        "n_permutations", "null_max_median", "null_max_q95", "null_max_q99", "n_null_ge_observed", "empirical_fwer_p",
    ]
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    for metric in metrics:
        values = sorted(float(row["max_value"]) for row in by_metric[metric])
        n = len(values)
        n_ge = sum(value >= observed[metric] for value in values)
        quantile = lambda p: values[round((n - 1) * p)]
        chrom, start, end = locations[metric]
        writer.writerow({
            "metric": metric,
            "observed_genome_max": f"{observed[metric]:.10g}",
            "observed_chrom": chrom,
            "observed_start": start,
            "observed_end": end,
            "n_permutations": n,
            "null_max_median": f"{quantile(0.5):.10g}",
            "null_max_q95": f"{quantile(0.95):.10g}",
            "null_max_q99": f"{quantile(0.99):.10g}",
            "n_null_ge_observed": n_ge,
            "empirical_fwer_p": f"{(1 + n_ge) / (1 + n):.10g}",
        })
