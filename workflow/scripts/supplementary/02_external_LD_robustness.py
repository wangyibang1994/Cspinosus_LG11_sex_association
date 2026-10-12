#!/usr/bin/env python3
"""Balanced-sample and sex-label-permutation robustness tests for external LD.

The analysis uses a fixed random sample of up to 2,000 SNP pairs per physical
distance bin. External females are repeatedly downsampled to the male sample
size. Sex labels are then permuted at the individual level while preserving
the observed 23:37 group sizes. The biological unit is the individual fish;
SNP pairs are not treated as independent replicates.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np


EDGES = np.array([0, 10_000, 25_000, 50_000, 100_000, 250_000, 500_000, 700_001])


def read_sex(path: Path) -> dict[str, str]:
    out = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = row["sex"].strip().upper()
            out[row["sample_id"]] = "F" if value.startswith("F") else "M"
    return out


def read_vcf(path: Path):
    samples = []
    positions = []
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            fields = line.rstrip().split("\t")
            fmt = fields[8].split(":")
            gt_i = fmt.index("GT")
            dosage = []
            for cell in fields[9:]:
                gt = cell.split(":")[gt_i].replace("|", "/")
                if "." in gt:
                    dosage.append(np.nan)
                else:
                    dosage.append(sum(int(x) for x in gt.split("/")))
            positions.append(int(fields[1]))
            rows.append(dosage)
    return samples, np.asarray(positions), np.asarray(rows, dtype=float).T


def r2_for_pairs(genotypes: np.ndarray, individuals: np.ndarray, left: np.ndarray, right: np.ndarray, chunk=5000):
    result = np.full(left.size, np.nan)
    for start in range(0, left.size, chunk):
        stop = min(start + chunk, left.size)
        x = genotypes[np.ix_(individuals, left[start:stop])].T
        y = genotypes[np.ix_(individuals, right[start:stop])].T
        valid = np.isfinite(x) & np.isfinite(y)
        n = valid.sum(axis=1)
        x0 = np.where(valid, x, 0.0); y0 = np.where(valid, y, 0.0)
        sx = x0.sum(axis=1); sy = y0.sum(axis=1)
        sxx = (x0 * x0).sum(axis=1); syy = (y0 * y0).sum(axis=1)
        sxy = (x0 * y0).sum(axis=1)
        cov_num = sxy - sx * sy / np.maximum(n, 1)
        vx = sxx - sx * sx / np.maximum(n, 1)
        vy = syy - sy * sy / np.maximum(n, 1)
        denominator = vx * vy
        ok = (n >= 4) & (denominator > 0)
        local = np.full(stop - start, np.nan)
        local[ok] = (cov_num[ok] * cov_num[ok]) / denominator[ok]
        result[start:stop] = local
    return result


def bin_means(genotypes, individuals, pair_sets):
    means = []
    for left, right in pair_sets:
        values = r2_for_pairs(genotypes, individuals, left, right)
        means.append(float(np.nanmean(values)))
    return np.asarray(means)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--sex", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--pairs-per-bin", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    sex = read_sex(args.sex)
    samples, positions, genotypes = read_vcf(args.vcf)
    female = np.asarray([i for i, s in enumerate(samples) if sex[s] == "F"])
    male = np.asarray([i for i, s in enumerate(samples) if sex[s] == "M"])
    if not (len(male) == 23 and len(female) == 37):
        raise SystemExit(f"Expected 23 males and 37 females, observed {len(male)} and {len(female)}")

    left_all, right_all = np.triu_indices(len(positions), 1)
    distance = np.abs(positions[right_all] - positions[left_all])
    pair_sets = []
    pair_manifest = []
    for bin_i in range(len(EDGES) - 1):
        eligible = np.flatnonzero((distance >= EDGES[bin_i]) & (distance < EDGES[bin_i + 1]))
        chosen = rng.choice(eligible, size=min(args.pairs_per_bin, len(eligible)), replace=False)
        pair_sets.append((left_all[chosen], right_all[chosen]))
        pair_manifest.append((bin_i, EDGES[bin_i], EDGES[bin_i + 1] - 1, len(eligible), len(chosen)))

    male_observed = bin_means(genotypes, male, pair_sets)
    female_observed = bin_means(genotypes, female, pair_sets)
    female_downsample = np.empty((args.replicates, len(pair_sets)))
    permuted_diff = np.empty((args.replicates, len(pair_sets)))
    all_individuals = np.arange(len(samples))
    for rep in range(args.replicates):
        female_subset = rng.choice(female, size=len(male), replace=False)
        female_downsample[rep, :] = bin_means(genotypes, female_subset, pair_sets)
        permuted_male = rng.choice(all_individuals, size=len(male), replace=False)
        is_male = np.zeros(len(samples), dtype=bool); is_male[permuted_male] = True
        permuted_female = all_individuals[~is_male]
        permuted_diff[rep, :] = (
            bin_means(genotypes, permuted_male, pair_sets)
            - bin_means(genotypes, permuted_female, pair_sets)
        )

    observed_diff = male_observed - female_observed
    observed_global = float(np.mean(observed_diff))
    permuted_global = np.mean(permuted_diff, axis=1)
    global_p = (1 + np.sum(permuted_global >= observed_global)) / (args.replicates + 1)

    with (args.outdir / "LD_pair_sampling_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["bin", "distance_start_bp", "distance_end_bp", "eligible_pairs", "sampled_pairs"])
        writer.writerows(pair_manifest)
    with (args.outdir / "LD_downsample_replicates.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["replicate"] + [f"bin_{i+1}" for i in range(len(pair_sets))])
        for i, values in enumerate(female_downsample, 1):
            writer.writerow([i] + list(values))
    with (args.outdir / "LD_permutation_replicates.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["replicate"] + [f"bin_{i+1}" for i in range(len(pair_sets))] + ["global_mean_difference"])
        for i, values in enumerate(permuted_diff, 1):
            writer.writerow([i] + list(values) + [float(np.mean(values))])
    with (args.outdir / "LD_robustness_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["distance_start_bp", "distance_end_bp", "male_mean_r2", "female_full_mean_r2", "female_downsample_median", "female_downsample_q025", "female_downsample_q975", "downsample_exceedance_p", "sex_label_permutation_p"])
        for i in range(len(pair_sets)):
            writer.writerow([
                EDGES[i], EDGES[i+1] - 1, male_observed[i], female_observed[i],
                np.median(female_downsample[:, i]), np.quantile(female_downsample[:, i], 0.025),
                np.quantile(female_downsample[:, i], 0.975),
                (1 + np.sum(female_downsample[:, i] >= male_observed[i])) / (args.replicates + 1),
                (1 + np.sum(permuted_diff[:, i] >= observed_diff[i])) / (args.replicates + 1),
            ])
        writer.writerow(["GLOBAL", "GLOBAL", float(np.mean(male_observed)), float(np.mean(female_observed)), "NA", "NA", "NA", "NA", global_p])

    print(f"variants={len(positions)} females={len(female)} males={len(male)} replicates={args.replicates}")
    print(f"observed_global_mean_difference={observed_global:.6f}; permutation_p={global_p:.6g}")
    for i in range(len(pair_sets)):
        print(f"{EDGES[i]}-{EDGES[i+1]-1}\tM={male_observed[i]:.5f}\tF={female_observed[i]:.5f}\tF23_95={np.quantile(female_downsample[:,i],0.025):.5f}-{np.quantile(female_downsample[:,i],0.975):.5f}")


if __name__ == "__main__":
    main()
