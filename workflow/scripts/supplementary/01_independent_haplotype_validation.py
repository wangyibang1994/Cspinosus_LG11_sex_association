#!/usr/bin/env python3
"""Discovery-frozen haplotype validation using separately phased cohorts.

The discovery and external cohorts are read from independent Beagle outputs.
Marker selection, male-allele orientation, K=2 score centres and the class
boundary are estimated from the discovery cohort only. The frozen rule is
then applied unchanged to the external cohort.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


def read_sex(path: Path, sample_col: str, sex_col: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = row[sex_col].strip().upper()
            out[row[sample_col]] = "F" if value.startswith("F") else "M"
    return out


def read_phased_vcf(path: Path):
    samples: list[str] = []
    records: dict[tuple[str, int, str, str], list[tuple[int | None, int | None]]] = {}
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
            genotypes = []
            for cell in fields[9:]:
                gt = cell.split(":")[gt_i]
                if "|" not in gt:
                    genotypes.append((None, None))
                    continue
                left, right = gt.split("|")
                genotypes.append(
                    (None if left == "." else int(left), None if right == "." else int(right))
                )
            records[(fields[0], int(fields[1]), fields[3], fields[4])] = genotypes
    return samples, records


def kmeans_two(values: list[float]) -> tuple[float, float]:
    centres = [min(values), max(values)]
    for _ in range(200):
        groups = [[], []]
        for value in values:
            group = 0 if abs(value - centres[0]) <= abs(value - centres[1]) else 1
            groups[group].append(value)
        updated = [sum(group) / len(group) if group else centres[i] for i, group in enumerate(groups)]
        if max(abs(updated[i] - centres[i]) for i in (0, 1)) < 1e-14:
            centres = updated
            break
        centres = updated
    return tuple(sorted(centres))


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = success / total
    denom = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return centre - half, centre + half


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery-vcf", type=Path, required=True)
    parser.add_argument("--external-vcf", type=Path, required=True)
    parser.add_argument("--discovery-sex", type=Path, required=True)
    parser.add_argument("--external-sex", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    discovery_sex = read_sex(args.discovery_sex, "Sample", "Group")
    external_sex = read_sex(args.external_sex, "sample_id", "sex")
    discovery_samples, discovery = read_phased_vcf(args.discovery_vcf)
    external_samples, external = read_phased_vcf(args.external_vcf)
    shared = sorted(set(discovery).intersection(external), key=lambda x: (x[0], x[1], x[2], x[3]))

    marker_rows = []
    selected = []
    for key in shared:
        calls = discovery[key]
        alleles = {"F": [], "M": []}
        for index, sample in enumerate(discovery_samples):
            sex = discovery_sex[sample]
            alleles[sex].extend(a for a in calls[index] if a is not None)
        if not alleles["F"] or not alleles["M"]:
            continue
        af_f = sum(alleles["F"]) / len(alleles["F"])
        af_m = sum(alleles["M"]) / len(alleles["M"])
        delta = af_m - af_f
        male_allele = 1 if delta > 0 else 0
        use = abs(delta) >= 0.30
        marker_rows.append((*key, af_f, af_m, delta, male_allele, use))
        if use:
            selected.append((key, male_allele))

    if len(selected) < 5:
        raise SystemExit(f"Only {len(selected)} discovery-defined markers passed |delta AF| >= 0.30")

    def score_cohort(samples, records, cohort, sex_map):
        rows = []
        for index, sample in enumerate(samples):
            for haplotype in (0, 1):
                matches = []
                for key, male_allele in selected:
                    allele = records[key][index][haplotype]
                    if allele is not None:
                        matches.append(allele == male_allele)
                rows.append({
                    "sample_id": sample,
                    "cohort": cohort,
                    "gonadal_sex": sex_map[sample],
                    "haplotype": haplotype + 1,
                    "male_allele_score": sum(matches) / len(matches) if matches else math.nan,
                    "n_informative_markers": len(matches),
                })
        return rows

    discovery_rows = score_cohort(discovery_samples, discovery, "discovery", discovery_sex)
    external_rows = score_cohort(external_samples, external, "external", external_sex)
    discovery_values = [r["male_allele_score"] for r in discovery_rows if not math.isnan(r["male_allele_score"])]
    low_centre, high_centre = kmeans_two(discovery_values)
    threshold = (low_centre + high_centre) / 2

    all_rows = discovery_rows + external_rows
    for row in all_rows:
        row["haplotype_class"] = (
            "male-associated" if row["male_allele_score"] >= threshold else "other"
        )

    # Export the discovery-defined marker matrix for a fully reproducible
    # neighbour-joining analysis. Cohorts remain independently phased; rows
    # are chromosome copies and columns are the same 219 frozen markers.
    matrix_rows = []
    for cohort, samples, records, sex_map in (
        ("discovery", discovery_samples, discovery, discovery_sex),
        ("external", external_samples, external, external_sex),
    ):
        for index, sample in enumerate(samples):
            for haplotype in (0, 1):
                score_row = next(
                    r for r in all_rows
                    if r["cohort"] == cohort and r["sample_id"] == sample and r["haplotype"] == haplotype + 1
                )
                matrix_rows.append({
                    "tip_id": f"{cohort}|{sample}|H{haplotype + 1}",
                    "sample_id": sample,
                    "cohort": cohort,
                    "gonadal_sex": sex_map[sample],
                    "haplotype": haplotype + 1,
                    "haplotype_class": score_row["haplotype_class"],
                    **{
                        f"{key[0]}:{key[1]}:{key[2]}:{key[3]}": records[key][index][haplotype]
                        for key, _male_allele in selected
                    },
                })

    individual_rows = []
    for cohort, samples in (("discovery", discovery_samples), ("external", external_samples)):
        for sample in samples:
            rows = [r for r in all_rows if r["cohort"] == cohort and r["sample_id"] == sample]
            n_high = sum(r["haplotype_class"] == "male-associated" for r in rows)
            individual_rows.append({
                "sample_id": sample,
                "cohort": cohort,
                "gonadal_sex": rows[0]["gonadal_sex"],
                "max_male_allele_score": max(r["male_allele_score"] for r in rows),
                "n_male_associated_haplotypes": n_high,
                "male_associated_carrier": n_high > 0,
            })

    with (args.outdir / "independent_haplotype_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(all_rows)
    with (args.outdir / "independent_haplotype_individuals.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=individual_rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(individual_rows)
    with (args.outdir / "discovery_defined_markers.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "AF_female_discovery", "AF_male_discovery", "delta_AF", "male_allele", "selected"])
        writer.writerows(marker_rows)
    with (args.outdir / "independent_haplotype_matrix.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=matrix_rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(matrix_rows)

    summary_rows = []
    for cohort in ("discovery", "external"):
        for sex in ("F", "M"):
            subset = [r for r in individual_rows if r["cohort"] == cohort and r["gonadal_sex"] == sex]
            carriers = sum(r["male_associated_carrier"] for r in subset)
            summary_rows.append((cohort, sex, len(subset), carriers, carriers / len(subset)))
    ext = [r for r in individual_rows if r["cohort"] == "external"]
    correct = sum((r["gonadal_sex"] == "M") == r["male_associated_carrier"] for r in ext)
    ci_low, ci_high = wilson(correct, len(ext))
    with (args.outdir / "independent_haplotype_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cohort", "sex", "n", "male_associated_carriers", "carrier_frequency"])
        writer.writerows(summary_rows)
    with (args.outdir / "frozen_rule_and_accuracy.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows([
            ("shared_exact_markers", len(shared)),
            ("discovery_selected_markers_abs_delta_AF_ge_0.30", len(selected)),
            ("discovery_low_score_centre", low_centre),
            ("discovery_high_score_centre", high_centre),
            ("frozen_score_threshold", threshold),
            ("external_correct", correct),
            ("external_total", len(ext)),
            ("external_accuracy", correct / len(ext)),
            ("external_wilson95_low", ci_low),
            ("external_wilson95_high", ci_high),
        ])

    print(f"shared={len(shared)} selected={len(selected)} centres={low_centre:.6f},{high_centre:.6f} threshold={threshold:.6f}")
    print(f"external={correct}/{len(ext)} ({correct/len(ext):.4%}); Wilson95={ci_low:.4%}-{ci_high:.4%}")
    for row in summary_rows:
        print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
