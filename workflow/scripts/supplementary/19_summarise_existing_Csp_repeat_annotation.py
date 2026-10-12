#!/usr/bin/env python3
"""Summarise the original project RepeatMasker GFF3 across LG11.

Coordinates in the GFF3 are 1-based inclusive. Reported covered bases are the
union of clipped intervals, so overlapping fragments are not double-counted in
the overall coverage estimate. Class-specific unions are calculated separately.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REGIONS = {
    "LG11_candidate_4.80_5.50Mb": (4_800_000, 5_500_000),
    "core_kmer_cluster_5.168_5.336Mb": (5_168_000, 5_336_000),
    "upstream_flank_4.10_4.80Mb": (4_100_000, 4_800_000),
    "downstream_flank_5.50_6.20Mb": (5_500_000, 6_200_000),
}


def merge_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    cur_start, cur_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    return total + cur_end - cur_start + 1


def repeat_class(attributes: str) -> str:
    match = re.search(r"(?:^|;)Class=([^;]+)", attributes)
    return match.group(1) if match else "Unclassified"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--seqid", default="LG11")
    parser.add_argument("--window", type=int, default=100_000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    records: list[tuple[int, int, str]] = []
    with args.gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != args.seqid:
                continue
            records.append((int(fields[3]), int(fields[4]), repeat_class(fields[8])))

    summary_rows = []
    for name, (region_start, region_end) in REGIONS.items():
        clipped_all: list[tuple[int, int]] = []
        by_class: dict[str, list[tuple[int, int]]] = defaultdict(list)
        n_records = 0
        for start, end, cls in records:
            left, right = max(start, region_start), min(end, region_end)
            if left <= right:
                clipped_all.append((left, right))
                by_class[cls].append((left, right))
                n_records += 1
        region_bp = region_end - region_start + 1
        covered = merge_length(clipped_all)
        summary_rows.append(
            [name, region_start, region_end, region_bp, "All_repeats", n_records,
             covered, covered / region_bp]
        )
        for cls, intervals in sorted(by_class.items()):
            class_bp = merge_length(intervals)
            summary_rows.append(
                [name, region_start, region_end, region_bp, cls, len(intervals),
                 class_bp, class_bp / region_bp]
            )

    summary_path = args.outdir / "Csp_original_repeat_annotation.region_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "region", "start_1based", "end_1based", "region_bp", "repeat_class",
            "n_gff_records", "union_covered_bp", "coverage_fraction",
        ])
        writer.writerows(summary_rows)

    max_end = max((end for _, end, _ in records), default=0)
    windows_path = args.outdir / "Csp_original_repeat_annotation.LG11_100kb_windows.tsv"
    with windows_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["seqid", "start_1based", "end_1based", "window_bp",
                         "repeat_covered_bp", "coverage_fraction", "n_gff_records"])
        for start in range(1, max_end + 1, args.window):
            end = min(start + args.window - 1, max_end)
            clipped = []
            count = 0
            for rec_start, rec_end, _ in records:
                left, right = max(rec_start, start), min(rec_end, end)
                if left <= right:
                    clipped.append((left, right))
                    count += 1
            window_bp = end - start + 1
            covered = merge_length(clipped)
            writer.writerow([args.seqid, start, end, window_bp, covered,
                             covered / window_bp, count])

    subset_path = args.outdir / "Csp_original_repeat_annotation.candidate_interval.gff3"
    with args.gff.open() as source, subset_path.open("w") as target:
        target.write("##gff-version 3\n")
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if (len(fields) >= 9 and fields[0] == args.seqid and
                    int(fields[3]) <= 5_500_000 and int(fields[4]) >= 4_800_000):
                target.write(line)

    print(f"LG11 records: {len(records):,}")
    print(summary_path)
    print(windows_path)
    print(subset_path)


if __name__ == "__main__":
    main()
