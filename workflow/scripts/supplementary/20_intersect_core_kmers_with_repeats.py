#!/usr/bin/env python3
"""Intersect exact core-31-mer alignments with the original repeat GFF3."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


def get_attribute(text: str, key: str, default: str = "") -> str:
    match = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", text)
    return match.group(1) if match else default


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignments", required=True, type=Path)
    parser.add_argument("--repeat-gff", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--total-panel-markers", type=int, default=83)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    repeats = []
    with args.repeat_gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != "LG11":
                continue
            start, end = int(fields[3]), int(fields[4])
            if start <= 5_500_000 and end >= 4_800_000:
                repeats.append({
                    "start": start,
                    "end": end,
                    "id": get_attribute(fields[8], "ID"),
                    "target": get_attribute(fields[8], "Target"),
                    "class": get_attribute(fields[8], "Class", "Unclassified"),
                })

    rows = []
    marker_classes: dict[str, set[str]] = defaultdict(set)
    mapped_markers = set()
    repeat_markers = set()
    repeat_alignment_count = 0
    overlapping_repeat_ids = set()
    overlapping_repeat_targets = set()
    class_alignment_counts = Counter()
    with args.alignments.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            marker, subject = fields[0], fields[1]
            region = re.fullmatch(r"([^:]+):(\d+)-(\d+)", subject)
            if not region:
                raise ValueError(f"Unexpected subject name: {subject}")
            seqid, region_start = region.group(1), int(region.group(2))
            rel_start, rel_end = int(fields[6]), int(fields[7])
            genomic_start = region_start + min(rel_start, rel_end) - 1
            genomic_end = region_start + max(rel_start, rel_end) - 1
            hits = [r for r in repeats
                    if r["start"] <= genomic_end and r["end"] >= genomic_start]
            mapped_markers.add(marker)
            if hits:
                repeat_alignment_count += 1
                repeat_markers.add(marker)
                for hit in hits:
                    marker_classes[marker].add(hit["class"])
                    class_alignment_counts[hit["class"]] += 1
                    if hit["id"]:
                        overlapping_repeat_ids.add(hit["id"])
                    if hit["target"]:
                        overlapping_repeat_targets.add(hit["target"])
            rows.append({
                "kmer_id": marker,
                "seqid": seqid,
                "genomic_start_1based": genomic_start,
                "genomic_end_1based": genomic_end,
                "strand": "+" if rel_start <= rel_end else "-",
                "n_overlapping_repeat_records": len(hits),
                "repeat_classes": ";".join(sorted({h["class"] for h in hits})) or "None",
                "repeat_ids": ";".join(sorted({h["id"] for h in hits if h["id"]})) or "None",
                "repeat_targets": ";".join(sorted({h["target"] for h in hits if h["target"]})) or "None",
            })

    detail_path = args.outdir / "Csp_core83_exact_alignments.repeat_overlap.tsv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    marker_path = args.outdir / "Csp_core83_marker_level.repeat_overlap.tsv"
    with marker_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["kmer_id", "has_exact_reference_alignment", "overlaps_repeat",
                         "repeat_classes"])
        for index in range(1, args.total_panel_markers + 1):
            marker = f"strict_{index}"
            writer.writerow([
                marker,
                int(marker in mapped_markers),
                int(marker in repeat_markers),
                ";".join(sorted(marker_classes[marker])) or "None",
            ])

    summary_path = args.outdir / "Csp_core83_repeat_overlap.summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["panel_markers", args.total_panel_markers])
        writer.writerow(["markers_with_exact_reference_alignment", len(mapped_markers)])
        writer.writerow(["exact_reference_alignments", len(rows)])
        writer.writerow(["markers_with_any_repeat_overlap", len(repeat_markers)])
        writer.writerow(["exact_alignments_with_any_repeat_overlap", repeat_alignment_count])
        writer.writerow(["distinct_overlapping_repeat_records", len(overlapping_repeat_ids)])
        writer.writerow(["distinct_overlapping_repeat_targets", len(overlapping_repeat_targets)])
        for cls, count in sorted(class_alignment_counts.items()):
            writer.writerow([f"repeat_overlap_alignments__{cls}", count])

    print(summary_path)


if __name__ == "__main__":
    main()
