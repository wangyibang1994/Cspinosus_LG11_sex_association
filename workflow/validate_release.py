#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
required = [
    "README.md",
    "docs/WORKFLOW_CN.md",
    "docs/WORKFLOW_EN.md",
    "docs/DATA_DEPOSITION_PLAN.md",
    "docs/REPRODUCIBILITY_AUDIT.md",
    "config/config.example.yaml",
    "metadata/discovery_samples.template.tsv",
    "metadata/external_samples.template.tsv",
]
missing = [p for p in required if not (root / p).is_file()]
if missing:
    raise SystemExit("Missing required files: " + ", ".join(missing))

source = sorted(
    p for p in (root / "source_data").glob("*")
    if p.is_file() and p.name != "README.md"
)
tables = sorted(
    p for p in (root / "supplementary_tables").glob("*")
    if p.is_file() and p.name != "README.md"
)
if not source:
    raise SystemExit("source_data/ is empty")
if not tables:
    raise SystemExit("supplementary_tables/ is empty")

manifest = root / "metadata" / "release_file_checksums.tsv"
with manifest.open("w", encoding="utf-8") as out:
    out.write("relative_path\tbytes\tsha256\n")
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p != manifest):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        out.write(f"{path.relative_to(root)}\t{path.stat().st_size}\t{digest}\n")

placeholder_hits = []
for path in root.rglob("*"):
    if path.is_file() and path.suffix.lower() in {".md", ".yaml", ".yml", ".tsv", ".cff"}:
        if "AUTHOR_INPUT_REQUIRED" in path.read_text(encoding="utf-8", errors="ignore"):
            placeholder_hits.append(str(path.relative_to(root)))
print(f"PASS: {len(source)} source-data files and {len(tables)} supplementary-table files found")
print(f"WROTE: {manifest.relative_to(root)}")
if placeholder_hits:
    print("RELEASE BLOCKED: unresolved AUTHOR_INPUT_REQUIRED fields in:")
    for item in placeholder_hits:
        print(f"  - {item}")
    raise SystemExit(2)
