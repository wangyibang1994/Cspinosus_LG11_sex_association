# Reproducibility audit and release blockers

## Overall judgement

The current materials are sufficient to reconstruct the central computational logic and to package the numerical source data. They are **not yet sufficient for an external researcher to reproduce every reported result from the deposited quality-filtered reads without author input**. Public release should occur only after the critical items below are resolved.

## Critical blockers

| Priority | Missing or inconsistent item | Why it matters | Required action |
|---|---|---|---|
| Critical | Discovery reference accession and annotation accession/version | Coordinates and annotations are otherwise not portable | Supply accession, source URL, release date and checksums |
| Critical | Discovery cohort read-cleaning software/parameters | Upstream cleaning cannot be reproduced from unprocessed instrument output | State clearly that only provider-delivered clean FASTQs are available, deposit those exact files, and preserve the provider QC summary |
| Critical | Ethics approval and any remaining mandatory archive metadata | Required Methods and archive metadata must be complete | Retain the confirmed sample manifests; resolve ethics approval and archive-assigned identifiers before submission |
| Critical | Exact mapping between discovery LG11 and external CM056831.1 reference coordinates | “Exact coordinate match” is uninterpretable across builds without provenance | Document whether assemblies are coordinate-identical, aligned, lifted over or separately defined; deposit crosswalk/alignment |
| Critical | Final scripts still contain hard-coded personal/server paths | Reviewers cannot run them | Parameterise every final script via command-line arguments/config |

## Major analytical issues to resolve or disclose

| Item | Current evidence | Resolution |
|---|---|---|
| Discovery genotype-level QC | Historical BCFtools workflow filters sites but does not mask low-DP/GQ genotypes | Add a pre-specified sensitivity analysis with genotype DP/GQ masking; preserve primary historical result |
| Coding consequence annotation | Custom script may not account correctly for CDS phase and transcript orientation | Reannotate with VEP/SnpEff or validate with phase-aware unit tests |
| Local heatmap marker count | Historical structure script used top200; final Figure 6 uses top120 | Freeze 120 in all public scripts and documentation |
| External k-mer cutoff | Midpoint 37.5 was derived from the external score gap | Describe as data-adaptive technical classification; report 58/60 coverage and do not present as a pre-registered diagnostic cutoff |
| BY1/BY5 | Both have dissected female sex but low k-mer evidence | Retain in full table; classify as inconclusive in strict k-mer analysis; do not delete post hoc |
| Global structure settings | Exact PLINK/ADMIXTURE versions/options incompletely recorded | Recover logs or label unresolved settings honestly |
| Duplicate handling | Discovery did not mark duplicates; external cohort did | State cohort-specific pipelines separately |
| Permutation exchangeability | Geographic/relatedness strata may violate unrestricted exchangeability | Assess structure/kinship; if needed, add stratified permutations or limitation statement |

## Frozen validation assets already packaged

- `metadata/frozen_top50_discovery_snps.csv`: first 50 SNPs in the final descending `|Delta AF|` ranking, with discovery-build chromosome, coordinate and alleles.
- `metadata/frozen_core83_kmers.csv`: all 83 core 31-mers with discovery-cohort sex-specific support and exact-reference alignment counts.
- `metadata/external60_dissection_validation.csv`: final 60-sample validation table with dissected sex, SNP prediction and k-mer result.

These files must receive checksums in the final release manifest and must not be regenerated using external phenotype labels.

## Files that must not be published as final analysis scripts

- Historical scripts that describe external sex records as phenotype-unknown.
- Scripts that load a personal workstation or institutional-server path without configuration.
- Scripts that exclude BY1/BY5 from the biological cohort or relabel them without preserving their dissected sex.
- Old Figure S3 leave-two-out analyses that were based on superseded discovery labels, unless retained under a clearly marked `legacy/` provenance directory and excluded from the manuscript evidence chain.

Earlier 58-sample external-validation source data and summary tables have been moved to `source_data/legacy_superseded/` and `supplementary_tables/legacy_superseded/`. They are retained for audit history but are not authoritative manuscript inputs.

## Release acceptance criteria

- `workflow/run_source_data_checks.sh` exits with code 0.
- No `AUTHOR_INPUT_REQUIRED` strings remain.
- No final script contains an absolute personal/server path.
- All figure/table inputs are present and checksummed.
- The manuscript's sample counts, interval coordinates, SNP counts, marker counts and external accuracy denominators agree with repository outputs.
- Both reference FASTA checksums and the GFF3 checksum are recorded.
- Raw-read accessions resolve publicly or through a reviewer link.
- A clean-environment smoke test reproduces summary tables and all final figures.
