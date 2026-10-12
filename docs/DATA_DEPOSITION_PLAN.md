# Data and code deposition plan

## Recommended public architecture

| Object | Repository | Public identifier | Release timing |
|---|---|---|---|
| Discovery quality-filtered FASTQ | NCBI SRA, ENA or NGDC GSA | BioProject + BioSample + Run accessions | Private/reviewer link at submission; public no later than publication |
| External 60-sample quality-filtered FASTQ | Same archive and preferably same BioProject | BioSample + Run accessions | As above |
| Reference assembly used for discovery | Existing primary accession if already public; otherwise GenBank/ENA/GSA assembly archive | Assembly accession | Before review if required for reanalysis |
| External reference | NCBI Assembly | GCA_029853015.1 | Already identifiable; verify exact downloaded file checksum |
| Filtered VCFs and frozen marker sets | Zenodo/Figshare/Dryad | DOI | At submission or acceptance according to journal policy |
| Figure source data and supplementary tables | Journal source-data files plus Zenodo release | Article link + DOI | Submission |
| Analysis code, configs and environments | GitHub/GitLab, archived by Zenodo | Version tag + concept/version DOI | Freeze `v1.0.0` at submission |

The deposited sequence files are provider-delivered, quality-filtered Illumina FASTQs rather than unprocessed instrument output. They belong in a recognised sequence archive; GitHub is not an appropriate home for FASTQ files. NCBI SRA accepts high-throughput sequence reads and links them to BioProject and BioSample records. A versioned GitHub release can be archived through Zenodo to obtain a stable DOI.

## Files to deposit with the code DOI

- All scripts used for final analyses and figures.
- Exact configuration file with no local absolute paths.
- Conda environment YAML plus an exported explicit package list or container digest.
- Final discovery sex metadata and external dissection-sex metadata.
- Discovery reference FASTA accession and SHA-256; external assembly accession and checksum.
- Candidate-interval and genome-wide filtered VCFs, subject to repository size limits.
- Frozen top-50 SNP panel with reference build, chromosome, coordinate, REF, ALT and selection rank.
- Frozen 83-k-mer FASTA/TSV with marker identifiers and SHA-256.
- Machine-readable source data for every main and supplementary figure.
- Supplementary tables as TSV/CSV, not only XLSX.
- A file-level checksum manifest.

## Metadata preparation

Create one BioProject for the study and one BioSample per biological individual. Do not use sequencing-library identifiers as biological sample identifiers. The same stable sample ID must join FASTQ, BAM/VCF, phenotype, k-mer score and figure source data. Record organism taxonomy, collection locality/date, tissue, biological sex, sex determination method, library strategy, library source, library selection, platform, layout and file checksums.

## Suggested release procedure

1. Resolve every `AUTHOR_INPUT_REQUIRED` item in this repository.
2. Upload the quality-filtered FASTQ files and metadata to SRA/GSA; obtain accession numbers and a reviewer/private access mechanism if available.
3. Add accessions to the two metadata manifests and manuscript.
4. Push this repository to GitHub as a private repository.
5. Run `bash workflow/run_source_data_checks.sh`; a public release must exit successfully with no unresolved placeholder.
6. Create an immutable Git tag, preferably `v1.0.0`.
7. Enable the repository in Zenodo and create the archived release DOI.
8. Insert the version DOI and repository URL into `CITATION.cff`, README and manuscript Code Availability statement.
9. Re-run checksums and retain the release manifest with the submission archive.

## Licence decisions

This release uses the MIT License for code (`LICENSE`) and CC BY 4.0 for author-owned
processed data, metadata and numerical source data (`LICENSE_DATA.md`). Third-party
reference sequences and annotations are cited by accession/DOI and are not
redistributed unless their licences permit it.

## Draft manuscript statements

### Data Availability

> The quality-filtered whole-genome resequencing reads generated for the discovery and external validation cohorts have been deposited in the [NCBI Sequence Read Archive/Genome Sequence Archive] under BioProject accession **AUTHOR_INPUT_REQUIRED**. Sample-level BioSample and Run accessions are provided in Supplementary Table **AUTHOR_INPUT_REQUIRED**. The discovery reference genome and annotation used for coordinate-based analyses are available under accession **AUTHOR_INPUT_REQUIRED**; the external cohort was aligned to NCBI assembly GCA_029853015.1. Filtered variant datasets, frozen SNP and k-mer marker panels, numerical source data underlying all figures, and supplementary machine-readable tables are available at Zenodo under DOI **AUTHOR_INPUT_REQUIRED**.

### Code Availability

> All custom scripts, workflow configuration files and software-environment specifications required to reproduce the analyses and figures are available at **AUTHOR_INPUT_REQUIRED_GITHUB_URL** and archived as version **v1.0.0** at Zenodo under DOI **AUTHOR_INPUT_REQUIRED**. The archived release contains a file-level SHA-256 manifest and a mapping from each manuscript figure and table to its generating script and source data.
