# Public code and processed-data release guide

## Purpose and scope

This package supports reproducibility of the analysis of a sex-associated genomic interval on LG11 in *Chelidonichthys spinosus*. It contains publication-level source data, curated metadata, frozen marker panels, software environments, and the key discovery and independent-validation scripts.

It does **not** contain raw FASTQ files, BAM/CRAM files, private NCBI preload credentials, NCBI private upload locations, personal disk paths, or institutional server paths. Raw sequence data must be deposited separately in NCBI SRA under the BioProject accession cited in the manuscript.

## Repository structure

| Directory or file | Public purpose |
|---|---|
| `README.md` | study summary, evidence design and quick-start commands |
| `workflow/` | discovery, external-validation, supplementary and figure workflow scripts |
| `config/config.example.yaml` | public configuration template; make a private `config.yaml` locally |
| `environment/` | Conda environments for figure/source-data checking and full reanalysis |
| `metadata/` | curated sample metadata, frozen SNP/k-mer panels and software versions |
| `source_data/` | machine-readable data underlying manuscript figures |
| `supplementary_tables/` | machine-readable supporting tables |
| `docs/` | workflow narrative, reference provenance, reproducibility audit and file mapping |

## Pre-publication checks

1. Obtain final public or reviewer-accessible SRA/BioProject accession(s) for both cohorts.
2. Replace all `AUTHOR_INPUT_REQUIRED` placeholders identified in `docs/REPRODUCIBILITY_AUDIT.md`.
3. Confirm that the public external-validation metadata retain the final gonadal-dissection sexes (37 female and 23 male).
4. Run `bash workflow/run_source_data_checks.sh` from the repository root and archive the terminal output as a release log.
5. Search for private information before pushing:

```bash
rg -n --hidden --glob '!.git/**' \
  'aspera\.openssh|uploads/[^ ]+|/Users/|/Volumes/|/public/home/|password|token|secret' .
```

The search must return no active private credential or personal/local path. Do not treat a comment in a generic template as permission to publish a real credential.

## GitHub release

1. Create a new GitHub repository named `Cspinosus_LG11_sex_association` under `wangyibang1994`.
2. Upload this package as the initial repository content. Do not upload raw reads, FASTQ MD5 files containing local paths, NCBI key files, browser screenshots or unpublished correspondence.
3. Ensure the repository has the included `MIT` code licence and `CC BY 4.0` data licence statement.
4. Create a tagged GitHub Release (recommended tag: `v1.0.0`) only after the manuscript files and source data are frozen.

## Zenodo archive

1. Connect Zenodo to the GitHub repository in the Zenodo GitHub integration settings.
2. Publish the GitHub Release; Zenodo will mint a version DOI and a concept DOI.
3. Use the concept DOI in the manuscript where a persistent code repository DOI is requested, and cite the specific version DOI for the submitted release.
4. Add both DOI(s) to `CITATION.cff`, the manuscript Data Availability statement, and the final repository README.

## Data availability mapping

| Evidence type | Deposit location |
|---|---|
| Raw discovery and validation WGS reads | NCBI SRA, linked to the study BioProject |
| Curated metadata, frozen marker panels, source data and supplementary tables | GitHub release archived by Zenodo |
| Scripts, environments and configuration template | GitHub release archived by Zenodo |
| Reference genome | Original public NCBI assembly cited in `docs/REFERENCE_PROVENANCE.md` |

Do not invent SRA accessions or DOI values in the manuscript. Insert only identifiers issued by the relevant repository.
