# Sequencing and processing provenance recovered from the project report

## Discovery cohort

The supplied provider report documents 60 discovery individuals, comprising 30
female-labelled and 30 male-labelled samples. The final public metadata must retain
the author-confirmed, gonadal-dissection sex labels.

## Library preparation and sequencing

The provider report states that DNA was sonically fragmented to approximately
300–350 bp, followed by end repair, A-tailing, adapter ligation, purification and
PCR amplification, before paired-end sequencing. A read-only inspection of the
analysis-ready FASTQ input on the analysis server confirmed 150-bp reads and paired
files (`*_1.fq.gz` / `*_2.fq.gz`).

The available report does **not** specify the instrument model. The FASTQ header
prefix must not be used to infer an instrument model. Obtain this field from the
sequencing provider's order/report before creating the SRA metadata.

## Quality control and alignment summary

| Metric | Value |
|---|---:|
| Samples | 60 |
| Total clean bases | 663.24 Gb |
| Mean clean bases per sample | 11.05 Gb |
| Range of clean bases per sample | 6.25–24.33 Gb |
| Mean mapping rate | 98.32% |
| Range of mapping rate | 94.39–99.42% |
| Mean mapped depth | 16.87× |
| Range of mapped depth | 9.45–36.99× |
| Mean genome coverage | 97.76% |

The provider report records read filtering with fastp using a 50-bp minimum retained
read length and four threads. Its stated filtering logic trims low-quality read ends,
removes adapters, excludes reads with more than 40% bases below Q15 or more than five
ambiguous bases, and excludes reads shorter than 50 bp after processing.

## Provider software record

The report lists fastp (version described as `latest`; command includes `-l 50 -w 4`),
BWA 0.7.12-r1039 (`mem -M`), Picard 1.124, GATK 4.4.0.0, VCFtools 0.1.16,
ANNOVAR 2016-02-01, VCF2Dis 1.47, VCF2PCACluster 1.40 and ADMIXTURE 1.3.0.

These provider-report versions must be distinguished from any later scripts used for
the manuscript's reanalysis, frozen panels and MBE figure generation. Do not replace
reanalysis versions with the provider versions without an explicit provenance check.

## Source records retained by the authors

- `60个小眼绿鳍鱼的分离性状关联分析结题报告/report/index.html`
- `result/01.Data_QC/FilterSummary.xls`
- `result/01.Data_QC/FilterSummary2.xls`
- `result/02.Reads_Mapping/Mapping.stat.xls`
- `result/02.Reads_Mapping/Mapping.stat2.xls`

The public code release should retain this concise derived provenance record and
machine-readable QC summaries, rather than publish the provider's complete HTML
report unless the authors hold redistribution permission.
