# Reproducible computational workflow

## Study design

The analysis combines four evidence layers: a balanced discovery cohort (30 females and 30 males), interval-restricted genetic structure, a reference-free k-mer screen, and an external cohort of 60 fish whose sex was established by direct gonadal dissection (37 females and 23 males). The discovery scan nominated LG11:4.80–5.50 Mb. Discovery-derived SNP and k-mer panels were frozen before being applied to the external cohort.

## Discovery read mapping and variant discovery

Paired-end reads were aligned to the 624.71-Mb, 24-linkage-group discovery reference with BWA-MEM v0.7.18. Alignments were coordinate-sorted and indexed with Samtools; `quickcheck` and `flagstat` were used for integrity and mapping summaries. The historical discovery workflow did not mark duplicates. The exact discovery-reference accession and upstream read-trimming parameters remain `AUTHOR_INPUT_REQUIRED`.

Variants were called with BCFtools/HTSlib v1.24. `bcftools mpileup` used minimum mapping and base qualities of 20 and emitted depth and allelic-depth annotations. Multiallelic calling used `bcftools call -m -v`, followed by normalisation and restriction to biallelic SNPs. Sites were retained at QUAL >=30, total depth 300–3,000, missingness <=10% and MAF >=0.05. This yielded 3,310 SNPs in LG11:4.80–5.50 Mb. The historical discovery workflow did not mask genotypes by DP, GQ or allelic balance; a stricter sensitivity analysis is therefore required before claiming that such filtering was applied.

## Genome-wide structure and differentiation scan

The genome-wide dataset comprised 3,891,433 filtered SNPs. LD pruning with PLINK (`--indep-pairwise 50 10 0.1`) retained 1,383,748 SNPs for PCA, neighbour-joining analysis and ADMIXTURE K=1–10. Exact historical PLINK/ADMIXTURE versions and remaining options must be recovered before release.

At each SNP we calculated absolute allele-frequency difference, `|Delta AF|=|pF-pM|`; Euclidean distance, `ED=sqrt(2)|Delta AF|`; and Hudson-type FST components. Statistics were summarised in 200-kb windows advanced by 50 kb. Window-level FST was the sum of Hudson numerators divided by the sum of denominators. The joint profiles of weighted FST, mean ED and the number of SNPs with `|Delta AF|>=0.3` nominated LG11:4.80–5.50 Mb as an operational candidate interval.

Genome-wide significance was calibrated by 500 label permutations. Each permutation retained 30 labels per group, recomputed the complete genome scan and recorded the maximum value of each statistic. Seeds were 1–500 and empirical family-wise P values were `(1 + number of permuted maxima >= observed)/(501)`.

## Candidate-interval structure

The 3,310 interval SNPs were analysed after mean imputation for PCA and distance calculations. Greedy local pruning removed the later SNP of any pair with r-squared >0.8 within 50 kb, retaining 2,533 SNPs. PCA used standardised genotype dosages. An IBS-like distance was calculated as the Manhattan distance between dosage vectors divided by twice the number of non-missing markers; a neighbour-joining tree was generated with `ape::nj`. Sex separation was evaluated with 10,000 label permutations (seed 20250826). The publication heatmap uses the 120 SNPs with the largest `|Delta AF|`.

## Functional annotation, heterozygosity and depth

High-`|Delta AF|` SNPs were intersected with a version-locked GFF3 annotation and classified by genomic context. The provisional custom coding-consequence script requires replacement or verification with a phase-aware tool such as VEP or SnpEff before coding changes are interpreted.

Individual heterozygosity was calculated across all 3,310 interval SNPs and the top 50, 120 and 200 differentiated markers. Sex contrasts used Welch tests and 10,000 label permutations (seed 20260831). Depth across LG11:4.20–6.10 Mb was obtained using `samtools depth -aa`, summarised in 10-kb bins and normalised to median flanking depth. A sensitivity analysis excluded 5.20–5.22 Mb.

## Reference-free k-mer analysis

Canonical 31-mers were counted with KMC v3.2.4. Sex-pooled databases used a minimum count of five. Male-only 31-mers with pooled male count >=100 formed 95,859 high-support candidates. Individual databases used a minimum count of two, and a marker was called present in an individual at count >=5. Eighty-three markers occurring in at least 29/30 males and 0/30 females were frozen as the core male-associated panel.

Core markers were aligned to the discovery reference with BWA-MEM v0.7.18 (`-a -k 15`). Only exact 31M, NM=0 matches were retained. Eighty-one of 83 markers mapped exactly, producing 101 matches concentrated at approximately LG11:5,168,315–5,336,051.

## External validation cohort

External reads were mapped independently to NCBI assembly GCA_029853015.1 with BWA-MEM v0.7.18. Duplicates were marked with GATK MarkDuplicates. Per-sample GVCFs were generated with GATK v4.2.0.0 HaplotypeCaller and jointly genotyped through GenomicsDBImport and GenotypeGVCFs. SNPs failing QD 2.0, FS 60.0, SOR 3.0, MQ 40.0, MQRankSum -12.5 or ReadPosRankSum -8.0 criteria were removed. Filtered genotypes with DP<4 or GQ<20 were converted to no-calls; retained sites had <=10% missingness and MAF>=0.01. The external interval was CM056831.1:4.80–5.50 Mb.

The discovery top-50 `|Delta AF|` SNPs were frozen before validation. Twenty-four exact coordinate-and-allele matches were represented externally. Dosages were mean-imputed and standardised. Unsupervised K-means used the first five PCs (K=2, nstart=100, seed 20260828), with Ward.D2 clustering as a corroborating analysis. Cluster labels were assigned only after clustering by comparison with dissected sex. Accuracy, sex-specific sensitivity/specificity, balanced accuracy, Matthews correlation coefficient and Wilson 95% confidence intervals should be reported.

For k-mer validation, the frozen 83-marker set was counted without reading external sex labels. A marker was detected at an individual count >=3. The largest observed score gap was between 7 and 68 detected markers, giving a midpoint of 37.5. Scores >37.5 were high-confidence positive, scores of zero were negative, and intermediate scores were inconclusive. BY1 and BY5 each scored seven and were retained as inconclusive. Binary concordance uses the 58 conclusive samples; coverage must additionally be reported as 58/60.

## Reproducibility boundary

The source-data route can be released after scripts are path-parameterised. Full raw-data reproducibility remains blocked until the discovery accession, read-processing parameters, sample collection/ethics metadata, GFF3 provenance, reference coordinate crosswalk and exact historical structure-analysis software options are supplied. These omissions must not be replaced with assumed values.
