# Reference and coordinate provenance

## Discovery reference

- Assembly description: 624.71 Mb, 24 linkage groups.
- Local FASTA contigs: 293.
- MD5: `479576ef44db39a532ce7b598cf01021`.
- SHA-256: `5cab70fe6ead189149a7a4d37e5cd78abde0e956ec7188d992abb40be26d392b`.
- Public assembly accession: `GCA_029853015.1` (ASM2985301v1); GenBank WGS master
  accession `JARXKY000000000.1`.
- FASTA source: NCBI Assembly/Datasets record for `GCA_029853015.1`; the release
  must record the date downloaded and the supplied FASTA checksum.
- GFF3 source: Figshare, DOI `10.6084/m9.figshare.22239901.v1`; the release must
  retain the exact downloaded-file checksum.
- Candidate interval: LG11:4,800,000–5,500,000 (coordinate convention must be stated per file format).

## External-cohort reference

- Assembly accession: GCA_029853015.1 (ASM2985301v1).
- Historical filename: `GCA_029853015.1_ASM2985301v1_genomic.fna.gz`.
- Candidate contig and interval: CM056831.1:4,800,000–5,500,000.
- Exact downloaded FASTA checksum: `AUTHOR_INPUT_REQUIRED`.

## Cross-build mapping

The discovery reference's LG11 and the external reference's CM056831.1 were
checked as sequence-identical for the 28,268,041-bp chromosome; their shared
SHA-256 is `4c1429cfae16e64c8dd2ad253dcd81548c65de9a3acd4c2bf271c0db31841908`.
No lift-over is required for the candidate interval. The release retains the
machine-readable chromosome crosswalk in `metadata/chromosome_crosswalk.tsv`.
