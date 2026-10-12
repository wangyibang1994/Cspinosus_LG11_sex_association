args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Usage: 11_haplotype_tree.R matrix.tsv output_prefix')
suppressPackageStartupMessages(library(ape))
d <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
x <- as.matrix(d[,7:ncol(d)])
storage.mode(x) <- 'numeric'
rownames(x) <- d$tip_id
# Normalised Hamming distance across discovery-oriented high-delta markers.
if (anyNA(x)) stop('Missing phased alleles in haplotype matrix')
dm <- dist(x, method='manhattan') / ncol(x)
tr <- nj(dm)
write.tree(tr, paste0(args[2], '.nj.newick'))
write.table(d[,1:6], paste0(args[2], '.tip_metadata.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
write.table(as.matrix(dm), paste0(args[2], '.distance_matrix.tsv'), sep='\t', quote=FALSE)
