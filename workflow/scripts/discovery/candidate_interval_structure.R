#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript candidate_interval_structure.R <multisample.vcf[.gz]> <sample.xls> <output_dir>")
vcf_file <- args[1]; sample_file <- args[2]; outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(ape))

region_chr <- "LG11"; region_start <- 4800000L; region_end <- 5500000L
max_missing <- 0.10; min_maf <- 0.05; ld_window_bp <- 50000L; ld_r2 <- 0.80
set.seed(20250826)

open_text <- function(path) if (grepl("\\.gz$", path, ignore.case=TRUE)) gzfile(path, "rt") else file(path, "rt")
con <- open_text(vcf_file); on.exit(close(con), add=TRUE)
samples <- NULL; records <- list(); meta_n <- 0L
repeat {
  z <- readLines(con, n=10000, warn=FALSE)
  if (!length(z)) break
  for (line in z) {
    if (startsWith(line, "##")) { meta_n <- meta_n + 1L; next }
    if (startsWith(line, "#CHROM")) { samples <- strsplit(line, "\t", fixed=TRUE)[[1]][10:length(strsplit(line, "\t", fixed=TRUE)[[1]])]; next }
    if (startsWith(line, "#")) next
    f <- strsplit(line, "\t", fixed=TRUE)[[1]]
    if (length(f) < 10 || f[1] != region_chr) next
    pos <- as.integer(f[2]); if (is.na(pos) || pos < region_start || pos > region_end) next
    if (grepl(",", f[5], fixed=TRUE) || nchar(f[4]) != 1L || nchar(f[5]) != 1L) next
    fmt <- strsplit(f[9], ":", fixed=TRUE)[[1]]; gt_i <- match("GT", fmt)
    if (is.na(gt_i)) next
    dosage <- vapply(f[10:length(f)], function(x) {
      q <- strsplit(x, ":", fixed=TRUE)[[1]]; g <- if (length(q) >= gt_i) q[gt_i] else "."
      if (g %in% c("0/0","0|0")) 0 else if (g %in% c("0/1","1/0","0|1","1|0")) 1 else if (g %in% c("1/1","1|1")) 2 else NA_real_
    }, numeric(1))
    records[[length(records)+1L]] <- list(chr=f[1], pos=pos, ref=f[4], alt=f[5], dosage=dosage)
  }
}
if (is.null(samples)) stop("VCF缺少#CHROM表头")
if (length(records) == 0L) stop("VCF中未找到LG11:4,800,000-5,500,000的双等位SNP")
G <- do.call(rbind, lapply(records, `[[`, "dosage")); colnames(G) <- samples
info <- data.frame(CHROM=region_chr, POS=vapply(records, `[[`, integer(1), "pos"), REF=vapply(records, `[[`, character(1), "ref"), ALT=vapply(records, `[[`, character(1), "alt"), stringsAsFactors=FALSE)

pheno <- read.delim(sample_file, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("Sample","Group") %in% names(pheno))) stop("样本表必须包含Sample和Group列")
if (!all(samples %in% pheno$Sample)) stop("VCF样本在分组表中不完整")
pheno <- pheno[match(samples, pheno$Sample), ]; sex <- factor(pheno$Group, levels=c("F","M"))
if (any(is.na(sex))) stop("Group列只能使用F和M")

miss <- rowMeans(is.na(G)); af <- rowMeans(G, na.rm=TRUE)/2; maf <- pmin(af,1-af)
keep <- is.finite(maf) & miss <= max_missing & maf >= min_maf
qc <- data.frame(metric=c("Raw biallelic SNPs in interval","SNPs after missingness and MAF filters","Samples"), value=c(nrow(G),sum(keep),ncol(G)))
write.csv(qc, file.path(outdir,"candidate_interval_qc.csv"), row.names=FALSE)
G <- G[keep,,drop=FALSE]; info <- info[keep,,drop=FALSE]
if (nrow(G) < 2L) stop("过滤后候选区间SNP不足")

impute_rows <- function(x) { for(i in seq_len(nrow(x))) x[i,is.na(x[i,])] <- mean(x[i,],na.rm=TRUE); x }
Gi <- impute_rows(G)

# Greedy local LD thinning, retaining the earlier marker when r2 exceeds threshold.
ord <- order(info$POS); G <- G[ord,,drop=FALSE]; Gi <- Gi[ord,,drop=FALSE]; info <- info[ord,,drop=FALSE]
ld_keep <- logical(nrow(Gi)); retained <- integer(0)
for(i in seq_len(nrow(Gi))) {
  near <- retained[info$POS[i] - info$POS[retained] <= ld_window_bp]
  reject <- FALSE
  if(length(near)) for(j in near) {
    r <- suppressWarnings(cor(Gi[i,],Gi[j,])); if(is.finite(r) && r*r > ld_r2) { reject <- TRUE; break }
  }
  if(!reject) { ld_keep[i] <- TRUE; retained <- c(retained,i) }
}
Gp <- Gi[ld_keep,,drop=FALSE]; if(nrow(Gp)<2) stop("LD简化后标记不足")

# PCA on sample-by-marker standardized dosages.
X <- t(Gp); nz <- apply(X,2,sd)>0; pc <- prcomp(X[,nz,drop=FALSE], center=TRUE, scale.=TRUE)
pve <- 100*pc$sdev^2/sum(pc$sdev^2)
pcs <- data.frame(Sample=samples, Group=as.character(sex), PC1=pc$x[,1], PC2=pc$x[,2], PC3=if(ncol(pc$x)>=3) pc$x[,3] else NA_real_)
write.csv(pcs,file.path(outdir,"candidate_interval_PCA_scores.csv"),row.names=FALSE)
cols <- c(F="#D95F8D",M="#2878B5")
plot_pca <- function(path, device=c("png","pdf")) { device<-match.arg(device); if(device=="png") png(path,1800,1500,res=220) else pdf(path,7.2,6); plot(pcs$PC1,pcs$PC2,pch=21,bg=cols[sex],col="white",cex=1.5,xlab=sprintf("PC1 (%.2f%%)",pve[1]),ylab=sprintf("PC2 (%.2f%%)",pve[2]),main="LG11:4.80–5.50 Mb PCA"); legend("topright",levels(sex),pt.bg=cols,pch=21,bty="n"); dev.off() }
plot_pca(file.path(outdir,"candidate_interval_PCA.png"),"png"); plot_pca(file.path(outdir,"candidate_interval_PCA.pdf"),"pdf")

# IBS-like dosage distance and NJ tree.
D <- as.matrix(dist(X, method="manhattan"))/ncol(X)/2; rownames(D)<-colnames(D)<-samples
tree <- nj(as.dist(D)); write.tree(tree,file.path(outdir,"candidate_interval_NJ_tree.nwk"))
plot_tree <- function(path,device=c("png","pdf")) { device<-match.arg(device); if(device=="png") png(path,2000,2000,res=220) else pdf(path,8,8); plot(tree,tip.color=cols[sex[match(tree$tip.label,samples)]],cex=.75,no.margin=TRUE); legend("topleft",levels(sex),col=cols,lwd=2,bty="n"); dev.off() }
plot_tree(file.path(outdir,"candidate_interval_NJ_tree.png"),"png"); plot_tree(file.path(outdir,"candidate_interval_NJ_tree.pdf"),"pdf")

# Sex clustering statistics with label permutations.
diag(D)<-Inf; nn <- apply(D,1,which.min); nn_same <- mean(sex == sex[nn]); diag(D)<-0
within <- mean(D[outer(sex,sex,"==") & row(D)!=col(D)]); between <- mean(D[outer(sex,sex,"!=")]); ratio <- between/within
B<-10000L; perm_nn<-perm_ratio<-numeric(B)
for(b in seq_len(B)){ s<-sample(sex); diag(D)<-Inf; perm_nn[b]<-mean(s==s[apply(D,1,which.min)]); diag(D)<-0; perm_ratio[b]<-mean(D[outer(s,s,"!=")])/mean(D[outer(s,s,"==") & row(D)!=col(D)]) }
stats <- data.frame(metric=c("Interval SNPs after QC","LD-thinned interval SNPs","Nearest-neighbour same-sex proportion","Permutation P (nearest-neighbour)","Mean within-sex distance","Mean between-sex distance","Between/within distance ratio","Permutation P (distance ratio)"),value=c(nrow(G),nrow(Gp),nn_same,(1+sum(perm_nn>=nn_same))/(B+1),within,between,ratio,(1+sum(perm_ratio>=ratio))/(B+1)))
write.csv(stats,file.path(outdir,"candidate_interval_sex_clustering_statistics.csv"),row.names=FALSE)

# Heatmap of top absolute allele-frequency differences.
f_idx<-which(sex=="F"); m_idx<-which(sex=="M"); delta<-abs(rowMeans(G[,f_idx,drop=FALSE],na.rm=TRUE)-rowMeans(G[,m_idx,drop=FALSE],na.rm=TRUE))/2
n_top<-min(200L,length(delta)); top<-order(delta,decreasing=TRUE)[seq_len(n_top)]
top_table<-cbind(info[top,],Delta_AF=delta[top],Missing_rate=rowMeans(is.na(G[top,,drop=FALSE]))); write.csv(top_table,file.path(outdir,"candidate_interval_top_markers.csv"),row.names=FALSE)
H<-G[top,,drop=FALSE]; rownames(H)<-paste0(info$CHROM[top],":",info$POS[top]); sample_ord<-order(sex,hclust(as.dist(D))$order)
draw_heatmap <- function(path,device=c("png","pdf")){device<-match.arg(device); if(device=="png") png(path,2600,2200,res=240) else pdf(path,11,9); par(mar=c(7,7,3,1)); image(seq_len(ncol(H)),seq_len(nrow(H)),t(H[,sample_ord,drop=FALSE][nrow(H):1,]),col=c("#2166AC","#F7F7F7","#B2182B"),breaks=c(-.5,.5,1.5,2.5),axes=FALSE,xlab="",ylab="",main="Top sex-differentiated genotypes in LG11 candidate interval"); axis(1,at=seq_len(ncol(H)),labels=samples[sample_ord],las=2,cex.axis=.45); abline(v=sum(sex[sample_ord]=="F")+.5,lwd=2); mtext("Females",side=3,at=mean(which(sex[sample_ord]=="F"))); mtext("Males",side=3,at=mean(which(sex[sample_ord]=="M"))); dev.off()}
draw_heatmap(file.path(outdir,"candidate_interval_genotype_heatmap.png"),"png"); draw_heatmap(file.path(outdir,"candidate_interval_genotype_heatmap.pdf"),"pdf")

cat("Completed candidate interval analysis with",nrow(G),"QC SNPs and",nrow(Gp),"LD-thinned SNPs.\n")
