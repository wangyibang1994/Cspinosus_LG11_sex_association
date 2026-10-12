#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(svglite)
  library(ragg)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
res <- file.path(root, "results")
legacy_res <- Sys.getenv("REFERENCE_CONTEXT_RESULTS_DIR", unset = "")
if (!nzchar(legacy_res) || !dir.exists(legacy_res)) {
  stop("Set REFERENCE_CONTEXT_RESULTS_DIR to the directory containing reference-context results")
}
out <- file.path(root, "figures", "Supplementary_reference_context")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

blue <- "#397AAE"
red <- "#D65F5F"
gold <- "#D69B2D"
dark <- "#26343B"

theme_paper <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      axis.title = element_text(size = base_size + 0.5, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      plot.title = element_text(size = base_size + 1.2, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4B555A"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.3),
      plot.margin = margin(5, 6, 5, 6)
    )
}

# a, independently phased haplotype NJ tree -------------------------------
hap <- read_tsv(file.path(res, "01_independent_haplotypes", "independent_haplotype_matrix.tsv"), show_col_types = FALSE)
meta <- hap %>%
  select(tip_id, sample_id, cohort, gonadal_sex, haplotype, haplotype_class) %>%
  mutate(haplotype_class = factor(haplotype_class, levels = c("other", "male-associated")),
         cohort = factor(cohort, levels = c("discovery", "external")))
marker_matrix <- hap %>% select(-tip_id, -sample_id, -cohort, -gonadal_sex, -haplotype, -haplotype_class) %>% as.matrix()
rownames(marker_matrix) <- hap$tip_id
tree <- nj(dist(marker_matrix, method = "manhattan") / ncol(marker_matrix))
write_tsv(meta, file.path(out, "SourceData_SuppFigure_a_haplotype_tree_metadata.tsv"))

p_a <- ggtree(tree, layout = "circular", size = 0.22, colour = "#6A7479") %<+% meta +
  geom_tippoint(aes(fill = haplotype_class, shape = cohort), size = 1.45, stroke = 0.25, colour = "white") +
  scale_fill_manual(values = c(other = "#B8C1C6", `male-associated` = gold),
                    guide = "none") +
  scale_shape_manual(values = c(discovery = 21, external = 24), guide = "none") +
  labs(title = "Independent-phase NJ tree",
       subtitle = "219 discovery-frozen markers\nGold: male-associated; grey: other") +
  theme_void(base_family = "sans", base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9.2), plot.subtitle = element_text(size = 7.8, colour = "#4B555A"),
        legend.position = "none", plot.margin = margin(5, 6, 5, 6))

# b, reference-individual SNP projection ----------------------------------
cent <- read_tsv(file.path(legacy_res, "04_reference_retrosex", "reference_juvenile.top50_centroid_markers.tsv"), show_col_types = FALSE) %>%
  arrange(POS) %>% mutate(marker_rank = row_number(), marker = paste0(CHROM, ":", POS))
heat <- cent %>%
  select(marker_rank, marker, female_mean_dosage, reference_alt_dosage, male_mean_dosage, closer_sex) %>%
  pivot_longer(c(female_mean_dosage, reference_alt_dosage, male_mean_dosage), names_to = "profile", values_to = "alt_dosage") %>%
  mutate(profile = factor(profile,
                          levels = c("female_mean_dosage", "reference_alt_dosage", "male_mean_dosage"),
                          labels = c("Female centroid", "Reference juvenile", "Male centroid")))
write_tsv(heat, file.path(out, "SourceData_SuppFigure_b_reference_SNP_projection.tsv"))

p_b <- ggplot(heat, aes(profile, marker_rank, fill = alt_dosage)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#2171B5"), limits = c(0, 2), name = "ALT dosage") +
  scale_y_reverse(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0, 0)) +
  labs(title = "Reference-juvenile SNP projection", subtitle = "49/50 markers were closer to the male centroid",
       x = NULL, y = "Marker rank along LG11") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "bottom")

# c, raw core-k-mer counts --------------------------------------------------
kmer_panel <- read_csv(file.path(legacy_res, "04_reference_retrosex", "frozen_core83_kmers.csv"), show_col_types = FALSE)
kmer_hits <- read_tsv(file.path(legacy_res, "04_reference_retrosex", "reference_core83_hits.tsv"),
                      col_names = c("kmer", "read_count"), show_col_types = FALSE) %>%
  left_join(select(kmer_panel, kmer_id, kmer), by = "kmer") %>% arrange(read_count) %>% mutate(rank = row_number())
write_tsv(kmer_hits, file.path(out, "SourceData_SuppFigure_c_reference_kmer_counts.tsv"))

p_c <- ggplot(kmer_hits, aes(rank, read_count)) +
  geom_col(fill = blue, width = 0.82) +
  annotate("text", x = 4, y = max(kmer_hits$read_count) * 0.96, label = "83/83 detected", hjust = 0,
           size = 3.0, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), breaks = c(1, 20, 40, 60, 83)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Frozen core 31-mers in reference reads", subtitle = "Counts from the original reference-individual libraries",
       x = "Core markers ordered by count", y = "31-mer read count") +
  theme_paper()

# d, original repeat-annotation context ------------------------------------
repwin <- read_tsv(file.path(legacy_res, "03_comparative", "original_Csp_repeat_annotation",
                             "Csp_original_repeat_annotation.LG11_100kb_windows.tsv"), show_col_types = FALSE) %>%
  mutate(mid_mb = (start_1based + end_1based) / 2e6) %>% filter(mid_mb >= 4.0, mid_mb <= 6.3)
write_tsv(repwin, file.path(out, "SourceData_SuppFigure_d_repeat_windows.tsv"))

p_d <- ggplot(repwin, aes(mid_mb, coverage_fraction)) +
  annotate("rect", xmin = 4.80, xmax = 5.50, ymin = -Inf, ymax = Inf, fill = "#F4E4B7", alpha = 0.45) +
  annotate("rect", xmin = 5.168, xmax = 5.336, ymin = -Inf, ymax = Inf, fill = gold, alpha = 0.28) +
  geom_line(linewidth = 0.65, colour = dark) +
  geom_point(size = 1.7, colour = dark) +
  scale_x_continuous(limits = c(4.0, 6.3), breaks = seq(4.0, 6.0, 0.5)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Original LG11 repeat-annotation context",
       subtitle = "Gold: candidate interval; darker band: core-marker cluster",
       x = "LG11 position (Mb)", y = "Repeat-union coverage per 100 kb") +
  theme_paper()

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11)))

width_in <- 183 / 25.4
height_in <- 154 / 25.4
grDevices::pdf(file.path(out, "FigureS7_reference_and_haplotype_context.pdf"), width = width_in, height = height_in,
               family = "Helvetica", useDingbats = FALSE)
print(fig)
dev.off()
svglite::svglite(file.path(out, "FigureS7_reference_and_haplotype_context.svg"), width = width_in, height = height_in,
                 bg = "white")
print(fig)
dev.off()
agg_tiff(file.path(out, "FigureS7_reference_and_haplotype_context_600dpi.tiff"), width = width_in, height = height_in,
         units = "in", res = 600, compression = "lzw", background = "white")
print(fig)
dev.off()
agg_png(file.path(out, "FigureS7_reference_and_haplotype_context_preview.png"), width = width_in, height = height_in,
        units = "in", res = 300, background = "white")
print(fig)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
message("Candidate supplementary figure written to: ", out)
