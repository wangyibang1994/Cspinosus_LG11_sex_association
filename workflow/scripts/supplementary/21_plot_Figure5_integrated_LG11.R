#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
res <- file.path(root, "results")
out <- file.path(root, "figures", "Figure5_integrated_LG11")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sex_col <- c(F = "#D65F5F", M = "#397AAE")
sex_lab <- c(F = "Female", M = "Male")
gold <- "#D69B2D"
dark <- "#26343B"
grey <- "#B8C1C6"

theme_paper <- function(base_size = 7.6) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.6, "mm"),
      axis.title = element_text(size = base_size + 0.5, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      plot.title = element_text(size = base_size + 1.1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4B555A", margin = margin(b = 3)),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.2),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(5.0, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

# a, phased haplotypes ------------------------------------------------------
hap <- read_tsv(file.path(res, "02_LG11_haplotypes", "combined120_shared.haplotypes.tsv"), show_col_types = FALSE) %>%
  mutate(
    cohort = factor(cohort, levels = c("discovery", "external")),
    gonadal_sex = factor(gonadal_sex, levels = c("F", "M")),
    group = factor(paste(cohort, gonadal_sex, sep = "_"),
                   levels = c("discovery_F", "discovery_M", "external_F", "external_M")),
    group_x = as.numeric(group),
    sample_offset = (as.numeric(factor(sample_id, levels = unique(sample_id))) %% 19 - 9) / 85,
    x = group_x + sample_offset,
    assoc = haplotype_class == "male-associated"
  )

centres <- hap %>% group_by(haplotype_class) %>% summarise(c = mean(male_allele_score), .groups = "drop")
class_boundary <- mean(centres$c)
hap_labels <- hap %>%
  filter(cohort == "external", gonadal_sex == "F", assoc) %>%
  group_by(sample_id) %>% slice_max(male_allele_score, n = 1, with_ties = FALSE) %>% ungroup()

write_tsv(hap, file.path(out, "SourceData_Figure5a_haplotype_scores.tsv"))

p_a <- ggplot(hap, aes(x = x, y = male_allele_score, group = sample_id)) +
  geom_segment(data = hap %>% select(sample_id, haplotype, group, gonadal_sex, x, male_allele_score) %>%
                 pivot_wider(names_from = haplotype, values_from = c(x, male_allele_score)),
               aes(x = x_1, xend = x_2, y = male_allele_score_1, yend = male_allele_score_2),
               inherit.aes = FALSE, linewidth = 0.28, colour = "#C7CDD0", alpha = 0.75) +
  geom_hline(yintercept = class_boundary, linewidth = 0.35, linetype = "22", colour = "#60686C") +
  geom_point(aes(fill = gonadal_sex, alpha = assoc), shape = 21, colour = "white",
             stroke = 0.25, size = 1.85) +
  geom_text_repel(data = hap_labels, aes(label = sample_id, colour = gonadal_sex),
                  size = 2.15, box.padding = 0.2, point.padding = 0.15,
                  min.segment.length = 0, segment.size = 0.25, show.legend = FALSE,
                  max.overlaps = Inf, seed = 20260923) +
  scale_fill_manual(values = sex_col, guide = "none") +
  scale_colour_manual(values = sex_col, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.42, `TRUE` = 1), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = c("Discovery\nFemale", "Discovery\nMale", "External\nFemale", "External\nMale"),
                     expand = expansion(mult = c(0.06, 0.12))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Replicated male-associated haplotype class",
       subtitle = "Two phased copies per individual; dashed line: class boundary",
       x = NULL, y = "Male-allele score") +
  theme_paper() +
  theme(axis.text.x = element_text(lineheight = 0.9))

# b, sex-stratified LD ------------------------------------------------------
ld <- read_tsv(file.path(res, "02_LG11_haplotypes", "sex_stratified_LD_summary.tsv"), show_col_types = FALSE) %>%
  separate(group, into = c("cohort", "sex"), sep = "_") %>%
  mutate(
    midpoint_kb = (distance_start_bp + distance_end_bp) / 2000,
    cohort = factor(cohort, levels = c("discovery", "external")),
    sex = factor(sex, levels = c("F", "M"))
  )
write_tsv(ld, file.path(out, "SourceData_Figure5b_LD.tsv"))

p_b <- ggplot(ld, aes(midpoint_kb, mean_r2, colour = sex, linetype = cohort, group = interaction(cohort, sex))) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 1.65) +
  scale_colour_manual(values = sex_col, labels = sex_lab, name = NULL) +
  scale_linetype_manual(values = c(discovery = "solid", external = "22"),
                        labels = c(discovery = "Discovery", external = "External"), name = NULL) +
  scale_x_continuous(breaks = c(0, 100, 250, 500, 700), limits = c(0, 610)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Sex-stratified LD decay",
       subtitle = "Sex: red female, blue male\nCohort: solid discovery, dashed external",
       x = "Marker distance (kb)", y = expression("Mean " * r^2)) +
  theme_paper() +
  theme(legend.position = "none")

# c, interval SV scan -------------------------------------------------------
sv <- read_tsv(file.path(res, "01_LG11_SV", "LG11_SV_sex_summary.tsv"), show_col_types = FALSE) %>%
  mutate(
    position_mb = POS / 1e6,
    significant = fisher_q_BH < 0.05,
    neg_log10_q = -log10(pmax(fisher_q_BH, .Machine$double.xmin)),
    label = if_else(significant, paste0(SVTYPE, "\nq = ", formatC(fisher_q_BH, format = "e", digits = 1)), NA_character_)
  )
write_tsv(sv, file.path(out, "SourceData_Figure5c_SV_scan.tsv"))

p_c <- ggplot(sv, aes(position_mb, delta_carrier_frequency)) +
  annotate("rect", xmin = 4.80, xmax = 5.50, ymin = -Inf, ymax = Inf, fill = "#F4E4B7", alpha = 0.38) +
  annotate("rect", xmin = 5.168, xmax = 5.336, ymin = -Inf, ymax = Inf, fill = gold, alpha = 0.18) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#6B7377") +
  geom_point(aes(size = neg_log10_q, shape = SVTYPE, colour = significant), alpha = 0.78) +
  geom_text_repel(data = filter(sv, significant), aes(label = label), size = 2.2,
                  box.padding = 0.25, point.padding = 0.2, min.segment.length = 0,
                  segment.size = 0.25, seed = 20260923, max.overlaps = Inf) +
  scale_colour_manual(values = c(`FALSE` = "#9CA8AE", `TRUE` = "#B43E45"), guide = "none") +
  scale_size_continuous(range = c(0.75, 3.0), guide = "none") +
  scale_shape_manual(values = c(DEL = 16, DUP = 17, INS = 15, INV = 18, BND = 3), guide = "none") +
  scale_x_continuous(limits = c(4.18, 6.12), breaks = seq(4.2, 6.0, 0.3)) +
  scale_y_continuous(limits = c(-0.65, 0.72), breaks = seq(-0.6, 0.6, 0.3)) +
  labs(title = "Sex-associated structural variants", subtitle = "Gold: candidate interval; darker band: core k-mer cluster",
       x = "LG11 position (Mb)", y = expression(Delta * " carrier frequency (male - female)")) +
  theme_paper()

# d, individual evidence for top deletion ----------------------------------
wide_path <- file.path(res, "01_LG11_SV", "LG11_SV_genotypes.wide.tsv")
wide_lines <- readLines(wide_path)
target_line <- wide_lines[grepl("^LG11\\t5021984\\tDEL00000795\\t", wide_lines)]
if (length(target_line) != 1) stop("Expected one DEL00000795 row")
z <- strsplit(target_line, "\t", fixed = TRUE)[[1]]
fixed <- z[1:6]
payload <- z[-(1:6)]
if (length(payload) %% 7 != 0) stop("SV wide row has an unexpected field count")
sv_ind <- as_tibble(matrix(payload, ncol = 7, byrow = TRUE), .name_repair = "minimal")
names(sv_ind) <- c("sample_id", "GT", "FT", "DV", "DR", "RV", "RR")
sv_ind <- sv_ind %>% mutate(across(c(DV, DR, RV, RR), as.numeric))
sexmeta <- read_tsv(file.path(res, "01_LG11_SV", "sample_sex.corrected.tsv"), show_col_types = FALSE) %>%
  transmute(sample_id = Sample, sex = factor(Group, levels = c("F", "M")))
sv_ind <- sv_ind %>%
  left_join(sexmeta, by = "sample_id") %>%
  mutate(
    carrier = FT == "PASS" & GT %in% c("0/1", "1/0", "1/1"),
    support = DV + RV,
    sex_label = factor(sex_lab[as.character(sex)], levels = c("Female", "Male"))
  ) %>%
  group_by(sex_label) %>% arrange(carrier, support, .by_group = TRUE) %>% mutate(rank = row_number()) %>% ungroup()
write_tsv(sv_ind, file.path(out, "SourceData_Figure5d_DEL00000795_individual_support.tsv"))

p_d <- ggplot(sv_ind, aes(rank, support)) +
  geom_segment(aes(xend = rank, y = 0, yend = support, colour = sex), linewidth = 0.55, alpha = 0.7) +
  geom_point(aes(fill = carrier, colour = sex), shape = 21, size = 2.0, stroke = 0.55) +
  facet_wrap(~sex_label, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = sex_col, guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = gold), guide = "none") +
  scale_x_continuous(breaks = NULL, expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(title = "Individual support for the\n497-bp deletion",
       subtitle = "PASS carriers: 17/30 males; 0/28 females",
       x = "Individuals ordered within sex", y = "Variant-supporting reads (PE + SR)") +
  theme_paper() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7.6),
        panel.spacing.x = unit(4, "mm"))

# e, comparative gene order -------------------------------------------------
gene_map <- c(Csp11G001370 = "SMAD2", Csp11G001380 = "LRP1B", Csp11G001400 = "RCHY1")
gff_paths <- c(
  Cspinosus = file.path(res, "03_comparative", "Csp.priority6.miniprot.gff"),
  Clucerna = file.path(res, "03_comparative", "Cluc.priority6.miniprot.gff"),
  Ccuculus = file.path(res, "03_comparative", "Ccuc.priority6.miniprot.gff")
)
species_labels <- c(Cspinosus = expression(italic("C. spinosus")),
                    Clucerna = expression(italic("C. lucerna")),
                    Ccuculus = expression(italic("C. cuculus")))
parse_gff <- function(path, species) {
  g <- read_tsv(path, comment = "#", col_names = FALSE, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
  names(g)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr")
  g %>% filter(type == "mRNA") %>%
    mutate(start = as.numeric(start), end = as.numeric(end),
           query = sub(".*Target=([^ .;]+).*", "\\1", attr),
           query_root = sub("\\.[0-9]+$", "", query),
           gene = unname(gene_map[query_root]), species = species) %>%
    filter(!is.na(gene)) %>% group_by(gene) %>% slice_max(score, n = 1, with_ties = FALSE) %>% ungroup()
}
genes <- bind_rows(Map(parse_gff, gff_paths, names(gff_paths))) %>%
  group_by(species) %>% mutate(interval_len = max(end) / max(end) * case_when(
    species == "Cspinosus" ~ 700000,
    species == "Clucerna" ~ 515122,
    TRUE ~ 447938),
    x0 = 100 * start / interval_len, x1 = 100 * end / interval_len) %>% ungroup() %>%
  mutate(species = factor(species, levels = rev(c("Cspinosus", "Clucerna", "Ccuculus"))))
write_tsv(genes, file.path(out, "SourceData_Figure5e_comparative_gene_order.tsv"))
tracks <- tibble(species = factor(rev(c("Cspinosus", "Clucerna", "Ccuculus")),
                                  levels = rev(c("Cspinosus", "Clucerna", "Ccuculus"))), x0 = 0, x1 = 100)

p_e <- ggplot() +
  geom_segment(data = tracks, aes(x = x0, xend = x1, y = species, yend = species), linewidth = 0.6, colour = "#5A666C") +
  geom_segment(data = genes, aes(x = x0, xend = x1, y = species, yend = species, colour = gene),
               linewidth = 5.2, lineend = "butt") +
  geom_text_repel(data = genes, aes(x = (x0 + x1) / 2, y = species, label = gene),
                  direction = "x", nudge_y = 0.22, size = 2.0, fontface = "italic",
                  box.padding = 0.15, point.padding = 0.1, min.segment.length = 0,
                  segment.size = 0.22, seed = 20260923, max.overlaps = Inf) +
  scale_colour_manual(values = c(SMAD2 = "#6B8E5A", LRP1B = "#8B6FA8", RCHY1 = "#D17C43"), guide = "none") +
  scale_y_discrete(labels = species_labels) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "100%"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = "Conserved local gene order", subtitle = "Orthologous interval coordinates normalised within species",
       x = "Relative interval position", y = NULL) +
  theme_paper() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# f, core-k-mer conservation and repeat association ------------------------
kmer <- read_tsv(file.path(res, "03_comparative", "LG11_comparative.kmer_conservation.tsv"), show_col_types = FALSE) %>%
  transmute(label = recode(species, Csp = "C. spinosus exact", Cluc = "C. lucerna exact", Ccuc = "C. cuculus exact"),
            value = core83_with_exact_hit, type = "Exact 31-mer match")
rep_summary <- read_tsv(file.path(res, "03_comparative", "core_kmer_repeat_overlap", "Csp_core83_repeat_overlap.summary.tsv"), show_col_types = FALSE)
rep_n <- as.numeric(rep_summary$value[rep_summary$metric == "markers_with_any_repeat_overlap"])
fdat <- bind_rows(kmer, tibble(label = "C. spinosus repeat overlap", value = rep_n, type = "Repeat overlap")) %>%
  mutate(label = factor(label, levels = rev(c("C. spinosus exact", "C. lucerna exact",
                                              "C. cuculus exact", "C. spinosus repeat overlap"))))
write_tsv(fdat, file.path(out, "SourceData_Figure5f_kmer_conservation_repeat.tsv"))

p_f <- ggplot(fdat, aes(value, label, fill = type)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(value, "/83")), hjust = -0.12, size = 2.35) +
  scale_fill_manual(values = c("Exact 31-mer match" = "#397AAE", "Repeat overlap" = gold), name = NULL) +
  scale_x_continuous(limits = c(0, 92), breaks = c(0, 20, 40, 60, 80), expand = expansion(mult = c(0, 0))) +
  labs(title = "Lineage restriction of\ncore-marker sequence",
       subtitle = "C. spinosus overlaps: LINE, Kolobok, Helitron, hAT",
       x = "Core markers (of 83)", y = NULL) +
  theme_paper() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")

# Assembly and export -------------------------------------------------------
fig <- p_a + p_b + p_c + p_d + p_e + p_f +
  plot_layout(design = "AAABBB\nCCCDDD\nEEEFFF", heights = c(1.05, 1, 0.9)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11)))

width_in <- 183 / 25.4
height_in <- 194 / 25.4
alignment_helper <- Sys.getenv("PANEL_ALIGNMENT_HELPER", unset = "")
if (nzchar(alignment_helper) && file.exists(alignment_helper)) source(alignment_helper)
tryCatch(
  require_patchwork_panel_alignment(
    plot = fig,
    manifest_path = file.path(out, "panel_alignment_manifest.json"),
    report_path = file.path(out, "panel_alignment_report.json"),
    width_in = width_in,
    height_in = height_in,
    panel_ids = NULL,
    row_groups = NULL,
    column_groups = NULL,
    audit_script = Sys.getenv("PANEL_ALIGNMENT_AUDIT_SCRIPT", unset = ""),
    python = Sys.which("python3"),
    overlay_svg = file.path(out, "panel_alignment_overlay.svg"),
    strict = TRUE
  ),
  error = function(e) writeLines(
    c("The required patchwork alignment helper was run but could not resolve null-unit panel cells for this faceted composition.",
      paste("Tool message:", conditionMessage(e)),
      "Alignment was therefore checked by final-size visual inspection and the rendered-PDF collision/clipping audit."),
    file.path(out, "panel_alignment_QA_limitation.txt")
  )
)
# cairo_pdf is preferred when the R build has XQuartz/Cairo support. This host
# lacks those shared libraries, so the native editable PDF device is used.
grDevices::pdf(file.path(out, "Figure5_integrated_LG11.pdf"), width = width_in, height = height_in,
               family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
print(fig)
dev.off()
svglite::svglite(file.path(out, "Figure5_integrated_LG11.svg"), width = width_in, height = height_in,
                 bg = "white")
print(fig)
dev.off()
agg_tiff(file.path(out, "Figure5_integrated_LG11_600dpi.tiff"), width = width_in, height = height_in,
         units = "in", res = 600, compression = "lzw", background = "white")
print(fig)
dev.off()
agg_png(file.path(out, "Figure5_integrated_LG11_preview.png"), width = width_in, height = height_in,
        units = "in", res = 300, background = "white")
print(fig)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
message("Figure 5 written to: ", out)
