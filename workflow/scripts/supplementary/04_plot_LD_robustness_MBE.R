#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
res <- file.path(root, "results", "02_LD_robustness")
out <- file.path(root, "figures", "Supplementary_LD_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

red <- "#D65F5F"
blue <- "#397AAE"
grey <- "#8C969B"

theme_paper <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.38, colour = "black"),
      axis.ticks = element_line(linewidth = 0.32),
      axis.title = element_text(size = base_size + 0.5, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      plot.title = element_text(size = base_size + 1.2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.1, colour = "#4B555A", margin = margin(b = 4)),
      plot.margin = margin(5, 7, 5, 7),
      legend.position = "none"
    )
}

summary <- read_tsv(file.path(res, "LD_robustness_summary.tsv"), show_col_types = FALSE) %>%
  filter(distance_start_bp != "GLOBAL") %>%
  mutate(
    distance_start_bp = as.numeric(distance_start_bp),
    distance_end_bp = as.numeric(distance_end_bp),
    midpoint_kb = (distance_start_bp + distance_end_bp) / 2000
  )

permutation <- read_tsv(file.path(res, "LD_permutation_replicates.tsv"), show_col_types = FALSE)
observed_global <- mean(summary$male_mean_r2 - summary$female_full_mean_r2)
empirical_p <- (1 + sum(permutation$global_mean_difference >= observed_global)) / (nrow(permutation) + 1)

write_tsv(summary, file.path(out, "SourceData_SuppFigure_LD_a.tsv"))
write_tsv(permutation, file.path(out, "SourceData_SuppFigure_LD_b.tsv"))

p_a <- ggplot(summary, aes(midpoint_kb)) +
  geom_ribbon(aes(ymin = female_downsample_q025, ymax = female_downsample_q975),
              fill = red, alpha = 0.16) +
  geom_line(aes(y = female_downsample_median), colour = red, linewidth = 0.75, linetype = "22") +
  geom_point(aes(y = female_downsample_median), colour = red, size = 1.8) +
  geom_line(aes(y = male_mean_r2), colour = blue, linewidth = 0.85) +
  geom_point(aes(y = male_mean_r2), colour = blue, size = 1.9) +
  annotate("text", x = 590, y = summary$male_mean_r2[nrow(summary)] + 0.010,
           label = "Male, n = 23", colour = blue, size = 2.7, hjust = 1) +
  annotate("text", x = 590, y = summary$female_downsample_median[nrow(summary)] - 0.008,
           label = "Female downsampling, n = 23", colour = red, size = 2.7, hjust = 1) +
  scale_x_continuous(breaks = c(0, 100, 250, 500, 700), limits = c(0, 610)) +
  scale_y_continuous(limits = c(0, 0.26), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Balanced-sample LD sensitivity analysis",
       subtitle = "Ribbon: 95% empirical interval from 1,000\nfemale downsampling replicates",
       x = "Marker distance (kb)", y = expression("Mean " * r^2)) +
  theme_paper()

p_b <- ggplot(permutation, aes(global_mean_difference)) +
  geom_histogram(bins = 35, fill = grey, colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = observed_global, colour = blue, linewidth = 0.9) +
  annotate("text", x = observed_global, y = Inf,
           label = sprintf("Observed = %.3f\nEmpirical P = %.3f", observed_global, empirical_p),
           vjust = 1.4, hjust = 1.05, colour = blue, size = 2.8) +
  labs(title = "Individual-level sex-label permutation",
       subtitle = "1,000 permutations preserving\n23 male and 37 female labels",
       x = expression("Mean male-female " * Delta * r^2 * " across distance bins"),
       y = "Permutations") +
  theme_paper()

fig <- p_a + p_b +
  plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11)))
width_in <- 183 / 25.4
height_in <- 88 / 25.4

grDevices::pdf(file.path(out, "Supplementary_LD_robustness.pdf"), width = width_in, height = height_in,
               family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
print(fig); dev.off()
svglite::svglite(file.path(out, "Supplementary_LD_robustness.svg"), width = width_in, height = height_in, bg = "white")
print(fig); dev.off()
agg_tiff(file.path(out, "Supplementary_LD_robustness_600dpi.tiff"), width = width_in, height = height_in,
         units = "in", res = 600, compression = "lzw", background = "white")
print(fig); dev.off()
agg_png(file.path(out, "Supplementary_LD_robustness_preview.png"), width = width_in, height = height_in,
        units = "in", res = 300, background = "white")
print(fig); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
