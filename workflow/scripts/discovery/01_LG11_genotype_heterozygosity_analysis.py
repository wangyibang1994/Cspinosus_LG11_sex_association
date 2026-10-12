#!/usr/bin/env python3
"""Assess genotype evidence that can distinguish male- from female-heterogamety.

This script uses the already filtered 60-individual LG11 VCF and corrected sex
metadata. It does not infer a sex system from allele-frequency differences alone.
It quantifies sex-specific heterozygosity, identifies whether high-|dAF| SNPs
are enriched for heterozygotes in males or females, and produces an evidence
summary for the subsequent coverage/k-mer/SV analyses.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, binomtest


ROOT = Path(__file__).resolve().parents[2]
VCF = ROOT / "candidate_interval_analysis/server_results/LG11_4.8_5.5Mb.filtered.snps.vcf.gz"
META = ROOT / "candidate_interval_analysis/corrected_results/sample_sex.corrected.tsv"
OUT = ROOT / "candidate_interval_analysis/sex_system_inference/results_genotype"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260831)
B = 10000


def parse_gt(field: str, gt_index: int) -> float:
    parts = field.split(":")
    gt = parts[gt_index] if gt_index < len(parts) else "."
    return {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}.get(gt, np.nan)


def load_vcf():
    samples, rows, geno = None, [], []
    with gzip.open(VCF, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or "," in fields[4] or len(fields[3]) != 1 or len(fields[4]) != 1:
                continue
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                continue
            gt_index = fmt.index("GT")
            g = [parse_gt(x, gt_index) for x in fields[9:]]
            rows.append((fields[0], int(fields[1]), fields[3], fields[4]))
            geno.append(g)
    if samples is None:
        raise RuntimeError("No #CHROM header in VCF")
    return samples, pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"]), np.asarray(geno, dtype=float)


def empirical_difference(values: np.ndarray, sexes: np.ndarray) -> tuple[float, float]:
    f = values[sexes == "F"]
    m = values[sexes == "M"]
    obs = float(np.nanmean(m) - np.nanmean(f))
    perm = np.empty(B)
    for b in range(B):
        shuffled = RNG.permutation(sexes)
        perm[b] = np.nanmean(values[shuffled == "M"]) - np.nanmean(values[shuffled == "F"])
    p_two_sided = (1 + np.sum(np.abs(perm) >= abs(obs))) / (B + 1)
    return obs, p_two_sided


def main():
    samples, info, G = load_vcf()
    meta = pd.read_csv(META, sep="\t")
    meta = meta.set_index("Sample").loc[samples].reset_index()
    sexes = meta["Group"].to_numpy()
    if set(sexes) != {"F", "M"}:
        raise ValueError("Corrected metadata must contain F and M")

    called = np.sum(~np.isnan(G), axis=0)
    n_het = np.sum(G == 1, axis=0)
    individual = meta.copy()
    individual["called_SNPs"] = called
    individual["heterozygous_SNPs"] = n_het
    individual["heterozygosity"] = n_het / called
    individual.to_csv(OUT / "LG11_individual_heterozygosity_by_sex.csv", index=False)

    fmask, mmask = sexes == "F", sexes == "M"
    f_af = np.nanmean(G[:, fmask], axis=1) / 2
    m_af = np.nanmean(G[:, mmask], axis=1) / 2
    f_het = np.nanmean(G[:, fmask] == 1, axis=1)
    m_het = np.nanmean(G[:, mmask] == 1, axis=1)
    called_f = np.sum(~np.isnan(G[:, fmask]), axis=1)
    called_m = np.sum(~np.isnan(G[:, mmask]), axis=1)
    sites = info.copy()
    sites["AF_F"] = f_af
    sites["AF_M"] = m_af
    sites["abs_delta_AF"] = np.abs(f_af - m_af)
    sites["heterozygosity_F"] = f_het
    sites["heterozygosity_M"] = m_het
    sites["delta_heterozygosity_M_minus_F"] = m_het - f_het
    sites["called_F"] = called_f
    sites["called_M"] = called_m
    sites = sites.sort_values("abs_delta_AF", ascending=False).reset_index(drop=True)
    sites.to_csv(OUT / "LG11_sitewise_AF_and_heterozygosity_by_sex.csv", index=False)

    panel_rows = []
    for panel, sub in [("all_QC_SNPs", sites), ("top_50_abs_delta_AF", sites.head(50)), ("top_120_abs_delta_AF", sites.head(120)), ("top_200_abs_delta_AF", sites.head(200))]:
        # For panels other than all sites, calculate per-individual heterozygosity specifically within panel.
        if panel != "all_QC_SNPs":
            idx = np.isin(info.POS.to_numpy(), sub.POS.to_numpy())
            vals = np.nanmean(G[idx] == 1, axis=0)
        else:
            vals = individual["heterozygosity"].to_numpy()
        score, perm_p = empirical_difference(vals, sexes)
        higher_m = int(np.sum(sub["delta_heterozygosity_M_minus_F"] > 0))
        higher_f = int(np.sum(sub["delta_heterozygosity_M_minus_F"] < 0))
        ties = int(len(sub) - higher_m - higher_f)
        non_ties = higher_m + higher_f
        sign_p = binomtest(higher_m, non_ties, 0.5).pvalue if non_ties else np.nan
        panel_rows.append({
            "panel": panel,
            "n_SNPs": len(sub),
            "mean_heterozygosity_F": float(np.nanmean(vals[fmask])),
            "mean_heterozygosity_M": float(np.nanmean(vals[mmask])),
            "M_minus_F": score,
            "empirical_two_sided_P": perm_p,
            "SNPs_M_het_excess": higher_m,
            "SNPs_F_het_excess": higher_f,
            "SNPs_tied": ties,
            "binomial_P_M_vs_F_het_excess": sign_p,
        })
    summary = pd.DataFrame(panel_rows)
    summary.to_csv(OUT / "LG11_heterozygosity_sex_system_evidence_summary.csv", index=False)

    # Standard parametric/non-parametric sample-level comparisons for reporting only.
    t = ttest_ind(individual.loc[fmask, "heterozygosity"], individual.loc[mmask, "heterozygosity"], equal_var=False)
    u = mannwhitneyu(individual.loc[fmask, "heterozygosity"], individual.loc[mmask, "heterozygosity"], alternative="two-sided")
    pd.DataFrame([{
        "Welch_t": t.statistic, "Welch_two_sided_P": t.pvalue,
        "Mann_Whitney_U": u.statistic, "Mann_Whitney_two_sided_P": u.pvalue,
        "interpretation_rule": "Positive M_minus_F supports a male-heterogametic-compatible genotype pattern; negative supports a female-heterogametic-compatible pattern. Neither alone identifies XY or ZW."
    }]).to_csv(OUT / "LG11_individual_heterozygosity_tests.csv", index=False)

    colors = {"F": "#CC6677", "M": "#4477AA"}
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), gridspec_kw={"width_ratios": [1.0, 1.25, 1.1]})

    # A: individual heterozygosity.
    for group, x in [("F", 0), ("M", 1)]:
        y = individual.loc[individual.Group == group, "heterozygosity"].to_numpy()
        jitter = RNG.normal(x, 0.045, len(y))
        axes[0].scatter(jitter, y, s=28, color=colors[group], alpha=0.85, edgecolor="white", linewidth=0.4)
        axes[0].boxplot(y, positions=[x], widths=0.42, patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor="none", color=colors[group]), medianprops=dict(color="black"), whiskerprops=dict(color=colors[group]), capprops=dict(color=colors[group]))
    axes[0].set_xticks([0, 1], ["Female", "Male"])
    axes[0].set_ylabel("Heterozygous genotype proportion")
    axes[0].set_title("A  Individual heterozygosity across LG11")

    # B: sex-specific heterozygosity at sites ranked by dAF.
    show = sites.head(200)
    axes[1].scatter(show.abs_delta_AF, show.delta_heterozygosity_M_minus_F, c=np.where(show.delta_heterozygosity_M_minus_F >= 0, colors["M"], colors["F"]), alpha=0.65, s=18, linewidth=0)
    axes[1].axhline(0, color="#666666", lw=0.8)
    axes[1].set_xlabel("|ΔAF| (top 200 SNPs)")
    axes[1].set_ylabel("Male − female heterozygosity")
    axes[1].set_title("B  SNP-level heterozygosity contrast")

    # C: top panel scores.
    plot_df = summary.loc[summary.panel != "all_QC_SNPs"].copy()
    labels = [x.replace("top_", "Top ").replace("_abs_delta_AF", "") for x in plot_df.panel]
    vals = plot_df["M_minus_F"].to_numpy()
    axes[2].bar(range(len(vals)), vals, color=[colors["M"] if x >= 0 else colors["F"] for x in vals], width=0.65)
    axes[2].axhline(0, color="#333333", lw=0.8)
    axes[2].set_xticks(range(len(vals)), labels)
    axes[2].set_ylabel("Male − female heterozygosity")
    axes[2].set_title("C  High-|ΔAF| SNP panels")
    axes[2].text(0.02, 0.98, "Positive: male-heterogametic-compatible\nNegative: female-heterogametic-compatible", transform=axes[2].transAxes, va="top", fontsize=8)

    fig.suptitle("LG11 genotype evidence relevant to XY versus ZW inference", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "FigureS_LG11_genotype_heterozygosity_XY_ZW_evidence.png", dpi=400, bbox_inches="tight")
    fig.savefig(OUT / "FigureS_LG11_genotype_heterozygosity_XY_ZW_evidence.pdf", bbox_inches="tight")
    plt.close(fig)

    with open(OUT / "README_interpretation.txt", "w") as handle:
        handle.write(
            "This analysis tests only whether called SNP genotypes show excess heterozygosity in males or females. "
            "Male excess is compatible with XY; female excess is compatible with ZW. A result is not sufficient to assign XY/ZW because reference-assembly composition, incomplete linkage, genotype error and structural variation can alter the pattern. Coverage, k-mer and SV/CNV analyses must be integrated before inference.\n"
        )


if __name__ == "__main__":
    main()
