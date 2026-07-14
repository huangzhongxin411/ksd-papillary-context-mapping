#!/usr/bin/env python3
from pathlib import Path
import csv
import textwrap

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin  # noqa: F401


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "supplementary_figures"
TAB_DIR = ROOT / "supplementary_tables"
MAT_DIR = ROOT / "supplementary_materials"
AUDIT_DIR = ROOT / "results/tables"

TEAL = "#245B63"
BLUE_GREY = "#7F9BA5"
GOLD = "#B89A55"
TERRACOTTA = "#A15D4B"
LIGHT = "#E9EFF0"
DARK = "#243238"


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def combine_with_record_type(specs):
    frames = []
    for record_type, relative_path in specs:
        path = ROOT / relative_path
        frame = read_tsv(path)
        frame.insert(0, "packaging_source_file", relative_path)
        frame.insert(0, "record_type", record_type)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False).fillna("")


def make_s1():
    gwas = pd.read_csv(
        ROOT / "source_data/figures/phase1_step4_gwas_plot_source_data.tsv.gz",
        sep="\t", usecols=["P"]
    )
    p = pd.to_numeric(gwas["P"], errors="coerce")
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)].to_numpy()
    ranked = read_tsv(ROOT / "results/tables/phase1_step3_magma_ranked_canonical.tsv")
    ranked["magma_p"] = pd.to_numeric(ranked["magma_p"])
    ranked["magma_rank"] = pd.to_numeric(ranked["magma_rank"])
    top = ranked.nsmallest(20, "magma_rank").sort_values("magma_rank", ascending=False)
    summary = read_tsv(ROOT / "results/tables/phase1_step3_magma_gene_set_summary.tsv")
    summary["n_genes"] = pd.to_numeric(summary["n_genes"])

    canvas = Image.new("RGB", (3600, 2600), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((120, 70), "Supplementary Figure S1 | GWAS and MAGMA diagnostics", fill=TEAL, font=font(62, bold=True))
    draw.text((120, 165), "QC and gene-priority summaries only; not fine mapping or causal-gene evidence.", fill=DARK, font=font(32))

    panels = [(120, 300, 1640, 1230), (1900, 300, 3420, 1230), (120, 1450, 1640, 2380), (1900, 1450, 3420, 2380)]

    def panel_frame(box, title):
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=16, outline="#C9D4D7", width=3, fill="white")
        draw.text((x0 + 35, y0 + 25), title, fill=DARK, font=font(36, bold=True))
        return x0 + 120, y0 + 130, x1 - 60, y1 - 110

    def axes_frame(plot_box, x_label, y_label):
        x0, y0, x1, y1 = plot_box
        draw.line((x0, y1, x1, y1), fill=DARK, width=3)
        draw.line((x0, y0, x0, y1), fill=DARK, width=3)
        draw.text(((x0 + x1) // 2 - 120, y1 + 35), x_label, fill=DARK, font=font(26))
        draw.text((x0 - 85, y0 - 55), y_label, fill=DARK, font=font(24))

    # A: QQ diagnostic, downsampled only for raster display.
    plot = panel_frame(panels[0], "A  GWAS QQ diagnostic")
    axes_frame(plot, "Expected -log10(P)", "Observed")
    expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    observed = -np.log10(np.sort(p))
    display_idx = np.linspace(0, len(p) - 1, min(6000, len(p))).astype(int)
    expected_max = float(expected.max())
    observed_max = float(observed.max())
    x0, y0, x1, y1 = plot
    reference_y = y1 - int(expected_max / observed_max * (y1 - y0))
    draw.line((x0, y1, x1, reference_y), fill=TERRACOTTA, width=4)
    for idx in display_idx:
        x = x0 + int(expected[idx] / expected_max * (x1 - x0))
        y = y1 - int(observed[idx] / observed_max * (y1 - y0))
        draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill=TEAL)

    # B: P-value histogram.
    plot = panel_frame(panels[1], "B  GWAS P-value distribution")
    axes_frame(plot, "P value", "Count")
    counts, _ = np.histogram(p, bins=80, range=(0, 1))
    x0, y0, x1, y1 = plot
    max_count = max(counts)
    bar_width = (x1 - x0) / len(counts)
    for index, count in enumerate(counts):
        left = x0 + index * bar_width
        right = x0 + (index + 1) * bar_width
        top_y = y1 - (count / max_count) * (y1 - y0)
        draw.rectangle((left, top_y, right, y1), fill=BLUE_GREY)

    # C: top 20 locked MAGMA ranks.
    plot = panel_frame(panels[2], "C  Top 20 MAGMA genes by locked rank")
    x0, y0, x1, y1 = plot
    values = -np.log10(top["magma_p"].to_numpy())
    max_value = values.max()
    row_height = (y1 - y0) / len(top)
    for index, ((_, row), value) in enumerate(zip(top.iterrows(), values)):
        y = y0 + index * row_height
        draw.text((x0, y), str(row["gene_symbol"]), fill=DARK, font=font(22))
        left = x0 + 220
        right = left + (value / max_value) * (x1 - left)
        draw.rectangle((left, y + 4, right, y + row_height - 6), fill=TEAL)
    draw.text((x0 + 520, y1 + 35), "-log10(MAGMA P)", fill=DARK, font=font(26))

    # D: frozen module sizes.
    plot = panel_frame(panels[3], "D  Frozen MAGMA module sizes")
    x0, y0, x1, y1 = plot
    labels = summary["gene_set_name"].str.replace("MAGMA_", "", regex=False).str.replace("_p1e4", "", regex=False).tolist()
    values = summary["n_genes"].tolist()
    colors = [TEAL, BLUE_GREY, TERRACOTTA, GOLD, "#879197"]
    slot = (x1 - x0) / len(values)
    max_value = max(values)
    for index, (label, value, color) in enumerate(zip(labels, values, colors)):
        left = x0 + index * slot + 30
        right = x0 + (index + 1) * slot - 30
        top_y = y1 - (value / max_value) * (y1 - y0 - 70)
        draw.rectangle((left, top_y, right, y1), fill=color)
        draw.text((left + 15, top_y - 42), str(value), fill=DARK, font=font(25, bold=True))
        draw.text((left, y1 + 25), label, fill=DARK, font=font(21))

    png = FIG_DIR / "Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.png"
    pdf = FIG_DIR / "Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.pdf"
    canvas.save(png, dpi=(300, 300), optimize=True)
    canvas.save(pdf, resolution=300)
    return png, pdf


def font(size, bold=False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def stack_components(title, subtitle, components, out_stem):
    target_width = 3600
    margin = 100
    header_height = 280
    gap = 80
    resized = []
    for label, path in components:
        image = Image.open(path).convert("RGB")
        ratio = (target_width - 2 * margin) / image.width
        image = image.resize((target_width - 2 * margin, round(image.height * ratio)), Image.Resampling.LANCZOS)
        resized.append((label, image))
    total_height = header_height + sum(image.height + 80 for _, image in resized) + gap * (len(resized) - 1) + margin
    canvas = Image.new("RGB", (target_width, total_height), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((margin, 60), title, fill=TEAL, font=font(54, bold=True))
    draw.text((margin, 145), subtitle, fill=DARK, font=font(30))
    y = header_height
    for label, image in resized:
        draw.text((margin, y), label, fill=DARK, font=font(36, bold=True))
        y += 65
        canvas.paste(image, (margin, y))
        y += image.height + gap
    png = FIG_DIR / f"{out_stem}.png"
    pdf = FIG_DIR / f"{out_stem}.pdf"
    canvas.save(png, dpi=(300, 300), optimize=True)
    canvas.save(pdf, resolution=300)
    return png, pdf


def make_s2():
    return stack_components(
        "Supplementary Figure S2 | Descriptive spatial projection",
        "10 complete sections: 4 GSE206306 + 6 GSE231630 | broad-compartment context only; no lesion ROI",
        [
            ("S2A  Spatial resource and label-transfer quality control", ROOT / "results/figures/phase3_step4_SuppFig_spatial_label_transfer_QC.png"),
            ("S2B  Descriptive MAGMA module overlays", ROOT / "results/figures/phase3_step4_SuppFig_spatial_MAGMA_module_overlays.png"),
        ],
        "Supplementary_Figure_S2_spatial_projection",
    )


def make_s3():
    return stack_components(
        "Supplementary Figure S3 | Kidney_Cortex TWAS proxy and candidate evidence",
        "Proxy annotation and reporting categories only; not papilla-specific regulation or causal tiers",
        [
            ("S3A  TWAS model-quality and MAGMA-overlap audit", ROOT / "results/figures/phase4_step3_SuppFig_TWAS_proxy_quality.png"),
            ("S3B  Full 235-gene candidate reporting matrix", ROOT / "supplementary_figures/Supplementary_Figure_S3_full_candidate_matrix_bulk_reviewed.png"),
        ],
        "Supplementary_Figure_S3_TWAS_candidate_proxy",
    )


def make_tables():
    outputs = {}
    table1 = combine_with_record_type([
        ("reconstructed_locus", "results/tables/phase1_step4_reconstructed_loci.tsv"),
        ("distance_pruned_lead_snp", "results/tables/phase1_step4_reconstructed_lead_snps.tsv"),
    ])
    outputs[1] = ("Supplementary_Table_1_locus_reconstruction.tsv", table1)

    table2 = combine_with_record_type([
        ("ranked_MAGMA_gene", "results/tables/phase1_step3_magma_ranked_canonical.tsv"),
        ("frozen_module_summary", "results/tables/phase1_step3_magma_gene_set_summary.tsv"),
    ])
    outputs[2] = ("Supplementary_Table_2_MAGMA_results_and_modules.tsv", table2)

    table3 = combine_with_record_type([
        ("donor_compartment_count", "results/tables/phase2_step2_donor_compartment_counts.tsv"),
        ("donor_LoopTAL_contrast", "results/tables/phase2_step2_loop_tal_vs_other_summary.tsv"),
        ("matched_random_benchmark", "results/tables/phase2_step3_matched_random_benchmark_summary.tsv"),
        ("low_cell_count_sensitivity", "results/tables/phase2_step3_low_cell_count_sensitivity.tsv"),
        ("original_vs_driver_removed", "results/tables/phase2_step4_original_vs_driver_removed_summary.tsv"),
        ("driver_removed_matched_random", "results/tables/phase2_step4_driver_removed_matched_random_benchmark_summary.tsv"),
    ])
    outputs[3] = ("Supplementary_Table_3_snRNA_context_and_sensitivity.tsv", table3)

    table4 = read_tsv(ROOT / "results/tables/phase4_step1_twas_proxy_gene_table.tsv")
    table4.insert(0, "source_file", "results/tables/phase4_step1_twas_proxy_gene_table.tsv")
    outputs[4] = ("Supplementary_Table_4_TWAS_proxy_results.tsv", table4)

    table5 = read_tsv(ROOT / "results/tables/phase4_step2_candidate_gene_evidence_matrix.tsv")
    table5.insert(0, "source_file", "results/tables/phase4_step2_candidate_gene_evidence_matrix.tsv")
    table5["bulk_status_if_available"] = "disease_context_reviewed_non_upgrading"
    table5["notes"] = table5["notes"].astype(str) + " Bulk status updated during Phase 6 packaging; no reporting-group reassignment."
    outputs[5] = ("Supplementary_Table_5_candidate_reporting_matrix.tsv", table5)

    table6 = combine_with_record_type([
        ("canonical_module_base_model", "results/tables/phase5_step2_bulk_paired_model_results.tsv"),
        ("tissue_state_paired_model", "results/tables/phase5_step3_tissue_state_paired_model_results.tsv"),
        ("tissue_state_adjusted_module_model", "results/tables/phase5_step3_tissue_state_adjusted_module_models.tsv"),
        ("module_tissue_state_correlation", "results/tables/phase5_step3_module_tissue_state_correlations.tsv"),
        ("interpretation_summary", "results/tables/phase5_step3_tissue_state_sensitivity_interpretation_summary.tsv"),
    ])
    outputs[6] = ("Supplementary_Table_6_bulk_disease_context_and_tissue_state.tsv", table6)

    for number, (filename, frame) in outputs.items():
        write_tsv(frame, TAB_DIR / filename)
    return outputs


def write_crosswalk(outputs):
    rows = [
        ["Supplementary Table 1", "Lead SNP and reconstructed locus table", "phase1_step4_reconstructed_loci.tsv; phase1_step4_reconstructed_lead_snps.tsv", "Older v0.2 S2 plus Phase 1-Step 4 outputs", "57 reconstructed loci and 60 lead SNPs", "Distance-based reconstruction; not fine mapping."],
        ["Supplementary Table 2", "MAGMA results and frozen modules", "phase1_step3_magma_ranked_canonical.tsv; phase1_step3_magma_gene_set_summary.tsv", "Older v0.2 S3 plus Phase 1-Step 3 outputs", "17,316 ranked genes and five frozen module definitions", "EUR LD reference boundary applies."],
        ["Supplementary Table 3", "snRNA context and sensitivity", "Phase 2-Step 2 to Step 4 donor, matched-random, low-cell and driver-removal tables", "Older v0.2 S4-S6", "Donor-level context and sensitivity outputs", "Donors are support units; no causal cell-type assignment."],
        ["Supplementary Table 4", "Kidney_Cortex TWAS proxy results", "phase4_step1_twas_proxy_gene_table.tsv", "Phase 4-Step 1 and older TWAS supplementary tables", "All 5,989 tested models and proxy-quality fields", "Proxy annotation only."],
        ["Supplementary Table 5", "Candidate reporting matrix", "phase4_step2_candidate_gene_evidence_matrix.tsv; Phase 6 bulk-status audit", "Phase 4-Step 2 candidate matrix", "Full 235-gene R1-R6 reporting matrix", "Bulk review is non-upgrading; groups are not causal tiers."],
        ["Supplementary Table 6", "Bulk disease context and tissue state", "Phase 5-Step 2 base models and Phase 5-Step 3 tissue-state sensitivity tables", "Older v0.2 S8-S11 plus Phase 5 outputs", "Paired module, proxy and adjusted-model results", "Filename-derived labels; proxies are not cell fractions."],
    ]
    path = TAB_DIR / "Supplementary_Table_crosswalk_old_to_final.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["final_table", "final_table_title", "source_files", "older_table_or_phase_outputs", "main_content", "notes"])
        writer.writerows(rows)


def write_legends():
    MAT_DIR.mkdir(parents=True, exist_ok=True)
    (MAT_DIR / "Supplementary_Figure_Legends_v3.2.md").write_text(
        """# Supplementary Figure Legends v3.2

## Supplementary Figure S1. GWAS and MAGMA diagnostic summaries

Panels show (A) the GWAS quantile-quantile diagnostic, (B) the sampled GWAS P-value distribution, (C) the top 20 genes by locked MAGMA rank and (D) the sizes of the five frozen MAGMA modules. GWAS panels use the cleaned diagnostic sample retained for visualization; MAGMA panels are derived from the locked canonical gene results without rerunning MAGMA. These summaries document input behavior and EUR-LD-reference-based gene priority, not locus fine mapping or causal-gene evidence.

## Supplementary Figure S2. Descriptive spatial broad-compartment projection

Panels combine (A) spatial-resource and label-transfer quality control with (B) descriptive locked-MAGMA-module overlays across ten complete sections, comprising four GSE206306 and six GSE231630 sections. Loop/TAL predictions were zero or sparse, and no plaque, mineral, fibrosis or lesion region-of-interest annotation was available. The figure provides supplementary anatomical/tissue-context projection and does not support plaque localization, Loop/TAL enrichment, disease-control comparison or causal spatial inference.

## Supplementary Figure S3. Kidney_Cortex TWAS proxy and candidate evidence display

Panels combine (A) TWAS model-quality and MAGMA-overlap summaries with (B) the full 235-gene candidate reporting matrix. Kidney_Cortex S-PrediXcan tested 5,989 models and yielded 51 FDR-supported results, including 42 one-SNP and nine multi-SNP models. R1-R6 are mutually exclusive reporting groups, not causal tiers. The figure presents cortex-proxy annotation and evidence organization, not papilla-specific regulation or causal-gene validation.
""",
        encoding="utf-8",
    )
    (MAT_DIR / "Supplementary_Table_Captions_v3.2.md").write_text(
        """# Supplementary Table Captions v3.2

## Supplementary Table 1. Lead SNP and reconstructed locus table
Distance-pruned lead SNPs and overlapping +/-1 Mb windows used to reconstruct 57 loci from the locked cleaned GWAS. This reproducible reconstruction does not replace the 59 loci reported by the source publication.

## Supplementary Table 2. MAGMA gene results and frozen modules
All ranked canonical MAGMA gene results plus definitions and membership counts for the top-50, top-100, Bonferroni, FDR and suggestive modules. MAGMA used the 1000 Genomes European LD reference and provides genetic priority rather than causal-gene evidence.

## Supplementary Table 3. Donor-level snRNA context and sensitivity analyses
Donor-by-compartment counts and Loop/TAL contrasts, expression- and detection-matched random benchmarks, low-cell-count sensitivity and known-driver-removal sensitivity. Donors are the biological support units.

## Supplementary Table 4. Kidney_Cortex TWAS proxy results
All 5,989 tested S-PrediXcan models with FDR status, contributing-SNP count, model class and MAGMA overlap fields. These are Kidney_Cortex proxy annotations, not papilla-specific regulatory evidence.

## Supplementary Table 5. Full candidate reporting matrix
Complete 235-gene evidence matrix and deterministic R1-R6 assignments. Bulk disease-context review is recorded as non-upgrading. Reporting groups are not causal tiers.

## Supplementary Table 6. Paired bulk disease-context and tissue-state sensitivity
Canonical-module paired models, tissue-state proxy models, adjusted sensitivity models, module-proxy correlations and interpretation summaries for GSE73680. Group labels are filename-derived and manually reviewable; tissue-state scores are expression proxies, not cell fractions.
""",
        encoding="utf-8",
    )


def write_audits(outputs, figure_outputs):
    figure_audit = AUDIT_DIR / "phase6_step6_supplementary_figure_materialization_audit.tsv"
    with figure_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["supplementary_figure", "expected_panels", "created_file", "component_files_used", "missing_components", "legend_matches", "manual_polishing_needed", "notes"])
        writer.writerow(["S1", "QQ; P-value distribution; MAGMA gene ranks; frozen module sizes", str(figure_outputs[1][1].relative_to(ROOT)), "GWAS plot source data; phase1_step3_magma_ranked_canonical.tsv; phase1_step3_magma_gene_set_summary.tsv", "none", "yes", "no", "Replotted from locked source data; no GWAS QC or MAGMA rerun."])
        writer.writerow(["S2", "label-transfer QC; no-ROI/Loop-TAL sparsity boundary; descriptive overlays", str(figure_outputs[2][1].relative_to(ROOT)), "phase3_step4_SuppFig_spatial_label_transfer_QC.png; phase3_step4_SuppFig_spatial_MAGMA_module_overlays.png", "none", "yes", "yes", "Composite header states 4 GSE206306 + 6 GSE231630; component headings retain historical internal labels."])
        writer.writerow(["S3", "TWAS proxy quality; one-/multi-SNP burden; MAGMA overlap; full candidate matrix", str(figure_outputs[3][1].relative_to(ROOT)), "phase4_step3_SuppFig_TWAS_proxy_quality.png; Supplementary_Figure_S3_full_candidate_matrix_bulk_reviewed.png", "none", "yes", "no", "Full matrix was re-exported from the locked 235-gene table with bulk status reviewed/non-upgrading; R1-R6 assignments were unchanged."])

    table_audit = AUDIT_DIR / "phase6_step6_supplementary_table_materialization_audit.tsv"
    with table_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["supplementary_table", "created_file", "source_files_used", "row_count", "column_count", "legend_or_caption_available", "manual_review_needed", "notes"])
        for number, (filename, frame) in outputs.items():
            source_column = "packaging_source_file" if "packaging_source_file" in frame.columns else "source_file"
            source_values = sorted(set(frame[source_column])) if source_column in frame.columns else []
            writer.writerow([
                f"Supplementary Table {number}", str((TAB_DIR / filename).relative_to(ROOT)), ";".join(source_values),
                len(frame), len(frame.columns), "yes", "yes" if number in (1, 2, 3, 6) else "no",
                "Multi-component TSV uses record_type/source_file columns." if number in (1, 2, 3, 6) else "Direct locked table materialization."
            ])


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    figure_outputs = {1: make_s1(), 2: make_s2(), 3: make_s3()}
    outputs = make_tables()
    write_crosswalk(outputs)
    write_legends()
    write_audits(outputs, figure_outputs)


if __name__ == "__main__":
    main()
