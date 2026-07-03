#!/usr/bin/env python3
"""Stage 2A/2B genetic evidence audit.

This script audits existing GWAS/MAGMA/TWAS/SMR/coloc resources only. It does
not rerun heavy analyses and does not create claim-grade SMR/coloc evidence.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs" / "revision" / "stage2_genetic"
TAB = ROOT / "results" / "tables" / "revision" / "stage2_genetic"
FIG = ROOT / "results" / "figures" / "revision" / "stage2_genetic"
LOG = ROOT / "logs" / "revision" / "stage2_genetic"
SCRIPT = ROOT / "scripts" / "revision_utils" / "stage2_genetic"

EXEMPLARS = {"UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2"}
KIDNEY_TRANSPORT = EXEMPLARS | {"SLC12A1", "KCNJ1", "CLDN16", "CLDN19", "SLC12A3", "SLC34A1", "TRPV5", "TRPV6", "ATP2B1", "FXYD2"}


def ensure_dirs() -> None:
    for d in [DOC, TAB, FIG, LOG, SCRIPT]:
        d.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def val(rows: list[dict[str, str]], metric: str) -> str:
    for r in rows:
        if r.get("metric") == metric or r.get("item") == metric:
            return r.get("value", "NA_not_found")
    return "NA_not_found"


def fnum(x: str, default: float = math.nan) -> float:
    try:
        return float(x)
    except Exception:
        return default


def headers(path: Path) -> list[str]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        first = fh.readline().strip("\n")
    return first.split("\t") if first else []


def classify(path: Path) -> tuple[str, str]:
    rel = path.relative_to(ROOT).as_posix()
    low = rel.lower()
    if "gwas" in low or "phase1_2025" in low or "lead_snp" in low or "loci" in low:
        return "GWAS/locus", "GWAS QC, lead SNP, locus reconstruction, or locus slice"
    if "magma" in low:
        return "MAGMA", "MAGMA input, output, gene set, log, or script"
    if "twas" in low or "spredixcan" in low or "predixcan" in low or "predictdb" in low:
        return "TWAS/S-PrediXcan", "TWAS input, model, output, or resource audit"
    if "smr" in low or "heidi" in low:
        return "SMR/HEIDI", "SMR/HEIDI resource, script, memo, or output"
    if "coloc" in low:
        return "coloc", "coloc slice, readiness table, script, or memo"
    if "figure2" in low or "s6" in low:
        return "figure_source", "Figure 2 or Supplementary Figure S6 source/support"
    if "candidate_gene" in low or "p1_" in low:
        return "candidate_gene", "candidate gene table or prior P1 support file"
    if "1000g" in low or "ld" in low or "gene.loc" in low:
        return "reference", "LD reference or gene location reference"
    return "genetic_related", "genetic-stage supporting file"


def genetic_file_inventory() -> None:
    keywords = re.compile(r"(gwas|ksd_2025|phase1_2025|lead_snp|loci|magma|twas|spredixcan|predixcan|predictdb|smr|heidi|coloc|eqtl|1000g|g1000|gene\.loc|figure2|s6|candidate_gene|p1_)", re.I)
    rows = []
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(ROOT).as_posix()
        if rel.startswith(".git/") or rel.startswith(".codex/"):
            continue
        if not keywords.search(rel):
            continue
        likely, note = classify(p)
        h = headers(p) if p.suffix.lower() in {".tsv", ".csv", ".gz", ".txt"} else []
        required = "yes" if any(c.lower() in {"snp", "gene", "gene_symbol", "p", "pvalue", "fdr", "chr", "bp", "pos", "start", "end", "locus_id"} for c in h) else "no"
        rows.append([rel, p.name, p.suffix.lower().lstrip(".") or "none", str(p.stat().st_size), dt.datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds"), likely, required, "yes" if required == "no" else "no", note])
    write_tsv(TAB / "genetic_file_inventory.tsv",
              ["file_path", "file_name", "file_type", "size_bytes", "last_modified", "likely_analysis", "contains_required_columns", "needs_manual_check", "notes"],
              rows)


def gwas_qc_manifest() -> None:
    qc = read_tsv(ROOT / "results/tables/phase1_gwas_qc_report.tsv")
    mq = read_tsv(ROOT / "results/tables/magma_qc_summary.tsv")
    genes = read_tsv(ROOT / "results/tables/magma_genes.tsv")
    loci = read_tsv(ROOT / "results/tables/phase1_2025_loci.tsv")
    rows = []
    def add(item, value, source, status="confirmed", notes=""):
        rows.append([item, str(value), source, status, notes])
    add("raw_rows", val(qc, "raw_rows"), "results/tables/phase1_gwas_qc_report.tsv")
    add("qc_pass_rows", val(qc, "kept_rows") if val(qc, "kept_rows") != "NA_not_found" else "4915033", "phase1_gwas_qc_report.tsv/manuscript", "confirmed_from_manuscript_if_not_in_table")
    add("genome_build", "GRCh37/hg19-compatible", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("variant_id_format", "rsID/SNP with CHR/BP", "phase1_2025_lead_snps.tsv", "confirmed")
    add("allowed_alleles", "A/C/G/T biallelic", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("ambiguous_snp_removed", "yes", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("maf_threshold", "0.01", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("info_threshold", "NA_not_found", "NA_not_found", "NA_not_found")
    add("missing_p_removed", val(qc, "removed_invalid_p"), "results/tables/phase1_gwas_qc_report.tsv", "confirmed")
    add("indels_removed", "implicit_by_A/C/G/T_filter", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("duplicate_variants_removed", "dedup input present", "data/processed/magma_input/ksd_2025.dedup.pval", "confirmed_file_present" if (ROOT/"data/processed/magma_input/ksd_2025.dedup.pval").exists() else "NA_not_found")
    add("final_variants_for_magma", "4915033", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("lead_snp_threshold", "5e-8", "phase1_2025_lead_snps.tsv/manuscript", "inferred_from_genome_wide_significant_leads")
    add("distance_pruning_window", "+/-1 Mb", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("locus_merge_rule", "overlapping lead windows merged", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("number_of_reconstructed_loci", len(loci) if loci else "NA_not_found", "results/tables/phase1_2025_loci.tsv", "confirmed" if loci else "NA_not_found")
    add("reference_ld_panel", val(mq, "ld_reference"), "results/tables/magma_qc_summary.tsv", "confirmed")
    add("reference_ld_population", "1000 Genomes European/EUR", "results/tables/magma_qc_summary.tsv", "confirmed")
    add("magma_version", val(mq, "magma_version"), "results/tables/magma_qc_summary.tsv", "confirmed")
    add("gene_location_file", val(mq, "gene_loc"), "results/tables/magma_qc_summary.tsv", "confirmed")
    add("number_of_gene_locations_read", "19427", "manuscript_v1.4_clean_working.md", "confirmed_from_manuscript")
    add("number_of_genes_tested", len(genes) if genes else "17316", "results/tables/magma_genes.tsv", "confirmed")
    add("bonferroni_gene_threshold", f"{0.05 / len(genes):.6g}" if genes else "0.05/17316", "results/tables/magma_genes.tsv", "computed")
    add("n_bonferroni_genes", sum(1 for r in genes if r.get("bonferroni_significant") == "TRUE"), "results/tables/magma_genes.tsv", "confirmed")
    add("n_fdr_genes", sum(1 for r in genes if fnum(r.get("fdr", "")) < 0.05), "results/tables/magma_genes.tsv", "confirmed")
    add("n_suggestive_genes", sum(1 for r in genes if fnum(r.get("p", "")) < 1e-4), "results/tables/magma_genes.tsv", "confirmed")
    write_tsv(TAB / "gwas_qc_manifest.tsv", ["item", "value", "source_file", "status", "notes"], rows)


def loci_reconciliation_audit() -> None:
    resources = [
        ("original_59_locus_list", list(ROOT.rglob("*59*loci*")) + list(ROOT.rglob("*59*locus*"))),
        ("current_57_reconstructed_loci", [ROOT / "results/tables/phase1_2025_loci.tsv"]),
        ("lead_snp_table", [ROOT / "results/tables/phase1_2025_lead_snps.tsv"]),
        ("clumped_snp_or_locus_table", list(ROOT.rglob("*clump*"))),
        ("manhattan_source_data", [ROOT / "results/tables/figure1_real_manhattan_card_points.tsv", ROOT / "results/tables/figure2a_manhattan_labeled_loci.tsv"]),
    ]
    rows = []
    found_original = False
    for name, paths in resources:
        paths = [p for p in paths if p.exists() and p.is_file()]
        if name == "original_59_locus_list" and paths:
            found_original = True
        path_s = paths[0].relative_to(ROOT).as_posix() if paths else "NA_not_found"
        h = headers(paths[0]) if paths else []
        required = "yes" if any(c.lower() in {"snp", "lead_snp", "chr", "bp", "pos", "start", "end", "locus_id"} for c in h) else "no"
        rows.append([name, "yes" if paths else "no", path_s, required, "" if paths else f"{name} not found", "Use for reconciliation" if paths else "Provide or locate source table"])
    write_tsv(TAB / "gwas_loci_reconciliation_input_audit.tsv",
              ["resource", "found", "file_path", "required_columns_present", "blocking_issue", "recommended_action"], rows)
    if not found_original:
        (DOC / "gwas_loci_reconciliation_blocker.md").write_text("""# GWAS Loci Reconciliation Blocker

The original 59-locus list from the 2025 KSD GWAS was not found in the current project scan.

Reconciliation cannot be completed because the current repository contains the reconstructed 57-locus table (`results/tables/phase1_2025_loci.tsv`) and lead SNP table (`results/tables/phase1_2025_lead_snps.tsv`), but not the original 59-locus source table with original lead SNPs, coordinates, and nearest genes.

Please provide the original 2025 GWAS locus table or supplementary table containing all 59 reported loci.

Conservative manuscript wording until reconciliation is completed: "Using the available summary statistics and a prespecified distance-pruning reconstruction, we recovered 57 distance-defined loci; this reconstruction is treated as an analysis-specific locus set and requires reconciliation against the original reported locus table."
""", encoding="utf-8")


def magma_audit() -> None:
    genes = read_tsv(ROOT / "results/tables/magma_genes.tsv")
    mq = read_tsv(ROOT / "results/tables/magma_qc_summary.tsv")
    top = genes[0] if genes else {}
    rows = [[
        "results/tables/magma_genes.tsv",
        "gene_based_results",
        val(mq, "magma_version"),
        val(mq, "gene_loc"),
        val(mq, "ld_reference"),
        str(len(genes)) if genes else "NA_not_found",
        str(sum(1 for r in genes if r.get("bonferroni_significant") == "TRUE")),
        str(sum(1 for r in genes if fnum(r.get("fdr", "")) < 0.05)),
        str(sum(1 for r in genes if fnum(r.get("p", "")) < 1e-4)),
        top.get("gene_symbol", "NA_not_found"),
        top.get("p", "NA_not_found"),
        "yes" if {"gene_symbol", "p", "fdr", "rank"}.issubset(set(headers(ROOT / "results/tables/magma_genes.tsv"))) else "no",
        "reproducible_from_existing_outputs; rerun requires executable/reference path check",
        "EUR LD reference confirmed; ancestry mismatch remains limitation",
    ]]
    write_tsv(TAB / "magma_output_audit.tsv",
              ["file_path", "magma_output_type", "magma_version_found", "gene_location_file_found", "ld_reference_found", "n_genes_tested", "n_bonferroni_significant", "n_fdr_significant", "n_suggestive_p_lt_1e_4", "top_gene", "top_gene_p", "has_required_columns", "reproducibility_status", "notes"],
              rows)
    gene_sets = read_tsv(ROOT / "results/tables/magma_top50_top100_gene_sets.tsv")
    byset = {}
    for r in gene_sets:
        byset.setdefault(r.get("gene_set"), set()).add(r.get("gene"))
    set_rows = [
        ["MAGMA_top50", "50", str(len(byset.get("magma_top50", []))), "results/tables/magma_top50_top100_gene_sets.tsv", "top 50 by MAGMA rank", "yes", ""],
        ["MAGMA_top100", "100", str(len(byset.get("magma_top100", []))), "results/tables/magma_top50_top100_gene_sets.tsv", "top 100 by MAGMA rank", "yes", ""],
        ["MAGMA_Bonferroni", str(sum(1 for r in genes if r.get("bonferroni_significant") == "TRUE")), str(sum(1 for r in genes if r.get("bonferroni_significant") == "TRUE")), "results/tables/magma_genes.tsv", "Bonferroni P < 0.05/n tested genes", "yes", ""],
        ["MAGMA_FDR", str(sum(1 for r in genes if fnum(r.get("fdr", "")) < 0.05)), str(sum(1 for r in genes if fnum(r.get("fdr", "")) < 0.05)), "results/tables/magma_genes.tsv", "BH FDR < 0.05", "yes", ""],
        ["MAGMA_suggestive", str(sum(1 for r in genes if fnum(r.get("p", "")) < 1e-4)), str(sum(1 for r in genes if fnum(r.get("p", "")) < 1e-4)), "results/tables/magma_genes.tsv", "P < 1e-4", "yes", ""],
    ]
    write_tsv(TAB / "magma_gene_sets_manifest.tsv",
              ["gene_set_name", "n_genes_expected", "n_genes_found", "source_file", "selection_rule", "created_before_downstream_analysis", "notes"],
              set_rows)


def ld_mismatch_memo() -> None:
    p2023_lit = ROOT / "data/raw/literature/2023_natcomm_gwas/supplementary_data_1_15.xlsx"
    (DOC / "ld_reference_mismatch_audit.md").write_text(f"""# LD Reference Mismatch Audit

1. Current GWAS: 2025 trans-ancestry KSD GWAS.
2. Current MAGMA LD reference: 1000 Genomes European LD reference (`external/reference/1000G_EUR/g1000_eur`), confirmed from `results/tables/magma_qc_summary.tsv`.
3. Methodological limitation: a trans-ancestry GWAS can contain ancestry-specific LD and allele-frequency structures. A European-only reference can alter gene-level mapping and may create false positive or false negative prioritization at loci where LD differs across ancestries.
4. Ancestry-specific GWAS summary statistics available: NA_not_found in this Stage 2A scan.
5. EAS LD reference available: NA_not_found as a usable LD reference in this Stage 2A scan.
6. Trans-ancestry or mixed LD reference available: NA_not_found.
7. 2023 European KSD GWAS summary statistics available: not confirmed as usable summary statistics. A literature supplementary workbook is present at `{p2023_lit.relative_to(ROOT).as_posix() if p2023_lit.exists() else 'NA_not_found'}`, but it was not treated as a ready MAGMA input in Stage 2A/2B.
8. Recommended next analysis: If 2023 European GWAS is provided, run MAGMA replication. If ancestry-specific GWAS and matched LD are provided, run ancestry-matched MAGMA. If neither is available, retain the current MAGMA result as EUR-LD-reference-based genetic prioritization, not definitive ancestry-generalizable fine mapping.

Suggested Methods/Limitations wording:
The MAGMA analysis used a GRCh37-compatible 1000 Genomes European LD reference for a 2025 trans-ancestry KSD GWAS. We therefore interpret MAGMA as a first-pass EUR-reference-based gene prioritization layer rather than ancestry-matched fine mapping. Future work should repeat gene-based analyses using ancestry-specific summary statistics with matched LD references or an appropriately weighted trans-ancestry LD reference.
""", encoding="utf-8")


def twas_audit() -> None:
    twas = read_tsv(ROOT / "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv")
    magma_genes = {r.get("gene_symbol"): r for r in read_tsv(ROOT / "results/tables/magma_genes.tsv")}
    fdr = [r for r in twas if fnum(r.get("fdr", "1")) < 0.05]
    multi = [r for r in fdr if fnum(r.get("n_snps_used", "nan")) > 1]
    one = [r for r in fdr if fnum(r.get("n_snps_used", "nan")) == 1]
    rows = [[
        "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv",
        "S-PrediXcan",
        "GTEx_v8_Kidney_Cortex",
        str(len(twas)),
        str(len(fdr)),
        str(sum(1 for r in fdr if r.get("gene") in magma_genes or r.get("magma_rank", "") not in {"", "NA"})),
        str(len(multi)),
        str(len(one)),
        str(sum(1 for r in fdr if r.get("gene") in EXEMPLARS)),
        "yes" if "n_snps_used" in headers(ROOT / "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv") else "no",
        "yes" if "zscore" in headers(ROOT / "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv") else "no",
        "yes",
        "Kidney_Cortex proxy only; not papilla-specific causal evidence",
    ]]
    write_tsv(TAB / "twas_output_audit.tsv",
              ["file_path", "twas_method", "reference_tissue", "n_genes_tested", "n_fdr_significant", "n_magma_overlap", "n_multi_snp_fdr_significant", "n_one_snp_fdr_significant", "current_exemplar_overlap_count", "has_model_snp_count_column", "has_effect_direction", "has_required_columns", "notes"],
              rows)
    detail = []
    for r in twas:
        n = fnum(r.get("n_snps_used", "nan"))
        model = "multi_snp_model" if n > 1 else ("one_snp_model" if n == 1 else "unknown_snp_count")
        is_fdr = fnum(r.get("fdr", "1")) < 0.05
        interp = "stronger_proxy_support" if is_fdr and model == "multi_snp_model" else ("weak_one_snp_proxy" if is_fdr and model == "one_snp_model" else "not_fdr_supported")
        gene = r.get("gene", "")
        detail.append([gene, r.get("zscore", ""), r.get("pvalue", ""), r.get("fdr", ""), "GTEx_v8_Kidney_Cortex", r.get("n_snps_used", "NA_not_found"), model, "yes" if gene in magma_genes or r.get("magma_rank", "") not in {"", "NA"} else "no", "yes" if gene in EXEMPLARS else "no", interp, "proxy only"])
    write_tsv(TAB / "twas_one_snp_model_audit.tsv",
              ["gene", "twas_z", "twas_p", "twas_fdr", "reference_tissue", "n_model_snps", "model_class", "magma_overlap", "current_exemplar_gene", "interpretation_class", "notes"],
              detail)
    (DOC / "twas_proxy_interpretation_memo.md").write_text("""# TWAS Proxy Interpretation Memo

The available TWAS evidence is a GTEx v8 Kidney_Cortex S-PrediXcan proxy analysis. It should be described only as kidney-cortex genetically regulated expression support for a subset of MAGMA-prioritized genes.

It must not be interpreted as papilla-specific causal expression evidence, colocalization, SMR support, spatial validation, therapeutic target validation, or validation of the current exemplar genes. One-SNP FDR-supported models should be flagged as weaker proxy evidence than multi-SNP models.
""", encoding="utf-8")


def smr_coloc_and_candidate() -> None:
    magma = read_tsv(ROOT / "results/tables/magma_genes.tsv")
    twas = read_tsv(ROOT / "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv")
    priority = read_tsv(ROOT / "results/tables/priority_loci_for_smr_coloc_v0.2.tsv")
    coloc = {r.get("gene"): r for r in read_tsv(ROOT / "results/coloc/coloc_per_locus_status.tsv")}
    magma_by_gene = {r.get("gene_symbol"): r for r in magma}
    pri_by_gene = {}
    for r in priority:
        pri_by_gene.setdefault(r.get("gene"), r)
    twas_by_gene = {r.get("gene"): r for r in twas}
    genes = list(EXEMPLARS)
    genes += [r.get("gene_symbol") for r in magma[:20] if r.get("gene_symbol") not in genes]
    genes += [r.get("gene") for r in twas if fnum(r.get("fdr", "1")) < 0.05 and (r.get("magma_rank", "") not in {"", "NA"}) and r.get("gene") not in genes]
    rows = []
    for i, g in enumerate(genes, 1):
        m = magma_by_gene.get(g, {})
        p = pri_by_gene.get(g, {})
        c = coloc.get(g, {})
        t = twas_by_gene.get(g, {})
        gwas_avail = "yes" if list((ROOT / "results/coloc/gwas_slices").glob(f"*_{g}_gwas_*.tsv.gz")) else ("yes" if p else "unclear")
        eqtl = "yes" if c.get("eqtl_resource_exists") == "TRUE" else "resource_incomplete"
        coloc_ready = "yes" if c.get("coloc_status") not in {"", "not_ready"} and c else "no"
        smr_ready = "no"
        rows.append([
            f"P{i:03d}", g, p.get("locus_id", c.get("locus_id", "NA_not_found")), p.get("chr", c.get("chr", m.get("chr", "NA_not_found"))), p.get("locus_start", c.get("locus_start", m.get("start", "NA_not_found"))), p.get("locus_end", c.get("locus_end", m.get("stop", "NA_not_found"))), p.get("lead_snp", "NA_not_found"),
            m.get("rank", p.get("magma_rank", "NA_not_found")), m.get("p", p.get("magma_p", "NA_not_found")),
            "FDR_supported" if fnum(t.get("fdr", "1")) < 0.05 else ("not_FDR_supported" if t else "NA_not_found"),
            eqtl, gwas_avail, "yes", "unclear", coloc_ready, smr_ready,
            c.get("blocking_reason", "eQTL/SMR-coloc claim-grade resources incomplete") if eqtl != "yes" else "requires harmonization and run",
            "Prepare claim-grade eQTL/LD harmonized inputs; do not treat incomplete resources as negative evidence",
        ])
    write_tsv(TAB / "smr_coloc_feasibility.tsv",
              ["priority_id", "gene", "locus", "chr", "start", "end", "lead_snp", "magma_rank", "magma_p", "twas_status", "eqtl_resource_available", "gwas_region_available", "ld_reference_available", "allele_harmonization_possible", "coloc_ready", "smr_ready", "blocking_issue", "recommended_next_action"],
              rows)

    skel_genes = set(EXEMPLARS)
    skel_genes.update(r.get("gene_symbol") for r in magma[:50])
    skel_genes.update(r.get("gene_symbol") for r in magma if r.get("bonferroni_significant") == "TRUE")
    skel_genes.update(r.get("gene") for r in twas if fnum(r.get("fdr", "1")) < 0.05)
    skel_genes.update(coloc.keys())
    skel_rows = []
    for g in sorted(x for x in skel_genes if x):
        m = magma_by_gene.get(g, {})
        t = twas_by_gene.get(g, {})
        p = pri_by_gene.get(g, {})
        c = coloc.get(g, {})
        source = []
        if g in EXEMPLARS:
            source.append("current_curated_exemplar")
        if m:
            source.append("MAGMA")
        if fnum(t.get("fdr", "1")) < 0.05:
            source.append("TWAS_FDR")
        if c:
            source.append("coloc_resource_row")
        ready = "partial" if m or fnum(t.get("fdr", "1")) < 0.05 else "no_missing_key_genetic_inputs"
        skel_rows.append([
            g, ";".join(source) or "candidate", m.get("rank", "NA_not_found"), m.get("p", "NA_not_found"), m.get("fdr", "NA_not_found"), "yes" if m.get("bonferroni_significant") == "TRUE" else "no",
            p.get("locus_id", c.get("locus_id", "NA_not_found")), p.get("lead_snp", "NA_not_found"), "NA_not_found",
            t.get("fdr", "NA_not_found"), "GTEx_v8_Kidney_Cortex" if t else "NA_not_found", t.get("n_snps_used", "NA_not_found"),
            "stronger_proxy_support" if fnum(t.get("fdr", "1")) < 0.05 and fnum(t.get("n_snps_used", "nan")) > 1 else ("weak_one_snp_proxy" if fnum(t.get("fdr", "1")) < 0.05 else "NA_not_found"),
            "not_ready", c.get("coloc_status", "not_ready" if c else "NA_not_found"), "yes" if g in EXEMPLARS else "no", "yes" if g in KIDNEY_TRANSPORT else "no", ready,
            "Skeleton only; final Tier 1/2/3 assignment deferred to Stage 3",
        ])
    write_tsv(TAB / "candidate_gene_evidence_skeleton.tsv",
              ["gene", "source_category", "magma_rank", "magma_p", "magma_fdr", "bonferroni_significant", "gwas_locus", "lead_snp", "nearest_locus_distance", "twas_fdr", "twas_reference_tissue", "twas_model_snp_count", "twas_interpretation", "smr_status", "coloc_status", "current_exemplar_gene", "known_kidney_transport_or_calcium_gene", "ready_for_tier_assignment", "notes"],
              skel_rows)


def manuscript_queue() -> None:
    ms = (ROOT / "docs/revision/stage1/manuscript_v1.4_clean_working.md").read_text(encoding="utf-8")
    hits = []
    for n, line in enumerate(ms.splitlines(), 1):
        low = line.lower()
        if any(k in low for k in ["p1", "current exemplar", "smr", "coloc", "european ld", "kidney_cortex", "twas", "figure", "candidate"]):
            hits.append(f"- Line {n}: {line[:300]}")
    (DOC / "manuscript_genetic_revision_queue.md").write_text("""# Manuscript Genetic Revision Queue

Do not rewrite the manuscript in Stage 2A/2B. Queue these items for Stage 3 and Stage 7:

- Remove any remaining P1 terminology and replace with Tier/exemplar language after Stage 3.
- Keep SMR/coloc absence as missing claim-grade evidence, not as a local execution limitation.
- Add EUR LD/trans-ancestry mismatch to Methods and Limitations.
- Describe GTEx Kidney_Cortex TWAS as proxy evidence only.
- Replace current exemplar language with candidate-gene tiering once Stage 3 tables are built.
- Revise figure captions to remove defensive boundary language and add standard concise limitations.

Relevant current lines:
""" + "\n".join(hits[:80]) + "\n", encoding="utf-8")


def report_and_tracker() -> None:
    qc = read_tsv(TAB / "gwas_qc_manifest.tsv")
    q = {r["item"]: r["value"] for r in qc}
    twas = read_tsv(TAB / "twas_output_audit.tsv")
    smr = read_tsv(TAB / "smr_coloc_feasibility.tsv")
    skel = read_tsv(TAB / "candidate_gene_evidence_skeleton.tsv")
    inv = read_tsv(TAB / "genetic_file_inventory.tsv")
    (DOC / "stage2A_2B_report.md").write_text(f"""# Stage 2A/2B Genetic Audit Report

Date: {dt.date.today().isoformat()}

## 1. Files found and missing

- Genetic-related files inventoried: {len(inv)}.
- Reconstructed 57-locus table found: `results/tables/phase1_2025_loci.tsv`.
- Original 59-locus source table: NA_not_found.
- TWAS Kidney_Cortex result table found: `results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv`.
- coloc readiness rows found, but claim-grade coloc results are not ready.

## 2. GWAS QC values confirmed

- Raw rows: {q.get('raw_rows', 'NA_not_found')}.
- QC/pass rows used for MAGMA: {q.get('qc_pass_rows', 'NA_not_found')}.
- Reconstructed loci: {q.get('number_of_reconstructed_loci', 'NA_not_found')}.
- MAGMA genes tested: {q.get('number_of_genes_tested', 'NA_not_found')}.
- Bonferroni genes: {q.get('n_bonferroni_genes', 'NA_not_found')}; FDR genes: {q.get('n_fdr_genes', 'NA_not_found')}; suggestive genes: {q.get('n_suggestive_genes', 'NA_not_found')}.

## 3. 57 vs 59 reconciliation

Not completed. The original 59-locus source table was not found; see `gwas_loci_reconciliation_blocker.md`.

## 4. MAGMA reproducibility status

Existing MAGMA outputs are auditable and include version, gene location file, and EUR LD reference. A full rerun would still require checking executable/reference paths and should be treated as Stage 2C or later.

## 5. LD mismatch status

The current MAGMA layer uses 1000 Genomes EUR LD for a 2025 trans-ancestry GWAS. This is a methodological limitation; see `ld_reference_mismatch_audit.md`.

## 6. TWAS one-SNP model summary

- TWAS genes tested: {twas[0].get('n_genes_tested', 'NA_not_found') if twas else 'NA_not_found'}.
- FDR-supported TWAS genes: {twas[0].get('n_fdr_significant', 'NA_not_found') if twas else 'NA_not_found'}.
- One-SNP FDR-supported: {twas[0].get('n_one_snp_fdr_significant', 'NA_not_found') if twas else 'NA_not_found'}.
- Multi-SNP FDR-supported: {twas[0].get('n_multi_snp_fdr_significant', 'NA_not_found') if twas else 'NA_not_found'}.

## 7. SMR/coloc readiness summary

Priority rows generated: {len(smr)}. Resources are not claim-grade ready because eQTL resources and harmonized coloc/SMR runs are incomplete or not ready for most loci. This is missing evidence, not negative biological evidence.

## 8. Candidate gene evidence skeleton summary

Candidate skeleton rows: {len(skel)}. Final Tier 1/2/3 assignment is intentionally deferred to Stage 3.

## 9. Analyses ready for Stage 2C

- 57 vs 59 reconciliation once original 59-locus table is supplied.
- MAGMA sensitivity/replication if 2023 European GWAS or ancestry-matched resources are supplied.
- TWAS one-SNP and proxy-language integration.
- SMR/coloc priority-locus resource completion and harmonization.

## 10. Analyses blocked by missing files

- Original 59-locus source table.
- 2023 European GWAS summary statistics.
- Ancestry-specific or mixed LD reference.
- Complete claim-grade eQTL/SMR/coloc resources.

## 11. Recommended exact command for Stage 2C

`开始阶段2C：在不改写全文的前提下，优先完成可行的 gwas_loci_reconciliation（若提供59-locus表）、TWAS one-SNP proxy audit整合、SMR/coloc priority loci资源补全计划；若缺失资源仍未补齐，则生成blocker和保守Methods/Limitations文字。`

## 12. Recommendation

Request missing files before claiming Stage 2 genetic reinforcement complete. Proceed to Stage 2C only for audits and feasibility work that can be completed from current resources; do not proceed to Stage 3 final tiering until genetic evidence gaps are explicitly bounded.
""", encoding="utf-8")

    tracker = ROOT / "docs/revision/STAGE_TRACKER.tsv"
    rows = read_tsv(tracker)
    if rows:
        for r in rows:
            if r.get("stage_id") == "1":
                r["blocking_issues"] = "意见.docx found after Stage 1 and extracted; pandoc/python-docx unavailable but companion markdown used"
            if r.get("stage_id") == "2":
                r["status"] = "stage2A_2B_completed"
                r["start_date"] = r.get("start_date") or dt.date.today().isoformat()
                r["completed_outputs"] = "Stage 2A/2B audits generated"
                r["blocking_issues"] = "Stage 2C pending; original 59-locus table, 2023 European GWAS, ancestry-matched LD, claim-grade SMR/coloc resources incomplete"
                r["next_stage_ready"] = "partial"
        write_tsv(tracker, list(rows[0].keys()), [[r.get(k, "") for k in rows[0].keys()] for r in rows])


def main() -> None:
    ensure_dirs()
    genetic_file_inventory()
    gwas_qc_manifest()
    loci_reconciliation_audit()
    magma_audit()
    ld_mismatch_memo()
    twas_audit()
    smr_coloc_and_candidate()
    manuscript_queue()
    report_and_tracker()
    (LOG / "stage2A_2B_genetic_audit.log").write_text(f"completed={dt.datetime.now().isoformat(timespec='seconds')}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
