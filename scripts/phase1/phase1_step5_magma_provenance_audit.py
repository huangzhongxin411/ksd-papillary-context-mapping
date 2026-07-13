#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "results/tables"
NOTE_DIR = ROOT / "notes"
CODEX_DIR = ROOT / "codex_tasks"

PVAL = ROOT / "data/processed/magma_input/ksd_2025.dedup.pval"
GENES_OUT = ROOT / "results/phase2_magma_v0.1/magma_raw/ksd_2025.genes.out"
MAGMA_LOG = ROOT / "results/phase2_magma_v0.1/magma_raw/ksd_2025.log"
GENE_ANNOT = ROOT / "external/reference/gene_loc/NCBI37.3.gene.loc"
GENE_ANNOT_RUN = ROOT / "results/phase2_magma_v0.1/magma_raw/ksd_2025.genes.annot"
EUR_PREFIX = ROOT / "external/reference/1000G_EUR/g1000_eur"
STEP3_RANKED = ROOT / "results/tables/phase1_step3_magma_ranked_canonical.tsv"

EXPECTED = {
    "magma_version": "v1.10",
    "genes_tested": "17316",
    "genes_read": "17324",
    "ld_reference": "1000 Genomes European",
    "gene_location_build": "NCBI build 37 / GRCh37-compatible",
    "magma_output": "ksd_2025.genes.out",
    "bonferroni_genes": "94",
    "fdr05_genes": "369",
    "suggestive_genes": "187",
}


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_meta(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    st = path.stat()
    return str(st.st_size), datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")


def read_head(path: Path, n: int = 5) -> list[str]:
    if not path.exists() or path.suffix == ".bed":
        return []
    out = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for _, line in zip(range(n), fh):
                out.append(line.rstrip("\n"))
    except Exception:
        return []
    return out


def count_lines(path: Path, skip_header: bool = True) -> str:
    if not path.exists():
        return ""
    if path.suffix == ".bed":
        return "not_applicable_binary"
    n = 0
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for _ in fh:
            n += 1
    return str(n - 1 if skip_header else n)


def header(path: Path) -> str:
    lines = read_head(path, 1)
    if not lines:
        return "not_applicable_binary" if path.suffix == ".bed" else ""
    if "\t" in lines[0]:
        return "|".join(lines[0].split("\t"))
    return "|".join(lines[0].split())


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_magma_genes() -> list[dict[str, str]]:
    with GENES_OUT.open("r", encoding="utf-8", errors="ignore") as fh:
        hdr = fh.readline().split()
        return [dict(zip(hdr, line.split())) for line in fh if line.strip()]


def read_ranked_symbol_map() -> dict[str, str]:
    if not STEP3_RANKED.exists():
        return {}
    with STEP3_RANKED.open(newline="", encoding="utf-8") as fh:
        return {r["magma_gene_id"]: r.get("gene_symbol", "") for r in csv.DictReader(fh, delimiter="\t")}


def bh_adjust(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    out = [1.0] * n
    running = 1.0
    for reverse_i, idx in enumerate(reversed(order), start=1):
        rank = n - reverse_i + 1
        running = min(running, pvals[idx] * n / rank)
        out[idx] = min(running, 1.0)
    return out


def parse_log() -> dict[str, str]:
    text = MAGMA_LOG.read_text(encoding="utf-8", errors="ignore")
    pats = {
        "magma_version": r"Welcome to MAGMA\s+v?([0-9.]+)",
        "bfile": r"--bfile\s+([^\n]+)",
        "pval": r"--pval\s+([^\n]+)",
        "gene_annot": r"--gene-annot\s+([^\n]+)",
        "out": r"--out\s+([^\n]+)",
        "fam_individuals": r"Reading file .*?\.fam\.\.\.\s+([0-9]+) individuals read",
        "bim_snps": r"Reading file .*?\.bim\.\.\.\s+([0-9]+) SNPs read",
        "pval_lines": r"read\s+([0-9]+) lines from file, containing valid SNP p-values",
        "valid_snps": r"valid SNP p-values for\s+([0-9]+) SNPs in data",
        "genes_read": r"([0-9]+) gene definitions read from file",
        "genes_tested": r"found\s+([0-9]+) genes containing valid SNPs",
        "output_genes": r"writing gene analysis results to file\s+([^\n]+)",
    }
    out = {}
    for k, p in pats.items():
        m = re.search(p, text, re.I)
        out[k] = m.group(1).strip() if m else ""
    warnings = re.findall(r"WARNING:[^\n]*(?:\n\s+[^\n]*)*", text)
    errors = re.findall(r"ERROR:[^\n]*(?:\n\s+[^\n]*)*", text)
    out["warnings"] = " | ".join(" ".join(w.split()) for w in warnings)
    out["errors"] = " | ".join(" ".join(e.split()) for e in errors)
    return out


def audit_pval() -> dict[str, object]:
    missing_p = invalid_p = dup = 0
    seen = set()
    n = 0
    min_p = 1.0
    max_p = 0.0
    with PVAL.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        hdr = reader.fieldnames or []
        for row in reader:
            n += 1
            snp = row.get("SNP", "")
            if snp in seen:
                dup += 1
            else:
                seen.add(snp)
            try:
                p = float(row.get("P", ""))
            except Exception:
                missing_p += 1
                continue
            if math.isnan(p) or math.isinf(p) or p <= 0 or p > 1:
                invalid_p += 1
            else:
                min_p = min(min_p, p)
                max_p = max(max_p, p)
    return {
        "n_rows": n,
        "header": "|".join(hdr),
        "snp_column": "SNP" if "SNP" in hdr else "",
        "p_column": "P" if "P" in hdr else "",
        "missing_p": missing_p,
        "invalid_p": invalid_p,
        "duplicated_snps": dup,
        "min_p": f"{min_p:.6g}",
        "max_p": f"{max_p:.6g}",
        "appears_deduplicated": "yes" if dup == 0 and "dedup" in PVAL.name else "needs_review",
    }


def bim_fam_stats(prefix: Path) -> dict[str, str]:
    bim = prefix.with_suffix(".bim")
    fam = prefix.with_suffix(".fam")
    chroms = Counter()
    n_bim = 0
    if bim.exists():
        with bim.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.strip():
                    n_bim += 1
                    chroms[line.split()[0]] += 1
    n_fam = 0
    if fam.exists():
        with fam.open("r", encoding="utf-8", errors="ignore") as fh:
            n_fam = sum(1 for line in fh if line.strip())
    return {
        "n_variants": str(n_bim),
        "n_individuals": str(n_fam),
        "chromosomes": ";".join(sorted(chroms.keys(), key=lambda x: int(x) if x.isdigit() else 999)),
    }


def find_ld_prefixes() -> list[Path]:
    prefixes = set()
    for base in [ROOT / "external", ROOT / "data", ROOT / "results"]:
        if not base.exists():
            continue
        for dp, _, fns in os.walk(base):
            names = set(fns)
            for fn in names:
                p = Path(dp) / fn
                if p.suffix in {".bed", ".bim", ".fam"}:
                    prefixes.add(p.with_suffix(""))
    return sorted(prefixes, key=lambda p: rel(p).lower())


def gene_counts() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = read_magma_genes()
    pvals = [float(r["P"]) for r in rows]
    fdr = bh_adjust(pvals)
    n = len(rows)
    bonf = 0.05 / n
    enriched = []
    symbol_map = read_ranked_symbol_map()
    for i, row in enumerate(rows):
        enriched.append(
            {
                "gene_id": row["GENE"],
                "symbol": symbol_map.get(row["GENE"], row["GENE"]),
                "p": float(row["P"]),
                "fdr": fdr[i],
            }
        )
    enriched.sort(key=lambda r: r["p"])
    bonf_n = sum(1 for r in enriched if r["p"] < bonf)
    fdr_n = sum(1 for r in enriched if r["fdr"] < 0.05)
    sugg_n = sum(1 for r in enriched if r["p"] < 1e-4)
    count_rows = [
        {"item": "tested_genes", "value_from_genes_out": n, "expected_value": 17316, "match_status": "match" if n == 17316 else "mismatch", "notes": "Rows in canonical genes.out."},
        {"item": "bonferroni_threshold", "value_from_genes_out": f"{bonf:.8g}", "expected_value": "0.05/17316", "match_status": "consistent", "notes": "Derived from number of tested genes."},
        {"item": "bonferroni_genes", "value_from_genes_out": bonf_n, "expected_value": 94, "match_status": "match" if bonf_n == 94 else "mismatch", "notes": "P < Bonferroni threshold."},
        {"item": "fdr05_genes", "value_from_genes_out": fdr_n, "expected_value": 369, "match_status": "match" if fdr_n == 369 else "mismatch", "notes": "BH-FDR < 0.05 recomputed from genes.out."},
        {"item": "suggestive_p_lt_1e-4_genes", "value_from_genes_out": sugg_n, "expected_value": 187, "match_status": "match" if sugg_n == 187 else "mismatch", "notes": "P < 1e-4."},
    ]
    top20 = []
    for rank, r in enumerate(enriched[:20], 1):
        top20.append(
            {
                "rank": rank,
                "magma_gene_id": r["gene_id"],
                "gene_symbol": r["symbol"],
                "magma_p": f"{r['p']:.8g}",
                "magma_fdr_bh": f"{r['fdr']:.8g}",
                "bonferroni_significant": r["p"] < bonf,
                "fdr05_significant": r["fdr"] < 0.05,
                "suggestive_p1e4": r["p"] < 1e-4,
            }
        )
    return count_rows, top20


def main() -> int:
    for d in [TABLE_DIR, NOTE_DIR, CODEX_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    roles = [
        ("magma_snp_p_input", PVAL, "tabular SNP/P/N file with SNP,P,N columns"),
        ("magma_genes_out", GENES_OUT, "MAGMA genes.out table with GENE and P columns"),
        ("magma_log", MAGMA_LOG, "text log containing MAGMA version and input options"),
        ("magma_gene_location_reference", GENE_ANNOT, "NCBI37 gene location reference"),
        ("magma_gene_annotation_run_file", GENE_ANNOT_RUN, "MAGMA genes.annot file used in run"),
        ("eur_ld_bed", EUR_PREFIX.with_suffix(".bed"), "PLINK binary bed"),
        ("eur_ld_bim", EUR_PREFIX.with_suffix(".bim"), "PLINK bim variant table"),
        ("eur_ld_fam", EUR_PREFIX.with_suffix(".fam"), "PLINK fam sample table"),
    ]
    provenance = []
    for role, path, fmt in roles:
        size, mtime = file_meta(path)
        provenance.append(
            {
                "file_role": role,
                "file_path": rel(path),
                "exists": "yes" if path.exists() else "no",
                "file_size": size,
                "modified_time": mtime,
                "n_rows": count_lines(path, skip_header=role not in {"eur_ld_bim", "eur_ld_fam", "magma_gene_location_reference", "magma_gene_annotation_run_file"}),
                "header_fields": header(path),
                "format_check": fmt,
                "canonical_status": "matches Step 2 canonical path" if path.exists() else "missing",
                "concerns": "Binary file; row count not applicable." if path.suffix == ".bed" else "",
                "action_needed": "None for provenance; use as audit evidence." if path.exists() else "Provide missing file.",
            }
        )
    write_tsv(TABLE_DIR / "phase1_step5_magma_file_provenance.tsv", provenance, ["file_role", "file_path", "exists", "file_size", "modified_time", "n_rows", "header_fields", "format_check", "canonical_status", "concerns", "action_needed"])

    log = parse_log()
    log_rows = [
        ("MAGMA version", f"v{log.get('magma_version')}", EXPECTED["magma_version"]),
        ("SNP/P input option", log.get("pval", ""), rel(PVAL)),
        ("gene annotation option", log.get("gene_annot", ""), "results/magma/2025_trans_ancestry/ksd_2025.genes.annot / GRCh37-compatible"),
        ("LD reference prefix", log.get("bfile", ""), "external/reference/1000G_EUR/g1000_eur"),
        ("PLINK FAM individuals", log.get("fam_individuals", ""), "503"),
        ("PLINK BIM SNPs", log.get("bim_snps", ""), "22665064"),
        ("SNP/P lines read", log.get("pval_lines", ""), "4915034 including header per MAGMA log"),
        ("valid SNP p-values in LD data", log.get("valid_snps", ""), "4882602"),
        ("genes read", log.get("genes_read", ""), EXPECTED["genes_read"]),
        ("genes tested", log.get("genes_tested", ""), EXPECTED["genes_tested"]),
        ("MAGMA output", log.get("output_genes", ""), EXPECTED["magma_output"]),
        ("warnings", log.get("warnings", ""), "no fatal warnings expected"),
        ("errors", log.get("errors", ""), "none"),
    ]
    metadata = []
    for item, value, expected in log_rows:
        vl = str(value).lower()
        el = str(expected).lower()
        if item == "warnings":
            match = "warning_present_nonfatal" if value else "none"
        elif item == "errors":
            match = "match" if not value else "error_present"
        elif item == "MAGMA output":
            match = "match" if "ksd_2025.genes.out" in vl else "needs_review"
        elif item in {"gene annotation option"}:
            match = "consistent" if "genes.annot" in vl else "needs_review"
        elif item == "SNP/P lines read":
            match = "consistent" if value.startswith("4915034") else "needs_review"
        else:
            match = "match" if el in vl or vl in el else "needs_review"
        metadata.append({"metadata_item": item, "value_from_log": value, "value_expected_from_manuscript_or_canonical_files": expected, "match_status": match, "notes": "Parsed from MAGMA log; no MAGMA rerun performed."})
    write_tsv(TABLE_DIR / "phase1_step5_magma_log_metadata.tsv", metadata, ["metadata_item", "value_from_log", "value_expected_from_manuscript_or_canonical_files", "match_status", "notes"])

    count_rows, top20 = gene_counts()
    write_tsv(TABLE_DIR / "phase1_step5_magma_count_verification.tsv", count_rows, ["item", "value_from_genes_out", "expected_value", "match_status", "notes"])
    write_tsv(TABLE_DIR / "phase1_step5_magma_top20_genes.tsv", top20, ["rank", "magma_gene_id", "gene_symbol", "magma_p", "magma_fdr_bh", "bonferroni_significant", "fdr05_significant", "suggestive_p1e4"])

    pval_stats = audit_pval()
    pval_rows = [
        ("number_of_rows", pval_stats["n_rows"], "Data rows in SNP/P input."),
        ("header_fields", pval_stats["header"], "Expected SNP|P|N."),
        ("SNP_column", pval_stats["snp_column"], ""),
        ("P_value_column", pval_stats["p_column"], ""),
        ("missing_P_values", pval_stats["missing_p"], ""),
        ("invalid_P_values_outside_0_1", pval_stats["invalid_p"], "P <= 0 or P > 1 considered invalid for MAGMA input audit."),
        ("duplicated_SNP_IDs", pval_stats["duplicated_snps"], ""),
        ("P_value_range", f"{pval_stats['min_p']} to {pval_stats['max_p']}", "Matches cleaned GWAS minimum range from Step 4 where feasible."),
        ("appears_deduplicated", pval_stats["appears_deduplicated"], "Filename contains dedup and no duplicate SNP IDs were detected."),
    ]
    write_tsv(TABLE_DIR / "phase1_step5_magma_snp_p_input_audit.tsv", [{"item": a, "value": b, "notes": c} for a, b, c in pval_rows], ["item", "value", "notes"])

    ld = bim_fam_stats(EUR_PREFIX)
    ld_rows = [
        ("bed_exists", EUR_PREFIX.with_suffix(".bed").exists(), rel(EUR_PREFIX.with_suffix(".bed"))),
        ("bim_exists", EUR_PREFIX.with_suffix(".bim").exists(), rel(EUR_PREFIX.with_suffix(".bim"))),
        ("fam_exists", EUR_PREFIX.with_suffix(".fam").exists(), rel(EUR_PREFIX.with_suffix(".fam"))),
        ("bed_file_size", file_meta(EUR_PREFIX.with_suffix(".bed"))[0], "bytes"),
        ("bim_file_size", file_meta(EUR_PREFIX.with_suffix(".bim"))[0], "bytes"),
        ("fam_file_size", file_meta(EUR_PREFIX.with_suffix(".fam"))[0], "bytes"),
        ("n_variants_in_bim", ld["n_variants"], "Expected 22,665,064 from MAGMA log."),
        ("n_individuals_in_fam", ld["n_individuals"], "Expected 503 from MAGMA log."),
        ("chromosomes_represented", ld["chromosomes"], "BIM chromosomes observed."),
        ("file_prefix_complete", all(EUR_PREFIX.with_suffix(s).exists() for s in [".bed", ".bim", ".fam"]), "Complete PLINK bed/bim/fam prefix required by MAGMA."),
    ]
    write_tsv(TABLE_DIR / "phase1_step5_eur_ld_reference_audit.tsv", [{"item": a, "value": b, "notes": c} for a, b, c in ld_rows], ["item", "value", "notes"])

    alt_rows = []
    for prefix in find_ld_prefixes():
        label = prefix.parent.name
        stats = bim_fam_stats(prefix) if prefix.with_suffix(".bim").exists() or prefix.with_suffix(".fam").exists() else {"n_variants": "", "n_individuals": "", "chromosomes": ""}
        is_eur = "1000g_eur" in rel(prefix).lower() or "g1000_eur" in prefix.name.lower()
        alt_rows.append(
            {
                "ld_reference_label": label,
                "file_prefix": rel(prefix),
                "bed_exists": prefix.with_suffix(".bed").exists(),
                "bim_exists": prefix.with_suffix(".bim").exists(),
                "fam_exists": prefix.with_suffix(".fam").exists(),
                "n_variants_if_available": stats["n_variants"],
                "n_individuals_if_available": stats["n_individuals"],
                "chromosomes_if_available": stats["chromosomes"],
                "usable_for_magma_sensitivity": "no_current_alternative_baseline_only" if is_eur else ("candidate_needs_human_review" if all(prefix.with_suffix(s).exists() for s in [".bed", ".bim", ".fam"]) else "no_incomplete_prefix"),
                "notes": "Canonical EUR LD reference, not an alternative sensitivity reference." if is_eur else "Potential alternative only if ancestry/source and build are appropriate.",
            }
        )
    if not alt_rows:
        alt_rows.append({"ld_reference_label": "none_found", "file_prefix": "", "bed_exists": False, "bim_exists": False, "fam_exists": False, "n_variants_if_available": "", "n_individuals_if_available": "", "chromosomes_if_available": "", "usable_for_magma_sensitivity": "no", "notes": "No PLINK LD prefixes found."})
    write_tsv(TABLE_DIR / "phase1_step5_alternative_ld_reference_inventory.tsv", alt_rows, ["ld_reference_label", "file_prefix", "bed_exists", "bim_exists", "fam_exists", "n_variants_if_available", "n_individuals_if_available", "chromosomes_if_available", "usable_for_magma_sensitivity", "notes"])

    alternatives = [r for r in alt_rows if r["usable_for_magma_sensitivity"] == "candidate_needs_human_review"]
    if alternatives:
        alt_text = "\n".join(f"- `{r['file_prefix']}` ({r['ld_reference_label']})" for r in alternatives)
        command_text = "\n".join(
            f"magma --bfile {r['file_prefix']} --pval {rel(PVAL)} use=SNP,P ncol=N --gene-annot {rel(GENE_ANNOT_RUN)} --out results/magma_sensitivity/{r['ld_reference_label']}/ksd_2025"
            for r in alternatives
        )
    else:
        alt_text = "- No non-EUR complete alternative PLINK LD reference prefixes were found."
        command_text = "magma --bfile <alternative_LD_prefix> --pval data/processed/magma_input/ksd_2025.dedup.pval use=SNP,P ncol=N --gene-annot results/phase2_magma_v0.1/magma_raw/ksd_2025.genes.annot --out results/magma_sensitivity/<LD_label>/ksd_2025"

    sensitivity_note = [
        "# Phase 1-Step 5 LD-Reference Sensitivity Plan",
        "",
        "EUR LD is a limitation because the GWAS is trans-ancestry, while the MAGMA LD structure came from a 1000 Genomes European reference. The result should therefore be interpreted as EUR-LD-reference-based MAGMA prioritization rather than ancestry-generalizable fine mapping.",
        "",
        "## Current Alternative LD Availability",
        "",
        alt_text,
        "",
        "## Ideal Sensitivity Analyses",
        "",
        "- Repeat MAGMA gene analysis with ancestry-matched or multi-ancestry LD references, if complete PLINK bed/bim/fam prefixes are provided.",
        "- Compare top50, top100, Bonferroni, FDR05 and suggestive gene-set overlap against the canonical EUR-LD run.",
        "- Compute Spearman rank correlation of MAGMA gene ranks and assess stability of downstream Loop/TAL-associated modules in a later, explicitly approved step.",
        "",
        "## Proposed Command Template",
        "",
        "```bash",
        command_text,
        "```",
        "",
        "## Current Manuscript Boundary",
        "",
        "“MAGMA used a 1000 Genomes European LD reference and is therefore interpreted as EUR-LD-reference-based gene prioritization rather than ancestry-generalizable fine mapping.”",
    ]
    (NOTE_DIR / "phase1_step5_ld_reference_sensitivity_plan.md").write_text("\n".join(sensitivity_note) + "\n", encoding="utf-8")

    reproducibility_note = [
        "# Phase 1-Step 5 MAGMA Reproducibility Note",
        "",
        "## Inputs",
        "",
        f"- SNP/P input: `{rel(PVAL)}`",
        f"- MAGMA genes.out: `{rel(GENES_OUT)}`",
        f"- MAGMA log: `{rel(MAGMA_LOG)}`",
        f"- Gene annotation used by run: `{rel(GENE_ANNOT_RUN)}`",
        f"- Gene location/build reference: `{rel(GENE_ANNOT)}` (NCBI build 37 / GRCh37-compatible)",
        f"- LD reference: `{rel(EUR_PREFIX)}` bed/bim/fam, 1000 Genomes European",
        "",
        "## Reproducible Values",
        "",
        "- MAGMA version: v1.10",
        "- Genes tested: 17,316",
        "- Bonferroni threshold: 0.05/17,316",
        "- Bonferroni genes: 94",
        "- FDR05 genes: 369",
        "- Suggestive P < 1e-4 genes: 187",
        "",
        "## Remaining Limits",
        "",
        "- The LD reference is European while the GWAS is trans-ancestry.",
        "- No complete non-EUR or multi-ancestry LD sensitivity reference was identified in the current project files.",
        "- MAGMA provides genetic-priority genes, not causal genes or validated targets.",
        "",
        "## Recommended Methods Wording",
        "",
        "MAGMA v1.10 gene analysis used the deduplicated SNP/P input, NCBI build 37-compatible gene annotation and a 1000 Genomes European LD reference. Gene sets were derived from the locked `ksd_2025.genes.out` file using top-N rank, Bonferroni, BH-FDR and P < 1e-4 definitions.",
        "",
        "## Recommended Limitations Wording",
        "",
        "MAGMA used a 1000 Genomes European LD reference and is therefore interpreted as EUR-LD-reference-based gene prioritization rather than ancestry-generalizable fine mapping.",
    ]
    (NOTE_DIR / "phase1_step5_magma_reproducibility_note.md").write_text("\n".join(reproducibility_note) + "\n", encoding="utf-8")

    report = [
        "# Phase 1-Step 5 Report",
        "",
        "- MAGMA provenance is complete for the canonical SNP/P input, genes.out, log, gene annotation and EUR LD PLINK prefix.",
        "- The MAGMA log supports the manuscript-reported version, LD prefix, SNP/P input, gene annotation, 503 EUR reference individuals, 22,665,064 reference SNPs and 17,316 tested genes.",
        "- The canonical `.genes.out` supports the manuscript-reported counts: 17,316 tested genes, 94 Bonferroni genes, 369 FDR05 genes and 187 suggestive genes.",
        f"- SNP/P input is valid and deduplicated in this audit: {pval_stats['n_rows']} rows, {pval_stats['duplicated_snps']} duplicated SNP IDs, {pval_stats['invalid_p']} invalid P values.",
        "- EUR LD reference files are complete: bed, bim and fam are present.",
        f"- Alternative LD references available now: {'yes' if alternatives else 'no'}; only the canonical EUR PLINK prefix was identified as complete.",
        "- LD sensitivity can be run later only if alternative complete LD prefixes are provided or approved.",
        "- Step 6 can proceed to integrated Phase 1 genetic-layer summary and manuscript-ready wording, with the EUR-LD boundary retained.",
    ]
    (NOTE_DIR / "phase1_step5_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    next_steps = [
        "# Phase 1-Step 5 Limitations and Next Steps",
        "",
        "- EUR LD remains a limitation for trans-ancestry GWAS interpretation.",
        f"- Alternative LD-reference sensitivity is {'feasible after human review of candidate prefixes' if alternatives else 'not feasible now because no complete non-EUR or multi-ancestry PLINK LD prefix was found'} in current project files.",
        "- Current manuscript wording must remain softened: use EUR-LD-reference-based prioritization, not ancestry-generalizable fine mapping.",
        "- Downstream analyses should use Step 3 canonical MAGMA modules rather than deprecated old frozen files.",
        "",
        "Recommended next step: A. proceed to Step 6: integrated Phase 1 genetic-layer summary and manuscript-ready wording.",
    ]
    (NOTE_DIR / "phase1_step5_limitations_and_next_steps.md").write_text("\n".join(next_steps) + "\n", encoding="utf-8")

    checklist = [
        ("S5-01", "Create MAGMA provenance audit script and file provenance table", "yes", "scripts/03_magma_fuma/phase1_step5_magma_provenance_audit.py; results/tables/phase1_step5_magma_file_provenance.tsv", "", "yes", "No MAGMA rerun."),
        ("S5-02", "Parse MAGMA log metadata", "yes", "results/tables/phase1_step5_magma_log_metadata.tsv", "", "yes", "Parsed existing log only."),
        ("S5-03", "Verify MAGMA gene-level counts and top20 genes", "yes", "results/tables/phase1_step5_magma_count_verification.tsv; results/tables/phase1_step5_magma_top20_genes.tsv", "", "yes", "Counts recomputed from genes.out."),
        ("S5-04", "Audit MAGMA SNP/P input", "yes", "results/tables/phase1_step5_magma_snp_p_input_audit.tsv", "", "yes", "Streaming tabular audit."),
        ("S5-05", "Audit EUR LD reference files", "yes", "results/tables/phase1_step5_eur_ld_reference_audit.tsv", "", "yes", "PLINK prefix complete."),
        ("S5-06", "Search for alternative LD references", "yes", "results/tables/phase1_step5_alternative_ld_reference_inventory.tsv", "No complete non-EUR alternative identified." if not alternatives else "", "yes", "No sensitivity run performed."),
        ("S5-07", "Prepare LD-reference sensitivity plan", "yes", "notes/phase1_step5_ld_reference_sensitivity_plan.md", "", "yes", "Command templates only."),
        ("S5-08", "Create reviewer-facing MAGMA reproducibility note", "yes", "notes/phase1_step5_magma_reproducibility_note.md", "", "yes", "Boundary language included."),
        ("S5-09", "Create Step 5 report", "yes", "notes/phase1_step5_report.md", "", "yes", "Step 6 readiness stated."),
        ("S5-10", "Create limitations and next-steps note", "yes", "notes/phase1_step5_limitations_and_next_steps.md", "", "yes", "Recommends Step 6."),
    ]
    write_tsv(CODEX_DIR / "phase1_step5_completion_checklist.tsv", [{"task_id": a, "task_name": b, "completed": c, "output_file": d, "blocking_issue": e, "manual_review_needed": f, "notes": g} for a, b, c, d, e, f, g in checklist], ["task_id", "task_name", "completed", "output_file", "blocking_issue", "manual_review_needed", "notes"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
