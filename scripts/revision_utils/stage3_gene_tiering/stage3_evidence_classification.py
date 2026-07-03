#!/usr/bin/env python3
"""Stage 3 conservative evidence-class system.

Uses the frozen Stage 3 input and does not incorporate snRNA, spatial, GSE73680,
SMR/coloc support, or causal claims.
"""

from __future__ import annotations

import csv
import datetime as dt
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs/revision/stage3_gene_tiering"
TAB = ROOT / "results/tables/revision/stage3_gene_tiering"
FIG = ROOT / "results/figures/revision/stage3_gene_tiering"
LOG = ROOT / "logs/revision/stage3_gene_tiering"
SCRIPT = ROOT / "scripts/revision_utils/stage3_gene_tiering"

FROZEN = TAB / "stage3_input_candidate_gene_evidence_skeleton_frozen.tsv"
CHECKSUM = TAB / "stage3_input_checksum.md"
MANUSCRIPT = ROOT / "docs/revision/stage1/manuscript_v1.4_clean_working.md"
TRACKER = ROOT / "docs/revision/STAGE_TRACKER.tsv"

EXEMPLARS = {
    "UMOD": "TAL identity / uromodulin biology",
    "CLDN10": "epithelial transport / TAL tight junction",
    "CLDN14": "calcium/ion handling / claudin biology",
    "CASR": "calcium sensing",
    "HIBADH": "supporting MAGMA-associated context gene",
    "PKD2": "broader renal epithelial context",
}

DISALLOWED = "causal gene; validated disease gene; papilla-specific TWAS gene; SMR-supported gene; coloc-supported gene; therapeutic target"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def fnum(value: str, default: float | None = None) -> float | None:
    try:
        if value in {"", "NA", "NA_not_found", "not_available"}:
            return default
        return float(value)
    except Exception:
        return default


def has_magma(row: dict[str, str]) -> bool:
    if "MAGMA" in row.get("source_category", ""):
        return True
    return fnum(row.get("magma_p", ""), None) is not None or fnum(row.get("magma_rank", ""), None) is not None


def twas_fdr_supported(row: dict[str, str]) -> bool:
    fdr = fnum(row.get("twas_fdr", ""), None)
    return fdr is not None and fdr < 0.05


def twas_proxy_class(row: dict[str, str]) -> str:
    if not twas_fdr_supported(row):
        return "not_twas_fdr_supported" if row.get("twas_fdr") not in {"", "NA_not_found", "not_available"} else "not_available"
    n = fnum(row.get("twas_model_snp_count", ""), None)
    if n is not None and n > 1:
        return "multi_snp_proxy"
    if n == 1:
        return "one_snp_proxy"
    return "not_available"


def assign_class(row: dict[str, str]) -> tuple[str, str, str, str]:
    magma = has_magma(row)
    twas = twas_fdr_supported(row)
    n = fnum(row.get("twas_model_snp_count", ""), None)
    bonf = row.get("bonferroni_significant", "").lower() == "yes"
    if magma and twas and n is not None and n > 1:
        return (
            "Class B1",
            "multi_snp_kidney_cortex_twas_proxy",
            "MAGMA-prioritized genes with multi-SNP Kidney_Cortex TWAS proxy support",
            "MAGMA-prioritized gene with multi-SNP Kidney_Cortex TWAS proxy support",
        )
    if magma and twas and n == 1:
        return (
            "Class B2",
            "one_snp_kidney_cortex_twas_proxy",
            "MAGMA-prioritized genes with one-SNP Kidney_Cortex TWAS proxy signal",
            "MAGMA-prioritized gene with one-SNP Kidney_Cortex TWAS proxy signal",
        )
    if magma and bonf:
        return ("Class A", "magma_bonferroni", "MAGMA-prioritized genes", "MAGMA-prioritized gene")
    if (not magma) and twas:
        return (
            "Class C-TWAS-proxy-only",
            "twas_fdr_without_magma_support",
            "TWAS FDR-supported kidney-cortex proxy-only genes",
            "Kidney_Cortex proxy TWAS signal only",
        )
    return ("Class C", "contextual_or_module_gene", "Contextual or module genes", "contextual module gene")


def definitions() -> None:
    (DOC / "evidence_class_definition_v0.1.md").write_text("""# Evidence Class Definition v0.1

This Stage 3 evidence-class system replaces the subjective P1 framework. It is conservative and uses currently available genetic and Kidney_Cortex proxy TWAS evidence only. It does not assign causal status, SMR support, coloc support, papilla-specific TWAS support, therapeutic target status, or validated disease-gene status.

## Class A: MAGMA-prioritized genes

Criteria:
- MAGMA Bonferroni significant; or MAGMA top-ranked gene if explicitly listed in the frozen skeleton.
- No requirement for TWAS, SMR, or coloc support.

Interpretation:
- Can be described as genetically prioritized by MAGMA.
- Cannot be described as causal or validated.

## Class B1: MAGMA-prioritized genes with multi-SNP Kidney_Cortex TWAS proxy support

Criteria:
- MAGMA-supported.
- TWAS FDR-supported in GTEx Kidney_Cortex.
- TWAS model supported by more than one SNP.

Interpretation:
- Stronger kidney-cortex proxy support than one-SNP models.
- Not papilla-specific causal expression evidence.

## Class B2: MAGMA-prioritized genes with one-SNP Kidney_Cortex TWAS proxy signal

Criteria:
- MAGMA-supported.
- TWAS FDR-supported in GTEx Kidney_Cortex.
- TWAS model supported by one SNP.

Interpretation:
- Weaker proxy signal.
- Should be interpreted cautiously and not used as robust TWAS evidence.

## Class C: Contextual or module genes

Criteria:
- MAGMA FDR, MAGMA suggestive, TWAS FDR-only, or module-included genes that do not meet Class A/B criteria.

Interpretation:
- Useful for module and pathway context.
- Not individually prioritized disease genes.

## Class C-TWAS-proxy-only

Criteria:
- TWAS FDR-supported but no MAGMA support.

Interpretation:
- Supplementary kidney-cortex proxy signal only.
- Not a main candidate class.

## Curated biological exemplar panel

Genes: UMOD, CASR, CLDN14, CLDN10, HIBADH, PKD2.

Definition:
A literature- and biology-informed explanatory panel used to interpret TAL identity, epithelial transport, calcium/ion handling, supporting context, and broader epithelial biology. This panel is not an unbiased P1 output and is not a validated disease-gene panel.
""", encoding="utf-8")


def class_table(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows = []
    header = [
        "gene", "assigned_class", "class_subtype", "class_label", "source_category",
        "magma_rank", "magma_p", "magma_fdr", "bonferroni_significant",
        "gwas_locus", "lead_snp", "twas_fdr", "twas_reference_tissue",
        "twas_model_snp_count", "twas_proxy_class", "current_curated_exemplar",
        "known_kidney_transport_or_calcium_gene", "smr_status", "coloc_status",
        "allowed_claim", "disallowed_claim", "notes",
    ]
    table = []
    for row in rows:
        assigned, subtype, label, allowed = assign_class(row)
        gene = row["gene"]
        exemplar = "yes" if gene in EXEMPLARS else "no"
        if exemplar == "yes":
            allowed = allowed + "; curated biological exemplar gene"
        notes = "No SMR/coloc support assigned; no causal or validated-gene interpretation."
        table.append([
            gene, assigned, subtype, label, row.get("source_category", ""),
            row.get("magma_rank", "NA_not_found"), row.get("magma_p", "NA_not_found"),
            row.get("magma_fdr", "NA_not_found"), row.get("bonferroni_significant", "no"),
            row.get("gwas_locus", "NA_not_found"), row.get("lead_snp", "NA_not_found"),
            row.get("twas_fdr", "NA_not_found"), row.get("twas_reference_tissue", "NA_not_found"),
            row.get("twas_model_snp_count", "NA_not_found"), twas_proxy_class(row),
            exemplar, row.get("known_kidney_transport_or_calcium_gene", "no"),
            "not_ready", "not_ready", allowed, DISALLOWED, notes,
        ])
        out_rows.append(dict(zip(header, table[-1])))
    write_tsv(TAB / "candidate_gene_evidence_classes_v0.1.tsv", header, table)
    return out_rows


def exemplar_panel(class_rows: list[dict[str, str]]) -> None:
    by_gene = {r["gene"]: r for r in class_rows}
    rows = []
    for gene, role in EXEMPLARS.items():
        r = by_gene.get(gene, {})
        twas_status = "FDR_supported" if r.get("twas_proxy_class") in {"multi_snp_proxy", "one_snp_proxy"} else "not_FDR_supported_or_not_available"
        rows.append([
            gene,
            r.get("assigned_class", "NA_not_found"),
            role,
            "Included as a curated biological exemplar for role-spectrum interpretation, not as an unbiased P1 output.",
            r.get("magma_rank", "NA_not_found"),
            r.get("magma_p", "NA_not_found"),
            r.get("magma_fdr", "NA_not_found"),
            r.get("bonferroni_significant", "NA_not_found"),
            twas_status,
            r.get("twas_proxy_class", "not_available"),
            "not_ready; no claim-grade SMR/coloc support assigned",
            "curated biological exemplar for role-spectrum interpretation",
            "unbiased P1 output; validated disease gene; causal gene; SMR-supported gene; coloc-supported gene; therapeutic target",
            "curated biological exemplar for role-spectrum interpretation",
        ])
    write_tsv(
        TAB / "curated_exemplar_panel_v0.1.tsv",
        [
            "gene", "assigned_class", "biological_role_label", "why_included",
            "magma_rank", "magma_p", "magma_fdr", "bonferroni_significant",
            "twas_status", "twas_proxy_class", "smr_coloc_status",
            "interpretation_allowed", "interpretation_not_allowed",
            "recommended_manuscript_role",
        ],
        rows,
    )


def summary_counts(class_rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    order = ["Class A", "Class B1", "Class B2", "Class C", "Class C-TWAS-proxy-only"]
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in class_rows:
        c = r["assigned_class"]
        counts[c]["n_genes"] += 1
        counts[c]["n_bonferroni"] += int(r.get("bonferroni_significant") == "yes")
        counts[c]["n_twas_fdr"] += int(r.get("twas_proxy_class") in {"multi_snp_proxy", "one_snp_proxy"})
        counts[c]["n_multi_snp_twas"] += int(r.get("twas_proxy_class") == "multi_snp_proxy")
        counts[c]["n_one_snp_twas"] += int(r.get("twas_proxy_class") == "one_snp_proxy")
        counts[c]["n_curated_exemplar"] += int(r.get("current_curated_exemplar") == "yes")
    rows = []
    interp = {
        "Class A": "MAGMA-prioritized genes without FDR-supported Kidney_Cortex proxy TWAS in this frozen table.",
        "Class B1": "MAGMA-prioritized genes with stronger multi-SNP Kidney_Cortex proxy TWAS support.",
        "Class B2": "MAGMA-prioritized genes with weaker one-SNP Kidney_Cortex proxy TWAS signal.",
        "Class C": "Contextual/module genes or MAGMA-supported genes not meeting Class A/B definitions.",
        "Class C-TWAS-proxy-only": "TWAS FDR-supported kidney-cortex proxy signal without MAGMA support.",
        "Curated exemplar panel": "Biology-informed explanatory panel; not an unbiased P1 or validated disease-gene panel.",
    }
    for c in order:
        rows.append([
            c, str(counts[c]["n_genes"]), str(counts[c]["n_bonferroni"]),
            str(counts[c]["n_twas_fdr"]), str(counts[c]["n_multi_snp_twas"]),
            str(counts[c]["n_one_snp_twas"]), str(counts[c]["n_curated_exemplar"]),
            interp[c],
        ])
    rows.append([
        "Curated exemplar panel", "6",
        str(sum(1 for r in class_rows if r["gene"] in EXEMPLARS and r.get("bonferroni_significant") == "yes")),
        str(sum(1 for r in class_rows if r["gene"] in EXEMPLARS and r.get("twas_proxy_class") in {"multi_snp_proxy", "one_snp_proxy"})),
        str(sum(1 for r in class_rows if r["gene"] in EXEMPLARS and r.get("twas_proxy_class") == "multi_snp_proxy")),
        str(sum(1 for r in class_rows if r["gene"] in EXEMPLARS and r.get("twas_proxy_class") == "one_snp_proxy")),
        "6", interp["Curated exemplar panel"],
    ])
    write_tsv(
        TAB / "evidence_class_summary_counts.tsv",
        [
            "category", "n_genes", "n_bonferroni", "n_twas_fdr",
            "n_multi_snp_twas", "n_one_snp_twas", "n_curated_exemplar",
            "interpretation",
        ],
        rows,
    )
    return counts


def replacement_text(counts: dict[str, dict[str, int]]) -> None:
    def n(c: str, key: str = "n_genes") -> int:
        return counts[c][key]
    (DOC / "manuscript_replacement_text_stage3.md").write_text(f"""# Stage 3 Manuscript Replacement Text

## A. Methods Replacement

### Candidate gene evidence classification and curated exemplar panel

Candidate genes were classified using the frozen Stage 3 evidence skeleton derived from auditable MAGMA gene-based results and GTEx v8 Kidney_Cortex S-PrediXcan proxy TWAS. Class A genes were MAGMA-prioritized genes without requiring TWAS, SMR, or coloc support. Class B1 genes were MAGMA-supported genes with FDR-supported Kidney_Cortex proxy TWAS models supported by more than one SNP, whereas Class B2 genes were MAGMA-supported genes with FDR-supported one-SNP Kidney_Cortex proxy TWAS signals. Class C genes captured contextual or module genes that did not meet Class A/B criteria, and Class C-TWAS-proxy-only genes captured TWAS FDR-supported kidney-cortex proxy signals without MAGMA support.

SMR/coloc support was not assigned because no claim-grade SMR/coloc resources were ready for the audited priority loci. UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were treated as curated biological exemplar genes for role-spectrum interpretation, not as an unbiased P1 output or a validated disease-gene panel.

## B. Results Replacement

### Evidence-classified candidate genes and curated exemplars support a TAL-to-broader epithelial interpretation

The conservative evidence-class system assigned {n('Class A')} genes to Class A, {n('Class B1')} genes to Class B1, {n('Class B2')} genes to Class B2, {n('Class C')} genes to Class C and {n('Class C-TWAS-proxy-only')} genes to Class C-TWAS-proxy-only. The six curated exemplar genes were all retained as biology-informed exemplars and were MAGMA-prioritized in the frozen evidence skeleton, but they were not assigned TWAS-, SMR- or coloc-supported status in this version.

The curated exemplar panel illustrates biological role diversity across TAL identity, epithelial transport, calcium/ion handling, supporting MAGMA-associated context and broader renal epithelial biology. This panel is used for interpretation of the candidate-gene role spectrum and does not validate individual disease genes.

## C. Figure 3 Replacement Concept

Redesign Figure 3 as a conservative evidence-class and exemplar-role figure:

- Panel A: evidence-class summary bar/table.
- Panel B: curated exemplar role spectrum.
- Panel C: exemplar gene MAGMA/TWAS evidence strip.
- Panel D: allowed versus disallowed claim box.

Do not include new snRNA expression plots here unless already present; deeper donor-level snRNA analysis belongs to Stage 4.

## D. Discussion Replacement Paragraph

We removed the P1 terminology because it could imply a subjective high-confidence or causal candidate-gene set. The revised evidence-class framework is intentionally conservative: it separates MAGMA-prioritized genes, kidney-cortex proxy TWAS signals and curated biological exemplars without assigning causal validation. No evidence class should be interpreted as causal, papilla-specific TWAS-supported, SMR-supported, coloc-supported or therapeutically validated. Future papilla-relevant eQTL resources, claim-grade SMR/coloc analyses, ancestry-matched genetic sensitivity analyses and experimental perturbation will be required to refine these candidates into causal disease mechanisms.
""", encoding="utf-8")


def p1_audit() -> None:
    patterns = [
        "P1", "priority-1", "P1 candidate", "P1 gene", "P1 validation",
        "current exemplar genes", "current exemplar", "single-gene exemplar",
        "core TAL-associated KSD candidates", "high-confidence", "causal gene",
        "validated gene",
    ]
    regex = re.compile("|".join(re.escape(p) for p in patterns), re.I)
    rows = []
    for i, line in enumerate(MANUSCRIPT.read_text(encoding="utf-8").splitlines(), 1):
        if regex.search(line):
            rec = "Replace with evidence-class genes or curated biological exemplar genes, depending on context."
            pri = "high" if re.search(r"P1|priority-1|high-confidence|causal gene|validated gene", line, re.I) else "medium"
            rows.append([str(i), line.strip()[:700], "yes", rec, pri])
    write_tsv(
        DOC / "p1_terminology_replacement_audit.tsv",
        ["line_or_section", "current_text", "replacement_needed", "recommended_replacement", "priority"],
        rows,
    )


def report(counts: dict[str, dict[str, int]]) -> None:
    checksum_text = CHECKSUM.read_text(encoding="utf-8") if CHECKSUM.exists() else "NA_not_found"
    checksum = "NA_not_found"
    m = re.search(r"SHA256: `([^`]+)`", checksum_text)
    if m:
        checksum = m.group(1)
    (DOC / "stage3_report.md").write_text(f"""# Stage 3 Report

Date: {dt.date.today().isoformat()}

## Input

- Frozen input file: `results/tables/revision/stage3_gene_tiering/stage3_input_candidate_gene_evidence_skeleton_frozen.tsv`
- Input SHA256: `{checksum}`

## Evidence Class Counts

- Class A: {counts['Class A']['n_genes']} genes.
- Class B1: {counts['Class B1']['n_genes']} genes.
- Class B2: {counts['Class B2']['n_genes']} genes.
- Class C: {counts['Class C']['n_genes']} genes.
- Class C-TWAS-proxy-only: {counts['Class C-TWAS-proxy-only']['n_genes']} genes.

## Curated Exemplar Panel

Curated exemplar genes: 6 genes: UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2. These are retained as curated biological exemplars, not as P1 candidates, unbiased priority-1 genes, validated disease genes or causal genes.

## TWAS Proxy Summary

- Multi-SNP TWAS proxy genes in the class table: {sum(v['n_multi_snp_twas'] for v in counts.values())}.
- One-SNP TWAS proxy genes in the class table: {sum(v['n_one_snp_twas'] for v in counts.values())}.

## SMR/Coloc Boundary

No genes are SMR-supported or coloc-supported in this version. SMR/coloc resources remain not claim-grade ready and are treated as missing evidence, not negative biological evidence.

## Recommended Manuscript Section Changes

- Replace the P1 candidate methods text with "Candidate gene evidence classification and curated exemplar panel".
- Replace P1 results with evidence-class counts and curated exemplar role-spectrum wording.
- Replace any P1-validation language with allowed claims from `candidate_gene_evidence_classes_v0.1.tsv`.
- Keep TWAS as Kidney_Cortex proxy wording only.

## Recommended Figure 3 Redesign

- Panel A: evidence-class summary.
- Panel B: curated exemplar role spectrum.
- Panel C: exemplar MAGMA/TWAS evidence strip.
- Panel D: allowed versus disallowed claim box.

## Stage 4 Readiness

Stage 4 snRNA strengthening can begin after this Stage 3 evidence-class table is accepted. Stage 4 should add donor-level and pseudobulk support without upgrading any gene to causal or validated status.
""", encoding="utf-8")


def update_tracker() -> None:
    rows = read_tsv(TRACKER)
    if not rows:
        return
    fields = list(rows[0].keys())
    today = dt.date.today().isoformat()
    for r in rows:
        if r.get("stage_id") == "3":
            r["status"] = "completed"
            r["start_date"] = r.get("start_date") or today
            r["end_date"] = today
            r["completed_outputs"] = "evidence class definitions; candidate_gene_evidence_classes_v0.1.tsv; curated_exemplar_panel_v0.1.tsv; class summary; replacement text; terminology audit; stage3_report"
            r["blocking_issues"] = "No claim-grade SMR/coloc; no causal/validated wording allowed"
            r["next_stage_ready"] = "yes"
    write_tsv(TRACKER, fields, [[r.get(f, "") for f in fields] for r in rows])


def main() -> None:
    for d in [DOC, TAB, FIG, LOG, SCRIPT]:
        d.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(FROZEN)
    definitions()
    classes = class_table(rows)
    exemplar_panel(classes)
    counts = summary_counts(classes)
    replacement_text(counts)
    p1_audit()
    report(counts)
    update_tracker()
    (LOG / "stage3_evidence_classification.log").write_text(f"completed={dt.datetime.now().isoformat(timespec='seconds')}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
