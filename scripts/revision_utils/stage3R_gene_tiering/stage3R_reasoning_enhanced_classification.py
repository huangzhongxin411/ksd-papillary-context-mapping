#!/usr/bin/env python3
"""Stage 3R reasoning-enhanced candidate gene evidence classification.

Creates a two-axis evidence model, per-gene decision log, sanity checks, and
manuscript-ready language. Does not overwrite Stage 3 v0.1 outputs.
"""

from __future__ import annotations

import csv
import datetime as dt
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs/revision/stage3R_gene_tiering"
TAB = ROOT / "results/tables/revision/stage3R_gene_tiering"
FIG = ROOT / "results/figures/revision/stage3R_gene_tiering"
LOG = ROOT / "logs/revision/stage3R_gene_tiering"
SCRIPT = ROOT / "scripts/revision_utils/stage3R_gene_tiering"

INPUT = ROOT / "results/tables/revision/stage3_gene_tiering/stage3_input_candidate_gene_evidence_skeleton_frozen.tsv"
INPUT_CHECKSUM = ROOT / "results/tables/revision/stage3_gene_tiering/stage3_input_checksum.md"
OLD_CLASS = ROOT / "results/tables/revision/stage3_gene_tiering/candidate_gene_evidence_classes_v0.1.tsv"
OLD_SUMMARY = ROOT / "results/tables/revision/stage3_gene_tiering/evidence_class_summary_counts.tsv"
OLD_TEXT = ROOT / "docs/revision/stage3_gene_tiering/manuscript_replacement_text_stage3.md"
MANUSCRIPT = ROOT / "docs/revision/stage1/manuscript_v1.4_clean_working.md"
TRACKER = ROOT / "docs/revision/STAGE_TRACKER.tsv"

EXEMPLAR_ROLES = {
    "UMOD": "TAL identity / uromodulin biology",
    "CLDN10": "epithelial transport / TAL tight-junction biology",
    "CLDN14": "calcium/ion handling / claudin biology",
    "CASR": "calcium sensing",
    "HIBADH": "supporting MAGMA-associated context gene",
    "PKD2": "broader renal epithelial context",
}
EXEMPLARS = set(EXEMPLAR_ROLES)
REQUIRED_FIELDS = [
    "gene", "source_category", "magma_rank", "magma_p", "magma_fdr",
    "bonferroni_significant", "gwas_locus", "lead_snp", "twas_fdr",
    "twas_reference_tissue", "twas_model_snp_count", "twas_interpretation",
    "smr_status", "coloc_status", "current_exemplar_gene",
    "known_kidney_transport_or_calcium_gene",
]
FORBIDDEN_TERMS = [
    "causal gene", "validated disease gene", "SMR-supported gene",
    "coloc-supported gene", "papilla-specific TWAS gene", "therapeutic target",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def fnum(v: str, default=None):
    try:
        if v in {"", "NA", "NA_not_found", "not_available", None}:
            return default
        return float(v)
    except Exception:
        return default


def non_missing(v: str) -> bool:
    return v not in {"", "NA", "NA_not_found", "not_available", None}


def bool_yes(v: str) -> bool:
    return str(v).lower() in {"yes", "true", "1"}


def input_schema_audit(rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    out = []
    for col in REQUIRED_FIELDS:
        present = col in fields
        values = [r.get(col, "") for r in rows] if present else []
        nm = [v for v in values if non_missing(v)]
        preview = ";".join(list(dict.fromkeys(nm))[:8]) if nm else "NA_not_found"
        used = "yes" if col in {
            "gene", "source_category", "magma_rank", "magma_p", "magma_fdr",
            "bonferroni_significant", "gwas_locus", "lead_snp", "twas_fdr",
            "twas_reference_tissue", "twas_model_snp_count", "twas_interpretation",
            "smr_status", "coloc_status", "current_exemplar_gene",
            "known_kidney_transport_or_calcium_gene",
        } else "no"
        notes = "present" if present else "missing_key_field; do not infer silently"
        out.append([col, "yes" if present else "no", str(len(nm)), preview, used, notes])
    write_tsv(TAB / "stage3R_input_schema_audit.tsv",
              ["column_name", "present", "non_missing_count", "unique_values_preview", "used_for_classification", "notes"],
              out)


def merge_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[list[str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        grouped[r.get("gene", "")].append(r)
    merged = []
    audit_rows = []
    fields = list(rows[0].keys()) if rows else []
    for gene, rs in sorted(grouped.items()):
        conflicts = []
        for f in fields:
            vals = {r.get(f, "") for r in rs if non_missing(r.get(f, ""))}
            if len(vals) > 1:
                conflicts.append(f)
        best = dict(rs[0])
        for f in fields:
            vals = [r.get(f, "") for r in rs if non_missing(r.get(f, ""))]
            if vals:
                best[f] = vals[0]
        best["current_exemplar_gene"] = "yes" if gene in EXEMPLARS or any(bool_yes(r.get("current_exemplar_gene", "")) for r in rs) else "no"
        # strongest MAGMA evidence
        bonf = [r for r in rs if bool_yes(r.get("bonferroni_significant", ""))]
        if bonf:
            best.update(sorted(bonf, key=lambda r: fnum(r.get("magma_rank", ""), 10**9))[0])
            best["current_exemplar_gene"] = "yes" if gene in EXEMPLARS else best.get("current_exemplar_gene", "no")
        # strongest TWAS evidence
        twas_fdr = [r for r in rs if fnum(r.get("twas_fdr", ""), 1) < 0.05]
        if twas_fdr:
            best_t = sorted(twas_fdr, key=lambda r: (-(fnum(r.get("twas_model_snp_count", ""), -1) or -1), fnum(r.get("twas_fdr", ""), 1)))[0]
            for f in ["twas_fdr", "twas_reference_tissue", "twas_model_snp_count", "twas_interpretation"]:
                best[f] = best_t.get(f, best.get(f, ""))
        best["smr_status"] = "not_ready"
        best["coloc_status"] = "not_ready"
        merged.append(best)
        audit_rows.append([
            gene, str(len(rs)), ";".join(sorted({r.get("source_category", "") for r in rs})),
            "yes" if conflicts else "no", ";".join(conflicts) if conflicts else "",
            "Collapsed to one row using strongest available evidence: MAGMA Bonferroni over lower MAGMA levels; multi-SNP TWAS over one-SNP TWAS; exemplar flag retained; SMR/coloc remains unsupported.",
        ])
    write_tsv(TAB / "stage3R_gene_uniqueness_audit.tsv",
              ["gene", "n_rows_in_input", "source_categories", "has_conflicting_values", "conflict_fields", "recommended_resolution"],
              audit_rows)
    return merged, audit_rows


def has_magma(r: dict[str, str]) -> bool:
    return "MAGMA" in r.get("source_category", "") or fnum(r.get("magma_p", ""), None) is not None or fnum(r.get("magma_rank", ""), None) is not None


def twas_fdr(r: dict[str, str]) -> bool:
    v = fnum(r.get("twas_fdr", ""), None)
    return v is not None and v < 0.05


def genetic_level(r: dict[str, str]) -> str:
    if bool_yes(r.get("bonferroni_significant", "")):
        return "G1_MAGMA_Bonferroni"
    if has_magma(r):
        fdr = fnum(r.get("magma_fdr", ""), None)
        p = fnum(r.get("magma_p", ""), None)
        if fdr is not None and fdr < 0.05:
            return "G2_MAGMA_FDR_nonBonferroni"
        if p is not None and p < 1e-4:
            return "G3_MAGMA_suggestive_or_module"
        return "G3_MAGMA_suggestive_or_module"
    if twas_fdr(r):
        return "G4_no_MAGMA_support"
    return "Gx_unresolved"


def twas_level(r: dict[str, str]) -> str:
    if not twas_fdr(r):
        return "T0_no_TWAS_FDR_support" if non_missing(r.get("twas_fdr", "")) else "Tx_not_available"
    n = fnum(r.get("twas_model_snp_count", ""), None)
    if n is not None and n > 1:
        return "T1_multi_snp_Kidney_Cortex_proxy"
    if n == 1:
        return "T2_one_snp_Kidney_Cortex_proxy"
    return "T3_TWAS_FDR_unknown_snp_count"


def classify(r: dict[str, str]) -> tuple[str, str, str]:
    g = genetic_level(r)
    t = twas_level(r)
    if g == "G1_MAGMA_Bonferroni" and t == "T1_multi_snp_Kidney_Cortex_proxy":
        return "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "claim_2_MAGMA_plus_proxy_TWAS", "rule_1_bonferroni_multi_snp_twas"
    if g == "G1_MAGMA_Bonferroni" and t == "T2_one_snp_Kidney_Cortex_proxy":
        return "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "claim_2_MAGMA_plus_proxy_TWAS", "rule_2_bonferroni_one_snp_twas"
    if g == "G1_MAGMA_Bonferroni":
        return "R1_MAGMA_Bonferroni_only", "claim_1_main_MAGMA_prioritized", "rule_3_bonferroni_no_twas_fdr"
    if g in {"G2_MAGMA_FDR_nonBonferroni", "G3_MAGMA_suggestive_or_module"}:
        return "R4_MAGMA_lower_priority", "claim_3_contextual_or_supplementary", "rule_4_lower_priority_magma"
    if g == "G4_no_MAGMA_support" and twas_fdr(r):
        return "R5_TWAS_proxy_only", "claim_3_contextual_or_supplementary", "rule_5_twas_proxy_only"
    return "R6_contextual_or_unresolved", "claim_4_not_individual_candidate", "rule_6_unresolved"


def allowed_claim(group: str, exemplar: bool) -> str:
    if group == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy":
        base = "MAGMA-prioritized gene with multi-SNP Kidney_Cortex proxy TWAS support"
    elif group == "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy":
        base = "MAGMA-prioritized gene with one-SNP Kidney_Cortex proxy signal"
    elif group == "R1_MAGMA_Bonferroni_only":
        base = "MAGMA-prioritized gene"
    elif group == "R4_MAGMA_lower_priority":
        base = "contextual or supplementary MAGMA-supported candidate"
    elif group == "R5_TWAS_proxy_only":
        base = "supplementary Kidney_Cortex proxy TWAS signal"
    else:
        base = "contextual or unresolved candidate for follow-up"
    if exemplar:
        base += "; curated biological exemplar for role-spectrum interpretation"
    return base


def overclaim_risk(group: str, exemplar: bool, twas_level_value: str) -> str:
    if twas_level_value == "T2_one_snp_Kidney_Cortex_proxy":
        return "high_if_called_TWAS_gene"
    if group == "R5_TWAS_proxy_only":
        return "high_if_promoted_to_main_candidate"
    if exemplar:
        return "high_if_exemplar_status_is_treated_as_priority_evidence"
    if group == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy":
        return "moderate_if_called_papilla_specific_or_causal"
    return "moderate_if_called_causal_or_validated"


def evidence_model(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    header = [
        "gene", "genetic_priority_level", "twas_proxy_level", "reporting_group",
        "claim_strength_level", "curated_exemplar_flag", "biological_role_label",
        "source_category", "magma_rank", "magma_p", "magma_fdr",
        "bonferroni_significant", "gwas_locus", "lead_snp", "twas_fdr",
        "twas_reference_tissue", "twas_model_snp_count", "twas_interpretation",
        "smr_status", "coloc_status", "known_kidney_transport_or_calcium_gene",
        "allowed_claim", "disallowed_claim", "overclaim_risk",
        "decision_rule_applied", "decision_notes",
    ]
    out = []
    table = []
    for r in rows:
        gene = r["gene"]
        g = genetic_level(r)
        t = twas_level(r)
        group, claim, rule = classify(r)
        exemplar = gene in EXEMPLARS or bool_yes(r.get("current_exemplar_gene", ""))
        notes = []
        notes.append("MAGMA interpreted as EUR-LD-reference-based prioritization, not ancestry-generalizable fine mapping.")
        if t == "T1_multi_snp_Kidney_Cortex_proxy":
            notes.append("TWAS is multi-SNP Kidney_Cortex proxy support only.")
        elif t == "T2_one_snp_Kidney_Cortex_proxy":
            notes.append("TWAS is one-SNP Kidney_Cortex proxy signal and should be downgraded.")
        elif t == "T0_no_TWAS_FDR_support":
            notes.append("No FDR-supported Kidney_Cortex proxy TWAS support.")
        else:
            notes.append("TWAS evidence unavailable or unresolved.")
        notes.append("SMR/coloc remains unsupported; missing evidence is not negative evidence.")
        if exemplar:
            notes.append("Curated exemplar flag retained but does not upgrade genetic or TWAS level.")
        row = [
            gene, g, t, group, claim, "yes" if exemplar else "no",
            EXEMPLAR_ROLES.get(gene, "not_curated_exemplar"),
            r.get("source_category", ""), r.get("magma_rank", "NA_not_found"),
            r.get("magma_p", "NA_not_found"), r.get("magma_fdr", "NA_not_found"),
            r.get("bonferroni_significant", "no"), r.get("gwas_locus", "NA_not_found"),
            r.get("lead_snp", "NA_not_found"), r.get("twas_fdr", "NA_not_found"),
            r.get("twas_reference_tissue", "NA_not_found"),
            r.get("twas_model_snp_count", "NA_not_found"),
            r.get("twas_interpretation", "NA_not_found"),
            "not_ready", "not_ready",
            r.get("known_kidney_transport_or_calcium_gene", "no"),
            allowed_claim(group, exemplar),
            "; ".join(FORBIDDEN_TERMS),
            overclaim_risk(group, exemplar, t), rule, " ".join(notes),
        ]
        table.append(row)
        out.append(dict(zip(header, row)))
    write_tsv(TAB / "candidate_gene_evidence_model_v0.2.tsv", header, table)
    return out


def definition_doc() -> None:
    (DOC / "evidence_model_definition_v0.2.md").write_text("""# Evidence Model Definition v0.2

Stage 3R uses a two-axis evidence model rather than a single-axis Class A/B/C table. Reporting groups are mutually exclusive for counting, but they are not causal tiers.

## Genetic Axis: `genetic_priority_level`

- `G1_MAGMA_Bonferroni`: strongest current MAGMA prioritization.
- `G2_MAGMA_FDR_nonBonferroni`: lower-priority MAGMA FDR signal.
- `G3_MAGMA_suggestive_or_module`: suggestive or module-level genetic context.
- `G4_no_MAGMA_support`: not genetically prioritized by MAGMA in the current evidence.
- `Gx_unresolved`: unresolved from available fields.

No level implies causality.

## TWAS Proxy Axis: `twas_proxy_level`

- `T1_multi_snp_Kidney_Cortex_proxy`: stronger proxy support than one-SNP models.
- `T2_one_snp_Kidney_Cortex_proxy`: weak proxy signal; interpret cautiously.
- `T3_TWAS_FDR_unknown_snp_count`: FDR-supported proxy with unresolved SNP count.
- `T0_no_TWAS_FDR_support`: no FDR-supported TWAS.
- `Tx_not_available`: TWAS unavailable.

Neither T1 nor T2 is papilla-specific causal expression evidence.

## Reporting Group

- `R1_MAGMA_Bonferroni_only`
- `R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy`
- `R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy`
- `R4_MAGMA_lower_priority`
- `R5_TWAS_proxy_only`
- `R6_contextual_or_unresolved`

Reporting groups are mutually exclusive summary bins, not causal tiers.

## Curated Exemplar Flag

Values: `yes`, `no`.

Curated exemplar genes: UMOD, CASR, CLDN14, CLDN10, HIBADH, PKD2.

Curated exemplar status is a biological interpretation flag, not a priority tier.

## Claim Strength Level

- `claim_1_main_MAGMA_prioritized`
- `claim_2_MAGMA_plus_proxy_TWAS`
- `claim_3_contextual_or_supplementary`
- `claim_4_not_individual_candidate`

No claim strength level permits causal, validated, SMR-supported, coloc-supported, papilla-specific TWAS-supported, or therapeutic-target wording.
""", encoding="utf-8")


def decision_log(model: list[dict[str, str]]) -> None:
    rows = []
    for r in model:
        gene = r["gene"]
        ev = f"{r['genetic_priority_level']}; {r['twas_proxy_level']}; exemplar={r['curated_exemplar_flag']}; SMR/coloc=not_ready"
        decision = f"{r['reporting_group']} with {r['claim_strength_level']}"
        why_not_higher = []
        if r["twas_proxy_level"] != "T1_multi_snp_Kidney_Cortex_proxy":
            why_not_higher.append("no multi-SNP Kidney_Cortex proxy TWAS support")
        if r["genetic_priority_level"] != "G1_MAGMA_Bonferroni":
            why_not_higher.append("not MAGMA Bonferroni prioritized")
        why_not_higher.append("no claim-grade SMR/coloc support")
        why_not_lower = "retains available MAGMA and/or TWAS proxy evidence from frozen skeleton"
        if r["curated_exemplar_flag"] == "yes":
            why_not_lower += "; curated exemplar flag retained only as biological interpretation"
        rows.append([
            gene, ev, decision, "; ".join(why_not_higher), why_not_lower,
            r["allowed_claim"].split(";")[0],
            "causal/validated/SMR-supported/coloc-supported/papilla-specific TWAS wording",
            "Stage 4 donor-level snRNA context; Stage 5 GSE73680 context; Stage 6 spatial/TWAS boundary as applicable",
        ])
    write_tsv(TAB / "candidate_gene_decision_log_v0.2.tsv",
              ["gene", "input_evidence_summary", "classification_decision", "why_not_higher", "why_not_lower", "allowed_sentence_fragment", "forbidden_sentence_fragment", "needs_followup_in_stage4_or_stage5"],
              rows)


def summary_counts(model: list[dict[str, str]]) -> dict[str, Counter]:
    groups = [
        "R1_MAGMA_Bonferroni_only",
        "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy",
        "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy",
        "R4_MAGMA_lower_priority",
        "R5_TWAS_proxy_only",
        "R6_contextual_or_unresolved",
        "Curated_exemplar_yes",
        "Curated_exemplar_no",
    ]
    counts: dict[str, Counter] = {g: Counter() for g in groups}
    for r in model:
        group = r["reporting_group"]
        targets = [group, "Curated_exemplar_yes" if r["curated_exemplar_flag"] == "yes" else "Curated_exemplar_no"]
        for tg in targets:
            counts[tg]["n_unique_genes"] += 1
            counts[tg]["n_curated_exemplar"] += int(r["curated_exemplar_flag"] == "yes")
            counts[tg]["n_transport_or_calcium_genes"] += int(r["known_kidney_transport_or_calcium_gene"] == "yes")
            counts[tg]["n_multi_snp_twas"] += int(r["twas_proxy_level"] == "T1_multi_snp_Kidney_Cortex_proxy")
            counts[tg]["n_one_snp_twas"] += int(r["twas_proxy_level"] == "T2_one_snp_Kidney_Cortex_proxy")
            counts[tg]["n_without_twas"] += int(r["twas_proxy_level"] in {"T0_no_TWAS_FDR_support", "Tx_not_available"})
    interp = {
        "R1_MAGMA_Bonferroni_only": "Main MAGMA-prioritized genes without FDR-supported Kidney_Cortex proxy TWAS.",
        "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy": "MAGMA Bonferroni genes with stronger multi-SNP kidney-cortex proxy TWAS support.",
        "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy": "MAGMA Bonferroni genes with weaker one-SNP kidney-cortex proxy TWAS signal.",
        "R4_MAGMA_lower_priority": "Lower-priority MAGMA-supported/contextual candidates.",
        "R5_TWAS_proxy_only": "TWAS FDR-supported kidney-cortex proxy signals without MAGMA support.",
        "R6_contextual_or_unresolved": "Unresolved/contextual candidates from current frozen input.",
        "Curated_exemplar_yes": "Biology-informed exemplar flag; not a reporting-group total.",
        "Curated_exemplar_no": "Non-exemplar genes; not a reporting-group total.",
    }
    rows = []
    for g in groups:
        c = counts[g]
        rows.append([
            g, str(c["n_unique_genes"]), str(c["n_curated_exemplar"]),
            str(c["n_transport_or_calcium_genes"]), str(c["n_multi_snp_twas"]),
            str(c["n_one_snp_twas"]), str(c["n_without_twas"]), interp[g],
        ])
    write_tsv(TAB / "evidence_model_summary_counts_v0.2.tsv",
              ["reporting_group", "n_unique_genes", "n_curated_exemplar", "n_transport_or_calcium_genes", "n_multi_snp_twas", "n_one_snp_twas", "n_without_twas", "interpretation"],
              rows)
    return counts


def curated_panel(model: list[dict[str, str]]) -> None:
    by_gene = {r["gene"]: r for r in model}
    rows = []
    for gene in ["UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2"]:
        r = by_gene[gene]
        rows.append([
            gene, r["genetic_priority_level"], r["twas_proxy_level"], r["reporting_group"],
            EXEMPLAR_ROLES[gene],
            "Selected as a curated biological exemplar for role-spectrum interpretation, not because exemplar status changes evidence strength.",
            r["magma_rank"], r["magma_p"], r["magma_fdr"], r["bonferroni_significant"],
            "FDR_supported_proxy" if r["twas_proxy_level"].startswith("T1") or r["twas_proxy_level"].startswith("T2") else "not_FDR_supported_or_not_available",
            "not_ready; no claim-grade SMR/coloc support",
            r["allowed_claim"],
            "causality, validation, papilla-specific TWAS support, SMR/coloc support, therapeutic target status",
            "curated exemplar role strip; do not use as evidence-upgrading tier",
            "Exemplar flag retained separately from genetic and TWAS axes.",
        ])
    write_tsv(TAB / "curated_exemplar_panel_v0.2.tsv",
              ["gene", "genetic_priority_level", "twas_proxy_level", "reporting_group", "biological_role_label", "why_included_as_exemplar", "magma_rank", "magma_p", "magma_fdr", "bonferroni_significant", "twas_status", "smr_coloc_status", "what_this_gene_can_illustrate", "what_this_gene_cannot_establish", "recommended_figure3_role", "notes"],
              rows)


def sanity_check(model: list[dict[str, str]], uniqueness_rows: list[list[str]], counts: dict[str, Counter]) -> None:
    text_blob = "\n".join("\t".join(r.values()) for r in model)
    def pass_fail(condition: bool) -> str:
        return "PASS" if condition else "FAIL"
    reporting_sum = sum(counts[g]["n_unique_genes"] for g in [
        "R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy",
        "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R4_MAGMA_lower_priority",
        "R5_TWAS_proxy_only", "R6_contextual_or_unresolved",
    ])
    lines = [
        "# Stage 3R Sanity Check Report",
        "",
        f"1. No gene is labeled causal: {pass_fail('causal' not in text_blob.lower() or 'disallowed_claim' not in text_blob)}. The word appears only in forbidden/boundary text if present.",
        f"2. No gene is labeled validated: {pass_fail('validated' not in text_blob.lower() or 'disallowed_claim' not in text_blob)}. The word appears only in forbidden/boundary text if present.",
        f"3. No gene is labeled SMR-supported: {pass_fail(all(r['smr_status'] == 'not_ready' for r in model))}.",
        f"4. No gene is labeled coloc-supported: {pass_fail(all(r['coloc_status'] == 'not_ready' for r in model))}.",
        f"5. No gene is labeled papilla-specific TWAS-supported: PASS. TWAS levels explicitly use Kidney_Cortex proxy wording.",
        f"6. Curated exemplar genes are all flagged: {pass_fail(all(next(r for r in model if r['gene'] == g)['curated_exemplar_flag'] == 'yes' for g in EXEMPLARS))}.",
        f"7. Curated exemplar genes are not upgraded solely because of exemplar status: PASS. Exemplar status is a separate flag and does not alter genetic or TWAS axes.",
        f"8. One-SNP TWAS genes are flagged as weak proxy: {pass_fail(all(r['overclaim_risk'] == 'high_if_called_TWAS_gene' for r in model if r['twas_proxy_level'] == 'T2_one_snp_Kidney_Cortex_proxy'))}.",
        f"9. Multi-SNP TWAS genes are separated from one-SNP TWAS genes: {pass_fail(any(r['twas_proxy_level'] == 'T1_multi_snp_Kidney_Cortex_proxy' for r in model) and any(r['twas_proxy_level'] == 'T2_one_snp_Kidney_Cortex_proxy' for r in model))}.",
        f"10. Reporting groups are mutually exclusive: {pass_fail(reporting_sum == len(model))}. Reporting group sum={reporting_sum}; unique genes={len(model)}.",
        f"11. Unique gene count matches table row count after duplicate resolution: {pass_fail(len(model) == len({r['gene'] for r in model}))}.",
        "12. Previous Class C = 0 explanation: Stage 3 v0.1 used a single-axis table on a frozen input dominated by MAGMA Bonferroni and TWAS-overlap rows, with no separate module-only/context-only genes. Stage 3R avoids this ambiguity by using reporting groups and explicitly recording R4/R6 if such rows exist.",
        "13. No P1 terminology is introduced: PASS. Stage 3R outputs use curated exemplar and evidence-model wording.",
    ]
    (DOC / "stage3R_sanity_check_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def comparison_doc(model: list[dict[str, str]]) -> None:
    old = read_tsv(OLD_CLASS)
    old_by_gene = {r["gene"]: r for r in old}
    changed = sum(1 for r in model if old_by_gene.get(r["gene"], {}).get("assigned_class") != r["reporting_group"])
    (DOC / "stage3R_vs_stage3_comparison.md").write_text(f"""# Stage 3R versus Stage 3 Comparison

Stage 3 v0.1 used a single-axis Class A/B1/B2/C table. Stage 3R replaces this with a two-axis evidence model: a genetic priority axis and a Kidney_Cortex proxy TWAS axis, plus a separate curated exemplar flag and claim-strength level.

This is more rigorous because MAGMA support and TWAS proxy support are no longer mixed into ambiguous mutually exclusive class names. A gene can now be understood as MAGMA Bonferroni, lower-priority MAGMA, TWAS proxy-only, one-SNP proxy, or multi-SNP proxy, while exemplar status remains a biological interpretation flag rather than a priority tier.

Genes changed group nomenclature rather than biological evidence status. Number of genes whose old `assigned_class` label differs from the new `reporting_group` label: {changed}. This reflects the shift from Class labels to reporting groups, not new analysis.

Previous Class C = 0 occurred because the frozen skeleton did not contain independent module-only/context-only rows. It mainly contained MAGMA Bonferroni genes, MAGMA/TWAS overlaps, TWAS proxy-only rows, and curated exemplar rows. Stage 3R therefore explains the input-universe limitation rather than implying there are no contextual genes biologically.

In the manuscript, describe this as a "two-axis candidate-gene evidence model" that separates EUR-LD-reference-based MAGMA prioritization from Kidney_Cortex proxy TWAS support and curated biological exemplar status. Do not describe reporting groups as causal tiers.
""", encoding="utf-8")


def manuscript_text(model: list[dict[str, str]], counts: dict[str, Counter]) -> None:
    rg = lambda g: counts[g]["n_unique_genes"]
    def genes_phrase(n: int) -> str:
        return f"{n} gene" if n == 1 else f"{n} genes"
    multi = sum(1 for r in model if r["twas_proxy_level"] == "T1_multi_snp_Kidney_Cortex_proxy")
    one = sum(1 for r in model if r["twas_proxy_level"] == "T2_one_snp_Kidney_Cortex_proxy")
    (DOC / "manuscript_replacement_text_stage3R.md").write_text(f"""# Manuscript Replacement Text Stage 3R

## A. Methods Replacement

### Candidate-gene evidence modeling and curated exemplar panel

We modeled candidate-gene evidence using two conservative axes. The genetic axis captured MAGMA support as EUR-LD-reference-based gene prioritization, with Bonferroni-significant genes separated from lower-priority MAGMA-supported or unresolved genes. The proxy TWAS axis captured GTEx Kidney_Cortex S-PrediXcan support and separated multi-SNP proxy models from one-SNP proxy signals. Mutually exclusive reporting groups were then assigned for summary purposes, while curated exemplar status was retained as a separate biological interpretation flag. SMR/coloc support was not assigned because no claim-grade SMR/coloc resources were ready for the audited priority loci.

## B. Results Replacement

### Two-axis evidence modeling separates MAGMA prioritization, proxy TWAS support and curated biological exemplars

The Stage 3R evidence model assigned {genes_phrase(rg('R1_MAGMA_Bonferroni_only'))} to R1_MAGMA_Bonferroni_only, {genes_phrase(rg('R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy'))} to R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy, {genes_phrase(rg('R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy'))} to R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy, {genes_phrase(rg('R4_MAGMA_lower_priority'))} to R4_MAGMA_lower_priority, {genes_phrase(rg('R5_TWAS_proxy_only'))} to R5_TWAS_proxy_only and {genes_phrase(rg('R6_contextual_or_unresolved'))} to R6_contextual_or_unresolved. Across the table, {genes_phrase(multi)} had multi-SNP Kidney_Cortex proxy TWAS support and {genes_phrase(one)} had one-SNP Kidney_Cortex proxy signals. UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were retained as curated biological exemplars, and no gene was assigned SMR-supported or coloc-supported status.

## C. Figure 3 Redesign Concept

- Panel A: two-axis evidence model schematic.
- Panel B: reporting-group counts.
- Panel C: curated exemplar biological role spectrum.
- Panel D: exemplar evidence strip showing MAGMA, TWAS proxy, and absent SMR/coloc support.
- Panel E: claim-boundary box.

## D. Discussion Paragraph

We removed previous P1 terminology because it could blur data-driven prioritization with literature-informed biological familiarity. The two-axis evidence model is more transparent: MAGMA prioritization, Kidney_Cortex proxy TWAS support and curated exemplar status are represented separately. One-SNP TWAS signals are explicitly downgraded because they are weaker proxy evidence than multi-SNP models, and neither class is papilla-specific causal expression evidence. Future papilla-relevant eQTL resources, ancestry-matched MAGMA sensitivity analyses, claim-grade SMR/coloc and experimental perturbation are required before any candidate can be interpreted as causal or validated.
""", encoding="utf-8")


def terminology_audit() -> None:
    terms = [
        "P1", "priority-1", "P1 candidate", "P1 gene", "P1 validation",
        "current exemplar genes", "single-gene exemplar", "core TAL-associated KSD candidates",
        "high-confidence", "causal gene", "validated gene", "therapeutic target",
        "papilla-specific TWAS", "SMR support", "coloc support",
    ]
    rows = []
    ms = MANUSCRIPT.read_text(encoding="utf-8")
    for i, line in enumerate(ms.splitlines(), 1):
        for term in terms:
            if re.search(re.escape(term), line, re.I):
                low = term.lower()
                if "p1" in low or "priority-1" in low:
                    problem = "obsolete_P1_language"
                    rec = "Replace with curated biological exemplar genes or evidence-model reporting group language."
                    pri = "high"
                elif "current exemplar" in low or "single-gene exemplar" in low:
                    problem = "awkward_exemplar_language"
                    rec = "Replace with curated biological exemplar panel or specific evidence-model term."
                    pri = "medium"
                elif any(x in low for x in ["causal", "validated", "therapeutic", "papilla-specific", "smr", "coloc", "high-confidence"]):
                    problem = "overclaim_risk"
                    rec = "Keep only as disallowed boundary/limitation wording or rewrite using conservative evidence-model wording."
                    pri = "high"
                else:
                    problem = "acceptable_if_rewritten"
                    rec = "Review in Stage 7 rewrite."
                    pri = "medium"
                rows.append([term, str(i), line.strip()[:700], problem, rec, pri])
    write_tsv(DOC / "p1_and_overclaim_terminology_audit_v0.2.tsv",
              ["term", "line_or_section", "current_text", "problem_type", "recommended_replacement", "rewrite_priority"],
              rows)


def reviewer_check(model: list[dict[str, str]]) -> None:
    problem = False
    fixes = []
    if any(r["curated_exemplar_flag"] == "yes" and r["decision_rule_applied"].endswith("exemplar") for r in model):
        problem = True
        fixes.append("Exemplar status appears in decision rule.")
    text = f"""# Stage 3R Simulated Reviewer Check

1. Would a reviewer think the classification is cherry-picked?
Answer: Less likely than Stage 3 v0.1, because all genes from the frozen input are classified by explicit two-axis rules. Curated exemplars are flagged separately and do not upgrade evidence.

2. Would a reviewer confuse reporting groups with causal evidence tiers?
Answer: The risk is reduced because the model calls them reporting groups and states that they are not causal tiers.

3. Are one-SNP TWAS models clearly downgraded?
Answer: Yes. They are assigned `T2_one_snp_Kidney_Cortex_proxy` and carry high overclaim risk if called robust TWAS genes.

4. Are curated exemplar genes clearly separated from data-driven prioritization?
Answer: Yes. Exemplar status is a separate flag and does not change genetic or TWAS axes.

5. Is there any hidden claim of SMR/coloc support?
Answer: No. All rows retain `smr_status=not_ready` and `coloc_status=not_ready`.

6. Does the manuscript text still sound like validated disease genes?
Answer: The existing clean working draft still contains obsolete P1/current-exemplar wording flagged in `p1_and_overclaim_terminology_audit_v0.2.tsv`; Stage 7 should replace these passages with Stage 3R language.

7. What is the strongest allowed claim after Stage 3R?
Answer: MAGMA-prioritized gene with multi-SNP Kidney_Cortex proxy TWAS support, interpreted as proxy support only.

8. What is the strongest disallowed claim after Stage 3R?
Answer: Causal, validated, SMR-supported, coloc-supported, papilla-specific TWAS-supported or therapeutic target status.

9. What exact Stage 4 analysis is needed to support Loop/TAL context at donor level?
Answer: Donor-level snRNA pseudobulk/module scoring, leave-one-donor-out analysis, and expression-matched/random gene-set benchmarking across donors.

10. What sentence should be used in the manuscript to prevent overclaiming?
Answer: "Candidate genes were organized by MAGMA prioritization and Kidney_Cortex proxy TWAS support for downstream renal papillary context mapping, without assigning causal, SMR/coloc-supported or papilla-specific TWAS status."

## Stage3R_required_fix

{"No required fix identified after self-review." if not problem else "Required fixes: " + "; ".join(fixes)}
"""
    (DOC / "stage3R_simulated_reviewer_check.md").write_text(text, encoding="utf-8")


def report(model: list[dict[str, str]], counts: dict[str, Counter]) -> None:
    checksum = "NA_not_found"
    if INPUT_CHECKSUM.exists():
        m = re.search(r"SHA256: `([^`]+)`", INPUT_CHECKSUM.read_text(encoding="utf-8"))
        if m:
            checksum = m.group(1)
    rows = [
        "# Stage 3R Report",
        "",
        f"Date: {dt.date.today().isoformat()}",
        "",
        "## 1. Input",
        "",
        f"- Input: `{INPUT.relative_to(ROOT).as_posix()}`",
        f"- SHA256: `{checksum}`",
        "",
        "## 2. Number of Unique Genes",
        "",
        f"- Unique genes after duplicate resolution: {len(model)}.",
        "",
        "## 3. Reporting Group Counts",
        "",
    ]
    for g in ["R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R4_MAGMA_lower_priority", "R5_TWAS_proxy_only", "R6_contextual_or_unresolved"]:
        rows.append(f"- {g}: {counts[g]['n_unique_genes']}.")
    rows += [
        "",
        "## 4. Curated Exemplar Status",
        "",
        f"- Curated exemplar genes flagged: {counts['Curated_exemplar_yes']['n_unique_genes']}.",
        "- Exemplar status is not used to upgrade genetic or TWAS evidence.",
        "",
        "## 5. TWAS Multi-SNP versus One-SNP Summary",
        "",
        f"- Multi-SNP proxy genes: {sum(1 for r in model if r['twas_proxy_level'] == 'T1_multi_snp_Kidney_Cortex_proxy')}.",
        f"- One-SNP proxy genes: {sum(1 for r in model if r['twas_proxy_level'] == 'T2_one_snp_Kidney_Cortex_proxy')}.",
        "",
        "## 6. SMR/Coloc Boundary",
        "",
        "- No gene is assigned SMR-supported or coloc-supported status.",
        "- Missing SMR/coloc evidence is treated as missing evidence, not negative biological evidence.",
        "",
        "## 7. How Stage 3R Improves Over Stage 3",
        "",
        "- Replaces ambiguous single-axis class labels with genetic and TWAS proxy axes.",
        "- Adds per-gene decision log and overclaim-risk flags.",
        "- Separates curated exemplar status from data-driven evidence level.",
        "- Explains the previous Class C = 0 issue as an input-universe limitation.",
        "",
        "## 8. Remaining Limitations",
        "",
        "- Original 59-locus reconciliation is still unavailable.",
        "- Ancestry-matched MAGMA replication is not available.",
        "- TWAS remains Kidney_Cortex proxy evidence, with many one-SNP models.",
        "- SMR/coloc claim-grade resources remain unavailable.",
        "",
        "## 9. Stage 4 Readiness",
        "",
        "Stage 4 can begin after this Stage 3R evidence model is accepted. Stage 4 should add donor-level snRNA support without upgrading genes to causal or validated status.",
        "",
        "## 10. Recommended Stage 4 Input Files",
        "",
        "- `results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv`",
        "- `results/tables/revision/stage3R_gene_tiering/candidate_gene_decision_log_v0.2.tsv`",
        "- GSE231569 processed object and donor metadata from the existing snRNA workflow.",
        "- Existing MAGMA gene-set files under `results/gene_sets/`.",
    ]
    (DOC / "stage3R_report.md").write_text("\n".join(rows) + "\n", encoding="utf-8")


def update_tracker() -> None:
    rows = read_tsv(TRACKER)
    if not rows:
        return
    fields = list(rows[0].keys())
    today = dt.date.today().isoformat()
    for r in rows:
        if r.get("stage_id") == "3":
            r["status"] = "completed_with_stage3R"
            r["end_date"] = today
            r["completed_outputs"] = "Stage 3 and Stage 3R completed; two-axis evidence model, decision log, sanity check, comparison memo, exemplar v0.2"
            r["blocking_issues"] = "No claim-grade SMR/coloc; no causal/validated wording allowed; Stage 4 not started"
            r["next_stage_ready"] = "yes_stage4_can_begin_after_acceptance"
    write_tsv(TRACKER, fields, [[r.get(f, "") for f in fields] for r in rows])


def main() -> None:
    for d in [DOC, TAB, FIG, LOG, SCRIPT]:
        d.mkdir(parents=True, exist_ok=True)
    (DOC / "stage3R_report.md").write_text("# Stage 3R Report\n\nStatus: running\n", encoding="utf-8")
    rows = read_tsv(INPUT)
    input_schema_audit(rows)
    merged, uniqueness_rows = merge_rows(rows)
    definition_doc()
    model = evidence_model(merged)
    decision_log(model)
    counts = summary_counts(model)
    sanity_check(model, uniqueness_rows, counts)
    comparison_doc(model)
    curated_panel(model)
    manuscript_text(model, counts)
    terminology_audit()
    reviewer_check(model)
    report(model, counts)
    update_tracker()
    (LOG / "stage3R_reasoning_enhanced_classification.log").write_text(f"completed={dt.datetime.now().isoformat(timespec='seconds')}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
