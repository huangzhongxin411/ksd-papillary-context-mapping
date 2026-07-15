#!/usr/bin/env python3
"""Build synchronized v3.3 manifests, checksums, and the draft submission package."""

from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "submission_package_v3.3_draft"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


main_sources = {
    "figures/Figure_1.pdf": "figures_v3.3/Figure_1_GWAS_QC_Manhattan_polished.pdf",
    "figures/Figure_1.png": "figures_v3.3/Figure_1_GWAS_QC_Manhattan_polished.png",
    "figures/Figure_2.pdf": "figures_v3.3/Figure_2_snRNA_context_polished.pdf",
    "figures/Figure_2.png": "figures_v3.3/Figure_2_snRNA_context_polished.png",
    "figures/Figure_3.pdf": "figures_v3.3/Figure_3_candidate_reporting_model_polished.pdf",
    "figures/Figure_3.png": "figures_v3.3/Figure_3_candidate_reporting_model_polished.png",
    "figures/Figure_4.pdf": "figures_v3.3/Figure_4_bulk_disease_context_polished.pdf",
    "figures/Figure_4.png": "figures_v3.3/Figure_4_bulk_disease_context_polished.png",
}
supp_sources = {
    "supplementary_figures/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.pdf": "supplementary_figures_v3.3/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics_polished.pdf",
    "supplementary_figures/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.png": "supplementary_figures_v3.3/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics_polished.png",
    "supplementary_figures/Supplementary_Figure_S2_spatial_projection.pdf": "supplementary_figures_v3.3/Supplementary_Figure_S2_spatial_projection_polished.pdf",
    "supplementary_figures/Supplementary_Figure_S2_spatial_projection.png": "supplementary_figures_v3.3/Supplementary_Figure_S2_spatial_projection_polished.png",
    "supplementary_figures/Supplementary_Figure_S3_TWAS_candidate_proxy.pdf": "supplementary_figures_v3.3/Supplementary_Figure_S3_TWAS_candidate_proxy_polished.pdf",
    "supplementary_figures/Supplementary_Figure_S3_TWAS_candidate_proxy.png": "supplementary_figures_v3.3/Supplementary_Figure_S3_TWAS_candidate_proxy_polished.png",
}


def build_figure_manifest(path: Path, mapping: dict[str, str]) -> None:
    rows = []
    for target, source in mapping.items():
        target_path, source_path = ROOT / target, ROOT / source
        if not target_path.is_file() or not source_path.is_file():
            raise FileNotFoundError(target if not target_path.is_file() else source)
        target_hash, source_hash = sha256(target_path), sha256(source_path)
        rows.append({
            "file_path": target, "source_polished_file": source,
            "file_type": target_path.suffix.lower().lstrip("."), "size_bytes": target_path.stat().st_size,
            "sha256": target_hash, "matches_polished_source": "yes" if target_hash == source_hash else "no",
            "version_status": "v3.3 manuscript-synchronized working tree",
            "notes": "Prepared for future v1.0.1 release; v1.0.0 unchanged.",
        })
    write_tsv(path, ["file_path", "source_polished_file", "file_type", "size_bytes", "sha256", "matches_polished_source", "version_status", "notes"], rows)


build_figure_manifest(ROOT / "figures/Figure_file_manifest_v3.3.tsv", main_sources)
build_figure_manifest(ROOT / "supplementary_figures/Supplementary_Figure_file_manifest_v3.3.tsv", supp_sources)

package_mapping = {
    "manuscript_v3.3_submission_clean.docx": "manuscript/manuscript_v3.3_submission_clean.docx",
    "manuscript_v3.3_submission_clean_rendered.pdf": "manuscript/manuscript_v3.3_submission_clean_rendered.pdf",
    "manuscript_v3.3_submission_clean.md": "manuscript/manuscript_v3.3_submission_clean.md",
    **{Path(k).name: k for k in main_sources},
    "Supplementary_Figure_S1.pdf": "supplementary_figures/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.pdf",
    "Supplementary_Figure_S1.png": "supplementary_figures/Supplementary_Figure_S1_GWAS_MAGMA_diagnostics.png",
    "Supplementary_Figure_S2.pdf": "supplementary_figures/Supplementary_Figure_S2_spatial_projection.pdf",
    "Supplementary_Figure_S2.png": "supplementary_figures/Supplementary_Figure_S2_spatial_projection.png",
    "Supplementary_Figure_S3.pdf": "supplementary_figures/Supplementary_Figure_S3_TWAS_candidate_proxy.pdf",
    "Supplementary_Figure_S3.png": "supplementary_figures/Supplementary_Figure_S3_TWAS_candidate_proxy.png",
    "Supplementary_Tables_1-6.xlsx": "supplementary_tables/Supplementary_Tables_1-6.xlsx",
    "Supplementary_Figure_Legends_v3.3.md": "supplementary_materials/Supplementary_Figure_Legends_v3.3.md",
    "Supplementary_Table_Captions_v3.3.md": "supplementary_materials/Supplementary_Table_Captions_v3.3.md",
    "Source_Data_manifest_v3.3.tsv": "source_data/Source_Data_manifest_v3.3.tsv",
}
for number, stem in [
    (1, "locus_reconstruction"), (2, "MAGMA_results_and_modules"),
    (3, "snRNA_context_and_sensitivity"), (4, "TWAS_proxy_results"),
    (5, "candidate_reporting_matrix"), (6, "bulk_disease_context_and_tissue_state"),
]:
    name = f"Supplementary_Table_{number}_{stem}.tsv"
    package_mapping[name] = f"supplementary_tables/{name}"

if PACKAGE.exists():
    shutil.rmtree(PACKAGE)
PACKAGE.mkdir(parents=True)
for destination, source in package_mapping.items():
    source_path = ROOT / source
    if not source_path.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source_path, PACKAGE / destination)

checklist_rows = [
    ("Clean v3.3 DOCX", "complete", "Authors", "Review final text and formatting", "yes", "15-page render passed visual QC."),
    ("Clean v3.3 rendered PDF", "complete", "Authors", "Use for final human review", "yes", "No draft header/footer or local file-name notes."),
    ("Main Figures 1-4", "complete", "Authors", "Confirm aesthetic approval", "yes", "Files match Phase 7-Step 2 polished outputs."),
    ("Supplementary Figures S1-S3", "complete", "Authors", "Confirm multipage S2/S3 approval", "yes", "S2=5 pages; S3=13 pages."),
    ("Supplementary Tables 1-6", "complete", "Authors", "Confirm workbook and TSV content", "yes", "Workbook and six individual TSVs included."),
    ("Ethics wording", "human_pending", "Corresponding authors/institution", "Verify institutional wording and remove placeholder", "yes", "Must not be inferred by automation."),
    ("CRediT contributions", "human_pending", "All authors", "Complete author-role statement", "yes", "Placeholder intentionally retained."),
    ("Funding and competing interests", "human_confirmation", "All authors", "Confirm existing statements", "yes", "No automated change made."),
    ("All-author approval", "human_pending", "Corresponding authors", "Obtain documented approval", "yes", "Required before release or submission."),
    ("Repository Zenodo DOI", "pending_release", "Repository owner", "Mint only after approved v1.0.1 release", "yes", "No repository DOI exists or is claimed."),
    ("GitHub Release v1.0.1", "not_created", "Repository owner", "Create only after branch review and merge", "yes", "v1.0.0 remains unchanged."),
]
write_tsv(
    PACKAGE / "BMC_Genomics_readiness_checklist_v3.3.tsv",
    ["check_item", "status", "owner", "required_action", "done_before_submission", "notes"],
    [dict(zip(["check_item", "status", "owner", "required_action", "done_before_submission", "notes"], row)) for row in checklist_rows],
)

manifest_rows = []
for path in sorted(PACKAGE.iterdir(), key=lambda p: p.name.lower()):
    if not path.is_file() or path.name == "Submission_file_manifest_v3.3.tsv":
        continue
    manifest_rows.append({
        "file_name": path.name, "file_type": path.suffix.lower().lstrip("."),
        "size_bytes": path.stat().st_size, "sha256": sha256(path),
        "submission_role": "manuscript" if path.name.startswith("manuscript_") else
                           "main_figure" if path.name.startswith("Figure_") else
                           "supplementary_figure" if path.name.startswith("Supplementary_Figure_S") and path.suffix.lower() in {".pdf", ".png"} else
                           "supplementary_table" if path.name.startswith("Supplementary_Table") else
                           "package_metadata",
        "status": "current_v3.3", "notes": "Draft package for human review; not submitted and not archived to Zenodo.",
    })
write_tsv(PACKAGE / "Submission_file_manifest_v3.3.tsv", ["file_name", "file_type", "size_bytes", "sha256", "submission_role", "status", "notes"], manifest_rows)

canonical = [
    "manuscript/manuscript_v3.3_submission_clean.md",
    "manuscript/manuscript_v3.3_submission_clean.docx",
    "manuscript/manuscript_v3.3_submission_clean_rendered.pdf",
    "manuscript/manuscript_v3.3_revision_change_log.md",
    *main_sources.keys(), *supp_sources.keys(),
    "figures/Figure_file_manifest_v3.3.tsv",
    "supplementary_figures/Supplementary_Figure_file_manifest_v3.3.tsv",
    "supplementary_materials/Supplementary_Figure_Legends_v3.3.md",
    "supplementary_materials/Supplementary_Table_Captions_v3.3.md",
    "source_data/Source_Data_manifest_v3.3.tsv",
    "supplementary_tables/Supplementary_Tables_1-6.xlsx",
]
canonical += [f"supplementary_tables/Supplementary_Table_{n}_{stem}.tsv" for n, stem in [
    (1, "locus_reconstruction"), (2, "MAGMA_results_and_modules"),
    (3, "snRNA_context_and_sensitivity"), (4, "TWAS_proxy_results"),
    (5, "candidate_reporting_matrix"), (6, "bulk_disease_context_and_tissue_state"),
]]
manifest = []
for relative in canonical:
    path = ROOT / relative
    if not path.is_file():
        raise FileNotFoundError(relative)
    manifest.append({
        "file_path": relative, "file_type": path.suffix.lower().lstrip("."),
        "description": path.stem.replace("_", " "), "package_role": "v3.3 current-facing file",
        "used_in_manuscript": "yes" if relative.startswith(("manuscript/", "figures/", "supplementary_tables/")) else "supporting",
        "deprecated_status": "current", "size_bytes": path.stat().st_size,
    })
write_tsv(ROOT / "MANIFEST.tsv", ["file_path", "file_type", "description", "package_role", "used_in_manuscript", "deprecated_status", "size_bytes"], manifest)

checksum_paths = canonical + ["MANIFEST.tsv"]
with (ROOT / "CHECKSUMS.sha256").open("w", encoding="utf-8", newline="") as handle:
    for relative in sorted(checksum_paths):
        handle.write(f"{sha256(ROOT / relative)}  {relative}\n")

print(f"Created synchronized manifests and package with {len(manifest_rows)} payload files.")
