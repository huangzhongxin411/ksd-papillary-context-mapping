#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const outDir = path.join(root, "results/tables/revision/stage7F_designer_artwork_brief");
const logDir = path.join(root, "logs/revision/stage7F_designer_artwork_brief");

function row(headers, values) {
  return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
}

function escapeTsv(value) {
  return String(value ?? "").replace(/[\t\r\n]+/g, " ").trim();
}

function columnName(count) {
  let number = count;
  let name = "";
  while (number > 0) {
    number -= 1;
    name = String.fromCharCode(65 + (number % 26)) + name;
    number = Math.floor(number / 26);
  }
  return name;
}

async function writeAuditedTsv(filename, headers, rows) {
  const matrix = [headers, ...rows.map((item) => headers.map((header) => item[header] ?? ""))];
  const workbook = Workbook.create();
  const sheet = workbook.worksheets.add("ArtworkBrief");
  sheet.getRangeByIndexes(0, 0, matrix.length, headers.length).values = matrix;
  const audit = await workbook.inspect({
    kind: "table",
    range: `ArtworkBrief!A1:${columnName(headers.length)}${matrix.length}`,
    include: "values",
    tableMaxRows: Math.min(matrix.length, 8),
    tableMaxCols: headers.length,
    maxChars: 7000,
  });
  if (!audit?.ndjson) throw new Error(`artifact-tool inspection failed for ${filename}`);
  const text = matrix.map((cells) => cells.map(escapeTsv).join("\t")).join("\n") + "\n";
  await fs.writeFile(path.join(outDir, filename), text, "utf8");
  return `${filename}\trows=${rows.length}\tcolumns=${headers.length}\tartifact_tool_audit=pass`;
}

const textHeaders = [
  "figure", "panel", "current_text_or_label", "problem", "recommended_short_text",
  "move_detail_to_legend", "claim_boundary_preserved", "notes",
];

const textValues = [
  ["Figure 1", "title", "Evidence-stratified renal papillary context mapping of KSD genetic risk", "Long embedded title competes with the artwork", "Evidence-stratified KSD papillary context map", "partial", "yes", "Remove from artwork entirely if the journal places titles only in legends."],
  ["Figure 1", "subtitle", "Genetic prioritization is linked to parallel, graded papillary context evidence - not cumulative proof", "Explanatory sentence duplicates the visual structure", "Remove", "yes", "yes", "Parallel layout and evidence labels must carry this meaning."],
  ["Figure 1", "A", "GWAS to MAGMA", "Process phrasing encourages a pipeline reading", "Genetic prioritization", "no", "yes", "Use an evidence rail rather than a process arrow."],
  ["Figure 1", "A", "4,915,033 QC rows; 17,316 genes tested; 94 Bonferroni genes; 1000G EUR LD reference", "Key anchors are split across multiple boxes", "4,915,033 QC rows | 17,316 genes | 94 Bonferroni | 1000G EUR LD", "partial", "yes", "Keep all numbers visible and assign MAGMA counts to MAGMA."],
  ["Figure 1", "B", "Two-axis candidate evidence model", "Panel title is longer than needed", "Evidence model", "yes", "yes", "Axes retain the full meaning."],
  ["Figure 1", "B", "42/51 FDR TWAS models are one-SNP", "Long technical phrase", "42/51 one-SNP models", "yes", "yes", "Legend defines FDR-supported Kidney_Cortex models."],
  ["Figure 1", "B", "Curated exemplars; INTERPRETIVE FLAG; Biology-informed; no evidence upgrade", "Three lines repeat one role", "Curated exemplars: interpretive only", "yes", "yes", "Must not become a tier or validation badge."],
  ["Figure 1", "C", "Donor-level snRNA context; 4 donors; Loop/TAL-associated; partial matched-random support", "Four separate labels can be condensed", "snRNA: 4 donors | Loop/TAL-associated | Partial random support", "partial", "yes", "Keep MODERATE CONTEXT as a separate evidence label."],
  ["Figure 1", "D", "Paired bulk disease context; 26 paired patients; injury/remodeling-associated; attenuated after adjustment", "Dense branch text", "Paired bulk: 26 patients | Injury/remodeling-associated | Attenuated after adjustment", "partial", "yes", "Keep ATTENUATED BULK CONTEXT visible."],
  ["Figure 1", "E", "Spatial / TWAS context; 5 sections / 7,747 spots; Loop/TAL co-distribution; No lesion ROI; Kidney_Cortex TWAS proxy; supplementary only", "Six text units crowd the branch", "Spatial: 5 sections; 7,747 spots | Loop/TAL co-distribution | No lesion ROI | TWAS proxy", "partial", "yes", "Keep SUPPLEMENTARY CONTEXT and specify Kidney_Cortex in the proxy badge."],
  ["Figure 1", "F", "MAGMA-prioritized genes map to a Loop/TAL-associated papillary context and an injury/remodeling-associated bulk background", "Full claim is too long for the footer", "Loop/TAL papillary context + injury/remodeling bulk background", "yes", "yes", "Exact locked claim remains in the legend/manuscript."],
  ["Figure 2", "title", "Figure 2 draft v0.2. Donor-level snRNA context mapping of MAGMA-prioritized modules", "Version label and full sentence dominate the page", "Donor-level snRNA context", "yes", "yes", "Remove draft/version text from artwork."],
  ["Figure 2", "subtitle", "Conservative visual refinement: moderate context support; partial matched-random support; descriptive driver-removal robustness", "Repeats later panels and boundaries", "Remove", "yes", "yes", "Retain evidence roles in the relevant panels/footer."],
  ["Figure 2", "A", "Audited GSE231569 context; 43,878 nuclei; Loop/TAL = 540; no cell-level inference", "Audit prose crowds a basic resource panel", "GSE231569 | 43,878 nuclei | 540 Loop/TAL | 4 donors", "partial", "yes", "Explain no cell-level inference in the legend; keep donor unit visible."],
  ["Figure 2", "B", "Donor x compartment scores; donor x compartment summaries; descriptive", "Title and subtitle repeat", "Donor × compartment module scores", "yes", "yes", "Retain all donor and module labels."],
  ["Figure 2", "C", "Loop/TAL within-donor rank; descriptive ranks; no cell-level P values", "Method boundary is too long for the panel", "Loop/TAL rank by donor", "yes", "yes", "Legend states descriptive ranking and absence of cell-level inference."],
  ["Figure 2", "D", "Single-donor exclusion robustness; Tile label = Loop/TAL rank", "Full robustness matrix is secondary", "Leave-one-donor-out robustness", "yes", "yes", "Recommended for supplementary placement."],
  ["Figure 2", "E", "Matched random benchmark; partial matched-random support; delta above matched random; rank metrics saturated; empirical p_delta in source table", "In-plot methods paragraph is unreadable", "Matched-random benchmark | Partial support", "yes", "yes", "Keep values and partial classification; move method details to legend."],
  ["Figure 2", "F", "Known-driver removal; Module-level descriptive robustness", "Secondary robustness panel competes with the anchor", "Known-driver removal", "yes", "yes", "Recommended for supplementary placement."],
  ["Figure 2", "G", "Moderate context support; Donor-level Loop/TAL pattern; Partial matched-random support; Robust driver-removal; Not causal cell type; Not plaque nucleation site", "Large text box duplicates the legend", "Moderate donor-level context; not causal cell type or plaque nucleation site", "yes", "yes", "Use one footer rail, not a standalone panel."],
  ["Figure 3", "title", "Figure 3 draft v0.2. Candidate-gene evidence model and exemplar boundaries", "Version label and long title", "Candidate evidence model", "yes", "yes", "Remove draft/version text."],
  ["Figure 3", "subtitle", "Conservative display: mutually exclusive reporting groups, curated exemplars as flags, no causal or validated gene labels", "Repeats the boundary panel", "Remove", "yes", "yes", "Boundary footer retains the essential limitation."],
  ["Figure 3", "A", "Two-axis evidence model; Genetic priority x Kidney_Cortex proxy TWAS", "Axis concept is repeated", "Genetic priority × Kidney_Cortex proxy", "partial", "yes", "Move exact G/T code mapping to legend/Supplementary Table 5."],
  ["Figure 3", "A", "G1_MAGMA_Bonferroni; G2_MAGMA_FDR_nonBonferroni; G3_MAGMA_suggestive_or_module; G4_no_MAGMA_support", "Internal codes are not reader-facing", "Bonferroni; Lower MAGMA priority; Suggestive/module; No MAGMA support", "yes", "yes", "Exact codes remain in legend."],
  ["Figure 3", "A", "T1_multi_snp_Kidney_Cortex_proxy; T2_one_snp_Kidney_Cortex_proxy; T0_no_TWAS_FDR_support; Tx_not_available", "Internal codes and repeated tissue name are long", "Multi-SNP proxy; One-SNP proxy; No FDR proxy; Unavailable", "yes", "yes", "Axis title must retain Kidney_Cortex proxy."],
  ["Figure 3", "B", "Reporting group counts; Reporting groups, not causal tiers", "Boundary can be stated once", "Reporting groups", "partial", "yes", "Place `not causal tiers` as a small footer note."],
  ["Figure 3", "C", "Curated biological exemplars; Role spectrum; no evidence upgrade", "Three ideas compete with the gene list", "Curated exemplars", "yes", "yes", "Use concise role labels; retain interpretive-only status."],
  ["Figure 3", "D", "Exemplar evidence strip; TWAS proxy does not upgrade exemplars; SMR/coloc not claim-grade ready", "Subtitle contains the whole interpretation", "Exemplar evidence roles", "yes", "yes", "Retain `not claim-grade ready` in the legend or column label."],
  ["Figure 3", "E", "Claim boundary; No causal/validated/therapeutic labels", "Large text-only panel", "Evidence roles, not causal or validated-gene labels", "yes", "yes", "Move full allowed/disallowed list to legend/QC."],
  ["Figure 4", "title", "Figure 4 draft v0.1. GSE73680 provides attenuated paired bulk disease-context support", "Version label and Results-style sentence", "Paired bulk papillary context", "yes", "yes", "Attenuation remains prominent in Panel E/footer."],
  ["Figure 4", "subtitle", "MAGMA-prioritized module shifts track injury/remodeling context and attenuate in marker-signature sensitivity models", "Duplicates panels D/E", "Remove", "yes", "yes", "Legend carries the complete interpretation."],
  ["Figure 4", "A", "Paired GSE73680 design and QC; Control / adjacent papillary bulk; Plaque / stone papillary bulk; 26 paired patients; 55 samples total | 29 patients", "Design panel is word-heavy", "GSE73680 paired design | 26 pairs | 55 samples; 29 patients", "partial", "yes", "Sample labels remain visible; transformation moves to legend."],
  ["Figure 4", "B", "Primary MAGMA module paired deltas; Unadjusted paired bulk shift; each point is one patient", "Title/subtitle repeat the unit", "Paired module-score shifts", "yes", "yes", "Legend states one point per paired patient."],
  ["Figure 4", "C", "Patient-aware bulk sensitivity; Patient fixed-effect estimates before composition-aware adjustment", "Long methodological title", "Patient-aware module effects", "yes", "yes", "Keep 95% CI axis and exact values in Source Data."],
  ["Figure 4", "D", "Bulk marker-signature response; Proxy scores, not validated cell fractions", "Claim boundary and title compete", "Bulk marker-signature shifts | proxies", "yes", "yes", "Legend states not cell fractions."],
  ["Figure 4", "E", "Adjustment reveals attenuation; All four modules: attenuated but direction retained overall", "Title can carry the message more directly", "Adjustment attenuates support", "yes", "yes", "Retain retained/attenuated/lost categories."],
  ["Figure 4", "F", "Injury/remodeling coupling; Paired-delta correlation; coupling is not mechanism", "Detailed correlations are secondary and caveat is repeated", "Bulk coupling, not mechanism", "yes", "yes", "Recommended for supplementary placement."],
  ["Figure 4", "G", "Attenuated bulk-context support; Paired plaque/stone papilla module shift; Linked to injury/remodeling signatures; Not cell-type-specific response; Not genetic causality; Not plaque nucleation mechanism", "Six-box strip is visually heavy", "Attenuated bulk context; not cell-specific, causal or plaque mechanism", "yes", "yes", "Use one boundary footer."],
  ["Figure 5", "title", "Figure 5 draft v0.1. Spatial and TWAS layers provide supplementary papillary context", "Version label and long title", "Supplementary spatial/TWAS context", "yes", "yes", "Remove draft/version text."],
  ["Figure 5", "subtitle", "Loop/TAL co-distribution is retained or attenuated after complexity adjustment; injury/ECM are not retained; no plaque ROI is available", "Full Results sentence crowds the page", "Remove", "yes", "yes", "Null/mixed findings remain visible in Panels C/D and legend."],
  ["Figure 5", "A", "Five spatial sections; no ROI; 5 sections | 7,747 spots | No plaque/mineral/lesion ROI; GSM6250309 marked in amber: QC-selected for display, not correlation", "Resource and selection detail are overpacked", "5 sections | 7,747 spots | No plaque/mineral/lesion ROI", "yes", "yes", "Move selection rule and section counts to legend/Source Data."],
  ["Figure 5", "B", "Representative papillary spatial projection; GSM6250309; visualization only, not lesion localization", "Important boundary is mixed with metadata", "Representative papillary projection | Visualization only; no lesion ROI", "yes", "yes", "Keep sample ID in legend if space is limited."],
  ["Figure 5", "C", "Section-level adjusted co-distribution; CP = consistent; ATT = positive but attenuated; MIX = mixed; NS = not supported", "Abbreviation key and title compete", "Adjusted section-level co-distribution", "yes", "yes", "Use direct short words or symbols; exact definitions in legend."],
  ["Figure 5", "D", "Complexity adjustment changes the context signal; Thin lines = section x primary-module pairs; black line = median; no spot-level P values shown", "Method detail overwhelms the comparison", "Raw versus complexity-adjusted", "yes", "yes", "Legend defines lines, residualization and P-value boundary."],
  ["Figure 5", "E", "TWAS boundary; Supplementary GTEx Kidney_Cortex proxy only; 51 FDR-supported genes: 42 one-SNP | 9 multi-SNP; R2 descriptive; R5 not score-feasible", "Full-width proxy panel is too prominent", "51 FDR genes; 42 one-SNP | Kidney_Cortex proxy", "yes", "yes", "Move full bar and R2/R5 details to supplementary material."],
  ["Figure 5", "F", "Moderate supplementary spatial context; Loop/TAL-associated tissue-context projection; No plaque ROI; Not plaque-specific validation; Not causal spatial niche; TWAS = Kidney_Cortex proxy", "Six-box strip duplicates other boundaries", "Supplementary Loop/TAL context; no lesion localization; Kidney_Cortex TWAS proxy", "yes", "yes", "Use one footer rail."],
];

const textRows = textValues.map((values) => row(textHeaders, values));

const priorityHeaders = [
  "figure", "panel", "priority", "reason", "recommended_action", "requires_source_data_change",
  "requires_claim_change", "human_designer_needed", "notes",
];

const priorityValues = [
  ["Figure 1", "A", "high", "Primary numbers and EUR LD reference need a cleaner genetic-prioritization hierarchy", "Rebuild as an unframed GWAS/MAGMA rail with MAGMA counts grouped correctly", "no", "no", "yes", "Verify all numbers against locked panel sources."],
  ["Figure 1", "B", "high", "Current matrix is conceptually central but visually sparse and boxed", "Rebuild as a light two-axis evidence model with non-causal tags and an offset exemplar flag", "no", "no", "yes", "Do not create a score or tier ladder."],
  ["Figure 1", "C", "high", "snRNA branch must preserve donor unit and partial random support", "Use donor and Loop/TAL motifs with a concise moderate-context label", "no", "no", "yes", "No nuclei-as-replicates cue."],
  ["Figure 1", "D", "medium", "Bulk attenuation should read immediately without a success arrow", "Use paired samples plus a shortening/fading attenuation motif", "no", "no", "yes", "Keep 26 patients and attenuation wording."],
  ["Figure 1", "E", "medium", "Spatial/TWAS branch is text-heavy and vulnerable to lesion misreading", "Use neutral section motifs, no lesion graphic, and a proxy badge", "no", "no", "yes", "Keep no-ROI and Kidney_Cortex boundary."],
  ["Figure 1", "F", "high", "The global interpretation and stopping rules need an elegant footer", "Replace heavy boxes with one boundary rail and four visible limit labels", "no", "no", "yes", "Exact locked claim remains in legend."],
  ["Figure 2", "A", "low", "Resource overview is needed but not the scientific anchor", "Compress to a count strip with donor-level unit", "no", "no", "yes", "Retain 43,878 nuclei and 540 Loop/TAL."],
  ["Figure 2", "B", "high", "Primary donor-by-compartment evidence", "Make the heatmap the hero panel and consolidate its legend", "no", "no", "yes", "All cells and labels remain data-linked."],
  ["Figure 2", "C", "medium", "Within-donor ranks support consistency", "Reduce subtitle and align directly beneath/next to Panel B", "no", "no", "yes", "No cell-level P values."],
  ["Figure 2", "D", "low", "Leave-one-donor-out is robustness detail", "Move full matrix to supplementary material", "no", "no", "yes", "Do not call it external validation."],
  ["Figure 2", "E", "high", "Matched-random benchmark defines the moderate evidence level", "Retain in main figure and remove in-plot methods prose", "no", "no", "yes", "Partial support must remain explicit."],
  ["Figure 2", "F", "medium", "Driver-removal is useful but secondary", "Move full matrix to supplementary material and summarize in legend", "no", "no", "yes", "Preserve all removal outcomes."],
  ["Figure 2", "G", "high", "Large claim box consumes space and repeats prose", "Convert to a narrow boundary footer", "no", "no", "yes", "Keep moderate/not-causal boundaries visible."],
  ["Figure 3", "A", "high", "Two-axis evidence model is the figure's main conceptual contribution", "Enlarge and replace internal codes with reader-facing labels", "no", "no", "yes", "Exact codes remain in legend/supplement."],
  ["Figure 3", "B", "high", "Group counts communicate the reporting structure", "Shorten labels and align with Panel A categories", "no", "no", "yes", "Reporting groups are not causal tiers."],
  ["Figure 3", "C", "medium", "Exemplar roles are biologically useful but text-heavy", "Shorten role descriptions and combine with Panel D layout", "no", "no", "yes", "Exemplar list is locked."],
  ["Figure 3", "D", "high", "Evidence boundaries for exemplars prevent overinterpretation", "Combine with Panel C while preserving separate evidence columns", "no", "no", "yes", "SMR/coloc remains not claim-grade ready."],
  ["Figure 3", "E", "high", "Text-only boundary panel is oversized", "Move full list to legend/QC and use one footer sentence", "no", "no", "yes", "Forbidden labels remain prohibited."],
  ["Figure 4", "A", "medium", "Study design is important but oversized", "Compress to paired-sample header with cohort counts", "no", "no", "yes", "Keep 55 samples, 29 patients and 26 pairs."],
  ["Figure 4", "B", "high", "Patient-level paired shifts are primary evidence", "Keep as a dominant panel with simplified title", "no", "no", "yes", "All points and median markers are locked."],
  ["Figure 4", "C", "high", "Patient-aware estimates provide the inferential summary", "Keep confidence-interval panel and remove repeated annotations", "no", "no", "yes", "Do not alter estimates or CIs."],
  ["Figure 4", "D", "medium", "Marker signatures establish injury/remodeling context", "Retain as a smaller contextual panel with one shared legend", "no", "no", "yes", "Signatures remain bulk proxies."],
  ["Figure 4", "E", "high", "Adjustment attenuation is the central boundary", "Make this the co-anchor and simplify category legend", "no", "no", "yes", "Retain lost/explained outcomes."],
  ["Figure 4", "F", "medium", "Correlation heatmap can be mistaken for mechanism", "Move to supplementary material or clearly label `not mechanism`", "no", "no", "yes", "Exact correlations remain in Source Data."],
  ["Figure 4", "G", "high", "Six-box strip is visually heavy", "Replace with one concise boundary footer", "no", "no", "yes", "Keep not cell-specific/causal/plaque-mechanism limits."],
  ["Figure 5", "A", "high", "No-ROI resource boundary must precede spatial interpretation", "Convert to a compact section/spot/no-ROI header", "no", "no", "yes", "Representative selection rule moves to legend."],
  ["Figure 5", "B", "high", "Spatial image is visually salient and vulnerable to lesion-localization overreading", "Keep image and spots locked; add persistent visualization-only/no-ROI label", "no", "no", "yes", "Do not alter histology or coordinates."],
  ["Figure 5", "C", "high", "Adjusted section-level evidence is the main spatial summary", "Use direct labels, preserve rho values and make Loop/TAL pattern legible", "no", "no", "yes", "Null/mixed cells remain equally visible."],
  ["Figure 5", "D", "high", "Raw-to-adjusted attenuation supports the conservative interpretation", "Simplify method text and preserve all paired lines/medians", "no", "no", "yes", "No spot-level P-value claim."],
  ["Figure 5", "E", "medium", "TWAS proxy panel is too visually dominant for supplementary evidence", "Move full bar to supplement and retain a compact proxy badge in main", "no", "no", "yes", "Keep 51/42/9 counts."],
  ["Figure 5", "F", "high", "Boundary strip duplicates text and consumes space", "Replace with one footer rail", "no", "no", "yes", "No lesion localization or causal niche wording."],
];

const priorityRows = priorityValues.map((values) => row(priorityHeaders, values));

function countWords(value) {
  if (value === "Remove") return 0;
  return value.trim().split(/\s+/).filter(Boolean).length;
}

async function main() {
  await fs.mkdir(outDir, { recursive: true });
  await fs.mkdir(logDir, { recursive: true });
  const logs = [];
  logs.push(await writeAuditedTsv("figure_text_minimization_table_v1.0.tsv", textHeaders, textRows));
  logs.push(await writeAuditedTsv("figure_panel_redesign_priority_table_v1.0.tsv", priorityHeaders, priorityRows));

  const currentWords = textRows.reduce((sum, item) => sum + countWords(item.current_text_or_label), 0);
  const proposedWords = textRows.reduce((sum, item) => sum + countWords(item.recommended_short_text), 0);
  const reduction = currentWords === 0 ? 0 : (1 - proposedWords / currentWords) * 100;
  const summary = [
    `timestamp=${new Date().toISOString()}`,
    "artifact_tool=used_for_in_memory_table_authoring_and_inspection",
    ...logs,
    `current_word_units=${currentWords}`,
    `proposed_word_units=${proposedWords}`,
    `estimated_reduction_percent=${reduction.toFixed(1)}`,
    `all_claim_boundaries_preserved=${textRows.every((item) => item.claim_boundary_preserved === "yes")}`,
    `all_source_data_changes_no=${priorityRows.every((item) => item.requires_source_data_change === "no")}`,
    `all_claim_changes_no=${priorityRows.every((item) => item.requires_claim_change === "no")}`,
  ];
  await fs.writeFile(path.join(logDir, "stage7F_artwork_table_build.log"), summary.join("\n") + "\n", "utf8");
  console.log(summary.join("\n"));
}

await main();
