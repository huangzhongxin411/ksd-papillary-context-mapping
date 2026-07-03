#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const outDir = path.join(root, "results/tables/revision/stage7E_figure2_editorial_redesign_v0.6");
const logDir = path.join(root, "logs/revision/stage7E_figure2_editorial_redesign_v0.6");

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

async function writeAuditedTsv(filename, sheetName, headers, rows) {
  const matrix = [headers, ...rows.map((row) => headers.map((header) => row[header] ?? ""))];
  const workbook = Workbook.create();
  const sheet = workbook.worksheets.add(sheetName);
  sheet.getRangeByIndexes(0, 0, matrix.length, headers.length).values = matrix;
  const audit = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!A1:${columnName(headers.length)}${matrix.length}`,
    include: "values",
    tableMaxRows: matrix.length,
    tableMaxCols: headers.length,
    maxChars: 15000,
  });
  if (!audit?.ndjson) throw new Error(`artifact-tool inspection failed for ${filename}`);
  const tsv = matrix.map((cells) => cells.map(escapeTsv).join("\t")).join("\n") + "\n";
  await fs.writeFile(path.join(outDir, filename), tsv, "utf8");
  return `${filename}\trows=${rows.length}\tcolumns=${headers.length}\tartifact_tool_audit=pass`;
}

await fs.mkdir(outDir, { recursive: true });
await fs.mkdir(logDir, { recursive: true });

const manifestHeaders = [
  "panel", "source_data_file", "source_exists", "row_count_if_readable",
  "column_count_if_readable", "data_changed_from_v0.5", "claim_changed_from_v0.5",
  "main_figure_or_source_detail", "notes",
];
const manifestRows = [
  {
    panel: "A",
    source_data_file: "results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelA_source.tsv",
    source_exists: "yes", row_count_if_readable: "6", column_count_if_readable: "6",
    "data_changed_from_v0.5": "no", "claim_changed_from_v0.5": "no",
    main_figure_or_source_detail: "main_figure_summary",
    notes: "Resource counts and donor-by-compartment unit are unchanged. The tiny MAGMA-input glyph is schematic only: no axes, threshold, locus label, discovery claim, fine mapping or genetic-validation meaning.",
  },
  {
    panel: "B",
    source_data_file: "results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelB_source.tsv",
    source_exists: "yes", row_count_if_readable: "92", column_count_if_readable: "17",
    "data_changed_from_v0.5": "no", "claim_changed_from_v0.5": "no",
    main_figure_or_source_detail: "main_figure_quantitative_hero",
    notes: "All four heatmaps, donor labels, compartments, score values and scale limits are unchanged.",
  },
  {
    panel: "C",
    source_data_file: "results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelC_source.tsv; results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelD_source.tsv; results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelF_source.tsv",
    source_exists: "yes;yes;yes", row_count_if_readable: "16;16;20", column_count_if_readable: "14;12;17",
    "data_changed_from_v0.5": "no", "claim_changed_from_v0.5": "no",
    main_figure_or_source_detail: "main_figure_summary; detailed_matrices_in_source_data",
    notes: "Old rank, single-donor-exclusion and driver-panel-removal matrices were summarized visually as 4/4, 4/4 and 5/5 cards. All 16, 16 and 20 detailed rows remain unchanged in Source Data.",
  },
  {
    panel: "D",
    source_data_file: "results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelE_source.tsv",
    source_exists: "yes", row_count_if_readable: "4", column_count_if_readable: "19",
    "data_changed_from_v0.5": "no", "claim_changed_from_v0.5": "no",
    main_figure_or_source_detail: "main_figure_summary_only_benchmark",
    notes: "Source-backed observed-delta percentiles only. Raw random-set distributions remain unavailable; no boxplot, violin, null band or interval was fabricated.",
  },
  {
    panel: "E",
    source_data_file: "results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panelG_source.tsv",
    source_exists: "yes", row_count_if_readable: "7", column_count_if_readable: "9",
    "data_changed_from_v0.5": "no", "claim_changed_from_v0.5": "no",
    main_figure_or_source_detail: "main_figure_claim_boundary",
    notes: "Supported and not-claimed interpretations are unchanged; the strip is widened and simplified editorially.",
  },
];

const redesignHeaders = [
  "old_panel", "new_panel", "issue_in_v0.5", "redesign_action_v0.6",
  "data_changed", "claim_changed", "reason", "notes",
];
const redesignRows = [
  {
    old_panel: "A", new_panel: "A",
    "issue_in_v0.5": "Resource metrics were still presented as a compact dashboard and did not explicitly show the restrained genetic-input handoff.",
    "redesign_action_v0.6": "Converted the panel into a narrow editorial strip with a tiny axis-free MAGMA-module glyph, one-way arrow, renal papilla snRNA label, inline metrics and analysis-unit note.",
    reason: "Establish input and analysis unit without duplicating Figure 1 or competing with the heatmap.",
    notes: "No full Manhattan plot; no discovery, fine-mapping, causal-locus or genetic-validation implication.",
  },
  {
    old_panel: "B", new_panel: "B",
    "issue_in_v0.5": "The central heatmaps were still constrained by seven-panel page competition.",
    "redesign_action_v0.6": "Expanded Panel B to nine of twelve layout columns and retained all four heatmaps with a compact shared color bar.",
    reason: "Make the donor-by-compartment pattern the unmistakable main result.",
    notes: "Quantitative values and scale are unchanged.",
  },
  {
    old_panel: "C/D/F", new_panel: "C",
    "issue_in_v0.5": "Three repetitive rank/retained matrices made the figure read like an audit dashboard.",
    "redesign_action_v0.6": "Collapsed the matrices into three calm source-backed cards: top-ranked in 4/4 donors, retained after 4/4 donor exclusions and retained after 5/5 removal panels.",
    reason: "Keep robustness visible while deferring repetitive detail to Source Data.",
    notes: "Detailed 16-row rank, 16-row exclusion and 20-row removal matrices remain unchanged and traceable.",
  },
  {
    old_panel: "E", new_panel: "D",
    "issue_in_v0.5": "The matched-random panel was one of several small audit panels rather than a clear secondary evidence layer.",
    "redesign_action_v0.6": "Expanded the percentile plot and consolidated the interpretation into one non-overlapping annotation box.",
    reason: "Communicate partial support and rank saturation without strong-enrichment cues.",
    notes: "No p_delta headline, stars or fabricated random distribution.",
  },
  {
    old_panel: "G", new_panel: "E",
    "issue_in_v0.5": "The interpretation boundary shared a row with another panel and remained relatively compressed.",
    "redesign_action_v0.6": "Moved the boundary to a full-width bottom strip with three supported and two not-claimed boxes.",
    reason: "End the visual argument clearly while keeping the boundary subordinate to Panel B.",
    notes: "Claim strength is unchanged.",
  },
].map((row) => ({ ...row, data_changed: "no", claim_changed: "no" }));

const logs = [];
logs.push(await writeAuditedTsv(
  "figure2_source_data_manifest_v0.6.tsv", "SourceManifest", manifestHeaders, manifestRows,
));
logs.push(await writeAuditedTsv(
  "figure2_visual_redesign_log_v0.6.tsv", "VisualRedesign", redesignHeaders, redesignRows,
));
logs.push("source_data_changed_from_v0.5=no", "claim_changed_from_v0.5=no");
await fs.writeFile(
  path.join(logDir, "stage7E_build_figure2_v0.6_tables.log"),
  logs.join("\n") + "\n",
  "utf8",
);
process.stdout.write(logs.join("\n") + "\n");
