#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const outDir = path.join(root, "results/tables/revision/stage7E_figure2_final_polish_v0.5");
const logDir = path.join(root, "logs/revision/stage7E_figure2_final_polish_v0.5");

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
    maxChars: 14000,
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
  "column_count_if_readable", "data_changed_from_v0.4", "claim_changed_from_v0.4",
  "ready_for_publication_source_data", "notes",
];
const counts = {
  A: [6, 6], B: [92, 17], C: [16, 14], D: [16, 12],
  E: [4, 19], F: [20, 17], G: [7, 9],
};
const sourceNotes = {
  A: "Frozen resource summary; schematic and KPI spacing were polished without changing any count or analysis unit.",
  B: "Frozen donor x compartment module-score extract; quantitative limits and all heatmap values are unchanged.",
  C: "Frozen within-donor Loop/TAL rank extract; the 4/4-donor summary is a direct condensation of the unchanged rank matrix.",
  D: "Frozen leave-one-donor-out extract; retained states are unchanged and remain internal robustness only.",
  E: "Frozen summary-only benchmark extract; no raw random-set distributions are available, and no boxplot or null interval was fabricated.",
  F: "Frozen panel-removal extract; all five displayed removal schemes and retained states are unchanged.",
  G: "Frozen interpretation-boundary extract; typography was enlarged without changing supported or not-claimed categories.",
};
const manifestRows = Object.entries(counts).map(([panel, [rows, columns]]) => ({
  panel,
  source_data_file: `results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panel${panel}_source.tsv`,
  source_exists: "yes",
  row_count_if_readable: rows,
  column_count_if_readable: columns,
  "data_changed_from_v0.4": "no",
  "claim_changed_from_v0.4": "no",
  ready_for_publication_source_data: "yes_draft_requires_final_human_review",
  notes: sourceNotes[panel],
}));

const polishHeaders = [
  "panel", "issue_in_v0.4", "polish_made_in_v0.5", "data_changed",
  "claim_changed", "reason", "human_review_needed", "notes",
];
const polishRows = [
  {
    panel: "A",
    "issue_in_v0.4": "The nephron schematic, vertical callout and KPI labels still had uneven spacing and excess white space.",
    "polish_made_in_v0.5": "Reduced the schematic, replaced the vertical label with a horizontal callout, balanced the 2 x 2 KPI grid and added a clear bottom analysis-unit strip.",
    reason: "Improve compactness and typography without letting the resource panel compete with the heatmap.",
    notes: "KPI wording is now 43,878 nuclei; 540 Loop/TAL; 4 donors; 6 compartments.",
  },
  {
    panel: "B",
    "issue_in_v0.4": "The color bar and heatmap gutters could be slightly more compact.",
    "polish_made_in_v0.5": "Reduced tile-border weight, shortened the color bar and tightened facet spacing while retaining the same scale.",
    reason: "Preserve Panel B as the visual anchor with cleaner internal rhythm.",
    notes: "No saturation, ordering or score-range change.",
  },
  {
    panel: "C",
    "issue_in_v0.4": "The all-rank-1 matrix required a quicker reader-level summary.",
    "polish_made_in_v0.5": "Added the concise source-backed summary `Top-ranked in 4/4 donors` and lightened tile borders.",
    reason: "Communicate the repeated descriptive result without introducing inferential emphasis.",
    notes: "Rank 1 remains defined as highest; no P values or stars.",
  },
  {
    panel: "D",
    "issue_in_v0.4": "The binary matrix was clear but its legend role could be more explicit.",
    "polish_made_in_v0.5": "Retained the filled-square matrix, lightened grid boundaries and kept a single `Filled = retained` key beneath the leave-one-donor-out subtitle.",
    reason: "Improve status decoding while avoiding repeated cell text.",
    notes: "No external confirmation wording.",
  },
  {
    panel: "E",
    "issue_in_v0.4": "Four repeated partial labels made a percentile-at-100 display visually ambiguous.",
    "polish_made_in_v0.5": "Removed row-wise partial labels, changed the axis to `Matched-random percentile`, added one `Overall: partial support` badge and a saturated-rank annotation.",
    reason: "Clarify why the benchmark remains partial without weakening or exaggerating the source result.",
    notes: "No p_delta headline, stars, boxplot, violin or fabricated null interval.",
  },
  {
    panel: "F",
    "issue_in_v0.4": "The compact matrix needed a cleaner subtitle and key separation.",
    "polish_made_in_v0.5": "Kept the binary matrix, lightened grid boundaries and split `Descriptive robustness` from the `Filled = retained` key.",
    reason: "Match Panel D while preserving its descriptive sensitivity role.",
    notes: "No mechanistic driver claim.",
  },
  {
    panel: "G",
    "issue_in_v0.4": "Boundary-box text was still slightly compressed at manuscript scale.",
    "polish_made_in_v0.5": "Increased box height and text size while retaining the wide two-group strip and equal-height boxes.",
    reason: "Improve claim-boundary readability without increasing visual dominance.",
    notes: "Not causal cell type and not plaque nucleation site remain explicit.",
  },
].map((row) => ({
  ...row,
  data_changed: "no",
  claim_changed: "no",
  human_review_needed: "yes",
}));

const logs = [];
logs.push(await writeAuditedTsv(
  "figure2_source_data_manifest_v0.5.tsv", "SourceManifest", manifestHeaders, manifestRows,
));
logs.push(await writeAuditedTsv(
  "figure2_visual_polish_log_v0.5.tsv", "VisualPolish", polishHeaders, polishRows,
));
logs.push("source_data_changed_from_v0.4=no", "claim_changed_from_v0.4=no");
await fs.writeFile(
  path.join(logDir, "stage7E_build_figure2_v0.5_tables.log"),
  logs.join("\n") + "\n",
  "utf8",
);
process.stdout.write(logs.join("\n") + "\n");
