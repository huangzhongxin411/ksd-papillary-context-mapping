#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const outDir = path.join(root, "results/tables/revision/stage7E_figure2_refinement");
const logDir = path.join(root, "logs/revision/stage7E_figure2_refinement");

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
    maxChars: 12000,
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
  "column_count_if_readable", "data_changed_from_v0.2", "claim_changed_from_v0.2",
  "ready_for_publication_source_data", "notes",
];

const sourceCounts = {
  A: [6, 6], B: [92, 17], C: [16, 14], D: [16, 12],
  E: [4, 19], F: [20, 17], G: [7, 9],
};

const sourceNotes = {
  A: "Frozen resource and compartment summary; compact card is a display-only redesign.",
  B: "Frozen donor x compartment module-score extract; reader-facing labels are documented in figure2_display_label_mapping_v0.3.tsv.",
  C: "Frozen within-donor Loop/TAL rank summary; no cell-level P values are displayed.",
  D: "Frozen leave-one-donor-out summary; retained status is presented as internal robustness only.",
  E: "Frozen matched-random benchmark extract; partial-support classification is preserved and p_delta is not foregrounded.",
  F: "Frozen driver-removal sensitivity extract; main figure shows the five prespecified panel removals only.",
  G: "Frozen interpretation-boundary extract; wording is condensed without changing the allowed claim.",
};

const manifestRows = Object.entries(sourceCounts).map(([panel, [rows, columns]]) => ({
  panel,
  source_data_file: `results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panel${panel}_source.tsv`,
  source_exists: "yes",
  row_count_if_readable: rows,
  column_count_if_readable: columns,
  "data_changed_from_v0.2": "no",
  "claim_changed_from_v0.2": "no",
  ready_for_publication_source_data: "yes_draft_requires_final_human_review",
  notes: sourceNotes[panel],
}));

const changeHeaders = [
  "panel", "visual_issue_in_v0.2", "change_made_in_v0.3", "data_changed",
  "claim_changed", "reason", "human_review_needed", "notes",
];

const changeRows = [
  {
    panel: "A",
    "visual_issue_in_v0.2": "Detailed bar-chart treatment and audit prose made a resource panel visually dense.",
    "change_made_in_v0.3": "Replaced with a compact 2 x 2 resource card showing nuclei, Loop/TAL nuclei, donors and broad compartments.",
    reason: "Keep the sampling unit and resource scale legible without implying nuclei-level replication.",
    notes: "Donor x compartment summaries remain explicitly stated.",
  },
  {
    panel: "B",
    "visual_issue_in_v0.2": "The heatmap was visually dense and relied on source-style module labels.",
    "change_made_in_v0.3": "Made the heatmap the visual anchor, shortened labels through a documented display mapping and restrained the teal scale.",
    reason: "Improve donor-by-compartment comparison while preserving moderate visual emphasis.",
    notes: "All four primary modules and all source compartments are retained.",
  },
  {
    panel: "C",
    "visual_issue_in_v0.2": "Redundant encodings and explanatory text obscured the simple descriptive rank result.",
    "change_made_in_v0.3": "Used direct rank dots with rank 1 labels and a concise no-cell-level-P-values subtitle.",
    reason: "Make the descriptive rank interpretation immediate and avoid significance cues.",
    notes: "No stars, brackets or P values are shown.",
  },
  {
    panel: "D",
    "visual_issue_in_v0.2": "Bright pass-fail coloring risked suggesting external confirmation.",
    "change_made_in_v0.3": "Used muted teal retained-status tiles under a single-donor-exclusion title.",
    reason: "Present internal robustness without elevating it to external evidence.",
    notes: "The panel uses retained language only.",
  },
  {
    panel: "E",
    "visual_issue_in_v0.2": "Dense in-panel methods and significance emphasis could overstate the matched-random result.",
    "change_made_in_v0.3": "Used an amber lollipop display of observed deltas with direct Partial labels and a rank-metrics-saturated subtitle.",
    reason: "Keep the key overclaim boundary visible and subordinate p_delta to the source table.",
    notes: "The visual conclusion remains partial matched-random support.",
  },
  {
    panel: "F",
    "visual_issue_in_v0.2": "Bright status colors and dense removal labels made a secondary sensitivity panel too dominant.",
    "change_made_in_v0.3": "Used muted retained-status tiles and limited the main panel to five prespecified panel removals.",
    reason: "Show descriptive robustness compactly without implying mechanistic driver evidence.",
    notes: "Single-gene removals remain in source data rather than the main panel.",
  },
  {
    panel: "G",
    "visual_issue_in_v0.2": "A large prose box duplicated the legend and competed with quantitative panels.",
    "change_made_in_v0.3": "Replaced it with a low-dominance five-box interpretation strip separating supported interpretation from not claimed.",
    reason: "Make claim boundaries scannable while keeping Panels B-F dominant.",
    notes: "Not causal cell type and not plaque nucleation site are explicit.",
  },
].map((row) => ({
  ...row,
  data_changed: "no",
  claim_changed: "no",
  human_review_needed: "yes",
}));

const logLines = [];
logLines.push(await writeAuditedTsv(
  "figure2_source_data_manifest_v0.3.tsv", "SourceManifest", manifestHeaders, manifestRows,
));
logLines.push(await writeAuditedTsv(
  "figure2_visual_change_log_v0.3.tsv", "VisualChanges", changeHeaders, changeRows,
));
logLines.push("source_data_changed=no", "claim_changed=no");
await fs.writeFile(
  path.join(logDir, "stage7E_build_figure2_audit_tables_v0.3.log"),
  logLines.join("\n") + "\n",
  "utf8",
);
process.stdout.write(logLines.join("\n") + "\n");
