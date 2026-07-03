#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const outDir = path.join(root, "results/tables/revision/stage7E_figure2_redesign_v0.4");
const logDir = path.join(root, "logs/revision/stage7E_figure2_redesign_v0.4");

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
  "column_count_if_readable", "data_changed_from_v0.3", "claim_changed_from_v0.3",
  "ready_for_publication_source_data", "notes",
];
const counts = {
  A: [6, 6], B: [92, 17], C: [16, 14], D: [16, 12],
  E: [4, 19], F: [20, 17], G: [7, 9],
};
const notes = {
  A: "Frozen resource summary; the papillary nephron/Loop-TAL schematic is display-only and does not encode a new biological result.",
  B: "Frozen donor x compartment module-score extract; compact donor and module labels are documented in the v0.4 display-label mapping.",
  C: "Frozen within-donor Loop/TAL rank summary; rank values are unchanged and shown as a tile matrix.",
  D: "Frozen leave-one-donor-out summary; filled squares replace repeated retained text.",
  E: "Frozen benchmark summary only. Raw random-set distributions were not retained, so v0.4 displays observed-delta random percentiles rather than a fabricated boxplot or null interval.",
  F: "Frozen panel-removal summary; filled squares replace repeated retained text and only the five prespecified panel removals are shown.",
  G: "Frozen claim-boundary extract; wording and evidence strength are unchanged.",
};
const manifestRows = Object.entries(counts).map(([panel, [rows, columns]]) => ({
  panel,
  source_data_file: `results/tables/revision/stage4C2R_draft_figures_v0.2/figure2_panel${panel}_source.tsv`,
  source_exists: "yes",
  row_count_if_readable: rows,
  column_count_if_readable: columns,
  "data_changed_from_v0.3": "no",
  "claim_changed_from_v0.3": "no",
  ready_for_publication_source_data: "yes_draft_requires_final_human_review",
  notes: notes[panel],
}));

const redesignHeaders = [
  "panel", "visual_issue_in_v0.3", "redesign_in_v0.4", "data_changed",
  "claim_changed", "scientific_role", "human_review_needed", "notes",
];
const redesignRows = [
  {
    panel: "A",
    "visual_issue_in_v0.3": "Rigid four-box dashboard with limited visual context and crowded small labels.",
    "redesign_in_v0.4": "Added a restrained vector papillary nephron/Loop-TAL schematic and rebuilt the four KPIs with clearer number-label hierarchy.",
    scientific_role: "Resource scale and donor-by-compartment analysis unit.",
    notes: "Schematic is illustrative only; resource counts remain source-backed.",
  },
  {
    panel: "B",
    "visual_issue_in_v0.3": "Heatmap was informative but did not dominate the page enough relative to secondary panels.",
    "redesign_in_v0.4": "Expanded the heatmap area, shortened donor labels, increased facet spacing and simplified the color-bar title.",
    scientific_role: "Hero evidence for the donor-level Loop/TAL-associated pattern.",
    notes: "Restrained teal scale and unchanged quantitative range.",
  },
  {
    panel: "C",
    "visual_issue_in_v0.3": "Circle marks duplicated the rank value and used more visual ink than needed.",
    "redesign_in_v0.4": "Replaced bubbles with a compact 4 x 4 tile matrix containing direct rank values.",
    scientific_role: "Descriptive within-donor rank summary.",
    notes: "No P values, stars or inferential symbols.",
  },
  {
    panel: "D",
    "visual_issue_in_v0.3": "The word retained was repeated in all 16 cells, weakening pattern recognition.",
    "redesign_in_v0.4": "Replaced repeated text with filled teal squares on a neutral matrix and a single filled-equals-retained key.",
    scientific_role: "Internal single-donor-exclusion robustness.",
    notes: "No external confirmation language.",
  },
  {
    panel: "E",
    "visual_issue_in_v0.3": "The lollipop showed observed deltas but did not visually communicate the matched-random reference structure.",
    "redesign_in_v0.4": "Used a neutral 0-100 matched-random percentile track, a gold observed-delta percentile point and direct partial labels.",
    scientific_role: "Partial expression/detection-matched random-set support.",
    notes: "Summary-only design; no raw null distribution or numerical null band was invented.",
  },
  {
    panel: "F",
    "visual_issue_in_v0.3": "Repeated retained labels made the secondary sensitivity grid text-heavy.",
    "redesign_in_v0.4": "Replaced cell text with filled teal squares and standardized compact removal-scheme labels.",
    scientific_role: "Descriptive driver-panel-removal robustness.",
    notes: "Does not imply mechanistic driver evidence.",
  },
  {
    panel: "G",
    "visual_issue_in_v0.3": "The boundary strip remained visually dense at reduced size.",
    "redesign_in_v0.4": "Widened the strip, grouped supported versus not-claimed labels and rewrapped text into equal-height outline boxes.",
    scientific_role: "Explicit interpretation and overclaim boundary.",
    notes: "Causal cell type and plaque nucleation site remain explicitly not claimed.",
  },
].map((row) => ({
  ...row,
  data_changed: "no",
  claim_changed: "no",
  human_review_needed: "yes",
}));

const logs = [];
logs.push(await writeAuditedTsv(
  "figure2_source_data_manifest_v0.4.tsv", "SourceManifest", manifestHeaders, manifestRows,
));
logs.push(await writeAuditedTsv(
  "figure2_visual_redesign_log_v0.4.tsv", "VisualRedesign", redesignHeaders, redesignRows,
));
logs.push("source_data_changed=no", "claim_changed=no");
await fs.writeFile(
  path.join(logDir, "stage7E_build_figure2_v0.4_tables.log"),
  logs.join("\n") + "\n",
  "utf8",
);
process.stdout.write(logs.join("\n") + "\n");
