suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

source_ok <- function(expr) {
  out <- tryCatch({ force(expr); TRUE }, error = function(e) FALSE)
  out
}

theme_ok <- source_ok(source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R"))
template_files <- list.files("scripts/09_manuscript/templates", pattern = "[.]R$", full.names = TRUE)
templates_ok <- all(vapply(template_files, function(f) source_ok(source(f)), logical(1)))
generic_profile_ok <- file.exists(".codex/skills/scientific_figure_design/project_profiles/generic_bioinformatics.profile.yml")
ksd_profile_ok <- file.exists(".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml")

out <- data.table(
  check = c("skill_md_exists", "reusable_rules_exist", "project_profiles_exist", "theme_v0.3_source", "templates_source", "generic_toy_demo", "ksd_toy_demo", "visual_qc_v0.2"),
  status = c(
    file.exists(".codex/skills/scientific_figure_design/SKILL.md"),
    dir.exists(".codex/skills/scientific_figure_design/reusable_rules"),
    generic_profile_ok && ksd_profile_ok,
    theme_ok,
    templates_ok,
    file.exists("results/figures/phase16b_generic_toy_demo.pdf") && file.exists("results/figures/phase16b_generic_toy_demo.png"),
    file.exists("results/figures/phase16b_ksd_toy_demo.pdf") && file.exists("results/figures/phase16b_ksd_toy_demo.png"),
    file.exists("results/tables/main_figure_visual_qc_v0.6.tsv")
  )
)
out[, note := fifelse(status, "pass", "review_needed")]
fwrite(out, "results/tables/phase16b_skill_validation_v0.1.tsv", sep = "\t")
message("Wrote results/tables/phase16b_skill_validation_v0.1.tsv")
