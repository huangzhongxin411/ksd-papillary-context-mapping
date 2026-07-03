suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

out_dir <- "results/spatial/phase28c_projection"
raw <- read_tsv(file.path(out_dir, "spatial_score_correlations_by_sample_v0.1.tsv"), show_col_types = FALSE)
adj <- read_tsv(file.path(out_dir, "spatial_score_correlations_adjusted_v0.1.tsv"), show_col_types = FALSE)

summary <- raw %>%
  select(sample_id, risk_signature, context_signature, rho_raw) %>%
  full_join(adj %>% select(sample_id, risk_signature, context_signature, rho_adjusted),
            by = c("sample_id", "risk_signature", "context_signature")) %>%
  group_by(risk_signature, context_signature) %>%
  summarise(
    n_samples_available = sum(is.finite(rho_adjusted)),
    median_rho_raw = median(rho_raw, na.rm = TRUE),
    median_rho_adjusted = median(rho_adjusted, na.rm = TRUE),
    positive_fraction_raw = mean(rho_raw > 0, na.rm = TRUE),
    positive_fraction_adjusted = mean(rho_adjusted > 0, na.rm = TRUE),
    min_rho_adjusted = min(rho_adjusted, na.rm = TRUE),
    max_rho_adjusted = max(rho_adjusted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    support_class = case_when(
      n_samples_available >= 4 & median_rho_adjusted >= 0.25 & positive_fraction_adjusted >= 0.80 ~ "strong_spatial_context_support",
      n_samples_available >= 4 & median_rho_adjusted >= 0.15 & positive_fraction_adjusted >= 0.60 ~ "moderate_spatial_context_support",
      TRUE ~ "no_consistent_support"
    )
  ) %>%
  arrange(desc(median_rho_adjusted), risk_signature, context_signature)

summary <- summary %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))
write_tsv(summary, file.path(out_dir, "spatial_projection_summary_v0.1.tsv"))
