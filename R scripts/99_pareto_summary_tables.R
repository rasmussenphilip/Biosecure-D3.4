library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ============================================================
# 15_pareto_summary_tables.R
# Export Pareto package summary tables by herd type
# ============================================================

out_dir <- "exports"

pareto <- read_csv(
  file.path(out_dir, "package_pareto_frontier.csv"),
  show_col_types = FALSE
)

net_benefits <- read_csv(
  file.path(out_dir, "package_net_benefits.csv"),
  show_col_types = FALSE
)

package_risk_reductions <- read_csv(
  file.path(out_dir, "package_risk_reductions.csv"),
  show_col_types = FALSE
)

# ============================================================
# 1. Recalculate payoff IQRs from package_net_benefits.csv
# ============================================================

payoff_iqr <- net_benefits %>%
  group_by(
    package_id,
    package_size,
    package_contents,
    benefit_type
  ) %>%
  summarise(
    payoff_iqr_low  = quantile(net_benefit_ref, 0.25, na.rm = TRUE),
    payoff_iqr_high = quantile(net_benefit_ref, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(
    context_benefit_type = benefit_type
  )

# ============================================================
# 2. Recalculate risk-reduction summaries, including IQRs
# ============================================================

risk_meta_cols <- c("package_id", "package_size", "package_contents")
risk_value_cols <- setdiff(names(package_risk_reductions), risk_meta_cols)

risk_long <- package_risk_reductions %>%
  pivot_longer(
    cols = all_of(risk_value_cols),
    names_to = "risk_col",
    values_to = "risk_value"
  ) %>%
  mutate(risk_col = str_to_lower(risk_col))

parsed_cols <- str_match(
  risk_long$risk_col,
  "^post_(naive|infected|average)_risk_(bvd|ptb|sal)_(.*)$"
)

risk_long <- risk_long %>%
  mutate(
    risk_type = parsed_cols[, 2],
    disease = parsed_cols[, 3],
    suffix = parsed_cols[, 4]
  ) %>%
  filter(
    !is.na(risk_type),
    !is.na(disease),
    !is.na(suffix)
  )

risk_states <- risk_long %>%
  select(
    package_id,
    package_size,
    package_contents,
    suffix,
    risk_type,
    disease,
    risk_value
  ) %>%
  unite(
    col = "risk_name",
    risk_type,
    disease,
    sep = "_"
  ) %>%
  pivot_wider(
    names_from = risk_name,
    values_from = risk_value
  )

empty_risk_states <- risk_states %>%
  filter(package_contents == "empty_package") %>%
  select(
    suffix,
    baseline_naive_bvd = naive_bvd,
    baseline_naive_ptb = naive_ptb,
    baseline_naive_sal = naive_sal,
    baseline_infected_bvd = infected_bvd,
    baseline_infected_ptb = infected_ptb,
    baseline_infected_sal = infected_sal,
    baseline_average_bvd = average_bvd,
    baseline_average_ptb = average_ptb,
    baseline_average_sal = average_sal
  )

risk_reduction_summary <- risk_states %>%
  left_join(empty_risk_states, by = "suffix") %>%
  mutate(
    naive_bvd_reduction = baseline_naive_bvd - naive_bvd,
    naive_ptb_reduction = baseline_naive_ptb - naive_ptb,
    naive_sal_reduction = baseline_naive_sal - naive_sal,
    
    infected_bvd_reduction = baseline_infected_bvd - infected_bvd,
    infected_ptb_reduction = baseline_infected_ptb - infected_ptb,
    infected_sal_reduction = baseline_infected_sal - infected_sal,
    
    average_bvd_reduction = baseline_average_bvd - average_bvd,
    average_ptb_reduction = baseline_average_ptb - average_ptb,
    average_sal_reduction = baseline_average_sal - average_sal
  ) %>%
  group_by(
    package_id,
    package_size,
    package_contents
  ) %>%
  summarise(
    across(
      c(
        naive_bvd_reduction,
        naive_ptb_reduction,
        naive_sal_reduction,
        infected_bvd_reduction,
        infected_ptb_reduction,
        infected_sal_reduction,
        average_bvd_reduction,
        average_ptb_reduction,
        average_sal_reduction
      ),
      list(
        iqr_low  = ~quantile(.x, 0.25, na.rm = TRUE),
        median   = ~median(.x, na.rm = TRUE),
        iqr_high = ~quantile(.x, 0.75, na.rm = TRUE),
        p05      = ~quantile(.x, 0.05, na.rm = TRUE),
        p95      = ~quantile(.x, 0.95, na.rm = TRUE)
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )

# ============================================================
# 3. Join payoff and risk summaries to Pareto table
# ============================================================

pareto <- pareto %>%
  select(
    -contains("_reduction")
  ) %>%
  left_join(
    payoff_iqr,
    by = c(
      "package_id",
      "package_size",
      "package_contents",
      "context_benefit_type"
    )
  ) %>%
  left_join(
    risk_reduction_summary,
    by = c(
      "package_id",
      "package_size",
      "package_contents"
    )
  )

# ============================================================
# 4. Helper function
# ============================================================

make_herd_table <- function(df, herd_type, risk_prefix) {
  
  df %>%
    filter(context_benefit_type == herd_type) %>%
    transmute(
      analysis_level,
      context_region,
      context_benefit_type,
      package_id,
      package_size,
      package_contents,
      
      n_payoffs,
      
      payoff_median = median_payoff,
      payoff_iqr_low,
      payoff_iqr_high,
      payoff_p05 = p05_payoff,
      payoff_p95 = p95_payoff,
      
      bvd_risk_reduction_median =
        .data[[paste0("median_", risk_prefix, "_bvd_reduction")]],
      bvd_risk_reduction_iqr_low =
        .data[[paste0("iqr_low_", risk_prefix, "_bvd_reduction")]],
      bvd_risk_reduction_iqr_high =
        .data[[paste0("iqr_high_", risk_prefix, "_bvd_reduction")]],
      bvd_risk_reduction_p05 =
        .data[[paste0("p05_", risk_prefix, "_bvd_reduction")]],
      bvd_risk_reduction_p95 =
        .data[[paste0("p95_", risk_prefix, "_bvd_reduction")]],
      
      ptb_risk_reduction_median =
        .data[[paste0("median_", risk_prefix, "_ptb_reduction")]],
      ptb_risk_reduction_iqr_low =
        .data[[paste0("iqr_low_", risk_prefix, "_ptb_reduction")]],
      ptb_risk_reduction_iqr_high =
        .data[[paste0("iqr_high_", risk_prefix, "_ptb_reduction")]],
      ptb_risk_reduction_p05 =
        .data[[paste0("p05_", risk_prefix, "_ptb_reduction")]],
      ptb_risk_reduction_p95 =
        .data[[paste0("p95_", risk_prefix, "_ptb_reduction")]],
      
      sal_risk_reduction_median =
        .data[[paste0("median_", risk_prefix, "_sal_reduction")]],
      sal_risk_reduction_iqr_low =
        .data[[paste0("iqr_low_", risk_prefix, "_sal_reduction")]],
      sal_risk_reduction_iqr_high =
        .data[[paste0("iqr_high_", risk_prefix, "_sal_reduction")]],
      sal_risk_reduction_p05 =
        .data[[paste0("p05_", risk_prefix, "_sal_reduction")]],
      sal_risk_reduction_p95 =
        .data[[paste0("p95_", risk_prefix, "_sal_reduction")]],
      
      frontier_rank_expected_value,
      frontier_rank_maximin,
      frontier_rank_maximax
    ) %>%
    arrange(
      analysis_level,
      context_region,
      frontier_rank_expected_value,
      package_size,
      package_id
    )
}

# ============================================================
# 5. Create one table per herd type
# ============================================================

pareto_naive <- make_herd_table(
  df = pareto,
  herd_type = "naive_herd",
  risk_prefix = "naive"
)

pareto_infected <- make_herd_table(
  df = pareto,
  herd_type = "infected_herd",
  risk_prefix = "infected"
)

pareto_average <- make_herd_table(
  df = pareto,
  herd_type = "average_herd",
  risk_prefix = "average"
)

# ============================================================
# 6. Export tables
# ============================================================

write_csv(
  pareto_naive,
  file.path(out_dir, "pareto_packages_naive_herds_summary.csv")
)

write_csv(
  pareto_infected,
  file.path(out_dir, "pareto_packages_infected_herds_summary.csv")
)

write_csv(
  pareto_average,
  file.path(out_dir, "pareto_packages_average_herds_summary.csv")
)

# ============================================================
# 7. Export top 5 packages by median net benefit
# (All regions only, formatted for reporting)
# ============================================================

measures_final <- read_csv(
  file.path("inputs", "measures_final.csv"),
  show_col_types = FALSE
)

measure_lookup <- measures_final %>%
  select(measure, name) %>%
  distinct()

clean_package_contents <- function(package_contents, lookup = measure_lookup) {
  
  if (is.na(package_contents)) {
    return(package_contents)
  }
  
  if (package_contents == "empty_package") {
    return("None")
  }
  
  measures <- str_split(package_contents, "\\s*\\+\\s*")[[1]]
  
  clean_names <- lookup %>%
    filter(measure %in% measures) %>%
    mutate(
      measure = factor(
        measure,
        levels = measures
      )
    ) %>%
    arrange(measure) %>%
    pull(name)
  
  paste(clean_names, collapse = "; ")
}

format_median_iqr <- function(median, p25, p75, digits = 3) {
  
  paste0(
    round(median, digits),
    " (",
    round(p25, digits),
    "-",
    round(p75, digits),
    ")"
  )
}

make_top5_table <- function(df, herd_label) {
  
  df %>%
    filter(context_region == "All regions") %>%
    arrange(
      desc(payoff_median),
      package_size,
      package_id
    ) %>%
    slice_head(n = 5) %>%
    mutate(
      package_contents = vapply(
        package_contents,
        clean_package_contents,
        character(1)
      )
    ) %>%
    transmute(
      `Herd type` = herd_label,
      `Package ID` = package_id,
      Size = package_size,
      `Package contents` = package_contents,
      
      `BVD; Median (IQR)` = format_median_iqr(
        bvd_risk_reduction_median,
        bvd_risk_reduction_iqr_low,
        bvd_risk_reduction_iqr_high
      ),
      
      `PTB; Median (IQR)` = format_median_iqr(
        ptb_risk_reduction_median,
        ptb_risk_reduction_iqr_low,
        ptb_risk_reduction_iqr_high
      ),
      
      `SAL; Median (IQR)` = format_median_iqr(
        sal_risk_reduction_median,
        sal_risk_reduction_iqr_low,
        sal_risk_reduction_iqr_high
      ),
      
      `Net benefit (EUR/cow-year); Median (IQR)` = format_median_iqr(
        payoff_median,
        payoff_iqr_low,
        payoff_iqr_high,
        digits = 1
      )
    )
}

top5_pareto_packages_by_herd_type <- bind_rows(
  make_top5_table(pareto_naive, "Naive herds"),
  make_top5_table(pareto_infected, "Infected herds"),
  make_top5_table(pareto_average, "Average herds")
) %>%
  mutate(
    `Herd type` = factor(
      `Herd type`,
      levels = c(
        "Naive herds",
        "Infected herds",
        "Average herds"
      )
    )
  ) %>%
  arrange(`Herd type`)

write_csv(
  top5_pareto_packages_by_herd_type,
  file.path(out_dir, "top5_pareto_packages_by_herd_type.csv")
)

# ============================================================
# 8. Export supplementary table with all summary-table rows
# ============================================================

make_supplementary_table <- function(df, herd_label) {
  
  df %>%
    mutate(
      package_contents = vapply(
        package_contents,
        clean_package_contents,
        character(1)
      )
    ) %>%
    transmute(
      Region = context_region,
      `Herd type` = herd_label,
      `Package ID` = package_id,
      Size = package_size,
      `Package contents` = package_contents,
      
      `BVD; Median (IQR)` = format_median_iqr(
        bvd_risk_reduction_median,
        bvd_risk_reduction_iqr_low,
        bvd_risk_reduction_iqr_high
      ),
      
      `PTB; Median (IQR)` = format_median_iqr(
        ptb_risk_reduction_median,
        ptb_risk_reduction_iqr_low,
        ptb_risk_reduction_iqr_high
      ),
      
      `SAL; Median (IQR)` = format_median_iqr(
        sal_risk_reduction_median,
        sal_risk_reduction_iqr_low,
        sal_risk_reduction_iqr_high
      ),
      
      `Net benefit (EUR/cow-year); Median (IQR)` = format_median_iqr(
        payoff_median,
        payoff_iqr_low,
        payoff_iqr_high,
        digits = 1
      )
    )
}

supplementary_pareto_packages <- bind_rows(
  make_supplementary_table(pareto_naive, "Naive herds"),
  make_supplementary_table(pareto_infected, "Infected herds"),
  make_supplementary_table(pareto_average, "Average herds")
) %>%
  mutate(
    `Herd type` = factor(
      `Herd type`,
      levels = c(
        "Naive herds",
        "Infected herds",
        "Average herds"
      )
    )
  ) %>%
  arrange(
    `Herd type`,
    Region,
    desc(`Net benefit (EUR/cow-year); Median (IQR)`)
  )

write_csv(
  supplementary_pareto_packages,
  file.path(out_dir, "supplementary_pareto_packages.csv")
)