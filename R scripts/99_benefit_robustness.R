library(dplyr)
library(readr)
library(stringr)

# ============================================================
# 16_supplementary_benefit_robustness.R
# Export supplementary robustness table for package net benefits
# ============================================================

out_dir <- "exports"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

include_empty_package <- FALSE

net_benefits <- read_csv(
  file.path(out_dir, "package_net_benefits.csv"),
  show_col_types = FALSE
)

# ============================================================
# 1. Checks
# ============================================================

required_cols <- c(
  "package_id",
  "package_size",
  "package_contents",
  "region",
  "benefit_type",
  "aggregation_method",
  "joint_state",
  "payoff_state",
  "net_benefit_ref"
)

missing_cols <- setdiff(required_cols, names(net_benefits))

if (length(missing_cols) > 0) {
  stop(
    "package_net_benefits.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# 2. Clean herd type labels
# ============================================================

nb <- net_benefits %>%
  mutate(
    herd_type = case_when(
      str_detect(str_to_lower(benefit_type), "naive|naïve") ~ "Naive",
      str_detect(str_to_lower(benefit_type), "infected") ~ "Infected",
      str_detect(str_to_lower(benefit_type), "average") ~ "Average",
      TRUE ~ benefit_type
    )
  )

if (!include_empty_package) {
  nb <- nb %>%
    filter(package_contents != "empty_package")
}

# ============================================================
# 3. Herd-type summary pooled across all regions
# ============================================================

supplementary_all_regions <- nb %>%
  group_by(
    herd_type,
    package_id,
    package_size,
    package_contents
  ) %>%
  summarise(
    n_contexts = n(),
    n_negative = sum(net_benefit_ref < 0, na.rm = TRUE),
    prob_negative = mean(net_benefit_ref < 0, na.rm = TRUE),
    ever_negative = any(net_benefit_ref < 0, na.rm = TRUE),
    always_positive = all(net_benefit_ref > 0, na.rm = TRUE),
    worst_net_benefit = min(net_benefit_ref, na.rm = TRUE),
    median_net_benefit = median(net_benefit_ref, na.rm = TRUE),
    best_net_benefit = max(net_benefit_ref, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(herd_type) %>%
  summarise(
    packages_evaluated = n_distinct(package_id),
    
    packages_ever_negative = sum(ever_negative),
    packages_ever_negative_pct = 100 * mean(ever_negative),
    
    packages_always_positive = sum(always_positive),
    packages_always_positive_pct = 100 * mean(always_positive),
    
    median_probability_negative_pct = 100 * median(prob_negative, na.rm = TRUE),
    mean_probability_negative_pct = 100 * mean(prob_negative, na.rm = TRUE),
    
    median_package_minimum = median(worst_net_benefit, na.rm = TRUE),
    worst_package_minimum = min(worst_net_benefit, na.rm = TRUE),
    
    median_package_maximum = median(best_net_benefit, na.rm = TRUE),
    best_package_maximum = max(best_net_benefit, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    herd_type = case_when(
      herd_type == "Naive" ~ "Naive (all regions)",
      herd_type == "Average" ~ "Average (all regions)",
      herd_type == "Infected" ~ "Infected (all regions)",
      TRUE ~ paste0(herd_type, " (all regions)")
    )
  )

# ============================================================
# 4. Herd-type by region summary
# ============================================================

supplementary_by_region <- nb %>%
  group_by(
    herd_type,
    region,
    package_id,
    package_size,
    package_contents
  ) %>%
  summarise(
    n_contexts = n(),
    n_negative = sum(net_benefit_ref < 0, na.rm = TRUE),
    prob_negative = mean(net_benefit_ref < 0, na.rm = TRUE),
    ever_negative = any(net_benefit_ref < 0, na.rm = TRUE),
    always_positive = all(net_benefit_ref > 0, na.rm = TRUE),
    worst_net_benefit = min(net_benefit_ref, na.rm = TRUE),
    median_net_benefit = median(net_benefit_ref, na.rm = TRUE),
    best_net_benefit = max(net_benefit_ref, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(herd_type, region) %>%
  summarise(
    packages_evaluated = n_distinct(package_id),
    
    packages_ever_negative = sum(ever_negative),
    packages_ever_negative_pct = 100 * mean(ever_negative),
    
    packages_always_positive = sum(always_positive),
    packages_always_positive_pct = 100 * mean(always_positive),
    
    median_probability_negative_pct = 100 * median(prob_negative, na.rm = TRUE),
    mean_probability_negative_pct = 100 * mean(prob_negative, na.rm = TRUE),
    
    median_package_minimum = median(worst_net_benefit, na.rm = TRUE),
    worst_package_minimum = min(worst_net_benefit, na.rm = TRUE),
    
    median_package_maximum = median(best_net_benefit, na.rm = TRUE),
    best_package_maximum = max(best_net_benefit, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    herd_type = paste0(herd_type, " (", region, ")")
  ) %>%
  select(-region)

# ============================================================
# 5. Combine and export
# ============================================================

supplementary_benefit_robustness <- bind_rows(
  supplementary_all_regions,
  supplementary_by_region
) %>%
  mutate(
    herd_type = factor(
      herd_type,
      levels = c(
        "Naive (all regions)",
        "Average (all regions)",
        "Infected (all regions)",
        sort(setdiff(
          herd_type,
          c(
            "Naive (all regions)",
            "Average (all regions)",
            "Infected (all regions)"
          )
        ))
      )
    )
  ) %>%
  arrange(herd_type) %>%
  mutate(
    herd_type = as.character(herd_type)
  ) %>%
  mutate(
    across(
      -herd_type,
      as.numeric
    )
  )

write.csv(
  supplementary_benefit_robustness,
  file.path(out_dir, "supplementary_benefit_robustness.csv"),
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

print(supplementary_benefit_robustness)