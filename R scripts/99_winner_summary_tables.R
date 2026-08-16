# ============================================================
# 9. Export game-winner tables
# ============================================================

package_game_results <- read_csv(
  file.path("exports", "package_game_results.csv"),
  show_col_types = FALSE
)

game_results <- package_game_results %>%
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

game_winners_long <- game_results %>%
  pivot_longer(
    cols = c(
      winner_maximin,
      winner_expected_value,
      winner_maximax
    ),
    names_to = "decision_rule",
    values_to = "winner"
  ) %>%
  filter(winner) %>%
  mutate(
    `Decision rule` = recode(
      decision_rule,
      winner_maximin = "Maximin",
      winner_expected_value = "Expected value",
      winner_maximax = "Maximax"
    )
  )

make_game_winner_table <- function(df, herd_type, herd_label, risk_prefix) {
  
  df %>%
    filter(context_benefit_type == herd_type) %>%
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
      `Decision rule`,
      `Package ID` = package_id,
      Size = package_size,
      `Package contents` = package_contents,
      
      `BVD; Median (IQR)` = format_median_iqr(
        .data[[paste0("median_", risk_prefix, "_bvd_reduction")]],
        .data[[paste0("iqr_low_", risk_prefix, "_bvd_reduction")]],
        .data[[paste0("iqr_high_", risk_prefix, "_bvd_reduction")]]
      ),
      
      `PTB; Median (IQR)` = format_median_iqr(
        .data[[paste0("median_", risk_prefix, "_ptb_reduction")]],
        .data[[paste0("iqr_low_", risk_prefix, "_ptb_reduction")]],
        .data[[paste0("iqr_high_", risk_prefix, "_ptb_reduction")]]
      ),
      
      `SAL; Median (IQR)` = format_median_iqr(
        .data[[paste0("median_", risk_prefix, "_sal_reduction")]],
        .data[[paste0("iqr_low_", risk_prefix, "_sal_reduction")]],
        .data[[paste0("iqr_high_", risk_prefix, "_sal_reduction")]]
      ),
      
      `Net benefit (EUR/cow-year); Median (IQR)` = format_median_iqr(
        median_payoff,
        payoff_iqr_low,
        payoff_iqr_high,
        digits = 1
      )
    )
}

game_winner_tables <- bind_rows(
  make_game_winner_table(
    game_winners_long,
    herd_type = "naive_herd",
    herd_label = "Naive herds",
    risk_prefix = "naive"
  ),
  make_game_winner_table(
    game_winners_long,
    herd_type = "infected_herd",
    herd_label = "Infected herds",
    risk_prefix = "infected"
  ),
  make_game_winner_table(
    game_winners_long,
    herd_type = "average_herd",
    herd_label = "Average herds",
    risk_prefix = "average"
  )
) %>%
  mutate(
    `Herd type` = factor(
      `Herd type`,
      levels = c(
        "Naive herds",
        "Infected herds",
        "Average herds"
      )
    ),
    `Decision rule` = factor(
      `Decision rule`,
      levels = c(
        "Maximin",
        "Expected value",
        "Maximax"
      )
    )
  ) %>%
  arrange(
    `Herd type`,
    Region,
    `Decision rule`,
    `Package ID`
  )

# Main game-winner table: All regions only
game_winners_all_regions <- game_winner_tables %>%
  filter(Region == "All regions") %>%
  select(-Region)

write_csv(
  game_winners_all_regions,
  file.path(out_dir, "game_winners_all_regions.csv")
)

# Supplementary game-winner table: all regions/contexts
supplementary_game_winners <- game_winner_tables

write_csv(
  supplementary_game_winners,
  file.path(out_dir, "supplementary_game_winners.csv")
)