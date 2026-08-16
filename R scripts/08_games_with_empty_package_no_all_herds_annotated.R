# ============================================================
# Script 08: Game-based summaries of package net benefits
# ============================================================
#
# Purpose
# -------
# This script converts the package-level net benefit estimates from the
# previous step into decision-rule summaries. Each package has already been
# evaluated across combinations of region, benefit type, aggregation method,
# and disease-state/payoff-state assumptions. Here, those payoff distributions
# are summarised into three game-theoretic decision criteria:
#
#   1. Maximin:          ranks packages by their worst-case payoff.
#   2. Expected value:   ranks packages by their mean payoff.
#   3. Maximax:          ranks packages by their best-case payoff.
#
# The script also creates a simple consensus ranking that favours packages
# that perform well across all three decision criteria.
#
# Important modelling choice
# --------------------------
# The empty package is intentionally retained. This means that "doing nothing"
# remains a valid comparator and can be selected if all intervention packages
# perform poorly after costs are considered.
#
# Scope of summaries
# ------------------
# This version produces herd-/benefit-type-specific summaries only. It does
# not create an "all herds" context. Results are produced for:
#
#   - each benefit type pooled across all regions;
#   - each region × benefit type combination.
#
# Main input
# ----------
#   exports/package_net_benefits.csv
#
# Main output
# -----------
#   exports/package_game_results.csv
#
# Interpretation of output
# ------------------------
# Each row represents one package within one analysis context. The payoff
# columns describe the distribution of net benefits across the underlying
# uncertain states. The rank columns then show how the package performs under
# different decision attitudes.
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# ============================================================
# Load package-level net benefits
# ============================================================
# The input file is produced by the previous script. Each row represents the
# net benefit of a package under one specific combination of region, benefit
# type, aggregation method, and disease-state/payoff-state assumptions.

net_benefits <- read_csv(
  "exports/package_net_benefits.csv",
  show_col_types = FALSE
)

# ============================================================
# Check that all required columns are present
# ============================================================
# These checks make the script fail early if an upstream script changes the
# structure of package_net_benefits.csv. This is preferable to producing
# incomplete or misleading rankings later in the pipeline.

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
# Standardise the payoff data for the game summaries
# ============================================================
# The game calculations only need package identifiers, context variables, and
# one payoff column. The net benefit estimate is therefore renamed to `payoff`
# to make the later decision-rule code easier to read.

game_long <- net_benefits %>%
  transmute(
    package_id,
    package_size,
    package_contents,
    region,
    benefit_type,
    aggregation_method,
    joint_state,
    payoff_state,
    payoff = net_benefit_ref
  )

# Warn rather than stop if missing payoff values are present. The summary
# statistics below use na.rm = TRUE, but this warning makes it clear that the
# input data should be inspected if missing values occur.

if (any(is.na(game_long$payoff))) {
  warning("Some net_benefit_ref values are NA.")
}

# ============================================================
# Helper function: summarise payoffs within an analysis context
# ============================================================
# For each package, this function summarises the distribution of payoffs across
# all states included in the requested context. The same function is reused for
# benefit-type summaries and region × benefit-type summaries.
#
# Arguments:
#   data           Data frame containing package payoffs.
#   group_vars     Context variables used to define the analysis level.
#   analysis_level Label added to the output so downstream scripts know which
#                  context each row belongs to.
#
# Key outputs:
#   - worst_case_payoff: payoff used by the maximin rule.
#   - mean_payoff:       payoff used by the expected-value rule.
#   - best_case_payoff:  payoff used by the maximax rule.
#   - p05/p10/p90/p95:   reporting percentiles for uncertainty summaries.
#   - worst/best states: labels identifying where each package performs worst
#                        or best.

summarise_game_context <- function(data, group_vars, analysis_level) {
  
  data %>%
    group_by(
      across(all_of(group_vars)),
      package_id,
      package_size,
      package_contents
    ) %>%
    summarise(
      n_payoffs = n(),
      worst_case_payoff = min(payoff, na.rm = TRUE),
      p05_payoff = quantile(payoff, 0.05, na.rm = TRUE),
      p10_payoff = quantile(payoff, 0.10, na.rm = TRUE),
      mean_payoff = mean(payoff, na.rm = TRUE),
      median_payoff = median(payoff, na.rm = TRUE),
      p90_payoff = quantile(payoff, 0.90, na.rm = TRUE),
      p95_payoff = quantile(payoff, 0.95, na.rm = TRUE),
      best_case_payoff = max(payoff, na.rm = TRUE),
      
      # If multiple states tie for worst or best payoff, all tied states are
      # retained. This is useful for hover labels and diagnostic tables.
      worst_payoff_state = paste(
        sort(unique(payoff_state[payoff == min(payoff, na.rm = TRUE)])),
        collapse = "; "
      ),
      best_payoff_state = paste(
        sort(unique(payoff_state[payoff == max(payoff, na.rm = TRUE)])),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    mutate(
      analysis_level = analysis_level
    )
}

# ============================================================
# Summarise payoffs by benefit type
# ============================================================
# This gives one set of package rankings for each benefit type after pooling
# over regions. The context_region value is set to "All regions" to distinguish
# this pooled context from region-specific summaries.

by_benefit_type_results <- summarise_game_context(
  data = game_long,
  group_vars = c("benefit_type"),
  analysis_level = "by_benefit_type"
) %>%
  mutate(
    context_region = "All regions",
    context_benefit_type = benefit_type
  ) %>%
  select(-benefit_type)

# ============================================================
# Summarise payoffs by region and benefit type
# ============================================================
# This creates a separate decision-rule summary for each region × benefit type
# combination. These results are used when recommendations need to be tailored
# to both geography and herd/payoff context.

by_region_benefit_type_results <- summarise_game_context(
  data = game_long,
  group_vars = c("region", "benefit_type"),
  analysis_level = "by_region_benefit_type"
) %>%
  mutate(
    context_region = region,
    context_benefit_type = benefit_type
  ) %>%
  select(-region, -benefit_type)

# ============================================================
# Combine summaries and calculate decision-rule rankings
# ============================================================
# Within each analysis context, packages are ranked according to three decision
# rules:
#
#   - rank_maximin ranks the highest worst-case payoff first.
#   - rank_expected_value ranks the highest mean payoff first.
#   - rank_maximax ranks the highest best-case payoff first.
#
# The consensus rank is a pragmatic summary ranking. It primarily orders
# packages by their average rank across the three rules, with small tie-breakers
# that favour packages with a better worst rank and more top-10 appearances.

package_game_results <- bind_rows(
  by_benefit_type_results,
  by_region_benefit_type_results
) %>%
  group_by(
    analysis_level,
    context_region,
    context_benefit_type
  ) %>%
  mutate(
    rank_maximin = min_rank(desc(worst_case_payoff)),
    rank_expected_value = min_rank(desc(mean_payoff)),
    rank_maximax = min_rank(desc(best_case_payoff)),
    
    # Average performance across the three decision attitudes.
    mean_rank = rowMeans(
      cbind(
        rank_maximin,
        rank_expected_value,
        rank_maximax
      ),
      na.rm = TRUE
    ),
    
    # Worst rank achieved by the package across the three decision rules.
    worst_rank = pmax(
      rank_maximin,
      rank_expected_value,
      rank_maximax,
      na.rm = TRUE
    ),
    
    # Count how often the package appears among the top 10 packages across the
    # three decision rules. This helps identify packages that are consistently
    # near the top, even if they are not always ranked first.
    top10_count =
      as.integer(rank_maximin <= 10) +
      as.integer(rank_expected_value <= 10) +
      as.integer(rank_maximax <= 10),
    
    # Consensus ranking. The 0.001 weights are intentionally small so that
    # mean_rank remains the dominant criterion, while worst_rank and top10_count
    # only act as tie-breakers or near-tie adjustments.
    rank_consensus = min_rank(
      mean_rank +
        0.001 * worst_rank -
        0.001 * top10_count
    ),
    
    # Logical flags identify the winning package(s) under each rule. Ties are
    # possible because min_rank assigns the same rank to equal payoff values.
    winner_maximin = rank_maximin == 1,
    winner_expected_value = rank_expected_value == 1,
    winner_maximax = rank_maximax == 1,
    winner_consensus = rank_consensus == 1
  ) %>%
  ungroup() %>%
  select(
    analysis_level,
    context_region,
    context_benefit_type,
    package_id,
    package_size,
    package_contents,
    n_payoffs,
    worst_case_payoff,
    p05_payoff,
    p10_payoff,
    mean_payoff,
    median_payoff,
    p90_payoff,
    p95_payoff,
    best_case_payoff,
    worst_payoff_state,
    best_payoff_state,
    rank_maximin,
    rank_expected_value,
    rank_maximax,
    rank_consensus,
    mean_rank,
    worst_rank,
    top10_count,
    winner_maximin,
    winner_expected_value,
    winner_maximax,
    winner_consensus
  ) %>%
  arrange(
    analysis_level,
    context_region,
    context_benefit_type,
    rank_consensus,
    package_size,
    package_id
  )

# ============================================================
# Export game summaries
# ============================================================
# The exported table is used by downstream scripts that construct Pareto
# frontiers, policy views, hover labels, and manuscript-ready summaries.

out_dir <- "exports"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  package_game_results,
  file.path(out_dir, "package_game_results.csv")
)
