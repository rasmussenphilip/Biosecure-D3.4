library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(plotly)
library(htmlwidgets)

# ============================================================
# 09_pareto_frontier_packages_with_empty_no_all_herds.R
# ============================================================
#
# Purpose
# -------
# This script identifies Pareto-optimal biosecurity packages using
# package-level net-benefit summaries generated in the previous script.
# It keeps the empty package as an explicit comparator and restricts the
# analysis to herd-type-specific contexts:
#
#   1. naive_herd
#   2. infected_herd
#   3. average_herd
#
# For each relevant analysis context, packages are compared across three
# payoff dimensions:
#
#   - worst_case_payoff  : conservative / maximin outcome
#   - mean_payoff        : expected-value outcome
#   - best_case_payoff   : optimistic / maximax outcome
#
# A package is considered Pareto-optimal if no other package performs at
# least as well on all retained payoff dimensions and strictly better on
# at least one of them. These Pareto-optimal packages are exported as the
# main frontier result.
#
# The script also creates interactive 3D Plotly figures in which:
#
#   - dominated packages are shown as background points;
#   - Pareto-optimal packages are highlighted;
#   - hover labels report net-benefit summaries, decision-rule ranks,
#     and disease-specific risk reductions.
#
# Main inputs
# -----------
# exports/package_game_results.csv
#   Output from the game/decision-rule script. Contains package-level
#   net-benefit summaries and ranks by context.
#
# exports/package_risk_reductions.csv
#   Output from the risk-reduction script. Contains post-package risks
#   for naive and infected herd states by disease and risk scenario.
#
# Main outputs
# ------------
# exports/package_pareto_frontier.csv
#   Table containing only Pareto-optimal packages in each context.
#
# exports/pareto_3d_plots/*.html
#   Interactive 3D Pareto plots for each analysis context.
#
# Notes
# -----
# - The frontier itself is based only on the three net-benefit dimensions:
#   worst-case, mean, and best-case payoff.
# - Percentiles and risk reductions are reported for interpretation only;
#   they do not determine Pareto-frontier membership.
# - The empty package is retained so that all packages can be interpreted
#   relative to a no-additional-measure comparator.
# - This script does not include an "all herds" context; it uses only the
#   herd-type-specific benefit contexts listed above.
# ============================================================

# Load package-level decision-rule results and post-package risk estimates.
# The first file supplies the economic payoff dimensions used for the frontier.
# The second file is used only to add risk-reduction diagnostics to the output
# table and plot hover labels.
package_game_results <- read_csv("exports/package_game_results.csv", show_col_types = FALSE)
package_risk_reductions <- read_csv("exports/package_risk_reductions.csv", show_col_types = FALSE)

# Define output folders. The CSV result is written directly to exports/,
# while the interactive HTML plots are stored in a dedicated subfolder.
out_dir <- "exports"
plot_dir <- file.path(out_dir, "pareto_3d_plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================

# Identify non-dominated packages within one analysis context.
#
# The input data frame should contain one row per package for a single
# context, with payoff columns named:
#   - worst_case_payoff
#   - mean_payoff
#   - best_case_payoff
#
# The function sorts packages from high to low worst-case payoff and then
# checks whether previously encountered packages already dominate the current
# one in the mean and best-case dimensions. Because the sorting handles the
# worst-case dimension, the loop only needs to track the best mean and best
# best-case payoffs seen so far.
get_frontier_fast <- function(df) {
  
  df_sorted <- df %>%
    arrange(
      desc(worst_case_payoff),
      desc(mean_payoff),
      desc(best_case_payoff)
    )
  
  best_mean_so_far <- -Inf
  best_best_so_far <- -Inf
  on_frontier <- logical(nrow(df_sorted))
  
  for (i in seq_len(nrow(df_sorted))) {
    
    dominated <- best_mean_so_far >= df_sorted$mean_payoff[i] &&
      best_best_so_far >= df_sorted$best_case_payoff[i] &&
      (
        best_mean_so_far > df_sorted$mean_payoff[i] ||
          best_best_so_far > df_sorted$best_case_payoff[i]
      )
    
    on_frontier[i] <- !dominated
    
    best_mean_so_far <- max(best_mean_so_far, df_sorted$mean_payoff[i])
    best_best_so_far <- max(best_best_so_far, df_sorted$best_case_payoff[i])
  }
  
  df_sorted %>%
    mutate(on_frontier = on_frontier)
}

# Convert context labels into safe file-name components for the HTML plots.
safe_file_label <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_") %>%
    str_remove("_$")
}

# ============================================================
# Risk-reduction diagnostics
# ============================================================

# Reshape the wide risk-reduction output into a long format so that the
# encoded column names can be parsed into:
#   - herd state / risk type: naive or infected
#   - disease: BVD, PTB, or Salmonella
#   - suffix: the remaining risk-scenario identifier
risk_meta_cols <- c("package_id", "package_size", "package_contents")
risk_value_cols <- setdiff(names(package_risk_reductions), risk_meta_cols)

risk_long <- package_risk_reductions %>%
  pivot_longer(
    cols = all_of(risk_value_cols),
    names_to = "risk_col",
    values_to = "risk_value"
  ) %>%
  mutate(risk_col = str_to_lower(risk_col))

# Parse column names of the form:
#   post_naive_risk_bvd_<suffix>
#   post_infected_risk_ptb_<suffix>
#
# The suffix is retained because it identifies matched baseline and
# post-package risks from the same underlying risk scenario.
parsed_cols <- str_match(
  risk_long$risk_col,
  "^post_(naive|infected)_risk_(bvd|ptb|sal)_(.*)$"
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

# Return to a wide state-by-disease layout, with one row per package and
# suffix. This format makes it straightforward to compare each package to
# the empty-package baseline within the same risk scenario.
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

# Extract the empty-package risks. These represent the no-additional-measure
# comparator and are used as the baseline for calculating absolute risk
# reductions.
empty_risk_states <- risk_states %>%
  filter(package_contents == "empty_package") %>%
  select(
    suffix,
    baseline_naive_bvd = naive_bvd,
    baseline_naive_ptb = naive_ptb,
    baseline_naive_sal = naive_sal,
    baseline_infected_bvd = infected_bvd,
    baseline_infected_ptb = infected_ptb,
    baseline_infected_sal = infected_sal
  )

# Calculate risk reductions for each package relative to the empty package.
#
# Naive and infected reductions are simple baseline-minus-post-risk
# differences within the same disease and risk scenario.
#
# The "average" reductions combine naive and infected components by comparing
# the product of the two risk terms before and after the package is applied.
# These summaries are diagnostic/reporting outputs only; they do not affect
# the Pareto-frontier calculation.
risk_reduction_summary <- risk_states %>%
  left_join(empty_risk_states, by = "suffix") %>%
  mutate(
    naive_bvd_reduction = baseline_naive_bvd - naive_bvd,
    naive_ptb_reduction = baseline_naive_ptb - naive_ptb,
    naive_sal_reduction = baseline_naive_sal - naive_sal,
    
    infected_bvd_reduction = baseline_infected_bvd - infected_bvd,
    infected_ptb_reduction = baseline_infected_ptb - infected_ptb,
    infected_sal_reduction = baseline_infected_sal - infected_sal,
    
    average_bvd_reduction =
      baseline_naive_bvd * baseline_infected_bvd -
      naive_bvd * infected_bvd,
    
    average_ptb_reduction =
      baseline_naive_ptb * baseline_infected_ptb -
      naive_ptb * infected_ptb,
    
    average_sal_reduction =
      baseline_naive_sal * baseline_infected_sal -
      naive_sal * infected_sal
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
        worst = ~min(.x, na.rm = TRUE),
        p05 = ~quantile(.x, 0.05, na.rm = TRUE),
        p10 = ~quantile(.x, 0.10, na.rm = TRUE),
        mean = ~mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        p90 = ~quantile(.x, 0.90, na.rm = TRUE),
        p95 = ~quantile(.x, 0.95, na.rm = TRUE),
        best = ~max(.x, na.rm = TRUE)
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )

# ============================================================
# Plot function
# ============================================================

# Create an interactive 3D Pareto plot for one context.
#
# Axes:
#   x = worst-case net benefit
#   y = mean net benefit
#   z = best-case net benefit
#
# Points are split into dominated and Pareto-optimal packages so that the
# frontier packages can be highlighted visually.
make_3d_plot <- function(df, title, file_name) {
  
  bg_df <- df %>% filter(!on_frontier)
  fr_df <- df %>% filter(on_frontier)
  
  # Construct detailed hover labels. These labels are intentionally verbose
  # because the interactive plots are used as exploratory diagnostics: they
  # allow each package's economic performance, ranks, and risk reductions to
  # be inspected without opening the CSV output.
  hover_text <- function(data) {
    paste(
      "Package ID:", data$package_id,
      "<br>Package size:", data$package_size,
      "<br>", data$package_contents,
      
      "<br><br><b>Net benefit summary</b>",
      "<br>Worst-case:", round(data$worst_case_payoff, 2),
      "<br>Mean:", round(data$mean_payoff, 2),
      "<br>Best-case:", round(data$best_case_payoff, 2),
      "<br>5th percentile:", round(data$p05_payoff, 2),
      "<br>10th percentile:", round(data$p10_payoff, 2),
      "<br>Median:", round(data$median_payoff, 2),
      "<br>90th percentile:", round(data$p90_payoff, 2),
      "<br>95th percentile:", round(data$p95_payoff, 2),
      
      "<br><br><b>Ranks</b>",
      "<br>Rank maximin:", data$rank_maximin,
      "<br>Rank expected value:", data$rank_expected_value,
      "<br>Rank maximax:", data$rank_maximax,
      "<br>Number of payoffs:", data$n_payoffs,
      
      "<br><br><b>Risk reductions shown as (BVD, PTB, SAL)</b>",
      
      "<br><br><b>Naive risk reduction</b>",
      "<br>Worst: (",
      round(data$worst_naive_bvd_reduction, 4), ", ",
      round(data$worst_naive_ptb_reduction, 4), ", ",
      round(data$worst_naive_sal_reduction, 4), ")",
      "<br>5th percentile: (",
      round(data$p05_naive_bvd_reduction, 4), ", ",
      round(data$p05_naive_ptb_reduction, 4), ", ",
      round(data$p05_naive_sal_reduction, 4), ")",
      "<br>10th percentile: (",
      round(data$p10_naive_bvd_reduction, 4), ", ",
      round(data$p10_naive_ptb_reduction, 4), ", ",
      round(data$p10_naive_sal_reduction, 4), ")",
      "<br>Mean: (",
      round(data$mean_naive_bvd_reduction, 4), ", ",
      round(data$mean_naive_ptb_reduction, 4), ", ",
      round(data$mean_naive_sal_reduction, 4), ")",
      "<br>Median: (",
      round(data$median_naive_bvd_reduction, 4), ", ",
      round(data$median_naive_ptb_reduction, 4), ", ",
      round(data$median_naive_sal_reduction, 4), ")",
      "<br>90th percentile: (",
      round(data$p90_naive_bvd_reduction, 4), ", ",
      round(data$p90_naive_ptb_reduction, 4), ", ",
      round(data$p90_naive_sal_reduction, 4), ")",
      "<br>95th percentile: (",
      round(data$p95_naive_bvd_reduction, 4), ", ",
      round(data$p95_naive_ptb_reduction, 4), ", ",
      round(data$p95_naive_sal_reduction, 4), ")",
      "<br>Best: (",
      round(data$best_naive_bvd_reduction, 4), ", ",
      round(data$best_naive_ptb_reduction, 4), ", ",
      round(data$best_naive_sal_reduction, 4), ")",
      
      "<br><br><b>Infected risk reduction</b>",
      "<br>Worst: (",
      round(data$worst_infected_bvd_reduction, 4), ", ",
      round(data$worst_infected_ptb_reduction, 4), ", ",
      round(data$worst_infected_sal_reduction, 4), ")",
      "<br>5th percentile: (",
      round(data$p05_infected_bvd_reduction, 4), ", ",
      round(data$p05_infected_ptb_reduction, 4), ", ",
      round(data$p05_infected_sal_reduction, 4), ")",
      "<br>10th percentile: (",
      round(data$p10_infected_bvd_reduction, 4), ", ",
      round(data$p10_infected_ptb_reduction, 4), ", ",
      round(data$p10_infected_sal_reduction, 4), ")",
      "<br>Mean: (",
      round(data$mean_infected_bvd_reduction, 4), ", ",
      round(data$mean_infected_ptb_reduction, 4), ", ",
      round(data$mean_infected_sal_reduction, 4), ")",
      "<br>Median: (",
      round(data$median_infected_bvd_reduction, 4), ", ",
      round(data$median_infected_ptb_reduction, 4), ", ",
      round(data$median_infected_sal_reduction, 4), ")",
      "<br>90th percentile: (",
      round(data$p90_infected_bvd_reduction, 4), ", ",
      round(data$p90_infected_ptb_reduction, 4), ", ",
      round(data$p90_infected_sal_reduction, 4), ")",
      "<br>95th percentile: (",
      round(data$p95_infected_bvd_reduction, 4), ", ",
      round(data$p95_infected_ptb_reduction, 4), ", ",
      round(data$p95_infected_sal_reduction, 4), ")",
      "<br>Best: (",
      round(data$best_infected_bvd_reduction, 4), ", ",
      round(data$best_infected_ptb_reduction, 4), ", ",
      round(data$best_infected_sal_reduction, 4), ")",
      
      "<br><br><b>Average risk reduction</b>",
      "<br>Worst: (",
      round(data$worst_average_bvd_reduction, 4), ", ",
      round(data$worst_average_ptb_reduction, 4), ", ",
      round(data$worst_average_sal_reduction, 4), ")",
      "<br>5th percentile: (",
      round(data$p05_average_bvd_reduction, 4), ", ",
      round(data$p05_average_ptb_reduction, 4), ", ",
      round(data$p05_average_sal_reduction, 4), ")",
      "<br>10th percentile: (",
      round(data$p10_average_bvd_reduction, 4), ", ",
      round(data$p10_average_ptb_reduction, 4), ", ",
      round(data$p10_average_sal_reduction, 4), ")",
      "<br>Mean: (",
      round(data$mean_average_bvd_reduction, 4), ", ",
      round(data$mean_average_ptb_reduction, 4), ", ",
      round(data$mean_average_sal_reduction, 4), ")",
      "<br>Median: (",
      round(data$median_average_bvd_reduction, 4), ", ",
      round(data$median_average_ptb_reduction, 4), ", ",
      round(data$median_average_sal_reduction, 4), ")",
      "<br>90th percentile: (",
      round(data$p90_average_bvd_reduction, 4), ", ",
      round(data$p90_average_ptb_reduction, 4), ", ",
      round(data$p90_average_sal_reduction, 4), ")",
      "<br>95th percentile: (",
      round(data$p95_average_bvd_reduction, 4), ", ",
      round(data$p95_average_ptb_reduction, 4), ", ",
      round(data$p95_average_sal_reduction, 4), ")",
      "<br>Best: (",
      round(data$best_average_bvd_reduction, 4), ", ",
      round(data$best_average_ptb_reduction, 4), ", ",
      round(data$best_average_sal_reduction, 4), ")"
    )
  }
  
  p <- plot_ly() %>%
    add_trace(
      data = bg_df,
      x = ~worst_case_payoff,
      y = ~mean_payoff,
      z = ~best_case_payoff,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, color = "steelblue", opacity = 0.25),
      text = hover_text(bg_df),
      hoverinfo = "text",
      name = "Dominated packages"
    ) %>%
    add_trace(
      data = fr_df,
      x = ~worst_case_payoff,
      y = ~mean_payoff,
      z = ~best_case_payoff,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 6, color = "red", opacity = 1),
      text = paste("PARETO FRONTIER<br>", hover_text(fr_df)),
      hoverinfo = "text",
      name = "Pareto-optimal packages"
    ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "Worst-case net benefit"),
        yaxis = list(title = "Mean net benefit"),
        zaxis = list(title = "Best-case net benefit")
      )
    )
  
  htmlwidgets::saveWidget(
    as_widget(p),
    file.path(plot_dir, file_name),
    selfcontained = TRUE
  )
}

# ============================================================
# Pareto input
# Keep only herd-type-specific contexts
# Empty package retained
# ============================================================

# Build the data set used for Pareto-frontier identification.
#
# Only herd-type-specific benefit contexts are retained:
#   - naive_herd
#   - infected_herd
#   - average_herd
#
# The empty package is not removed. Keeping it in the candidate set makes the
# no-additional-measure comparator visible in both the frontier table and the
# 3D plots.
pareto_input <- package_game_results %>%
  filter(
    analysis_level %in% c(
      "by_benefit_type",
      "by_region_benefit_type"
    )
  ) %>%
  filter(
    context_benefit_type %in% c(
      "naive_herd",
      "infected_herd",
      "average_herd"
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
# Identify Pareto frontier packages
# ============================================================

# Apply the Pareto-frontier algorithm separately within each context.
# Contexts are defined by:
#   - analysis level
#   - region
#   - benefit type / herd state
#
# The script then records how many packages were evaluated, how many are on
# the frontier, and what share of all packages the frontier represents.
pareto_all_packages <- pareto_input %>%
  group_by(
    analysis_level,
    context_region,
    context_benefit_type
  ) %>%
  group_modify(~ get_frontier_fast(.x)) %>%
  ungroup() %>%
  group_by(
    analysis_level,
    context_region,
    context_benefit_type
  ) %>%
  mutate(
    n_packages_in_context = n(),
    n_frontier_packages = sum(on_frontier),
    frontier_share = n_frontier_packages / n_packages_in_context,
    frontier_rank_maximin = ifelse(on_frontier, rank_maximin, NA_integer_),
    frontier_rank_expected_value = ifelse(on_frontier, rank_expected_value, NA_integer_),
    frontier_rank_maximax = ifelse(on_frontier, rank_maximax, NA_integer_)
  ) %>%
  ungroup()

# ============================================================
# Export Pareto-optimal packages only
# ============================================================

# Keep only Pareto-optimal packages for the exported frontier table.
#
# Percentile payoffs and risk-reduction summaries are retained here to support
# interpretation and manuscript/table preparation, but they were not used to
# decide whether a package is Pareto-optimal.
package_pareto_frontier <- pareto_all_packages %>%
  filter(on_frontier) %>%
  select(
    analysis_level,
    context_region,
    context_benefit_type,
    package_id,
    package_size,
    package_contents,
    n_payoffs,
    n_packages_in_context,
    n_frontier_packages,
    frontier_share,
    
    worst_case_payoff,
    p05_payoff,
    p10_payoff,
    mean_payoff,
    median_payoff,
    p90_payoff,
    p95_payoff,
    best_case_payoff,
    
    frontier_rank_maximin,
    frontier_rank_expected_value,
    frontier_rank_maximax,
    
    contains("naive_bvd_reduction"),
    contains("naive_ptb_reduction"),
    contains("naive_sal_reduction"),
    contains("infected_bvd_reduction"),
    contains("infected_ptb_reduction"),
    contains("infected_sal_reduction"),
    contains("average_bvd_reduction"),
    contains("average_ptb_reduction"),
    contains("average_sal_reduction")
  ) %>%
  arrange(
    analysis_level,
    context_region,
    context_benefit_type,
    frontier_rank_expected_value,
    package_size,
    package_id
  )

write_csv(
  package_pareto_frontier,
  file.path(out_dir, "package_pareto_frontier.csv")
)

# ============================================================
# Create one 3D Pareto plot per context
# ============================================================

# Create one interactive HTML plot for each analysis context.
# File names are assembled from the context labels after converting them to
# safe file-name strings.
plot_contexts <- pareto_all_packages %>%
  distinct(
    analysis_level,
    context_region,
    context_benefit_type
  )

for (i in seq_len(nrow(plot_contexts))) {
  
  ctx <- plot_contexts[i, ]
  
  df_plot <- pareto_all_packages %>%
    filter(
      analysis_level == ctx$analysis_level,
      context_region == ctx$context_region,
      context_benefit_type == ctx$context_benefit_type
    )
  
  plot_title <- paste(
    "Pareto frontier:",
    ctx$analysis_level,
    "|",
    ctx$context_region,
    "|",
    ctx$context_benefit_type
  )
  
  plot_file <- paste0(
    "pareto_3d_",
    safe_file_label(ctx$analysis_level),
    "_",
    safe_file_label(ctx$context_region),
    "_",
    safe_file_label(ctx$context_benefit_type),
    ".html"
  )
  
  make_3d_plot(
    df = df_plot,
    title = plot_title,
    file_name = plot_file
  )
}