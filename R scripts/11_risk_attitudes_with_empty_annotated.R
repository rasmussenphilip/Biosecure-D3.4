# ============================================================
# Script 11: Risk-attitude analysis of winning biosecurity packages
# ============================================================
#
# Purpose
# -------
# This script compares which biosecurity packages are selected under
# alternative farmer or policy-maker risk attitudes. It uses the package-level
# game results produced upstream, where each context has three decision-rule
# indicators: maximin, expected value, and maximax.
#
# The script asks three related questions:
#   1. Do the same packages win under cautious, average-outcome, and
#      optimistic decision rules?
#   2. Do the same individual biosecurity measures appear inside those winning
#      packages?
#   3. Which measures are most characteristic of each risk attitude?
#
# Inputs
# ------
#   exports/package_game_results.csv
#     Package-level decision-rule output from the previous script. This file
#     identifies which package wins within each analysis context.
#
#   inputs/measures_final.csv
#     Lookup table linking measure identifiers to human-readable measure names.
#
# Outputs
# -------
#   exports/risk_attitude/risk_attitude_venn_panel.png
#     Two-panel Venn figure showing overlap in winning packages and overlap in
#     measures contained within winning packages.
#
#   exports/risk_attitude/risk_attitude_signature_measures.csv
#     Supplementary table reporting the percentage of contexts in which each
#     measure appears in a winning package under each decision rule.
#
#   exports/risk_attitude/risk_attitude_measure_heatmap.png
#     Heatmap of measure-selection frequency by risk attitude. The empty
#     package is retained in the CSV but removed from the heatmap so the figure
#     focuses on active biosecurity measures.
#
# Interpretation
# --------------
# Maximin represents a cautious decision rule because it selects the package
# with the best worst-case payoff. Expected value selects the package with the
# best average payoff across modelled contexts. Maximax represents an optimistic
# decision rule because it selects the package with the best best-case payoff.
# Comparing these three rules shows whether recommendations are robust to
# different attitudes toward uncertainty.
#
# The code below is descriptive: it does not re-estimate package payoffs or
# change the Pareto/game results. It only reshapes, summarizes, visualizes, and
# exports the already-computed winning-package results.
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggVennDiagram)
library(paletteer)
library(patchwork)

out_dir <- "exports/risk_attitude"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Set this to FALSE if the all-region pooled context should be excluded from
# the risk-attitude summaries. Keeping it TRUE includes both region-specific
# and pooled contexts in the descriptive comparison.
include_all_regions <- TRUE

# ============================================================
# Load input data
# ============================================================
# The package-game file contains the winning-package flags for each
# decision rule and context. The measure lookup file is used later to
# replace compact measure identifiers with readable labels.

package_game_results <- read_csv(
  "exports/package_game_results.csv",
  show_col_types = FALSE
)

measure_lookup <- read_csv(
  "inputs/measures_final.csv",
  show_col_types = FALSE
) %>%
  distinct(
    measure_id = measure,
    name
  )

# Optional sensitivity setting: remove the pooled "All regions" context if the
# analysis should focus only on individual regional contexts.
if (!include_all_regions) {
  package_game_results <- package_game_results %>%
    filter(context_region != "All regions")
}

# ============================================================
# Prepare winning packages in long format
# ============================================================
# The upstream file stores the three decision-rule winners in separate logical
# columns. For plotting and counting, convert these columns to one row per
# winning package per decision rule.

winner_long <- package_game_results %>%
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
    decision_rule = recode(
      decision_rule,
      winner_maximin = "Maximin",
      winner_expected_value = "Expected\nvalue",
      winner_maximax = "Maximax"
    ),
    decision_rule = factor(
      decision_rule,
      levels = c("Maximin", "Expected\nvalue", "Maximax")
    ),
    package_label = paste0("P", package_id)
  )

# ============================================================
# Venn diagram: winning packages
# ============================================================
# For each decision rule, collect the unique package labels that ever win.
# The Venn diagram then shows how much package choice overlaps across
# risk attitudes.

package_sets <- winner_long %>%
  distinct(decision_rule, package_label) %>%
  group_by(decision_rule) %>%
  summarise(
    items = list(package_label),
    .groups = "drop"
  ) %>%
  deframe()

p_package_venn <- ggVennDiagram(
  package_sets,
  label = "count",
  label_alpha = 0,
  set_size = 0
) +
  annotate(
    "text",
    x = -3.2,
    y = 3.8,
    label = "Maximin",
    size = 5
  ) +
  annotate(
    "text",
    x = 7.3,
    y = 3.8,
    label = "Expected\nvalue",
    size = 5
  ) +
  annotate(
    "text",
    x = 2,
    y = -8.5,
    label = "Maximax",
    size = 5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#C73D7B"
  ) +
  labs(
    title = "Overlap in winning packages across risk attitudes"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# ============================================================
# Extract individual measures from winning packages
# ============================================================
# Package contents are stored as a combined string. Split each winning package
# into its component measures so the analysis can compare individual measure
# frequencies across decision rules. The empty package is intentionally kept
# here because choosing no active measure is also a possible decision-rule
# outcome.

winner_measures <- winner_long %>%
  select(
    analysis_level,
    context_region,
    context_benefit_type,
    decision_rule,
    package_id,
    package_contents
  ) %>%
  separate_rows(
    package_contents,
    sep = "\\s*\\+\\s*|\\s*;\\s*"
  ) %>%
  mutate(
    measure_id = str_squish(package_contents),
    measure_id = if_else(
      str_detect(
        measure_id,
        regex("^empty_package$|^empty package$", ignore_case = TRUE)
      ),
      "empty_package",
      measure_id
    )
  ) %>%
  filter(
    !is.na(measure_id),
    measure_id != ""
  ) %>%
  select(-package_contents)

# ============================================================
# Venn diagram: measures in winning packages
# ============================================================
# This repeats the Venn analysis at the measure level. A large overlap means
# that the same measures tend to appear in winners regardless of risk attitude,
# even if the exact package combinations differ.

measure_sets <- winner_measures %>%
  distinct(decision_rule, measure_id) %>%
  group_by(decision_rule) %>%
  summarise(
    items = list(measure_id),
    .groups = "drop"
  ) %>%
  deframe()

p_measure_venn <- ggVennDiagram(
  measure_sets,
  label = "count",
  label_alpha = 0,
  set_size = 0
) +
  annotate(
    "text",
    x = -3.2,
    y = 3.8,
    label = "Maximin",
    size = 5
  ) +
  annotate(
    "text",
    x = 7.3,
    y = 3.8,
    label = "Expected\nvalue",
    size = 5
  ) +
  annotate(
    "text",
    x = 2,
    y = -8.5,
    label = "Maximax",
    size = 5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#C73D7B"
  ) +
  labs(
    title = "Overlap in measures selected by risk attitude"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# ============================================================
# Combine package- and measure-level Venn diagrams
# ============================================================
# Panel A summarizes overlap in whole packages. Panel B summarizes overlap in
# the measures contained within those packages.

p_package_venn_panel <- p_package_venn +
  labs(
    title = "Winning packages",
    subtitle = "A."
  ) +
  coord_fixed(clip = "off") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0,
      face = "bold",
      size = 16,
      margin = margin(t = -15, b = 0)
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

p_measure_venn_panel <- p_measure_venn +
  labs(
    title = "Measures within winners",
    subtitle = "B."
  ) +
  coord_fixed(clip = "off") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0,
      face = "bold",
      size = 16,
      margin = margin(t = -15, b = 0)
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

p_venn_panel <- p_package_venn_panel +
  p_measure_venn_panel +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  ) +
  plot_annotation(
    title = "Overlap across strategies",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
      plot.margin = margin(2, 2, 2, 2)
    )
  )

ggsave(
  file.path(out_dir, "risk_attitude_venn_panel.png"),
  p_venn_panel,
  width = 10,
  height = 6,
  dpi = 600
)

# ============================================================
# Measure frequencies by risk attitude
# ============================================================
# Build a complete measure-by-decision-rule table. Measures that never appear
# in a winning package are retained with zero frequency, which makes the CSV
# complete and prevents the heatmap from silently dropping available measures.

all_measures <- package_game_results %>%
  distinct(package_contents) %>%
  separate_rows(
    package_contents,
    sep = "\\s*\\+\\s*|\\s*;\\s*"
  ) %>%
  transmute(
    measure_id = str_squish(package_contents),
    measure_id = if_else(
      str_detect(
        measure_id,
        regex("^empty_package$|^empty package$", ignore_case = TRUE)
      ),
      "empty_package",
      measure_id
    )
  ) %>%
  filter(
    !is.na(measure_id),
    measure_id != ""
  ) %>%
  distinct(measure_id) %>%
  pull(measure_id)

measure_lookup <- bind_rows(
  measure_lookup,
  tibble(
    measure_id = "empty_package",
    name = "Empty package"
  )
) %>%
  distinct(measure_id, .keep_all = TRUE)

# Count the number of analysis contexts represented for each decision rule.
# These denominators are used to convert raw counts into percentages.
n_contexts_by_rule <- winner_long %>%
  distinct(
    decision_rule,
    analysis_level,
    context_region,
    context_benefit_type
  ) %>%
  count(
    decision_rule,
    name = "n_contexts"
  )

# Count the number and percentage of contexts in which each measure appears
# in at least one winning package for each decision rule.
measure_rule_summary <- winner_measures %>%
  distinct(
    decision_rule,
    analysis_level,
    context_region,
    context_benefit_type,
    measure_id
  ) %>%
  count(
    decision_rule,
    measure_id,
    name = "n_contexts_with_measure"
  ) %>%
  complete(
    decision_rule,
    measure_id = all_measures,
    fill = list(n_contexts_with_measure = 0)
  ) %>%
  left_join(
    n_contexts_by_rule,
    by = "decision_rule"
  ) %>%
  mutate(
    pct_contexts = 100 * n_contexts_with_measure / n_contexts
  ) %>%
  left_join(
    measure_lookup,
    by = "measure_id"
  ) %>%
  mutate(
    name = if_else(is.na(name), measure_id, name)
  )

# Derive a simple risk-attitude signature for each measure. The signature is
# the decision rule under which the measure has its highest selection
# percentage. Preference strength is the spread between the highest and lowest
# percentages, so larger values indicate more rule-specific selection.
signature_measures <- measure_rule_summary %>%
  select(
    measure_id,
    name,
    decision_rule,
    pct_contexts
  ) %>%
  pivot_wider(
    names_from = decision_rule,
    values_from = pct_contexts,
    values_fill = 0
  ) %>%
  mutate(
    maximin_pct = Maximin,
    expected_value_pct = `Expected\nvalue`,
    maximax_pct = Maximax,
    preference_strength = pmax(
      maximin_pct,
      expected_value_pct,
      maximax_pct
    ) - pmin(
      maximin_pct,
      expected_value_pct,
      maximax_pct
    ),
    risk_attitude_signature = case_when(
      maximin_pct == pmax(maximin_pct, expected_value_pct, maximax_pct) ~ "Maximin",
      expected_value_pct == pmax(maximin_pct, expected_value_pct, maximax_pct) ~ "Expected value",
      maximax_pct == pmax(maximin_pct, expected_value_pct, maximax_pct) ~ "Maximax"
    )
  ) %>%
  select(
    measure_id,
    name,
    maximin_pct,
    expected_value_pct,
    maximax_pct,
    preference_strength,
    risk_attitude_signature
  ) %>%
  arrange(name)

# ============================================================
# Export supplementary CSV
# ============================================================
# The CSV keeps the empty package and provides rounded percentages suitable
# for checking results or reporting in supplementary material.

supplementary_table <- signature_measures %>%
  mutate(
    across(
      c(
        maximin_pct,
        expected_value_pct,
        maximax_pct,
        preference_strength
      ),
      ~ round(.x, 1)
    )
  ) %>%
  select(
    name,
    maximin_pct,
    expected_value_pct,
    maximax_pct,
    preference_strength,
    risk_attitude_signature
  ) %>%
  rename(
    `Measure` = name,
    `Maximin (%)` = maximin_pct,
    `Expected value (%)` = expected_value_pct,
    `Maximax (%)` = maximax_pct,
    `Preference strength` = preference_strength,
    `Risk-attitude signature` = risk_attitude_signature
  )

write_csv(
  supplementary_table,
  file.path(out_dir, "risk_attitude_signature_measures.csv")
)

# ============================================================
# Heatmap: measures favored by risk attitude
# ============================================================
# The heatmap uses the same measure frequencies as the CSV, but excludes the
# empty package so the figure focuses on active biosecurity measures.

heatmap_measure_summary <- measure_rule_summary %>%
  filter(name != "Empty package")

# Order measures by their total selection frequency across all three decision
# rules. This places the most commonly selected measures together in the plot.
measure_order <- heatmap_measure_summary %>%
  group_by(name) %>%
  summarise(
    overall_share = sum(pct_contexts),
    .groups = "drop"
  ) %>%
  arrange(
    desc(overall_share),
    name
  ) %>%
  pull(name)

# Convert percentages to proportions for the legend and fix the decision-rule
# order so the plot reads from cautious to optimistic.
heatmap_data <- heatmap_measure_summary %>%
  mutate(
    name = factor(
      name,
      levels = rev(measure_order)
    ),
    decision_rule = factor(
      decision_rule,
      levels = c("Maximin", "Expected\nvalue", "Maximax")
    ),
    share = pct_contexts / 100
  )

# Draw the heatmap. Each tile is the share of contexts in which a measure
# appears in a winning package under a given decision rule.
p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = decision_rule,
    y = name,
    fill = share
  )
) +
  geom_tile(
    color = "black",
    linewidth = 0.4
  ) +
  scale_x_discrete(
    labels = c(
      "Maximin" = "Maximin",
      "Expected\nvalue" = "Expected value",
      "Maximax" = "Maximax"
    ),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_discrete(
    expand = expansion(mult = 0.01)
  ) +
  scale_fill_gradientn(
    colours = as.character(
      paletteer_c("grDevices::Plasma", 30)
    ),
    limits = c(0, 1),
    name = "Share",
    guide = guide_colourbar(
      frame.colour = "black",
      frame.linewidth = 0.4,
      ticks.colour = NA
    )
  ) +
  labs(
    title = "Measure frequency among winning packages",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "risk_attitude_measure_heatmap.png"),
  p_heatmap,
  width = 6,
  height = 7,
  dpi = 600
)