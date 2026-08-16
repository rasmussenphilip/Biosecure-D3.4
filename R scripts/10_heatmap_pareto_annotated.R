library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)
library(paletteer)

# ============================================================
# 10_heatmap_pareto.R
# ============================================================
#
# Purpose:
# This script summarizes how often each biosecurity measure appears
# among Pareto-optimal packages, and visualizes those frequencies as
# heatmaps. The heatmaps are used to identify which measures recur most
# often among non-dominated packages across regions, herd-benefit types,
# package sizes, and region-by-herd-type contexts.
#
# Inputs:
#   1. exports/package_pareto_frontier.csv
#      - Pareto-frontier package results from the upstream decision-rule
#        and frontier analysis.
#      - Expected to include package_id, package_size, package_contents,
#        analysis_level, context_region, and context_benefit_type.
#
#   2. exports/measures_inverted.csv
#      - Measure metadata from the upstream inversion step.
#      - Used here to map internal measure codes to readable labels.
#
# Outputs:
#   Written to exports/pareto_measure_heatmaps/
#   1. pareto_measure_frequency_vertical_panel.png
#      - Two-panel figure showing measure frequency by region and by
#        herd type.
#   2. pareto_measure_frequency_by_package_size.png
#      - Standalone heatmap showing measure frequency by package size.
#   3. pareto_measure_frequency_by_context.png
#      - Three-panel heatmap showing measure frequency by region within
#        naïve, infected, and average herd contexts.
#
# Interpretation:
# The fill value is the share of Pareto-optimal packages in a given
# context that contain each measure. A value of 1 means the measure is
# present in all Pareto-optimal packages in that context; a value of 0
# means it is absent from all such packages.
#
# Important assumptions:
#   - Package contents are stored as measure names separated by " + ".
#   - The empty package is excluded from measure-frequency calculations.
#   - Frequencies are descriptive summaries of Pareto membership, not
#     causal estimates of individual measure effectiveness.
#   - Region and herd-type summaries average across the relevant Pareto
#     contexts after package-level frontier construction.
#
# Notes:
#   - The computational logic is intentionally unchanged in this annotated
#     version.
#   - The original script header referred to "11_heatmap_pareto.R"; this
#     annotated file follows the uploaded filename while preserving all
#     downstream file names and calculations.
# ============================================================

# Load the Pareto-frontier package results and the measure metadata.
# The Pareto file contains the selected/non-dominated packages from the
# previous script, while the measure file provides readable labels for
# each coded measure.
pareto <- read_csv("exports/package_pareto_frontier.csv", show_col_types = FALSE)
measures <- read_csv("exports/measures_inverted.csv", show_col_types = FALSE)

# Create a dedicated output folder for all heatmap figures generated here.
out_dir <- "exports/pareto_measure_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Measure labels
# ============================================================

# Build a lookup table linking each internal measure code to a readable
# plotting label. If a `name` column exists, it is used as the preferred
# label; otherwise, the raw measure code is retained.
measure_lookup <- measures %>%
  mutate(
    measure = str_trim(as.character(measure)),
    measure_label = if ("name" %in% names(.)) {
      coalesce(str_trim(as.character(name)), measure)
    } else {
      measure
    }
  ) %>%
  select(measure, measure_label) %>%
  distinct() %>%
  filter(!is.na(measure), measure != "")

# Keep a complete list of measures so that the heatmaps include measures
# with zero frequency in a given context, rather than silently dropping
# them from that panel.
all_measures <- measure_lookup %>%
  arrange(measure_label)

# ============================================================
# Clean Pareto data
# ============================================================

# Standardize the Pareto data before plotting.
# - Missing context_region values are treated as pooled "All regions"
#   contexts.
# - Only the three herd-benefit contexts are retained.
# - The empty package is removed because it contains no real measures and
#   would otherwise distort measure-frequency summaries.
pareto <- pareto %>%
  mutate(
    context_region = coalesce(context_region, "All regions"),
    context_benefit_type = as.character(context_benefit_type)
  ) %>%
  filter(
    context_benefit_type %in% c("naive_herd", "infected_herd", "average_herd"),
    package_size > 0
  ) %>%
  mutate(
    herd_type = recode(
      context_benefit_type,
      naive_herd = "Naïve herds",
      infected_herd = "Infected herds",
      average_herd = "Average herds"
    ),
    herd_type = factor(
      herd_type,
      levels = c("Naïve herds", "Infected herds", "Average herds")
    ),
    context_label = paste(
      analysis_level,
      context_region,
      herd_type,
      sep = " | "
    )
  )

# Convert package-level records into measure-level records by splitting
# package_contents into one row per measure. This makes it possible to
# count how often each measure appears among Pareto-optimal packages.
pareto_measures <- pareto %>%
  mutate(measure = str_split(package_contents, " \\+ ")) %>%
  unnest(measure) %>%
  mutate(measure = str_trim(measure)) %>%
  filter(
    measure != "empty_package",
    !is.na(measure),
    measure != ""
  ) %>%
  left_join(measure_lookup, by = "measure") %>%
  mutate(measure_label = coalesce(measure_label, measure))

# ============================================================
# Frequencies by context
# ============================================================

# Count how many distinct Pareto-optimal packages exist in each context.
# These counts become the denominators for the share calculations below.
context_counts <- pareto %>%
  group_by(analysis_level, context_region, herd_type, context_label) %>%
  summarise(
    n_pareto_packages = n_distinct(package_id),
    .groups = "drop"
  )

# Count how many distinct Pareto-optimal packages in each context contain
# each individual measure.
measure_context_counts <- pareto_measures %>%
  group_by(
    analysis_level,
    context_region,
    herd_type,
    context_label,
    measure,
    measure_label
  ) %>%
  summarise(
    n_with_measure = n_distinct(package_id),
    .groups = "drop"
  )

# Combine every context with every measure, then join the observed counts.
# The crossing() step ensures that zero-frequency combinations are retained
# and shown as empty/low-frequency cells in the heatmaps.
context_frequency <- context_counts %>%
  crossing(all_measures) %>%
  left_join(
    measure_context_counts,
    by = c(
      "analysis_level",
      "context_region",
      "herd_type",
      "context_label",
      "measure",
      "measure_label"
    )
  ) %>%
  mutate(
    n_with_measure = coalesce(n_with_measure, 0L),
    share = n_with_measure / n_pareto_packages
  )

# Order measures by their average frequency across all contexts so that
# the y-axis places the most commonly recurring Pareto measures together.
measure_order <- context_frequency %>%
  group_by(measure_label) %>%
  summarise(
    avg_share = mean(share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_share), measure_label) %>%
  pull(measure_label)

# ============================================================
# Shared plot elements
# ============================================================

# Define a shared fill scale so all heatmaps use the same 0-1 frequency
# interpretation and colour mapping.
shared_fill_scale <- scale_fill_gradientn(
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
)

# Define a shared visual theme to keep the individual figures consistent.
shared_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

# ============================================================
# A. By region
# ============================================================

# Average measure frequencies across all non-pooled contexts within each
# region. The pooled "All regions" context is excluded so this plot shows
# country/region-specific patterns only.
by_region <- context_frequency %>%
  filter(context_region != "All regions") %>%
  group_by(context_region, measure, measure_label) %>%
  summarise(
    share = mean(share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    measure_label = factor(measure_label, levels = rev(measure_order))
  )

# Heatmap A: rows are measures and columns are regions.
p_region <- ggplot(
  by_region,
  aes(x = context_region, y = measure_label, fill = share)
) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_x_discrete(expand = expansion(mult = 0.01)) +
  scale_y_discrete(expand = expansion(mult = 0.01)) +
  shared_fill_scale +
  labs(
    title = "By region",
    x = NULL,
    y = NULL
  ) +
  shared_theme

# ============================================================
# B. By herd type
# ============================================================

# Average measure frequencies across contexts within each herd-benefit
# type: naïve, infected, and average herds.
by_herd <- context_frequency %>%
  group_by(herd_type, measure, measure_label) %>%
  summarise(
    share = mean(share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    measure_label = factor(measure_label, levels = rev(measure_order))
  )

# Heatmap B: rows are measures and columns are herd-benefit types.
p_herd <- ggplot(
  by_herd,
  aes(x = herd_type, y = measure_label, fill = share)
) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_x_discrete(
    labels = c(
      "Naïve herds" = "Naïve",
      "Infected herds" = "Infected",
      "Average herds" = "Average"
    ),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_discrete(expand = expansion(mult = 0.01)) +
  shared_fill_scale +
  labs(
    title = "By herd type",
    x = NULL,
    y = NULL
  ) +
  shared_theme +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

# ============================================================
# Combined vertical panel: region + herd type only
# ============================================================

# Combine the region and herd-type heatmaps into a vertical figure with a
# shared legend. Patchwork tags the panels as A. and B.
combined_panel <- (
  p_region /
    p_herd
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.tag.position = c(0.01, 0.99)
  )

# Export the two-panel summary figure.
ggsave(
  file.path(out_dir, "pareto_measure_frequency_vertical_panel.png"),
  combined_panel,
  width = 6,
  height = 8,
  dpi = 600
)

# ============================================================
# C. By package size — standalone plot
# ============================================================

# Count Pareto-optimal packages by package size. These counts are used as
# denominators for package-size-specific measure frequencies.
size_counts <- pareto %>%
  group_by(package_size) %>%
  summarise(
    n_pareto_packages = n_distinct(package_id),
    .groups = "drop"
  )

# Count how many packages of each size contain each measure.
measure_size_counts <- pareto_measures %>%
  group_by(package_size, measure, measure_label) %>%
  summarise(
    n_with_measure = n_distinct(package_id),
    .groups = "drop"
  )

# Create the complete package-size by measure grid, fill missing counts
# with zero, and calculate the within-size frequency of each measure.
by_size <- size_counts %>%
  crossing(all_measures) %>%
  left_join(
    measure_size_counts,
    by = c("package_size", "measure", "measure_label")
  ) %>%
  mutate(
    n_with_measure = coalesce(n_with_measure, 0L),
    share = n_with_measure / n_pareto_packages,
    measure_label = factor(measure_label, levels = rev(measure_order))
  )

# Heatmap C: rows are measures and columns are package sizes.
p_size <- ggplot(
  by_size,
  aes(x = factor(package_size), y = measure_label, fill = share)
) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_x_discrete(expand = expansion(mult = 0.01)) +
  scale_y_discrete(expand = expansion(mult = 0.01)) +
  shared_fill_scale +
  labs(
    title = "Measure frequency among Pareto-optimal packages",
    x = "Package size",
    y = NULL
  ) +
  shared_theme

# Export the standalone package-size heatmap.
ggsave(
  file.path(out_dir, "pareto_measure_frequency_by_package_size.png"),
  p_size,
  width = 9,
  height = 7,
  dpi = 600
)

# ============================================================
# Full context heatmap — three-panel horizontal plot
# ============================================================

# Prepare the full context heatmap, keeping separate panels for herd type
# and separate x-axis columns for regions. This is the most detailed view
# of measure frequencies across decision contexts.
context_panel_data <- context_frequency %>%
  filter(context_region != "All regions") %>%
  mutate(
    region_label = recode(
      context_region,
      US = "USA",
      .default = context_region
    ),
    herd_panel = factor(
      herd_type,
      levels = c(
        "Naïve herds",
        "Infected herds",
        "Average herds"
      ),
      labels = c(
        "Naïve",
        "Infected",
        "Average"
      )
    ),
    region_label = factor(
      region_label,
      levels = c("CAN", "DK", "IT", "USA")
    ),
    measure_label = factor(
      measure_label,
      levels = rev(measure_order)
    )
  )

# Full context heatmap: each facet is a herd-benefit type and each column
# within a facet is a region.
p_context <- ggplot(
  context_panel_data,
  aes(x = region_label, y = measure_label, fill = share)
) +
  geom_tile(color = "black", linewidth = 0.4) +
  facet_wrap(
    ~ herd_panel,
    nrow = 1
  ) +
  scale_x_discrete(expand = expansion(mult = 0.01)) +
  scale_y_discrete(expand = expansion(mult = 0.01)) +
  shared_fill_scale +
  labs(
    title = "Measure frequency among Pareto-optimal packages",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 8
    ),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(b = -3),
    legend.box.margin = margin(b = -5)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),
      barheight = unit(0.35, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.4,
      ticks.colour = NA
    )
  )

# Export the full context heatmap.
ggsave(
  file.path(out_dir, "pareto_measure_frequency_by_context.png"),
  p_context,
  width = 10,
  height = 7,
  dpi = 600
)