# ============================================================
# 07_net_gmpcy_across_packages.R
# ============================================================
#
# Detailed script summary
# ------------------------------------------------------------
# This script converts package-level gross benefits into package-level
# net benefits by subtracting estimated implementation costs.
#
# The previous script produces package-level changes in gross margin per
# cow-year (delta_gmpcy). Those values represent the economic benefit of
# each biosecurity package before accounting for the cost of implementing
# the package. This script adds the cost side of the analysis.
#
# The workflow is:
#   1. Load package-level GM benefits and labelled per-measure costs.
#   2. Check that both input files contain the columns required downstream.
#   3. Split each package into its component measures.
#   4. Match each measure to its region-specific low, reference, and high
#      cost estimate.
#   5. Sum measure costs to package-level costs within each region.
#   6. Assign zero cost to the empty package.
#   7. Join package-level costs onto package-level benefits.
#   8. Calculate net benefit under low, reference, and high cost scenarios.
#   9. Export the final package-level net benefit file.
#
# Conceptually, this script answers:
#   "After accounting for implementation costs, how much value does each
#    biosecurity package generate per cow-year in each region and disease
#    state scenario?"
#
# Important modelling assumptions
# ------------------------------------------------------------
# - Costs are additive across measures within a package.
# - Costs are region-specific, because implementation costs differ across
#   countries/regions.
# - The same package cost is applied across benefit types, aggregation
#   methods, and joint disease states within a region.
# - The empty package has zero implementation cost.
# - Net benefits are expressed in the same unit as delta_gmpcy and costs:
#   euros per cow-year.
#
# Inputs
# ------------------------------------------------------------
#   exports/package_gm_values.csv
#     Package-level gross margin benefits from the previous step.
#
#   exports/costs_euro_pcy_labeled.csv
#     Region-specific low, reference, and high cost estimates for each
#     individual biosecurity measure, expressed in euros per cow-year.
#
# Output
# ------------------------------------------------------------
#   exports/package_net_benefits.csv
#     One row per package × region × benefit type × aggregation method ×
#     joint state, with package costs and net benefits included.
#
# ============================================================

# Load packages used for data manipulation, string handling, and CSV I/O.
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ============================================================
# 1. Load inputs
# ============================================================

# Package-level gross margin benefits. These are the gross benefits before
# subtracting implementation costs.
gm_values <- read_csv(
  "exports/package_gm_values.csv",
  show_col_types = FALSE
)

# Labelled cost estimates for individual measures. These should already be
# harmonised to euros per cow-year and matched to the measure names used in
# the package definitions.
costs <- read_csv(
  "exports/costs_euro_pcy_labeled.csv",
  show_col_types = FALSE
)

# ============================================================
# 2. Check required columns
# ============================================================

# Define the columns needed from the gross-margin benefit file. These identify
# the package, the region, the scenario/state, and the gross benefit estimate.
required_gm_cols <- c(
  "package_id",
  "package_size",
  "package_contents",
  "region",
  "benefit_type",
  "aggregation_method",
  "joint_state",
  "payoff_state",
  "delta_gmpcy"
)

# Define the columns needed from the cost file. Costs are provided for each
# measure and region, with low/reference/high values used for sensitivity.
required_cost_cols <- c(
  "measure",
  "region",
  "cost_low",
  "cost_ref",
  "cost_high"
)

# Identify any missing required columns before doing joins or calculations.
# This makes errors easier to diagnose than allowing downstream failures.
missing_gm_cols <- setdiff(required_gm_cols, names(gm_values))
missing_cost_cols <- setdiff(required_cost_cols, names(costs))

if (length(missing_gm_cols) > 0) {
  stop(
    "package_gm_values.csv is missing required columns: ",
    paste(missing_gm_cols, collapse = ", ")
  )
}

if (length(missing_cost_cols) > 0) {
  stop(
    "costs_euro_pcy_labeled.csv is missing required columns: ",
    paste(missing_cost_cols, collapse = ", ")
  )
}

# ============================================================
# 3. Create package-to-measure map
# ============================================================

# Extract the unique package definitions from the gross-margin results. Each
# package is identified by an ID, its size, and a text string listing the
# measures it contains.
packages <- gm_values %>%
  distinct(
    package_id,
    package_size,
    package_contents
  )

# Convert the package contents from a single string into one row per
# package-measure pair. Package contents are assumed to be separated by
# " + ". The empty package is excluded here because it contains no measures
# and will be assigned a zero cost later.
package_measure_map <- packages %>%
  filter(package_contents != "empty_package") %>%
  mutate(
    measure = str_split(package_contents, " \\+ ")
  ) %>%
  unnest(measure) %>%
  mutate(
    measure = str_trim(measure)
  )

# ============================================================
# 4. Sum measure costs to package-level costs by region
# ============================================================

# Use the regions that appear in the gross-margin benefit results. This keeps
# the cost expansion aligned with the regions actually used in the model.
gm_regions <- sort(unique(gm_values$region))

# Expand each package-measure pair across all modelled regions, then attach
# the corresponding measure-level costs for each region.
package_cost_components <- package_measure_map %>%
  crossing(region = gm_regions) %>%
  left_join(
    costs,
    by = c("measure", "region")
  )

# Check whether any package-measure-region combinations failed to match to a
# complete set of cost estimates. Missing costs would otherwise produce NA net
# benefits, so the script stops with an explicit list of missing combinations.
missing_cost_rows <- package_cost_components %>%
  filter(
    is.na(cost_low) |
      is.na(cost_ref) |
      is.na(cost_high)
  ) %>%
  distinct(measure, region)

if (nrow(missing_cost_rows) > 0) {
  stop(
    "Missing cost estimates for some measure-region combinations:\n",
    paste(
      paste(missing_cost_rows$measure, missing_cost_rows$region, sep = " / "),
      collapse = "\n"
    )
  )
}

# Sum individual measure costs to the package level. This assumes that total
# package cost is the sum of all component measure costs within the package.
package_costs <- package_cost_components %>%
  group_by(
    package_id,
    package_size,
    package_contents,
    region
  ) %>%
  summarise(
    package_cost_low  = sum(cost_low),
    package_cost_ref  = sum(cost_ref),
    package_cost_high = sum(cost_high),
    .groups = "drop"
  )

# Add the empty package back into the package-cost table. It represents the
# no-intervention comparator and therefore has zero implementation cost in
# every region.
empty_package_costs <- packages %>%
  filter(package_contents == "empty_package") %>%
  crossing(region = gm_regions) %>%
  mutate(
    package_cost_low = 0,
    package_cost_ref = 0,
    package_cost_high = 0
  )

# Combine ordinary package costs with the zero-cost empty package and sort for
# readability.
package_costs <- bind_rows(
  package_costs,
  empty_package_costs
) %>%
  arrange(region, package_id)

# ============================================================
# 5. Join GM benefits and costs
# ============================================================

# Attach package-level costs to every gross-margin benefit row. Costs are joined
# by package and region, then subtracted from delta_gmpcy to obtain net benefits
# under low, reference, and high cost assumptions.
package_net_benefits <- gm_values %>%
  left_join(
    package_costs,
    by = c(
      "package_id",
      "package_size",
      "package_contents",
      "region"
    )
  ) %>%
  mutate(
    net_benefit_low  = delta_gmpcy - package_cost_low,
    net_benefit_ref  = delta_gmpcy - package_cost_ref,
    net_benefit_high = delta_gmpcy - package_cost_high
  ) %>%
  select(
    package_id,
    package_size,
    package_contents,
    region,
    benefit_type,
    aggregation_method,
    joint_state,
    payoff_state,
    delta_gmpcy,
    package_cost_low,
    package_cost_ref,
    package_cost_high,
    net_benefit_low,
    net_benefit_ref,
    net_benefit_high
  ) %>%
  arrange(
    region,
    package_id,
    benefit_type,
    aggregation_method,
    joint_state
  )

# ============================================================
# 6. Export final net-benefit file
# ============================================================

# Create the output directory if needed.
out_dir <- "exports"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save one final file containing both gross benefits, package costs, and net
# benefits. This file is used by downstream game-theoretic and frontier scripts.
write_csv(
  package_net_benefits,
  file.path(out_dir, "package_net_benefits.csv")
)
