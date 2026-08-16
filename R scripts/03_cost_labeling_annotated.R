# ============================================================
# Script 03: Label low, reference, and high cost scenarios
# ============================================================
#
# Purpose
# -------
# This script takes the regional per-cow-year cost estimates created in
# Script 02 (`02_price_extrapolation.R`) and converts multiple candidate
# estimates for each biosecurity measure into three labelled cost scenarios:
#
#   1. `cost_low`:  the lowest available estimate for a measure-region pair
#   2. `cost_ref`:  the reference or central estimate
#   3. `cost_high`: the highest available estimate for a measure-region pair
#
# These scenarios are later used in the package-level cost-benefit model to
# evaluate how sensitive net benefits are to uncertainty in implementation
# costs. Keeping low/reference/high cost scenarios separate makes it possible
# to test optimistic, central, and pessimistic cost assumptions without
# re-estimating the biological effects of the biosecurity packages.
#
# Input
# -----
# exports/costs_euro_pcy.csv
#
# Expected structure of the input file
# ------------------------------------
# The input file should contain one row per cost estimate or cost option and
# separate region columns for the annual cost per cow-year:
#
#   - DK
#   - IRL
#   - CAN
#   - US
#   - IT
#
# It should also contain identifying/descriptive columns such as:
#
#   - measure
#   - name
#   - short_description
#
# Output
# ------
# exports/costs_euro_pcy_labeled.csv
#
# The output file contains one row per unique combination of:
#
#   - measure
#   - name
#   - region
#
# and includes the low, reference, and high cost estimates, together with
# short descriptions of where those scenario values came from.
#
# Main assumptions
# ----------------
# For each measure-region pair, the script assumes that the available cost
# estimates can be ordered from lowest to highest and summarized as follows:
#
#   - If only one estimate is available, it is used for all central/reference
#     purposes.
#   - If two estimates are available, the reference value is the midpoint.
#   - If three estimates are available, the reference value is the median.
#
# The script currently returns `NA` for the reference value if more than three
# estimates are available, because no rule has been defined for that case.
#
# ============================================================

library(dplyr)
library(tidyr)
library(readr)

# ============================================================
# 1. Read extrapolated regional cost estimates
# ------------------------------------------------------------
# This file is produced by the previous price extrapolation script.
#
# At this stage, costs are still in a wide regional format, meaning that each
# region has its own column (DK, IRL, CAN, US, IT). This format is convenient
# for manual checking but not ideal for grouped summaries, so the next step
# reshapes it to long format.
# ============================================================

costs_raw <- read_csv("exports/costs_euro_pcy.csv", show_col_types = FALSE)

# ============================================================
# 2. Convert regional cost columns to long format
# ------------------------------------------------------------
# `pivot_longer()` turns the five region columns into two columns:
#
#   - `region`: identifies which region the cost applies to
#   - `cost_euro_pcy`: stores the annual cost in EUR per cow-year
#
# After this step, each row represents one measure-cost-option-region
# combination. This makes it possible to group by measure and region and then
# calculate low, reference, and high cost scenarios.
# ============================================================

costs_long <- costs_raw %>%
  pivot_longer(
    cols = c(DK, IRL, CAN, US, IT),
    names_to = "region",
    values_to = "cost_euro_pcy"
  )

# ============================================================
# 3. Assign low, reference, and high cost scenarios
# ------------------------------------------------------------
# For each measure, measure name, and region, the script:
#
#   1. Sorts available cost estimates from lowest to highest.
#   2. Stores the number of available estimates.
#   3. Defines the low-cost scenario as the first sorted value.
#   4. Defines the high-cost scenario as the last sorted value.
#   5. Defines the reference scenario using a simple rule based on the number
#      of available estimates.
#
# The corresponding `short_description` is retained for the low and high
# scenarios. For the reference scenario:
#
#   - one estimate: use that estimate's description
#   - two estimates: label the reference as the midpoint of two estimates
#   - three estimates: use the description attached to the middle estimate
#
# This provides transparent labels for later tables and sensitivity analyses.
# ============================================================

costs_summary <- costs_long %>%
  group_by(measure, name, region) %>%
  arrange(cost_euro_pcy, .by_group = TRUE) %>%
  summarise(
    # Number of cost estimates available for this measure-region pair.
    n_estimates = n(),

    # Lowest available annual cost estimate and its description.
    cost_low = first(cost_euro_pcy),
    cost_low_description = first(short_description),

    # Highest available annual cost estimate and its description.
    cost_high = last(cost_euro_pcy),
    cost_high_description = last(short_description),

    # Reference cost estimate.
    # The rule depends on whether one, two, or three estimates are available.
    cost_ref = case_when(
      n_estimates == 1 ~ first(cost_euro_pcy),
      n_estimates == 2 ~ mean(cost_euro_pcy, na.rm = TRUE),
      n_estimates == 3 ~ median(cost_euro_pcy, na.rm = TRUE),
      TRUE ~ NA_real_
    ),

    # Human-readable description of the reference estimate.
    cost_ref_description = case_when(
      n_estimates == 1 ~ first(short_description),
      n_estimates == 2 ~ "Midpoint of two estimates",
      n_estimates == 3 ~ nth(short_description, 2),
      TRUE ~ NA_character_
    ),

    .groups = "drop"
  ) %>%
  arrange(measure, region)

# ============================================================
# 4. Export labelled cost scenarios
# ------------------------------------------------------------
# The output is the cost input used by later scripts in the modelling
# pipeline. It provides one low, one reference, and one high annual cost
# estimate for each measure-region pair.
# ============================================================

write_csv(costs_summary, "exports/costs_euro_pcy_labeled.csv")
