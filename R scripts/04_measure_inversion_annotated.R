# ============================================================
# Script 04: Convert risk-factor estimates into protective effects
# ============================================================
#
# Purpose
# -------
# This script prepares the measure-level effect estimates for the
# downstream biosecurity package model.
#
# The input file, `inputs/measures_final.csv`, contains published effect
# estimates for individual biosecurity-related measures. Some estimates are
# framed as risk-increasing effects, where values above 1 indicate higher
# disease risk when the risk factor is present. However, the package model
# needs all measures to be expressed in a consistent protective direction,
# where smaller values indicate reduced risk after implementing the measure.
#
# This script therefore:
#   1. Loads the final measure-effect table.
#   2. Converts the mean, lower confidence limit, and upper confidence limit
#      to numeric values.
#   3. Identifies estimates with mean > 1 as risk-factor-framed estimates.
#   4. Inverts those estimates so they can be interpreted as protective
#      effects.
#   5. Keeps estimates with mean <= 1 unchanged, because they are already
#      expressed in the protective direction.
#   6. Exports the converted table for use in later scripts.
#
# Important interpretation
# ------------------------
# For an effect estimate such as an odds ratio or relative risk:
#
#   - mean > 1  means the exposure increases risk.
#   - mean < 1  means the exposure reduces risk.
#
# Because the modelling pipeline evaluates biosecurity measures as protective
# interventions, estimates with mean > 1 are transformed using the reciprocal:
#
#   converted mean = 1 / original mean
#
# The confidence interval is also inverted. When inverting a confidence
# interval, the lower and upper limits swap order:
#
#   converted lower = 1 / original upper
#   converted upper = 1 / original lower
#
# The script also keeps the original values in separate columns so that the
# conversion is transparent and can be checked later.
#
# Main input
# ----------
# inputs/measures_final.csv
#
# Main output
# -----------
# exports/measures_inverted.csv
#
# Required columns in the input file
# ----------------------------------
# mean   : mean effect estimate
# lower  : lower confidence interval limit
# upper  : upper confidence interval limit
#
# Notes
# -----
# This script does not decide whether an estimate is biologically plausible or
# whether it should be included in the analysis. It only standardizes the
# direction of the estimates already selected for the model.
# ============================================================

library(readr)
library(dplyr)
library(writexl)

# ============================================================
# 1. Load measure effect estimates
# ============================================================
#
# The input table contains one row per biosecurity-related measure-effect
# estimate. These estimates may be reported in either a protective direction
# or a risk-increasing direction, depending on how the source study framed the
# exposure.

measures <- read_csv("inputs/measures_final.csv")

# ============================================================
# 2. Invert only risk-increasing effects
# ============================================================
#
# The downstream model expects all effects to be expressed as protective
# effects. Therefore, estimates with mean > 1 are treated as risk-factor-
# framed estimates and are inverted.
#
# Examples:
#   - Original estimate = 2.00
#     Interpretation: the exposure doubles risk.
#     Converted estimate = 1 / 2.00 = 0.50
#     Interpretation: removing or preventing that exposure reduces risk by
#     approximately 50% on the ratio scale.
#
#   - Original estimate = 0.70
#     Interpretation: the measure is already protective.
#     Converted estimate = 0.70
#     No inversion is needed.
#
# The original mean, lower limit, and upper limit are retained so that later
# scripts and manuscript checks can identify which values were converted.

measures_converted <- measures %>%
  mutate(
    # Convert the effect estimate columns to numeric values.
    # This protects against problems if the CSV reader imports any of these
    # columns as character values because of formatting or missing values.
    mean_original  = as.numeric(mean),
    lower_original = as.numeric(lower),
    upper_original = as.numeric(upper),
    
    # Record whether each estimate was inverted.
    # This flag is useful for diagnostics and for documenting which published
    # estimates were originally framed as risk factors.
    inverted = if_else(mean_original > 1, "yes", "no"),
    
    # Invert the mean estimate only when it is greater than 1.
    # Estimates less than or equal to 1 are already in the protective direction
    # and are left unchanged.
    mean = if_else(
      mean_original > 1,
      1 / mean_original,
      mean_original
    ),
    
    # Invert the confidence interval for risk-increasing estimates.
    # The upper original limit becomes the lower converted limit after taking
    # the reciprocal.
    lower = if_else(
      mean_original > 1,
      1 / upper_original,
      lower_original
    ),
    
    # The lower original limit becomes the upper converted limit after taking
    # the reciprocal.
    upper = if_else(
      mean_original > 1,
      1 / lower_original,
      upper_original
    )
  )

# ============================================================
# 3. Export converted effects
# ============================================================
#
# The exported file is the standardized measure-effect table used by the
# downstream risk and package modelling scripts. It includes both the converted
# values and the original values, plus the `inverted` flag.

write_csv(measures_converted, "exports/measures_inverted.csv")
