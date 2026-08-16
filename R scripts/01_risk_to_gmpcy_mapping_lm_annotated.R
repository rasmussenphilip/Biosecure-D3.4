library(dplyr)

# ============================================================
# 01_risk_to_gmpcy_mapping_lm.R
# ============================================================
#
# Purpose
# -------
# This script converts raw SimHerd output into a standardized mapping between
# disease risk and gross margin per cow-year (GM/PCY). The resulting file is
# used downstream to translate simulated disease-risk reductions from
# biosecurity packages into economic benefits.
#
# In practical terms, the script does four things:
#   1. Loads raw SimHerd outputs for each country/region and disease-risk
#      scenario.
#   2. Renames and standardizes the disease-risk and gross-margin columns.
#   3. Converts SimHerd gross margins from DKK per cow-year to EUR per cow-year.
#   4. Calculates absolute and percentage changes in GM/PCY relative to the
#      disease-free baseline within each region.
#
# Key assumptions
# ---------------
# - SimHerd reports gross margin values in DKK per cow-year.
# - Downstream cost and benefit calculations are expressed in EUR per cow-year,
#   so all SimHerd gross-margin values are converted from DKK to EUR here.
# - Disease-risk values in the raw SimHerd file are expressed as percentages and
#   are converted to proportions by dividing by 100.
# - The disease-free baseline is defined separately for each region as the row
#   where BVD, PTB, and Salmonella risks are all zero.
# - Economic disease impact is measured as the change in GM/PCY relative to that
#   region-specific disease-free baseline.
#
# Main output
# -----------
# exports/risk_to_gmpcy.csv
#
# This output contains one row per region-risk scenario, including:
#   - region
#   - risk_bvd, risk_ptb, risk_sal
#   - original SimHerd GM/PCY in DKK
#   - converted GM/PCY in EUR
#   - region-specific disease-free baseline GM/PCY
#   - absolute change in GM/PCY relative to baseline
#   - percentage change in GM/PCY relative to baseline
#
# ============================================================
# 1. Inputs and currency conversion
# ============================================================

# All SimHerd GM values are reported in DKK per cow-year.
# The downstream package model expresses costs and benefits in EUR per cow-year,
# so the first step is to define a DKK-to-EUR conversion factor.
#
# The value below corresponds to 745.8895 DKK per 100 EUR, i.e. approximately
# 7.458895 DKK per EUR.
dkk_per_euro <- 745.8895 / 100

# Convert from DKK to EUR by taking the inverse of DKK per EUR.
eur_per_dkk <- 1 / dkk_per_euro

# ============================================================
# 2. Load raw SimHerd output
# ============================================================

# The raw SimHerd file contains simulated gross-margin outcomes under different
# combinations of disease risks for BVD, PTB, and Salmonella.
#
# Expected raw columns used in this script:
#   - Country: region/country identifier
#   - BDV_inc: BVD risk/incidence input used by SimHerd
#   - PTB_inc: PTB risk/incidence input used by SimHerd
#   - SAL_inc: Salmonella risk/incidence input used by SimHerd
#   - NetReturn_PrCowYr: gross margin per cow-year in DKK
#
# Note: the raw column is named BDV_inc in the SimHerd output, but it is treated
# as BVD risk in the cleaned dataset.
simherd_risk_gmpcy <- read.csv(
  "inputs/SimHerd/results_raw_PR.txt",
  header = TRUE,
  sep = ";"
)

# ============================================================
# 3. Clean, standardize, and convert GM to EUR
# ============================================================

risk_gmpcy <- simherd_risk_gmpcy %>%
  
  # Keep only the columns needed for the risk-to-economics mapping.
  select(Country, BDV_inc, PTB_inc, SAL_inc, NetReturn_PrCowYr) %>%
  
  # Rename columns to the names used throughout the downstream biosecurity
  # package scripts.
  rename(
    region      = Country,
    risk_bvd    = BDV_inc,
    risk_ptb    = PTB_inc,
    risk_sal    = SAL_inc,
    gmpcy_dkk   = NetReturn_PrCowYr
  ) %>%
  mutate(
    # Convert disease-risk inputs from percentages to proportions.
    # Example: 20 becomes 0.20.
    risk_bvd = risk_bvd / 100,
    risk_ptb = risk_ptb / 100,
    risk_sal = risk_sal / 100,
    
    # Convert absolute GM/PCY from DKK to EUR.
    gmpcy = gmpcy_dkk * eur_per_dkk
  ) %>%
  
  # Baseline gross margin needs to be calculated within each region because the
  # disease-free SimHerd gross margin may differ across countries/regions.
  group_by(region) %>%
  mutate(
    # Disease-free GM/PCY in EUR.
    # This is the reference point for all economic impacts in the same region.
    # The [1] is used in case there is more than one matching baseline row.
    baseline_gmpcy = gmpcy[
      risk_bvd == 0 &
        risk_ptb == 0 &
        risk_sal == 0
    ][1],
    
    # Absolute change in GM/PCY relative to the region-specific disease-free
    # baseline. Values below zero indicate economic losses relative to the
    # disease-free state.
    c_gmpcy = gmpcy - baseline_gmpcy,
    
    # Percentage change in GM/PCY relative to the disease-free baseline.
    # This is useful for reporting and diagnostics, while downstream economic
    # calculations primarily use the EUR per cow-year values.
    p_gmpcy = 100 * c_gmpcy / baseline_gmpcy
  ) %>%
  ungroup()

# ============================================================
# 4. Export cleaned risk-to-GM/PCY mapping
# ============================================================

# This file is used by downstream scripts to estimate the economic value of
# reducing disease risk through biosecurity packages.
write.csv(
  risk_gmpcy,
  "exports/risk_to_gmpcy.csv",
  row.names = FALSE
)
