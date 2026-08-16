# ============================================================
# Script 02: Extrapolate Danish biosecurity measure costs
#            to study regions and express them as €/cow-year
# ============================================================
#
# Purpose
# -------
# This script takes Danish reference cost estimates for individual
# biosecurity measures and converts them into comparable annual costs
# per cow-year for each study region: Denmark (DK), Ireland (IRL),
# Canada (CAN), the United States (US), and Italy (IT).
#
# The resulting table is used downstream when net benefits are calculated
# by subtracting region-specific implementation costs from region-specific
# disease-loss reductions / gross-margin gains.
#
# Main workflow
# -------------
# 1. Load Danish cost estimates and regional economic indexes.
# 2. Define economic constants used for inflation, currency conversion,
#    and annualisation.
# 3. Clean regional wage and GDP-per-capita information.
# 4. Construct regional scaling factors relative to Denmark.
# 5. Clean measure-level cost inputs and inflate older Danish prices to
#    the target year.
# 6. Expand every Danish cost row to every region.
# 7. Scale investment, labour, maintenance, and consumable costs by the
#    relevant regional economic index.
# 8. Annualise one-time investments and installation labour over the
#    stated lifespan using an annuity formula.
# 9. Convert all costs to EUR and divide by herd size to obtain
#    annual cost per cow-year.
# 10. Collapse multi-row grouped measures into a single cost estimate.
# 11. Export one wide table with one row per measure and one cost column
#     per region.
#
# Key assumptions
# ---------------
# - Denmark is the reference country for the original cost estimates.
# - Danish costs from 2015 are inflated to 2024 using Danish CPI.
# - Costs already in other years are assumed to be in target-year terms.
# - Labour costs scale with regional wage rates.
# - Investment, maintenance, and consumables scale with GDP per capita.
# - One-time costs are converted to annual costs using a 5% discount rate.
# - Final costs are expressed as EUR per cow-year.
#
# Required input files
# --------------------
# - inputs/economic_indexes.csv
#     Region-level herd size, GDP per capita, wage rate, and wage currency.
# - inputs/costs_final.csv
#     Danish reference costs and measure metadata.
#
# Output file
# -----------
# - exports/costs_euro_pcy.csv
#     Region-specific implementation costs in EUR per cow-year.
#
# Notes for interpretation
# ------------------------
# The output is not intended to represent newly observed cost data for each
# country. It is an extrapolated cost table derived from Danish reference
# costs, adjusted using simple regional economic scalers. The estimates are
# therefore best interpreted as harmonised, model-ready cost approximations
# for scenario analysis rather than precise country-specific accounting data.
# ============================================================

# Load packages
library(dplyr)
library(tidyr)
library(readr)

# ============================================================
# 1. Inputs
# ------------------------------------------------------------
# Define price-index, exchange-rate, and discount-rate inputs
# used to standardize Danish reference costs and extrapolate
# them across regions.
# ============================================================

# Manually set Danish CPI values
dk_cpi_2015 <- 107.1332917
dk_cpi_2024 <- 127.2922061

# Manually set exchange rates
dkk_per_euro <- 745.8895 / 100
euro_per_usd <- 0.92

# Manually set discount rate for annualization
r <- 0.05

# Load economic index values for each region
economic_indexes <- read.csv("inputs/economic_indexes.csv")

# Load Danish reference cost estimates for each measure
costs_estimates <- read.csv("inputs/costs_final.csv")


# ============================================================
# 2. Helper functions
# ------------------------------------------------------------
# annualize_cost():
#   Converts a one-time cost into an equivalent annual cost
#   using the annuity formula over the stated lifespan.
#
# inflate_to_target():
#   Inflates 2015 DKK values into 2024 DKK using Danish CPI.
#   Other years are assumed already to be in target-year terms.
# ============================================================

annualize_cost <- function(cost, life_span, r) {
  ifelse(
    is.na(cost) | is.na(life_span) | life_span <= 0,
    0,
    cost * r / (1 - (1 + r)^(-life_span))
  )
}

inflate_to_target <- function(x, year, cpi_2015, cpi_2024) {
  ifelse(
    is.na(x),
    0,
    dplyr::case_when(
      year == 2015 ~ x * (cpi_2024 / cpi_2015),
      TRUE ~ x
    )
  )
}


# ============================================================
# 3. Clean and prepare regional economic indexes
# ------------------------------------------------------------
# Convert wage rates to EUR and construct region-specific
# scaling factors relative to Denmark.
#
# Assumptions:
# - labor costs scale with wage rates
# - capital and consumable costs scale with GDP per capita
# ============================================================

economic_indexes <- economic_indexes %>%
  mutate(
    wage_currency = tolower(trimws(wage_currency)),
    wage_rate_eur = case_when(
      wage_currency %in% c("euro", "eur") ~ wage_rate,
      wage_currency == "usd" ~ wage_rate * euro_per_usd,
      TRUE ~ NA_real_
    )
  )

# Extract Danish reference values
dk_gdppc <- economic_indexes$gdppc[economic_indexes$region == "DK"]
dk_wage  <- economic_indexes$wage_rate_eur[economic_indexes$region == "DK"]

economic_indexes <- economic_indexes %>%
  mutate(
    labor_scaler = wage_rate_eur / dk_wage,
    capital_scaler = gdppc / dk_gdppc,
    consumable_scaler = gdppc / dk_gdppc
  )


# ============================================================
# 4. Clean and prepare measure-level cost data
# ------------------------------------------------------------
# Standardize text fields, convert numeric inputs, replace
# blanks with NA where appropriate, and set missing cost
# components to zero.
#
# Monetary values are then inflated from 2015 DKK to 2024 DKK
# where needed.
# ============================================================

costs_estimates <- costs_estimates %>%
  mutate(
    source_row_id = row_number(),
    
    measure = trimws(measure),
    name = trimws(name),
    short_description = trimws(short_description),
    
    label = trimws(label),
    label = ifelse(label == "investment capacity", "investment_capacity", label),
    
    group_id = trimws(group_id),
    group_id = ifelse(group_id == "", NA, group_id),
    
    consumable_type_01 = trimws(consumable_type_01),
    consumable_type_01 = ifelse(consumable_type_01 == "", NA, consumable_type_01),
    
    source = trimws(source),
    source = ifelse(source == "", NA, source),
    
    comment = trimws(comment),
    comment = ifelse(comment == "", NA, comment),
    
    investment_cost = coalesce(investment_cost, 0),
    investment_annual_maintenance = coalesce(investment_annual_maintenance, 0),
    life_span = suppressWarnings(as.numeric(life_span)),
    investment_units_per_cow = coalesce(investment_units_per_cow, 0),
    investment_units_per_farm = coalesce(investment_units_per_farm, 0),
    reference_farm_size = suppressWarnings(as.numeric(reference_farm_size)),
    labor_hours_per_cow_installation = coalesce(labor_hours_per_cow_installation, 0),
    labor_hours_per_cow_annual = coalesce(labor_hours_per_cow_annual, 0),
    consumable_type_01_per_cow = coalesce(consumable_type_01_per_cow, 0),
    consumable_type_01_cost_per_unit = coalesce(consumable_type_01_cost_per_unit, 0)
  ) %>%
  mutate(
    investment_cost_target_dkk =
      inflate_to_target(investment_cost, investment_cost_year, dk_cpi_2015, dk_cpi_2024),
    
    investment_annual_maintenance_target_dkk =
      inflate_to_target(investment_annual_maintenance, investment_cost_year, dk_cpi_2015, dk_cpi_2024),
    
    consumable_type_01_cost_per_unit_target_dkk =
      inflate_to_target(consumable_type_01_cost_per_unit, investment_cost_year, dk_cpi_2015, dk_cpi_2024)
  )


# ============================================================
# 5. Expand each Danish cost row to all regions
# ------------------------------------------------------------
# Each measure-row is replicated across regions and scaled using
# region-specific herd sizes, wage rates, and GDP-per-capita
# multipliers.
#
# Costs are separated into:
# - investment
# - annual maintenance
# - installation labor
# - annual labor
# - consumables
#
# The final measure is annual cost per cow-year in EUR.
# ============================================================

costs_rowlevel <- costs_estimates %>%
  crossing(
    economic_indexes %>%
      select(region, herd_size, wage_rate_eur, labor_scaler, capital_scaler, consumable_scaler)
  ) %>%
  mutate(
    # Total investment units needed for this region
    investment_units_total =
      investment_units_per_cow * herd_size +
      case_when(
        investment_units_per_farm > 0 & !is.na(reference_farm_size) & reference_farm_size > 0 ~
          investment_units_per_farm * (herd_size / reference_farm_size),
        TRUE ~ investment_units_per_farm
      ),
    
    # Capital / investment costs scaled by GDP per capita
    investment_total_dkk =
      investment_cost_target_dkk * investment_units_total * capital_scaler,
    
    annualized_investment_dkk =
      annualize_cost(investment_total_dkk, life_span, r),
    
    annual_maintenance_dkk =
      investment_annual_maintenance_target_dkk * investment_units_total * capital_scaler,
    
    # Installation labor treated as a one-time cost and annualized
    installation_labor_total_eur =
      labor_hours_per_cow_installation * herd_size * wage_rate_eur,
    
    annualized_installation_labor_eur =
      annualize_cost(installation_labor_total_eur, life_span, r),
    
    # Ongoing annual labor cost
    annual_labor_eur =
      labor_hours_per_cow_annual * herd_size * wage_rate_eur,
    
    # Consumables scaled by GDP per capita
    annual_consumable_01_dkk =
      consumable_type_01_per_cow * herd_size *
      consumable_type_01_cost_per_unit_target_dkk * consumable_scaler,
    
    # Convert DKK-denominated costs to EUR
    annualized_investment_eur = annualized_investment_dkk / dkk_per_euro,
    annual_maintenance_eur = annual_maintenance_dkk / dkk_per_euro,
    annual_consumable_01_eur = annual_consumable_01_dkk / dkk_per_euro,
    
    # Total annual cost in EUR
    total_annual_cost_eur =
      annualized_investment_eur +
      annual_maintenance_eur +
      annualized_installation_labor_eur +
      annual_labor_eur +
      annual_consumable_01_eur,
    
    # Express as annual cost per cow-year
    total_cost_eur_per_cow_year =
      total_annual_cost_eur / herd_size
  ) %>%
  select(
    source_row_id,
    region,
    measure,
    name,
    short_description,
    group_id,
    total_cost_eur_per_cow_year
  )


# ============================================================
# 6. Collapse grouped rows
# ------------------------------------------------------------
# Some measures are represented by multiple component rows in
# the source file. Where group_id is present, sum these rows
# within region so each grouped measure is represented once.
# ============================================================

grouped_rows <- costs_rowlevel %>%
  filter(!is.na(group_id)) %>%
  group_by(group_id, region) %>%
  summarise(
    order_id = min(source_row_id),
    measure = first(measure),
    name = first(name),
    short_description = first(short_description),
    total_cost_eur_per_cow_year = sum(total_cost_eur_per_cow_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    row_id = group_id
  )


# ============================================================
# 7. Keep ungrouped rows as separate measures
# ============================================================

ungrouped_rows <- costs_rowlevel %>%
  filter(is.na(group_id)) %>%
  transmute(
    row_id = paste0("row_", source_row_id),
    order_id = source_row_id,
    measure = measure,
    name = name,
    short_description = short_description,
    region = region,
    total_cost_eur_per_cow_year = total_cost_eur_per_cow_year
  )


# ============================================================
# 8. Combine grouped and ungrouped rows and reshape wide
# ------------------------------------------------------------
# Final output is one row per measure, with separate columns
# for each region's annual cost per cow-year in EUR.
# ============================================================

final_export <- bind_rows(ungrouped_rows, grouped_rows) %>%
  group_by(row_id, order_id, measure, name, short_description, region) %>%
  summarise(
    total_cost_eur_per_cow_year = sum(total_cost_eur_per_cow_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = region,
    values_from = total_cost_eur_per_cow_year,
    values_fill = NA_real_
  ) %>%
  arrange(order_id) %>%
  select(measure, name, short_description, DK, IRL, CAN, US, IT)


# ============================================================
# 9. Export final region-specific cost table
# ============================================================

write.csv(
  final_export,
  "exports/costs_euro_pcy.csv",
  row.names = FALSE
)