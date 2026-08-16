# ============================================================
# Script 06: Estimate package-level changes in gross margin
# ============================================================
#
# Purpose
# -------
# This script converts package-level disease-risk reductions into
# economic benefits, expressed as changes in gross margin per cow-year
# (delta_gmpcy). It combines two upstream outputs:
#
#   1. exports/risk_to_gmpcy.csv
#      - The response surface linking BVD, PTB, and Salmonella risk to
#        gross margin per cow-year for each region.
#
#   2. exports/package_risk_reductions.csv
#      - The post-package risk values generated for each candidate
#        biosecurity package, aggregation method, and joint disease state.
#
# The script fits a region-specific quadratic linear model to approximate
# the SimHerd-derived relationship between disease risk and gross margin.
# It then predicts the economic benefit of each package under three herd
# interpretations:
#
#   - infected_herd: benefits from reducing within-herd disease risk once
#     infection is already present.
#
#   - average_herd: benefits from reducing the combined herd-level and
#     within-herd disease burden, represented as their product.
#
#   - naive_herd: benefits from reducing introduction risk into a disease-free
#     herd; outbreak severity is held at the infected-herd baseline so that
#     the economic value reflects prevention of introduction rather than
#     within-herd mitigation after infection.
#
# Main output
# -----------
# The script writes:
#
#   exports/package_gm_values.csv
#
# This file contains one row per package, region, herd interpretation,
# aggregation method, and joint disease state, with the associated
# delta_gmpcy payoff. It is the main economic-benefit input for the
# downstream cost, net-benefit, and game-theoretic package-ranking scripts.
#
# Important modelling notes
# -------------------------
# - The quadratic response surface is fitted separately for each region.
# - The model is fit without an intercept because the baseline gross margin
#   is handled separately using the disease-free regional baseline.
# - Package benefits are calculated relative to the empty-package baseline
#   within the same state suffix, so comparisons are state-specific.
# - The script does not change package costs; it estimates benefits only.
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

# ============================================================
# Step 1: Load the regional risk-to-GM mapping data
# ------------------------------------------------------------
# This file contains the simulated/derived relationship between disease
# risk and gross margin per cow-year. The columns c_gmpcy and gmpcy are
# both retained because c_gmpcy is used for the response-surface model,
# while gmpcy at zero risk provides the disease-free regional baseline.
# ============================================================

risk_gmpcy <- read_csv(
  "exports/risk_to_gmpcy.csv",
  show_col_types = FALSE
)

# Check that the mapping file has all columns needed for the regional
# response-surface models and the later baseline merge.
required_gmpcy_cols <- c(
  "region",
  "risk_bvd",
  "risk_ptb",
  "risk_sal",
  "gmpcy",
  "c_gmpcy"
)

missing_gmpcy_cols <- setdiff(required_gmpcy_cols, names(risk_gmpcy))

if (length(missing_gmpcy_cols) > 0) {
  stop(
    "risk_to_gmpcy.csv is missing required columns: ",
    paste(missing_gmpcy_cols, collapse = ", ")
  )
}

# ============================================================
# Step 2: Fit one quadratic response-surface model per region
# ------------------------------------------------------------
# The model approximates how combined BVD, PTB, and Salmonella risk affects
# c_gmpcy. Squared terms allow non-linear single-disease effects, while
# interaction terms allow the marginal cost of one disease to depend on the
# risk level of another disease.
#
# The model has no intercept (-1) because the disease-free baseline is added
# back separately after prediction. This keeps predicted disease effects
# conceptually separate from the regional baseline GM level.
# ============================================================

model_data <- split(risk_gmpcy, risk_gmpcy$region)

models <- lapply(model_data, function(df) {
  lm(
    c_gmpcy ~
      risk_bvd + risk_ptb + risk_sal +
      I(risk_bvd^2) + I(risk_ptb^2) + I(risk_sal^2) +
      risk_bvd:risk_ptb + risk_bvd:risk_sal + risk_ptb:risk_sal - 1,
    data = df
  )
})

# Extract the disease-free GM level for each region. This is later added
# back to the predicted disease-related component to recover absolute GM.
region_baselines <- risk_gmpcy %>%
  filter(risk_bvd == 0, risk_ptb == 0, risk_sal == 0) %>%
  select(region, disease_free_gmpcy = gmpcy)

if (nrow(region_baselines) == 0) {
  stop("No disease-free baseline rows found in risk_to_gmpcy.csv.")
}

# ============================================================
# Step 3: Load the package-level risk reductions
# ------------------------------------------------------------
# This compact wide file has one row per package and many risk columns.
# The risk columns encode herd interpretation, disease, aggregation method,
# and joint state in their names. The next section reshapes and parses
# those encoded column names into explicit variables.
# ============================================================

risk_wide <- read_csv(
  "exports/package_risk_reductions.csv",
  show_col_types = FALSE
)

meta_cols <- c("package_id", "package_size", "package_contents")
missing_meta_cols <- setdiff(meta_cols, names(risk_wide))

if (length(missing_meta_cols) > 0) {
  stop(
    "package_risk_reductions.csv is missing required columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

risk_value_cols <- setdiff(names(risk_wide), meta_cols)

if (length(risk_value_cols) == 0) {
  stop("package_risk_reductions.csv contains no risk columns.")
}

# ============================================================
# Step 4: Reshape package risk values from wide to long and back
# ------------------------------------------------------------
# The original column names follow the pattern:
#
#   post_<naive|infected>_risk_<bvd|ptb|sal>_<suffix>
#
# The suffix identifies the aggregation method and joint disease state.
# Parsing this structure creates one clean row per package-state
# combination, with separate columns for naive and infected risks for each
# disease.
# ============================================================

risk_long <- risk_wide %>%
  pivot_longer(
    cols = all_of(risk_value_cols),
    names_to = "risk_col",
    values_to = "risk_value"
  ) %>%
  mutate(
    risk_col = str_to_lower(risk_col)
  )

# Parse risk-column names into three components: herd risk type, disease,
# and suffix. Columns that do not match the expected pattern are dropped.
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

package_states <- risk_long %>%
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
  ) %>%
  arrange(package_id, suffix)

# Ensure each package-state row has all six disease-risk components needed
# for the infected, average, and naive herd calculations.
required_package_state_cols <- c(
  "naive_bvd", "naive_ptb", "naive_sal",
  "infected_bvd", "infected_ptb", "infected_sal"
)

missing_package_state_cols <- setdiff(
  required_package_state_cols,
  names(package_states)
)

if (length(missing_package_state_cols) > 0) {
  stop(
    "Missing disease-specific risk columns after reshaping: ",
    paste(missing_package_state_cols, collapse = ", ")
  )
}

if (any(is.na(package_states$naive_bvd)) ||
    any(is.na(package_states$naive_ptb)) ||
    any(is.na(package_states$naive_sal)) ||
    any(is.na(package_states$infected_bvd)) ||
    any(is.na(package_states$infected_ptb)) ||
    any(is.na(package_states$infected_sal))) {
  stop("Missing naive or infected risk values found after reshaping.")
}

# ============================================================
# Step 5: Recover aggregation method and joint disease state
# ------------------------------------------------------------
# The suffix is split into:
#
#   - aggregation_method: precision_based, effect_based, or joint_weighted
#   - joint_state: the combined BVD/PTB/Salmonella disease-state label
#
# Average-herd risk is then calculated as naive risk multiplied by infected
# risk. In this framework, naive risk represents herd-level/introduction
# probability and infected risk represents within-herd severity conditional
# on infection.
# ============================================================

package_states <- package_states %>%
  mutate(
    aggregation_method = case_when(
      str_starts(suffix, "precision_based_") ~ "precision_based",
      str_starts(suffix, "effect_based_") ~ "effect_based",
      str_starts(suffix, "joint_weighted_") ~ "joint_weighted",
      TRUE ~ NA_character_
    ),
    joint_state = case_when(
      aggregation_method == "precision_based" ~ str_remove(suffix, "^precision_based_"),
      aggregation_method == "effect_based" ~ str_remove(suffix, "^effect_based_"),
      aggregation_method == "joint_weighted" ~ str_remove(suffix, "^joint_weighted_"),
      TRUE ~ NA_character_
    ),
    average_bvd = naive_bvd * infected_bvd,
    average_ptb = naive_ptb * infected_ptb,
    average_sal = naive_sal * infected_sal
  )

if (any(is.na(package_states$aggregation_method))) {
  stop("Some suffix values could not be matched to an aggregation method.")
}

if (any(is.na(package_states$joint_state))) {
  stop("Some suffix values could not be converted to joint_state.")
}

# ============================================================
# Step 6: Add empty-package baseline risks for state-matched comparisons
# ------------------------------------------------------------
# Benefits are not measured against zero disease risk. They are measured
# against the empty_package within the same aggregation-method/state suffix.
# This keeps each package comparison aligned with the same baseline disease
# scenario.
# ============================================================

empty_package_states <- package_states %>%
  filter(package_contents == "empty_package") %>%
  select(
    suffix,
    baseline_naive_bvd = naive_bvd,
    baseline_naive_ptb = naive_ptb,
    baseline_naive_sal = naive_sal,
    baseline_infected_bvd = infected_bvd,
    baseline_infected_ptb = infected_ptb,
    baseline_infected_sal = infected_sal,
    baseline_average_bvd = average_bvd,
    baseline_average_ptb = average_ptb,
    baseline_average_sal = average_sal
  )

if (nrow(empty_package_states) == 0) {
  stop("No empty_package row found. Cannot create baseline comparisons.")
}

package_states <- package_states %>%
  left_join(
    empty_package_states,
    by = "suffix"
  )

# ============================================================
# Step 7: Define prediction helper functions
# ------------------------------------------------------------
# predict_c_gmpcy() applies a region-specific response-surface model to a
# set of BVD, PTB, and Salmonella risks.
#
# expected_c_gmpcy_naive() handles disease-free herds by averaging over all
# possible introduction combinations. For each combination, it calculates:
#
#   probability of that disease combination
#     x predicted GM loss if those introduced diseases occur
#
# and sums over all eight possible BVD/PTB/SAL presence-absence states.
# ============================================================

predict_c_gmpcy <- function(model, risk_bvd, risk_ptb, risk_sal) {
  as.numeric(
    predict(
      model,
      newdata = tibble(
        risk_bvd = risk_bvd,
        risk_ptb = risk_ptb,
        risk_sal = risk_sal
      )
    )
  )
}

expected_c_gmpcy_naive <- function(
    model,
    p_bvd,
    p_ptb,
    p_sal,
    severity_bvd,
    severity_ptb,
    severity_sal
) {
  
  n <- length(p_bvd)
  expected_value <- rep(0, n)
  
    # Enumerate the 2^3 possible introduction outcomes across the three
  # diseases: absent/present for BVD, PTB, and Salmonella.
  disease_combinations <- crossing(
    bvd_present = c(0, 1),
    ptb_present = c(0, 1),
    sal_present = c(0, 1)
  )
  
    # For each disease combination, calculate its probability, assign the
  # corresponding severity risks, predict c_gmpcy, and add the probability-
  # weighted value to the expected outcome.
  for (i in seq_len(nrow(disease_combinations))) {
    
    combo <- disease_combinations[i, ]
    
    prob_bvd <- if (combo$bvd_present == 1) p_bvd else 1 - p_bvd
    prob_ptb <- if (combo$ptb_present == 1) p_ptb else 1 - p_ptb
    prob_sal <- if (combo$sal_present == 1) p_sal else 1 - p_sal
    
    combo_probability <- prob_bvd * prob_ptb * prob_sal
    
    risk_bvd <- if (combo$bvd_present == 1) severity_bvd else 0
    risk_ptb <- if (combo$ptb_present == 1) severity_ptb else 0
    risk_sal <- if (combo$sal_present == 1) severity_sal else 0
    
    combo_c_gmpcy <- predict_c_gmpcy(
      model = model,
      risk_bvd = risk_bvd,
      risk_ptb = risk_ptb,
      risk_sal = risk_sal
    )
    
    expected_value <- expected_value + combo_probability * combo_c_gmpcy
  }
  
  expected_value
}

# ============================================================
# Step 8: Predict GM benefits by region and herd interpretation
# ------------------------------------------------------------
# For every region, the script predicts package and baseline c_gmpcy under
# three interpretations.
#
# infected_herd:
#   Package-adjusted infected risks are used directly. This captures the
#   value of lowering within-herd disease burden in herds where infection is
#   already present.
#
# average_herd:
#   Uses package-adjusted herd-level/introduction risk multiplied by
#   package-adjusted infected risk. This represents the expected burden in
#   an average herd across both infection occurrence and within-herd severity.
#
# naive_herd:
#   Uses expected c_gmpcy over possible introductions. Package-adjusted naive
#   risks change the probability of introduction, while outbreak severity is
#   fixed at the empty-package infected baseline. This makes naive-herd
#   benefits reflect prevention of introduction.
# ============================================================

regions <- names(models)

gm_results <- map_dfr(regions, function(region_name) {
  
  model <- models[[region_name]]
  df <- package_states
  
    # Infected-herd prediction: use the package-adjusted within-herd risks.
  infected_pred <- predict_c_gmpcy(
    model = model,
    risk_bvd = df$infected_bvd,
    risk_ptb = df$infected_ptb,
    risk_sal = df$infected_sal
  )
  
  infected_baseline <- predict_c_gmpcy(
    model = model,
    risk_bvd = df$baseline_infected_bvd,
    risk_ptb = df$baseline_infected_ptb,
    risk_sal = df$baseline_infected_sal
  )
  
    # Average-herd prediction: combine occurrence/introduction and severity.
  average_pred <- predict_c_gmpcy(
    model = model,
    risk_bvd = df$average_bvd,
    risk_ptb = df$average_ptb,
    risk_sal = df$average_sal
  )
  
  average_baseline <- predict_c_gmpcy(
    model = model,
    risk_bvd = df$baseline_average_bvd,
    risk_ptb = df$baseline_average_ptb,
    risk_sal = df$baseline_average_sal
  )
  
    # Naive-herd prediction: calculate expected GM loss over all possible
  # introduction combinations, holding severity at the infected baseline.
  naive_pred <- expected_c_gmpcy_naive(
    model = model,
    p_bvd = df$naive_bvd,
    p_ptb = df$naive_ptb,
    p_sal = df$naive_sal,
    severity_bvd = df$baseline_infected_bvd,
    severity_ptb = df$baseline_infected_ptb,
    severity_sal = df$baseline_infected_sal
  )
  
  naive_baseline <- expected_c_gmpcy_naive(
    model = model,
    p_bvd = df$baseline_naive_bvd,
    p_ptb = df$baseline_naive_ptb,
    p_sal = df$baseline_naive_sal,
    severity_bvd = df$baseline_infected_bvd,
    severity_ptb = df$baseline_infected_ptb,
    severity_sal = df$baseline_infected_sal
  )
  
  bind_rows(
    df %>%
      mutate(
        region = region_name,
        benefit_type = "infected_herd",
        risk_bvd = infected_bvd,
        risk_ptb = infected_ptb,
        risk_sal = infected_sal,
        baseline_risk_bvd = baseline_infected_bvd,
        baseline_risk_ptb = baseline_infected_ptb,
        baseline_risk_sal = baseline_infected_sal,
        predicted_c_gmpcy = infected_pred,
        baseline_c_gmpcy = infected_baseline
      ),
    df %>%
      mutate(
        region = region_name,
        benefit_type = "average_herd",
        risk_bvd = average_bvd,
        risk_ptb = average_ptb,
        risk_sal = average_sal,
        baseline_risk_bvd = baseline_average_bvd,
        baseline_risk_ptb = baseline_average_ptb,
        baseline_risk_sal = baseline_average_sal,
        predicted_c_gmpcy = average_pred,
        baseline_c_gmpcy = average_baseline
      ),
    df %>%
      mutate(
        region = region_name,
        benefit_type = "naive_herd",
        risk_bvd = naive_bvd * baseline_infected_bvd,
        risk_ptb = naive_ptb * baseline_infected_ptb,
        risk_sal = naive_sal * baseline_infected_sal,
        baseline_risk_bvd = baseline_naive_bvd * baseline_infected_bvd,
        baseline_risk_ptb = baseline_naive_ptb * baseline_infected_ptb,
        baseline_risk_sal = baseline_naive_sal * baseline_infected_sal,
        predicted_c_gmpcy = naive_pred,
        baseline_c_gmpcy = naive_baseline
      )
  )
})

# ============================================================
# Step 9: Convert predicted disease components into delta GM values
# ------------------------------------------------------------
# The response-surface models predict c_gmpcy, the disease-related component
# of GM. The disease-free regional baseline is added back to obtain absolute
# predicted GM. The final delta_gmpcy is the difference between each package
# and the state-matched empty-package baseline.
#
# payoff_state concatenates benefit_type and suffix so downstream scripts
# can treat each region/herd/state/aggregation combination as a separate
# payoff scenario.
# ============================================================

package_gm_values <- gm_results %>%
  left_join(
    region_baselines,
    by = "region"
  ) %>%
  mutate(
    delta_c_gmpcy = predicted_c_gmpcy - baseline_c_gmpcy,
    predicted_gmpcy = disease_free_gmpcy + predicted_c_gmpcy,
    baseline_package_gmpcy = disease_free_gmpcy + baseline_c_gmpcy,
    delta_gmpcy = predicted_gmpcy - baseline_package_gmpcy,
    payoff_state = paste(
      benefit_type,
      suffix,
      sep = "_"
    )
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
    delta_gmpcy
  ) %>%
  arrange(
    region,
    package_id,
    benefit_type,
    aggregation_method,
    joint_state
  )

# ============================================================
# Step 10: Export package-level GM benefits
# ------------------------------------------------------------
# This output is intentionally narrow: it contains the economic benefit
# estimate only. Costs and net benefits are added downstream.
# ============================================================

out_dir <- "exports"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(
  package_gm_values,
  file.path(out_dir, "package_gm_values.csv")
)