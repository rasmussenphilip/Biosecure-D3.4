# ============================================================
# SCRIPT 13: Sensitivity analysis of package net benefits
# ============================================================
#
# Purpose
# -------
# This script provides a lightweight post-hoc sensitivity analysis for the
# biosecurity package model. It does not rerun the full disease-risk,
# production-loss, or cost pipeline. Instead, it uses the already exported
# package-level net benefit estimates and asks a simpler question:
#
#   Which parts of the modelling setup explain most of the variation in
#   estimated net benefits?
#
# The response variable is net_benefit_ref, interpreted as the net economic
# benefit of each biosecurity package in EUR per cow-year under a given
# combination of region, herd/benefit type, aggregation method, and joint
# disease state.
#
# Conceptual model
# ----------------
# Each row in package_net_benefits.csv represents one package evaluated under
# one modelling context. Those contexts differ by:
#
#   1. package_id
#      The specific package of biosecurity measures being evaluated.
#
#   2. region
#      The regional economic setting used for gross-margin mapping and/or
#      price/cost assumptions.
#
#   3. benefit_type
#      The herd-type or benefit scenario, for example naive, average, or
#      infected herds, depending on the current version of the pipeline.
#
#   4. aggregation_method
#      The method used to combine individual measure effects within packages,
#      such as precision-based, effect-based, or hybrid aggregation.
#
#   5. joint_state
#      The combined disease-risk state across the diseases included in the
#      package model.
#
# Rather than changing each input one at a time, this script fits mixed-effects
# models and decomposes the observed variation in net benefits into variance
# components. These components indicate how much variation is associated with
# each modelling dimension.
#
# Statistical approach
# --------------------
# The script uses linear mixed-effects models of the form:
#
#   y_ijk... = intercept + random effects + residual error
#
# All explanatory dimensions are treated as random intercepts. This is useful
# here because the goal is not to estimate a fixed coefficient for a specific
# region or disease state. Instead, the goal is to quantify how much variation
# is attributable to each class of assumptions.
#
# Two related models are fitted:
#
#   Analysis A: Package-centered deviations
#   --------------------------------------
#   Response: net_benefit_deviation
#
#   Each package's own mean net benefit is subtracted before modelling:
#
#      net_benefit_deviation = net_benefit_ref - mean(net_benefit_ref within package)
#
#   This removes the fact that some packages are generally better or worse than
#   others. The model therefore asks:
#
#      After centering each package around its own average, which contextual
#      assumptions make that package perform above or below its typical value?
#
#   This is useful for understanding sensitivity around package performance,
#   independent of the overall quality of the package.
#
#   Analysis B: Overall net benefits including package
#   --------------------------------------------------
#   Response: net_benefit_ref
#
#   This model includes package_id as an additional random effect. It asks:
#
#      Across all results, how much variation is due to packages themselves
#      compared with region, benefit type, aggregation method, disease state,
#      and residual variation?
#
#   This is useful for checking whether package choice is the dominant driver,
#   or whether contextual assumptions dominate the results.
#
# Assumptions and interpretation
# ------------------------------
# 1. Random-effect interpretation
#    Each modelling dimension is interpreted as a source of variation. A larger
#    variance component means that differences across that dimension explain
#    more of the spread in net benefits.
#
# 2. Additive variance structure
#    The models assume that region, benefit type, aggregation method, joint
#    state, and package effects contribute additively to variation. Interactions
#    are not explicitly modelled here. This keeps the analysis interpretable and
#    lightweight, but means that interaction effects are absorbed into the
#    residual variance.
#
# 3. Post-hoc diagnostic use
#    These models are not intended to replace the main package-selection
#    framework. They are diagnostic tools for explaining where variation in the
#    already-generated outputs comes from.
#
# 4. No causal interpretation
#    The variance components do not imply causal effects. For example, a large
#    region component means that results vary substantially across regions, not
#    that region itself causally changes biological effectiveness.
#
# 5. Cost sensitivity
#    Costs are already embedded in net_benefit_ref. Therefore, any cost-driven
#    variation appears through package_id, region, benefit type, and the residual,
#    depending on how costs vary in the upstream pipeline. This script does not
#    independently perturb costs or herd size.
#
# Main output
# -----------
# The script exports:
#
#   exports/sensitivity/net_benefit_variance_components_comparison.csv
#
# This table contains variance components from both models and can be used to
# report whether net benefits are mainly driven by package choice, herd/benefit
# type, region, aggregation method, disease-state assumptions, or unexplained
# residual variation.
#
# ============================================================

library(dplyr)
library(readr)
library(lme4)

# ============================================================
# 1. Define input and output folders
# ============================================================
# The script reads package-level net benefits from the main exports folder and
# writes sensitivity outputs to a dedicated subfolder. This keeps diagnostic
# outputs separate from the main package-ranking and Pareto-frontier results.

export_dir <- "exports"
sensitivity_dir <- file.path(export_dir, "sensitivity")

if (!dir.exists(sensitivity_dir)) {
  dir.create(sensitivity_dir, recursive = TRUE)
}

# ============================================================
# 2. Helper function: extract variance components from lmer models
# ============================================================
# lme4 stores random-effect variances in a compact model object. This function
# converts them into a simple data frame with one row per model component.
#
# The returned variance is in squared units of the response variable. Since the
# response is EUR per cow-year, these variances are in squared EUR per cow-year.
# For interpretation, the relative size of the components is usually more useful
# than the absolute value.

extract_variance_components <- function(model, model_label) {
  as.data.frame(VarCorr(model)) %>%
    transmute(
      model = model_label,
      component = grp,
      variance = vcov
    )
}

# ============================================================
# 3. Load package-level net benefits
# ============================================================
# Each row should represent one evaluated package under one combination of
# assumptions. The key outcome is net_benefit_ref, which should already include
# both gross benefits and package costs.

df <- read_csv(
  file.path(export_dir, "package_net_benefits.csv"),
  show_col_types = FALSE
)

# ============================================================
# 4. Check that the expected columns are available
# ============================================================
# The model cannot be fitted unless the output from the previous scripts includes
# the package identifier, contextual variables, aggregation method, joint disease
# state, and net benefit outcome. This check fails early with a clear message if
# the upstream file has changed.

required_cols <- c(
  "package_id",
  "region",
  "benefit_type",
  "aggregation_method",
  "joint_state",
  "net_benefit_ref"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    "package_net_benefits.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# 5. Prepare the analysis data
# ============================================================
# Rows with missing net benefits are removed because lmer cannot use them.
# Categorical modelling dimensions are converted to factors so that lme4 treats
# them as grouping variables in the random-effects models.
#
# The package-centered response is also created here. For each package, the
# script calculates that package's mean net benefit across all contexts and then
# subtracts this mean from every package-specific result. This makes the first
# model focus on contextual deviations rather than overall package quality.

analysis_df <- df %>%
  filter(!is.na(net_benefit_ref)) %>%
  mutate(
    package_id = factor(package_id),
    region = factor(region),
    benefit_type = factor(benefit_type),
    aggregation_method = factor(aggregation_method),
    joint_state = factor(joint_state)
  ) %>%
  group_by(package_id) %>%
  mutate(
    package_mean_net_benefit = mean(net_benefit_ref, na.rm = TRUE),
    net_benefit_deviation = net_benefit_ref - package_mean_net_benefit
  ) %>%
  ungroup()

# ============================================================
# 6. Set model-fitting controls
# ============================================================
# Mixed-effects models can occasionally have convergence issues when there are
# many groups or unbalanced combinations of packages and scenarios. The bobyqa
# optimizer and a higher maximum number of function evaluations make fitting
# more stable without changing the statistical model.

lmer_control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# ============================================================
# 7. Analysis A: variance in package-centered deviations
# ============================================================
# Model question:
#   Once each package is centered around its own average net benefit, which
#   assumptions explain whether that package performs better or worse than its
#   typical value?
#
# Model structure:
#   net_benefit_deviation ~ 1
#     + random intercept for region
#     + random intercept for benefit_type
#     + random intercept for aggregation_method
#     + random intercept for joint_state
#
# Interpretation:
#   A large variance component for benefit_type would mean that package results
#   move substantially above or below their own average depending on herd/benefit
#   scenario. A large aggregation_method component would mean that the choice of
#   effect-aggregation method materially changes package performance. A large
#   residual component means that remaining row-level variation is not captured
#   by these additive grouping terms, possibly because of interactions among
#   package, region, disease state, and aggregation method.

vc_model_deviation <- lmer(
  net_benefit_deviation ~ 1 +
    (1 | region) +
    (1 | benefit_type) +
    (1 | aggregation_method) +
    (1 | joint_state),
  data = analysis_df,
  REML = TRUE,
  control = lmer_control
)

vc_deviation <- extract_variance_components(
  vc_model_deviation,
  "Package-centered deviations"
)

# ============================================================
# 8. Analysis B: variance in overall net benefits
# ============================================================
# Model question:
#   Across all package results, how much of the total variation in net benefit is
#   attributable to package identity versus modelling context?
#
# This model includes package_id as a random effect. The package variance
# component captures stable differences between packages: for example, whether
# some packages are consistently more profitable than others across regions,
# herd types, disease states, and aggregation methods.
#
# Interpretation:
#   - Large package_id variance:
#       Package choice is a major driver of net benefits.
#   - Large region variance:
#       Economic setting or regional GM/cost mapping strongly affects results.
#   - Large benefit_type variance:
#       Results differ substantially between herd/benefit scenarios.
#   - Large aggregation_method variance:
#       The method used to combine measure effects strongly affects outputs.
#   - Large joint_state variance:
#       Disease-risk configuration strongly affects net benefits.
#   - Large residual variance:
#       Much variation remains at the row level, likely reflecting combinations
#       or interactions not represented by this simple additive model.

vc_model_overall <- lmer(
  net_benefit_ref ~ 1 +
    (1 | package_id) +
    (1 | region) +
    (1 | benefit_type) +
    (1 | aggregation_method) +
    (1 | joint_state),
  data = analysis_df,
  REML = TRUE,
  control = lmer_control
)

vc_overall <- extract_variance_components(
  vc_model_overall,
  "Overall net benefits including package"
)

# ============================================================
# 9. Combine and export results
# ============================================================
# The two models answer related but distinct questions, so their variance
# components are stacked into one comparison table. This output can be inspected
# directly or used to create a bar plot of variance components.

vc_combined <- bind_rows(
  vc_deviation,
  vc_overall
)

write_csv(
  vc_combined,
  file.path(sensitivity_dir, "net_benefit_variance_components_comparison.csv")
)

# Print the table to the console so the main result is visible immediately after
# running the script in RStudio.

print(vc_combined)
