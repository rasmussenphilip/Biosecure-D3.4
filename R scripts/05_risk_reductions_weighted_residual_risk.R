# ============================================================
# Script 05_risk_reductions_weighted_residual_risk.R
# ------------------------------------------------------------
# Package effects are combined multiplicatively on residual risk.
# Precision-, effect-, and joint-based interpretations determine
# how much of each literature-derived RR contributes to the
# package's uncertain true effectiveness.
# ============================================================

library(tidyverse)

# ============================================================
# Load inputs
# ============================================================

measures <- read_csv("exports/measures_inverted.csv", show_col_types = FALSE)
risk     <- read_csv("inputs/risk_final.csv", show_col_types = FALSE)

required_measures_cols <- c("measure", "domain", "disease", "mean", "lower", "upper")
required_risk_cols <- c(
  "disease",
  "within_low", "within_mid", "within_high",
  "herd_low", "herd_mid", "herd_high"
)

missing_measures <- setdiff(required_measures_cols, names(measures))
missing_risk <- setdiff(required_risk_cols, names(risk))

if (length(missing_measures) > 0) {
  stop("measures_inverted.csv is missing required columns: ",
       paste(missing_measures, collapse = ", "))
}

if (length(missing_risk) > 0) {
  stop("risk_final.csv is missing required columns: ",
       paste(missing_risk, collapse = ", "))
}

# ============================================================
# Helper functions
# ============================================================

clamp01 <- function(x) pmin(pmax(x, 0), 1)

or_to_rr <- function(or, p0) {
  clamp01(or / ((1 - p0) + (p0 * or)))
}

calc_se <- function(lower, upper) {
  bad <- is.na(lower) | is.na(upper) | lower <= 0 | upper <= 0 | upper < lower
  out <- rep(NA_real_, length(lower))
  ok <- !bad
  out[ok] <- (log(upper[ok]) - log(lower[ok])) / (2 * 1.96)
  out
}

scale_by_max_safe <- function(x) {
  # Convert positive raw weights to bounded attenuation factors
  # in [0, 1], preserving relative differences without forcing
  # the smallest positive weight to zero.
  out <- rep(0, length(x))
  ok <- is.finite(x) & x > 0

  if (!any(ok)) return(out)

  xmax <- max(x[ok], na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) return(out)

  out[ok] <- x[ok] / xmax
  clamp01(out)
}

generate_package_matrix <- function(measure_names, include_empty = TRUE) {
  n <- length(measure_names)
  ids <- if (include_empty) 0:(2^n - 1) else 1:(2^n - 1)
  
  mat <- sapply(seq_len(n), function(j) {
    as.integer(bitwAnd(ids, bitwShiftL(1L, j - 1)) != 0)
  })
  
  if (n == 1) {
    mat <- matrix(mat, ncol = 1)
  }
  
  colnames(mat) <- measure_names
  
  package_contents <- apply(mat, 1, function(x) {
    selected <- measure_names[as.logical(x)]
    if (length(selected) == 0) {
      "empty_package"
    } else {
      paste(selected, collapse = " + ")
    }
  })
  
  tibble(
    package_id = seq_along(ids),
    package_index = ids,
    package_size = rowSums(mat),
    package_contents = package_contents
  ) %>%
    bind_cols(as_tibble(mat))
}

aggregate_rr_all_packages <- function(package_matrix_subset, rr_vec, alpha_vec) {
  # Weighted multiplicative residual-risk aggregation:
  #
  #   RR_package = exp(sum_i alpha_i * log(RR_i))
  #              = product_i RR_i ^ alpha_i
  #
  # alpha_i is a bounded [0, 1] interpretation weight.
  # Because all included RRs are protective (0 < RR <= 1),
  # adding another positively weighted protective measure can
  # only maintain or reduce the package RR.

  n_packages <- nrow(package_matrix_subset)

  if (
    is.null(package_matrix_subset) ||
    length(rr_vec) == 0 ||
    ncol(package_matrix_subset) == 0
  ) {
    return(rep(1, n_packages))
  }

  rr_vec <- clamp01(rr_vec)
  alpha_vec <- clamp01(alpha_vec)

  valid <- is.finite(rr_vec) &
    rr_vec > 0 &
    rr_vec <= 1 &
    is.finite(alpha_vec) &
    alpha_vec > 0

  if (!any(valid)) {
    return(rep(1, n_packages))
  }

  X <- package_matrix_subset[, valid, drop = FALSE]
  rr_use <- rr_vec[valid]
  alpha_use <- alpha_vec[valid]

  weighted_log_rr <- as.vector(
    X %*% (alpha_use * log(rr_use))
  )

  out <- exp(weighted_log_rr)
  clamp01(out)
}

make_aggregation_blocks <- function(
    packages,
    package_matrix,
    idx,
    rr_mean,
    rr_lower,
    rr_upper
) {
  
  n_packages <- nrow(packages)
  
  if (length(idx) == 0) {
    return(tibble(
      package_id = packages$package_id,
      precision_based = rep(1, n_packages),
      effect_based = rep(1, n_packages),
      joint_weighted = rep(1, n_packages)
    ))
  }
  
  se <- calc_se(rr_lower, rr_upper)

  # Raw interpretation weights.
  w_precision <- ifelse(is.na(se) | se <= 0, 0, 1 / (se^2))
  w_effect    <- ifelse(is.na(rr_mean) | rr_mean <= 0, 0, abs(log(rr_mean)))

  # Convert raw weights to bounded attenuation factors in [0, 1].
  # Scaling is across all measures available for the current
  # disease-domain-state, not separately within each package.
  alpha_precision <- scale_by_max_safe(w_precision)
  alpha_effect    <- scale_by_max_safe(w_effect)
  alpha_joint     <- sqrt(alpha_precision * alpha_effect)

  X <- package_matrix[, idx, drop = FALSE]

  tibble(
    package_id = packages$package_id,
    precision_based = aggregate_rr_all_packages(X, rr_mean, alpha_precision),
    effect_based    = aggregate_rr_all_packages(X, rr_mean, alpha_effect),
    joint_weighted  = aggregate_rr_all_packages(X, rr_mean, alpha_joint)
  )
}

# ============================================================
# Clean inputs
# ============================================================

measures <- measures %>%
  mutate(
    measure = as.character(measure),
    disease = as.character(disease) %>% str_to_lower() %>% str_trim(),
    domain  = as.character(domain) %>% str_trim()
  ) %>%
  filter(domain %in% c("herd_level", "within_herd"))

risk <- risk %>%
  mutate(
    disease = as.character(disease) %>% str_to_lower() %>% str_trim()
  )

dup_check <- measures %>%
  count(measure, disease, domain) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  warning("Duplicate measure-disease-domain rows found. Keeping first row in each group.")
  
  measures <- measures %>%
    group_by(measure, disease, domain) %>%
    slice(1) %>%
    ungroup()
}

# ============================================================
# Create low/mid/high disease risk states
# ============================================================

risk_states <- risk %>%
  transmute(
    disease,
    low_within  = clamp01(within_low),
    mid_within  = clamp01(within_mid),
    high_within = clamp01(within_high),
    low_herd    = clamp01(herd_low),
    mid_herd    = clamp01(herd_mid),
    high_herd   = clamp01(herd_high)
  ) %>%
  pivot_longer(
    cols = -disease,
    names_to = c("state", ".value"),
    names_pattern = "(low|mid|high)_(within|herd)"
  ) %>%
  rename(
    baseline_within = within,
    baseline_herd   = herd
  ) %>%
  mutate(
    state = factor(state, levels = c("low", "mid", "high"))
  ) %>%
  arrange(disease, state)

expected_disease_states <- risk_states %>%
  count(disease)

if (any(expected_disease_states$n != 3)) {
  stop("Each disease must have exactly 3 paired states: low, mid, and high.")
}

# ============================================================
# Generate packages
# ============================================================

measure_names <- sort(unique(measures$measure))

packages <- generate_package_matrix(
  measure_names = measure_names,
  include_empty = TRUE
)

package_matrix <- as.matrix(packages[, measure_names, drop = FALSE])
storage.mode(package_matrix) <- "numeric"

n_packages <- nrow(package_matrix)

# ============================================================
# Create disease-domain lookup
# ============================================================

lookup <- list()

for (d in sort(unique(risk_states$disease))) {
  for (dom in c("herd_level", "within_herd")) {
    
    sub <- measures %>%
      filter(disease == d, domain == dom) %>%
      select(measure, mean, lower, upper)
    
    if (nrow(sub) == 0) {
      lookup[[paste(d, dom, sep = "::")]] <- NULL
    } else {
      lookup[[paste(d, dom, sep = "::")]] <- list(
        idx   = match(sub$measure, measure_names),
        mean  = sub$mean,
        lower = sub$lower,
        upper = sub$upper
      )
    }
  }
}

# ============================================================
# Estimate package-level risk by disease state
# ------------------------------------------------------------
# Within each disease and domain:
#
#   RR_package,m = product_i RR_i ^ alpha_i,m
#
# where m is the precision-, effect-, or joint-based
# interpretation of uncertain true effectiveness.
# ============================================================

disease_state_blocks <- vector("list", nrow(risk_states))

for (s in seq_len(nrow(risk_states))) {
  
  sc <- risk_states[s, ]
  
  disease_value   <- sc$disease
  state_value     <- as.character(sc$state)
  baseline_herd   <- sc$baseline_herd
  baseline_within <- sc$baseline_within
  
  herd_info <- lookup[[paste(disease_value, "herd_level", sep = "::")]]
  
  if (is.null(herd_info)) {
    herd_rr <- tibble(
      package_id = packages$package_id,
      precision_based = rep(1, n_packages),
      effect_based = rep(1, n_packages),
      joint_weighted = rep(1, n_packages)
    )
  } else {
    herd_rr <- make_aggregation_blocks(
      packages = packages,
      package_matrix = package_matrix,
      idx = herd_info$idx,
      rr_mean  = or_to_rr(herd_info$mean,  baseline_herd),
      rr_lower = or_to_rr(herd_info$lower, baseline_herd),
      rr_upper = or_to_rr(herd_info$upper, baseline_herd)
    )
  }
  
  within_info <- lookup[[paste(disease_value, "within_herd", sep = "::")]]
  
  if (is.null(within_info)) {
    within_rr <- tibble(
      package_id = packages$package_id,
      precision_based = rep(1, n_packages),
      effect_based = rep(1, n_packages),
      joint_weighted = rep(1, n_packages)
    )
  } else {
    within_rr <- make_aggregation_blocks(
      packages = packages,
      package_matrix = package_matrix,
      idx = within_info$idx,
      rr_mean  = or_to_rr(within_info$mean,  baseline_within),
      rr_lower = or_to_rr(within_info$lower, baseline_within),
      rr_upper = or_to_rr(within_info$upper, baseline_within)
    )
  }
  
  block <- packages %>%
    select(package_id, package_size, package_contents) %>%
    left_join(
      herd_rr %>%
        pivot_longer(
          cols = c(precision_based, effect_based, joint_weighted),
          names_to = "aggregation_method",
          values_to = "aggregated_rr_herd"
        ),
      by = "package_id"
    ) %>%
    left_join(
      within_rr %>%
        pivot_longer(
          cols = c(precision_based, effect_based, joint_weighted),
          names_to = "aggregation_method",
          values_to = "aggregated_rr_within"
        ),
      by = c("package_id", "aggregation_method")
    ) %>%
    mutate(
      disease = disease_value,
      disease_state = state_value,
      
      naive_secondary_rr =
        pmax(aggregated_rr_herd, aggregated_rr_within),
      
      infected_secondary_rr =
        pmax(aggregated_rr_within, aggregated_rr_herd),
      
      post_naive_risk = clamp01(
        baseline_herd *
          aggregated_rr_herd *
          naive_secondary_rr
      ),
      
      post_infected_risk = clamp01(
        baseline_within *
          aggregated_rr_within *
          infected_secondary_rr
      ),
      
      post_average_risk = clamp01(
        post_naive_risk * post_infected_risk
      )
    )
  
  disease_state_blocks[[s]] <- block
}

disease_state_results <- bind_rows(disease_state_blocks)

# ============================================================
# Build joint disease states
# ============================================================

required_diseases <- c("bvd", "ptb", "sal")
missing_required_diseases <- setdiff(required_diseases, unique(disease_state_results$disease))

if (length(missing_required_diseases) > 0) {
  stop(
    "Script 06 expects disease labels: ",
    paste(required_diseases, collapse = ", "),
    ". Missing from inputs: ",
    paste(missing_required_diseases, collapse = ", ")
  )
}

joint_states <- crossing(
  bvd_state = c("low", "mid", "high"),
  ptb_state = c("low", "mid", "high"),
  sal_state = c("low", "mid", "high")
) %>%
  mutate(
    joint_state = paste0(
      "bvd_", bvd_state,
      "__ptb_", ptb_state,
      "__sal_", sal_state
    )
  ) %>%
  arrange(bvd_state, ptb_state, sal_state)

# ============================================================
# Assemble joint-state long table
# ============================================================

aggregation_methods <- c(
  "precision_based",
  "effect_based",
  "joint_weighted"
)

base_joint_index <- crossing(
  packages %>% select(package_id, package_size, package_contents),
  aggregation_method = aggregation_methods,
  joint_states
)

joint_results_long <- bind_rows(
  base_joint_index %>%
    transmute(
      package_id,
      package_size,
      package_contents,
      aggregation_method,
      joint_state,
      bvd_state,
      ptb_state,
      sal_state,
      disease = "bvd",
      disease_state = bvd_state
    ),
  base_joint_index %>%
    transmute(
      package_id,
      package_size,
      package_contents,
      aggregation_method,
      joint_state,
      bvd_state,
      ptb_state,
      sal_state,
      disease = "ptb",
      disease_state = ptb_state
    ),
  base_joint_index %>%
    transmute(
      package_id,
      package_size,
      package_contents,
      aggregation_method,
      joint_state,
      bvd_state,
      ptb_state,
      sal_state,
      disease = "sal",
      disease_state = sal_state
    )
) %>%
  left_join(
    disease_state_results,
    by = c(
      "package_id",
      "package_size",
      "package_contents",
      "aggregation_method",
      "disease",
      "disease_state"
    )
  ) %>%
  arrange(package_id, aggregation_method, joint_state, disease)

if (nrow(joint_results_long) == 0) {
  stop("No rows produced in joint_results_long.")
}

if (any(is.na(joint_results_long$post_naive_risk))) {
  stop("Missing post_naive_risk values found after assembling joint states.")
}

if (any(is.na(joint_results_long$post_infected_risk))) {
  stop("Missing post_infected_risk values found after assembling joint states.")
}

if (any(is.na(joint_results_long$post_average_risk))) {
  stop("Missing post_average_risk values found after assembling joint states.")
}

# ============================================================
# Create single package-level wide export
# ------------------------------------------------------------
# One row per package.
#
# Naive risk:
#   herd-level risk after herd_level measures
#
# Infected risk:
#   within-herd risk after within_herd measures
#
# Average risk:
#   combined animal-level risk after both herd_level and
#   within_herd measures
# ============================================================

results_long <- joint_results_long %>%
  mutate(
    suffix = paste(aggregation_method, joint_state, sep = "_"),
    result_name = paste(disease, suffix, sep = "_")
  )

package_risk_reductions <- results_long %>%
  select(
    package_id,
    package_size,
    package_contents,
    result_name,
    post_naive_risk,
    post_infected_risk,
    post_average_risk
  ) %>%
  pivot_wider(
    names_from = result_name,
    values_from = c(
      post_naive_risk,
      post_infected_risk,
      post_average_risk
    )
  ) %>%
  arrange(package_id)

# ============================================================
# Export only one file
# ============================================================

out_dir <- "exports"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(
  package_risk_reductions,
  file.path(out_dir, "package_risk_reductions.csv")
)