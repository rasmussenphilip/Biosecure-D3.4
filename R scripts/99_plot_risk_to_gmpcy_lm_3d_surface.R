library(dplyr)
library(plotly)

# Load prepared data ####
# Values should now be in EUR per cow-year from 01_risk_to_gmpcy_mapping.R
risk_gmpcy <- read.csv("exports/risk_to_gmpcy.csv")

# Get baseline GM per region ####
baseline_df <- risk_gmpcy %>%
  filter(risk_bvd == 0, risk_ptb == 0, risk_sal == 0) %>%
  select(region, gmpcy)

baseline_lookup <- setNames(baseline_df$gmpcy, baseline_df$region)

# Refit LM models ####
models <- risk_gmpcy %>%
  group_split(region) %>%
  setNames(unique(risk_gmpcy$region)) %>%
  lapply(function(df) {
    lm(
      c_gmpcy ~ 
        risk_bvd + risk_ptb + risk_sal +
        I(risk_bvd^2) + I(risk_ptb^2) + I(risk_sal^2) +
        risk_bvd:risk_ptb + risk_bvd:risk_sal + risk_ptb:risk_sal - 1,
      data = df
    )
  })

# Prediction: change in GM relative to baseline, EUR per cow-year ####
predict_c_gmpcy <- function(region, r_bvd, r_ptb, r_sal) {
  model <- models[[region]]
  
  predict(model, newdata = data.frame(
    risk_bvd = r_bvd,
    risk_ptb = r_ptb,
    risk_sal = r_sal
  ))
}

# Prediction: gross margin per cow-year, EUR ####
predict_gmpcy <- function(region, r_bvd, r_ptb, r_sal) {
  baseline <- baseline_lookup[[region]]
  baseline + predict_c_gmpcy(region, r_bvd, r_ptb, r_sal)
}

# Plot 3D surfaces ####
plot_gm_surface <- function(region,
                            fixed_var = c("risk_sal", "risk_ptb", "risk_bvd"),
                            fixed_value = 0.10,
                            n = 50,
                            plot_type = c("loss", "gmpcy", "c_gmpcy", "p_gmpcy")) {
  
  fixed_var <- match.arg(fixed_var)
  plot_type <- match.arg(plot_type)
  vals <- seq(0, 0.5, length.out = n)
  vals_pct <- vals * 100
  
  if (fixed_var == "risk_sal") {
    grid <- expand.grid(
      risk_bvd = vals,
      risk_ptb = vals,
      risk_sal = fixed_value
    )
    xtitle <- "BVD risk (%)"
    ytitle <- "PTB risk (%)"
    subtitle <- paste0("SAL fixed at ", fixed_value * 100, "%")
  }
  
  if (fixed_var == "risk_ptb") {
    grid <- expand.grid(
      risk_bvd = vals,
      risk_ptb = fixed_value,
      risk_sal = vals
    )
    xtitle <- "BVD risk (%)"
    ytitle <- "SAL risk (%)"
    subtitle <- paste0("PTB fixed at ", fixed_value * 100, "%")
  }
  
  if (fixed_var == "risk_bvd") {
    grid <- expand.grid(
      risk_bvd = fixed_value,
      risk_ptb = vals,
      risk_sal = vals
    )
    xtitle <- "PTB risk (%)"
    ytitle <- "SAL risk (%)"
    subtitle <- paste0("BVD fixed at ", fixed_value * 100, "%")
  }
  
  grid$c_gmpcy <- predict_c_gmpcy(
    region,
    grid$risk_bvd,
    grid$risk_ptb,
    grid$risk_sal
  )
  
  grid$loss_gmpcy <- -grid$c_gmpcy
  
  grid$gmpcy_pred <- predict_gmpcy(
    region,
    grid$risk_bvd,
    grid$risk_ptb,
    grid$risk_sal
  )
  
  grid$p_gmpcy_pred <- 100 * grid$c_gmpcy / baseline_lookup[[region]]
  
  if (plot_type == "loss") {
    z_mat <- matrix(grid$loss_gmpcy, nrow = n, ncol = n)
    plot_title <- paste(region, ": Economic loss surface")
    z_title <- "Loss (EUR per cow-year)"
  } else if (plot_type == "gmpcy") {
    z_mat <- matrix(grid$gmpcy_pred, nrow = n, ncol = n)
    plot_title <- paste(region, ": Gross margin surface")
    z_title <- "GM (EUR per cow-year)"
  } else if (plot_type == "c_gmpcy") {
    z_mat <- matrix(grid$c_gmpcy, nrow = n, ncol = n)
    plot_title <- paste(region, ": Gross margin change surface")
    z_title <- "ΔGM (EUR per cow-year)"
  } else if (plot_type == "p_gmpcy") {
    z_mat <- matrix(grid$p_gmpcy_pred, nrow = n, ncol = n)
    plot_title <- paste(region, ": Percent change in gross margin surface")
    z_title <- "%ΔGM"
  }
  
  plot_ly(
    x = vals_pct,
    y = vals_pct,
    z = z_mat,
    type = "surface"
  ) %>%
    layout(
      title = list(
        text = paste0(plot_title, "<br><sup>", subtitle, "</sup>"),
        x = 0.5
      ),
      scene = list(
        xaxis = list(title = xtitle),
        yaxis = list(title = ytitle),
        zaxis = list(title = z_title)
      )
    )
}

# Example plots ####

# Loss surfaces
plot_gm_surface("DK", "risk_sal", 0.1, plot_type = "loss")
plot_gm_surface("DK", "risk_sal", 0.3, plot_type = "loss")
plot_gm_surface("US", "risk_bvd", 0.2, plot_type = "loss")

# Gross margin surfaces
plot_gm_surface("DK", "risk_sal", 0.1, plot_type = "gmpcy")
plot_gm_surface("DK", "risk_sal", 0.3, plot_type = "gmpcy")
plot_gm_surface("US", "risk_bvd", 0.2, plot_type = "gmpcy")

# Change in gross margin surfaces
plot_gm_surface("DK", "risk_sal", 0.1, plot_type = "c_gmpcy")
plot_gm_surface("DK", "risk_sal", 0.3, plot_type = "c_gmpcy")
plot_gm_surface("US", "risk_bvd", 0.2, plot_type = "c_gmpcy")

# Percent change in gross margin surfaces
plot_gm_surface("US", "risk_ptb", 0.3, plot_type = "p_gmpcy")
plot_gm_surface("DK", "risk_sal", 0.3, plot_type = "p_gmpcy")
plot_gm_surface("US", "risk_bvd", 0.2, plot_type = "p_gmpcy")
