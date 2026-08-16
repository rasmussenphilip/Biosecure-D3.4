library(dplyr)
library(readr)
library(stringr)
library(plotly)
library(htmlwidgets)

# ============================================================
# 15_pareto_plots.R
# 3D Pareto frontier plots, excluding all-regions herd-type plots
# ============================================================

out_dir <- "exports"
plot_dir <- file.path(out_dir, "pareto_3d_plots")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

pareto_all_packages <- read_csv(
  file.path(out_dir, "pareto_all_packages_for_plotting.csv"),
  show_col_types = FALSE
)

safe_file_label <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_") %>%
    str_remove("_$")
}

make_3d_plot <- function(df, title, file_name) {
  
  bg_df <- df %>% filter(!on_frontier)
  fr_df <- df %>% filter(on_frontier)
  
  hover_text <- function(data) {
    paste(
      "Package ID:", data$package_id,
      "<br>Package size:", data$package_size,
      "<br>", data$package_contents,
      
      "<br><br><b>Net benefit summary</b>",
      "<br>Worst-case:", round(data$worst_case_payoff, 2),
      "<br>Mean:", round(data$mean_payoff, 2),
      "<br>Best-case:", round(data$best_case_payoff, 2),
      "<br>5th percentile:", round(data$p05_payoff, 2),
      "<br>10th percentile:", round(data$p10_payoff, 2),
      "<br>Median:", round(data$median_payoff, 2),
      "<br>90th percentile:", round(data$p90_payoff, 2),
      "<br>95th percentile:", round(data$p95_payoff, 2),
      
      "<br><br><b>Ranks</b>",
      "<br>Rank maximin:", data$rank_maximin,
      "<br>Rank expected value:", data$rank_expected_value,
      "<br>Rank maximax:", data$rank_maximax,
      "<br>Number of payoffs:", data$n_payoffs
    )
  }
  
  p <- plot_ly() %>%
    add_trace(
      data = bg_df,
      x = ~worst_case_payoff,
      y = ~mean_payoff,
      z = ~best_case_payoff,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, color = "steelblue", opacity = 0.25),
      text = hover_text(bg_df),
      hoverinfo = "text",
      name = "Dominated packages"
    ) %>%
    add_trace(
      data = fr_df,
      x = ~worst_case_payoff,
      y = ~mean_payoff,
      z = ~best_case_payoff,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 6, color = "red", opacity = 1),
      text = paste("PARETO FRONTIER<br>", hover_text(fr_df)),
      hoverinfo = "text",
      name = "Pareto-optimal packages"
    ) %>%
    layout(
      title = title,
      width = 1400,
      height = 1400,
      font = list(
        family = "Arial",
        size = 12
      ),
      scene = list(
        xaxis = list(title = "Worst-case net benefit (maximin)"),
        yaxis = list(title = "Mean net benefit (expected value)"),
        zaxis = list(title = "Best-case net benefit (maximax)")
      )
    )
  
  htmlwidgets::saveWidget(
    as_widget(p),
    file.path(plot_dir, file_name),
    selfcontained = TRUE
  )
}

# ============================================================
# Create one 3D Pareto plot per context,
# excluding all-regions Naive/Average/Infected plots
# ============================================================

plot_contexts <- pareto_all_packages %>%
  distinct(
    analysis_level,
    context_region,
    context_benefit_type
  ) %>%
  filter(
    !(
      context_region == "All regions" &
        str_detect(
          str_to_lower(context_benefit_type),
          "naive|naïve|average|infected"
        )
    )
  )

for (i in seq_len(nrow(plot_contexts))) {
  
  ctx <- plot_contexts[i, ]
  
  df_plot <- pareto_all_packages %>%
    filter(
      analysis_level == ctx$analysis_level,
      context_region == ctx$context_region,
      context_benefit_type == ctx$context_benefit_type
    )
  
  plot_title <- paste(
    "Pareto frontier:",
    ctx$analysis_level,
    "|",
    ctx$context_region,
    "|",
    ctx$context_benefit_type
  )
  
  plot_file <- paste0(
    "pareto_3d_",
    safe_file_label(ctx$analysis_level),
    "_",
    safe_file_label(ctx$context_region),
    "_",
    safe_file_label(ctx$context_benefit_type),
    ".html"
  )
  
  make_3d_plot(
    df = df_plot,
    title = plot_title,
    file_name = plot_file
  )
}