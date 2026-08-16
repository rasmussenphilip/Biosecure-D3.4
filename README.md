<b>BIOSECURE PROJECT<br>DELIVERABLE 3.4<br>
PARTIAL BUDGETING-BASED COST-BENEFIT ANALYSIS OF BIOSECURITY MEASURES</b><br>

This repository contains the full R-based pipeline for estimating the net economic benefits of biosecurity adoption in dairy herds across bovine viral diarrhoea (BVD), paratuberculosis (PTB), and salmonellosis (SAL). It integrates published measure-effect estimates, the results of a stochastic dairy herd simulation model (SimHerd), and country-specific cost data to rank 16,384 biosecurity packages under three decision rules and identify Pareto-optimal combinations across herd disease statuses and countries.

Note: scripts 01–12 form a sequential pipeline and must be run in numerical order, each reading from the exports/ folder populated by the previous step. Scripts prefixed 99_ are independent utilities that can be run at any point after the main pipeline to generate specific summary tables or plots.

<b>Project Structure</b>

<b>Data Sources</b><br>
-SimHerd simulation output: gross margin (GM) per cow-year as a function of BVD, PTB, and SAL incidence, by country<br>
-Measure effect estimates (odds ratios) from the published literature, by disease and domain<br>
-Assumed baseline prevalence values, by disease and domain<br>
-Danish reference implementation costs, by measure<br>
-Country-level wage rates and GDP per capita (Canada, Denmark, Italy, United States)<br>

<b>Main Outputs</b><br>
-Cleaned, cost-standardized, and country-extrapolated measure and cost datasets<br>
-Package-level risk reductions, gross-margin changes, and net benefits across 16,384 packages, three herd disease statuses, 27 disease-risk states, and three weighting schemes<br>
-Decision-rule rankings (maximin, expected value, maximax) and Pareto-optimal package sets<br>
-Variance decomposition of net benefit by herd disease status, package, disease-risk state, country, and weighting scheme<br>
-Static and interactive visualizations of measure frequency, decision-rule agreement, and the Pareto frontier<br>

<b>FAIR Data Principles</b><br>
This repository follows the FAIR data principles to ensure that all resources are Findable, Accessible, Interoperable, and Reusable. The open-source R code is fully documented and version-controlled, and the harmonised datasets are structured to support reproducibility and integration with other tools and workflows. By aligning with FAIR standards, this project promotes transparency, encourages collaboration, and supports the broader research community in advancing economic decision analysis for livestock biosecurity.

<b>Pipeline Overview</b><br><br>
Data Preparation & Cost Standardization<br>
-Converts SimHerd gross-margin output from DKK to EUR<br>
-Extrapolates Danish reference costs to Canada, Denmark, Italy, and the United States using wage and GDP-per-capita scalars<br>
-Derives low, reference, and high cost scenarios for each measure and country<br><br>
Measure Effects & Risk Reduction<br>
-Inverts risk-factor-oriented odds ratios into protective estimates<br>
-Converts odds ratios to risk ratios and builds all 16,384 candidate packages<br>
-Aggregates package-level risk reductions under three weighting schemes (precision, effect, joint), with a bounded cross-domain interaction between herd-level and within-herd measures<br><br>
Economic Modeling<br>
-Fits country-specific risk-to-gross-margin response surfaces<br>
-Computes economic benefit by herd disease status (non-infected, infected, intermediate)<br>
-Subtracts annualized package cost to compute net benefit for every package, country, herd status, disease-risk state, and weighting scheme<br><br>
Decision Analysis<br>
-Ranks packages under the maximin, expected value, and maximax decision rules<br>
-Identifies the set of Pareto-optimal packages<br><br>
Sensitivity Analysis<br>
-Runs a one-factor variance decomposition of net benefit, both directly and on package-centred deviations<br><br>
Visualization & Summary Tables<br>
-Maps measure frequency across package size, herd disease status, and country using ggplot2<br>
-Builds decision-rule agreement diagrams (package- and measure-level overlap)<br>
-Produces interactive 3D Pareto-frontier and risk-surface plots<br>
-Compiles summary tables of winning and top Pareto-optimal packages by herd status and country<br>
-Checks the robustness of estimated net benefits (frequency of negative outcomes)<br>

<b>Dependencies</b>

Make sure the following R packages are installed:

library(dplyr) library(tidyverse) library(tidyr) library(readr) library(purrr) library(stringr) library(ggplot2) library(ggVennDiagram) library(paletteer) library(patchwork) library(plotly) library(htmlwidgets) library(writexl)

<b>Author</b>

This R code was developed and written by Philip Rasmussen for estimating the net economic benefits of biosecurity adoption in dairy cattle across BVD, paratuberculosis, and salmonellosis. Developed for the Biosecure project (EU Horizon 101083923).
