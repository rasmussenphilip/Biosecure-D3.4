Biosecure Project – Deliverable 3.4 
PARTIAL BUDGETING-BASED COST-BENEFIT ANALYSIS OF BIOSECURITY MEASURES

This repository contains the full R-based pipeline for estimating the net economic benefits of biosecurity adoption in dairy herds across bovine viral diarrhoea (BVD), paratuberculosis (PTB), and salmonellosis (SAL). It integrates published measure-effect estimates, the results of a stochastic dairy herd simulation model (SimHerd), and country-specific cost data to rank 16,384 biosecurity packages under three decision rules and identify Pareto-optimal combinations across herd disease statuses and countries.

Note: scripts 01–12 form a sequential pipeline and must be run in numerical order, each reading from the exports/ folder populated by the previous step. Scripts prefixed 99_ are independent utilities that can be run at any point after the main pipeline to generate specific summary tables or plots.

Project Structure

Data Sources:
SimHerd simulation output: gross margin (GM) per cow-year as a function of BVD, PTB, and SAL incidence, by country
Measure effect estimates (odds ratios) from the published literature, by disease and domain
Assumed baseline prevalence values, by disease and domain
Danish reference implementation costs, by measure
Country-level wage rates and GDP per capita (Canada, Denmark, Italy, United States)
Main Outputs:
Cleaned, cost-standardized, and country-extrapolated measure and cost datasets
Package-level risk reductions, gross-margin changes, and net benefits across 16,384 packages, three herd disease statuses, 27 disease-risk states, and three weighting schemes
Decision-rule rankings (maximin, expected value, maximax) and Pareto-optimal package sets
Variance decomposition of net benefit by herd disease status, package, disease-risk state, country, and weighting scheme
Static and interactive visualizations of measure frequency, decision-rule agreement, and the Pareto frontier

FAIR Data Principles

This repository follows the FAIR data principles to ensure that all resources are Findable, Accessible, Interoperable, and Reusable. The open-source R code is fully documented and version-controlled, and the harmonised datasets are structured to support reproducibility and integration with other tools and workflows. By aligning with FAIR standards, this project promotes transparency, encourages collaboration, and supports the broader research community in advancing economic decision analysis for livestock biosecurity.

Pipeline Overview

Data Preparation & Cost Standardization
Converts SimHerd gross-margin output from DKK to EUR
Extrapolates Danish reference costs to Canada, Denmark, Italy, and the United States using wage and GDP-per-capita scalars
Derives low, reference, and high cost scenarios for each measure and country
Measure Effects & Risk Reduction
Inverts risk-factor-oriented odds ratios into protective estimates
Converts odds ratios to risk ratios and builds all 16,384 candidate packages
Aggregates package-level risk reductions under three weighting schemes (precision, effect, joint), with a bounded cross-domain interaction between herd-level and within-herd measures
Economic Modeling
Fits country-specific risk-to-gross-margin response surfaces
Computes economic benefit by herd disease status (non-infected, infected, intermediate)
Subtracts annualized package cost to compute net benefit for every package, country, herd status, disease-risk state, and weighting scheme
Decision Analysis
Ranks packages under the maximin, expected value, and maximax decision rules
Identifies the set of Pareto-optimal packages
Sensitivity Analysis
Runs a one-factor variance decomposition of net benefit, both directly and on package-centred deviations
Visualization & Summary Tables
Maps measure frequency across package size, herd disease status, and country using ggplot2
Builds decision-rule agreement diagrams (package- and measure-level overlap)
Produces interactive 3D Pareto-frontier and risk-surface plots
Compiles summary tables of winning and top Pareto-optimal packages by herd status and country
Checks the robustness of estimated net benefits (frequency of negative outcomes)

Dependencies

Make sure the following R packages are installed:

cowplot, dplyr, ggforce, ggplot2, grid, janitor, paletteer, patchwork, readr, readxl, stringr, tibble, and tidyr. 

Author

This R code was developed and written by Philip Rasmussen, for estimating the net economic benefits of biosecurity adoption in dairy cattle across BVD, paratuberculosis, and salmonellosis. Developed for the Biosecure project (EU Horizon 101083923).
