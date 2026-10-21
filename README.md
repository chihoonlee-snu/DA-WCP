# DA-WCP

This repository contains the R code used in **"Predicting Current Outcomes From Historical Survey Data With Weighted Conformal Prediction"** by Chihoon Lee, Sungkyu Jung, and Hyokyoung G. Hong.

The code implements design-adapted weighted conformal prediction (DA-WCP) and reproduces the simulation studies and NHANES application reported in the manuscript and its online supplementary material.

## Directory structure

```text
DA-WCP/
|-- Code/       R functions and analysis files
|-- Data/       ACS population counts used in the NHANES analysis
|-- Figures/    Figures, summary tables, and saved R results
`-- README.md
```

All file paths are relative to the main repository directory. The analyses should be run from this directory. The R files read the supplied data from `Data/` and save numerical results and figures in `Figures/`.

## Software requirements

The analyses require R and the following packages:

```r
install.packages(c(
  "mvtnorm", "dplyr", "quantreg", "ggplot2", "tidyr", "nnet",
  "xgboost", "ranger", "caret", "car", "pls", "haven", "readxl",
  "readr"
))
```

The NHANES analyses also require an internet connection because the public NHANES XPT files are obtained directly from the CDC website.

## Description of the R files

| File | Description |
|---|---|
| `Code/KLIEP.R` | Implements the Kullback–Leibler Importance Estimation Procedure (KLIEP) with Gaussian kernels and estimates subgroup-specific density ratios. |
| `Code/DA_WCP.R` | Defines the main `DA_WCP()` function, including weighted conformal cutoffs, continuous AR and CQR scores, categorical scores, and marginal and group-conditional prediction sets for DA-WCP and GWCP. |
| `Code/DA_WCP_simulation.R` | Conducts simulations under the five main scenarios (S1)–(S5) for marginal and group-conditional DA-WCP and GWCP. |
| `Code/DA_WCP_NHANES.R` | Prepares the NHANES 1999–2004 and 2015–March 2020 data, uses ACS counts to construct the target subgroup proportions, applies the continuous and categorical methods, and evaluates coverage. |
| `Code/DA_WCP_NHANES_overlap.R` | Produces subgroup-specific mirrored histograms for assessing overlap in BMI, HDL-C, and total cholesterol between the historical and current NHANES surveys. |
| `Code/DA_WCP_supp_setting_1.R` | Conducts simulations for Setting 1 of the additional simulations. |
| `Code/DA_WCP_supp_setting_2.R` | Conducts simulations for Setting 2 of the additional simulations. |
| `Code/DA_WCP_supp_H.R` | Contains two sensitivity analyses: a comparison of one and four strata, and the PPSWOR-within-strata approximation under varying degrees of PSU subgroup composition heterogeneity. |
| `Code/DA_WCP_supp_ICC.R` | Examines sensitivity to violations of Assumption 1 under varying levels of within-PSU dependence. |
| `Code/DA_WCP_supp_smallgroup.R` | Contains the basic small-subgroup analysis and the analyses based on truncation and regularization of the estimated density ratios. |
| `Code/DA_WCP_supp_nonresponse.R` | Examines violations of the nonresponse assumption (Assumption 2) by allowing the response probability at $t=2$ to depend on the covariates. |
| `Code/DA_WCP_supp_temporalshift.R` | Examines violations of temporal invariance in the conditional outcome distribution (Assumption 5). |

<!--
The section, figure, and table numbers below refer to the current revised manuscript and its online supplementary material.

| File | Location in the paper | Related results |
|---|---|---|
| `Code/KLIEP.R` | Main Section 4; Supplementary Section S1 | This file supplies the density-ratio estimates used by the other analyses and does not produce a figure by itself. |
| `Code/DA_WCP.R` | Main Sections 4–5; Supplementary Section S7 | This file supplies the principal estimation procedure and does not produce a figure by itself. |
| `Code/DA_WCP_simulation.R` | Main Section 6; detailed design in Supplementary Section S11 | Main Figures 2 and 3. |
| `Code/DA_WCP_NHANES.R` | Main Section 7; Supplementary Sections S14.1, S14.2, S14.4, and S14.5 | Main Figure 4, main Table 4, and Supplementary Figures S18–S19. |
| `Code/DA_WCP_NHANES_overlap.R` | Supplementary Section S14.3 | Supplementary Figures S15–S17. |
| `Code/DA_WCP_supp_setting_1.R` | Supplementary Section S12.1 | Panel (d) of Supplementary Figures S1, S2, and S4. |
| `Code/DA_WCP_supp_setting_2.R` | Supplementary Section S12.2 | Panel (b) of Supplementary Figures S6–S8. |
| `Code/DA_WCP_supp_H.R` | Supplementary Section S13.1 | Supplementary Table S3 and Figure S9. |
| `Code/DA_WCP_supp_ICC.R` | Supplementary Section S13.2 | Supplementary Figure S10. |
| `Code/DA_WCP_supp_smallgroup.R` | Supplementary Section S13.3 | Supplementary Tables S4–S5 and Figures S11–S12. |
| `Code/DA_WCP_supp_nonresponse.R` | Supplementary Section S13.4 | Supplementary Figure S13. |
| `Code/DA_WCP_supp_temporalshift.R` | Supplementary Section S13.5 | Supplementary Figure S14. |
-->

<!--
## Correspondence between output files and manuscript figures

Some figures in the manuscript combine two or more files produced by the R analyses. Their correspondence is given below.

| R output in `Figures/` | Location in the paper |
|---|---|
| `simulation_coverage_m.pdf` | Main Figure 2 |
| `simulation_coverage_g.pdf` | Main Figure 3 |
| `LDL_continuous_coverage_m.pdf` and `LDL_continuous_length_m.pdf` | Left and right panels of main Figure 4 |
| `simulation_coverage_m_supp_1_d.pdf` | Supplementary Figure S1, panel (d) |
| `simulation_length_m_supp_1_d.pdf` | Supplementary Figure S2, panel (d) |
| `simulation_coverage_g_supp_1_d.pdf` | Supplementary Figure S4, panel (d) |
| `simulation_coverage_m_supp_2_b.pdf` | Supplementary Figure S6, panel (b) |
| `simulation_length_m_supp_2_b.pdf` | Supplementary Figure S7, panel (b) |
| `simulation_coverage_g_supp_2_b.pdf` | Supplementary Figure S8, panel (b) |
| `simulation_H_PSU_coverage_m.pdf` | Supplementary Figure S9 |
| `simulation_ICC_coverage_m.pdf` | Supplementary Figure S10 |
| `simulation_smallgroup_2_balanced_truncation_coverage.pdf` and `simulation_smallgroup_2_balanced_regularization_coverage.pdf` | Supplementary Figure S11 |
| `simulation_smallgroup_2_small_truncation_coverage.pdf` and `simulation_smallgroup_2_small_regularization_coverage.pdf` | Supplementary Figure S12 |
| `simulation_nonresponse_shift_coverage.pdf` | Supplementary Figure S13 |
| `simulation_temporal_outcome_shift_coverage.pdf` | Supplementary Figure S14 |
| `DA_WCP_NHANES_overlap_bmi.pdf` | Supplementary Figure S15 |
| `DA_WCP_NHANES_overlap_hdl.pdf` | Supplementary Figure S16 |
| `DA_WCP_NHANES_overlap_tc.pdf` | Supplementary Figure S17 |
| `LDL_continuous_coverage_g.pdf` | Supplementary Figure S18 |
| `LDL_categorical_coverage_g.pdf` | Supplementary Figure S19 |
-->


## Data

The `Data/` directory contains:

- `2015-2016_acs.xlsx`
- `2017-March2020_acs.xlsx`

These files contain the ACS-based population control totals used to construct the target demographic subgroup probabilities for the NHANES analysis. The population totals were obtained from the [NCHS NHANES Response Rates and Population Totals webpage](https://wwwn.cdc.gov/nchs/nhanes/ResponseRates.aspx). The NHANES examination, laboratory, and demographic data files are downloaded directly from the CDC by `DA_WCP_NHANES.R`.



## Reproducing the results

Run each analysis in a new R session with the main repository directory as the working directory. Each of the following files begins its analysis when sourced; therefore, source one file at a time:

```r
source("Code/DA_WCP_simulation.R")
source("Code/DA_WCP_NHANES.R")
source("Code/DA_WCP_NHANES_overlap.R")
source("Code/DA_WCP_supp_setting_1.R")
source("Code/DA_WCP_supp_setting_2.R")
```

The other supplementary files define the functions for each sensitivity analysis. After sourcing the relevant file, call the corresponding function as follows:

```r
source("Code/DA_WCP_supp_H.R")
H_results <- DA_WCP_simul_H(M = 100)

source("Code/DA_WCP_supp_H.R")
H_PSU_results <- DA_WCP_simul_H_PSU(M = 100)

source("Code/DA_WCP_supp_ICC.R")
ICC_results <- DA_WCP_simul_ICC(M = 100)

source("Code/DA_WCP_supp_smallgroup.R")
smallgroup_results <- DA_WCP_simul_smallgroup_1(M = 100)

source("Code/DA_WCP_supp_smallgroup.R")
smallgroup_stabilization_results <- DA_WCP_simul_smallgroup_2(M = 100)

source("Code/DA_WCP_supp_nonresponse.R")
nonresponse_results <- DA_WCP_simul_nonresponse(M = 100)

source("Code/DA_WCP_supp_temporalshift.R")
temporal_shift_results <- DA_WCP_simul_temporalshift(M = 100)
```

The results reported in the paper are based on 100 repetitions. A short test run may be conducted with `M = 1` and without saving the output, for example:

```r
source("Code/DA_WCP_supp_ICC.R")
ICC_test <- DA_WCP_simul_ICC(M = 1, save_outputs = FALSE)
```
