if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# NHANES covariate-overlap diagnostics.
# Each helper prepares or validates subgroup-specific data before drawing mirrored histograms.

library(mvtnorm)
library(dplyr)
library(quantreg)
library(ggplot2)
library(tidyr)
library(nnet)
library(xgboost)
library(ranger)
library(xgboost)
library(caret)
library(car)
library(pls)
library(haven)

source("KLIEP.R")
source("DA_WCP.R")


#####################################################################################
#
#                 Collecting NHANES data (1999-2004, 2015-2020/03)
#
#####################################################################################


yearletters = c("", "_B", "_C", "_D", "_E", "_F", "_G", "_H", "_I", "") # 99-00, 01-02, 03-04, 05-06, 07-08, 09-10, 11-12, 13-14, 15-16, 17-20/03
years = c("1999", "2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017")

yearsvec = c(1, 2, 3, 9, 10)

# variable file names(varcodes), variable names(varlist)
# suggested predictor variables : age(RIDAGEYR), sex(RIAGENDR), race(RIDRETH1), income(INDFMPIR), education level(DMDEDUC2), physical activity(PAQ665) 2,
# fiber intake(DR1TFIBE + DR2TFIBE), total saturated fatty acid(DR1TSFAT + DR2TSFAT), total dietary cholesterol(DR1TCHOL + DR2TCHOL),
# bmi(BMXBMI), waist circumference(BMXWAIST), blood pressure(BPXSY1), diabetes(DIQ010)

varcodes99 = c("DEMO", "PAQ", "DRXTOT", "BMX", "BPX", "DIQ", "LAB10AM", "LAB13AM", "LAB13")
varlist99 = list(c("WTINT2YR", "WTMEC2YR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "INDFMPIR", "DMDEDUC2"), c("PAD320"), c("DRXTFIBE", "DRXTSFAT", "DRXTCHOL"),
                 c("BMXBMI", "BMXWAIST"), c("BPXSY1"), c("DIQ010"), c("LBXGLU"), c("WTSAF2YR", "LBDLDL", "LBXTR"), c("LBXTC", "LBDHDL"))

varcodes01 = c("DEMO", "PAQ", "DRXTOT", "BMX", "BPX", "DIQ", "L10AM", "L13AM", "L13")
varlist01 = list(c("WTINT2YR", "WTMEC2YR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "INDFMPIR", "DMDEDUC2"), c("PAD320"), c("DRXTFIBE", "DRXTSFAT", "DRXTCHOL"),
                 c("BMXBMI", "BMXWAIST"), c("BPXSY1"), c("DIQ010"), c("LBXGLU"), c("WTSAF2YR", "LBDLDL", "LBXTR"), c("LBXTC", "LBDHDL"))

varcodes03 = c("DEMO", "PAQ", "DR1TOT", "DR2TOT", "BMX", "BPX", "DIQ", "L10AM", "L13AM", "L13")
varlist03 = list(c("WTINT2YR", "WTMEC2YR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "INDFMPIR", "DMDEDUC2"), c("PAD320"), c("DR1TFIBE", "DR1TSFAT", "DR1TCHOL"),
                 c("DR2TFIBE", "DR2TSFAT", "DR2TCHOL"), c("BMXBMI", "BMXWAIST"), c("BPXSY1"), c("DIQ010"), c("LBXGLU"), c("WTSAF2YR", "LBDLDL", "LBXTR"), c("LBXTC", "LBXHDD"))

varcodes15 = c("DEMO", "PAQ", "DR1TOT", "DR2TOT", "BMX", "BPX", "DIQ", "GLU", "TRIGLY", "TCHOL", "HDL")
varlist15 = list(c("WTINT2YR", "WTMEC2YR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "INDFMPIR", "DMDEDUC2"), c("PAQ665"), c("DR1TFIBE", "DR1TSFAT", "DR1TCHOL"),
                 c("DR2TFIBE", "DR2TSFAT", "DR2TCHOL"), c("BMXBMI", "BMXWAIST"), c("BPXSY1"), c("DIQ010"), c("LBXGLU"), c("WTSAF2YR", "LBDLDL", "LBXTR"), c("LBXTC"), c("LBDHDD"))

varcodes17 = c("P_DEMO", "P_PAQ", "P_DR1TOT", "P_DR2TOT", "P_BMX", "P_BPXO", "P_DIQ", "P_GLU", "P_TRIGLY", "P_TCHOL", "P_HDL")
varlist17 = list(c("WTINTPRP", "WTMECPRP", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "INDFMPIR", "DMDEDUC2"), c("PAQ665"), c("DR1TFIBE", "DR1TSFAT", "DR1TCHOL"),
                 c("DR2TFIBE", "DR2TSFAT", "DR2TCHOL"), c("BMXBMI", "BMXWAIST"), c("BPXOSY1"), c("DIQ010"), c("LBXGLU"), c("WTSAFPRP", "LBDLDL", "LBXTR"), c("LBXTC"), c("LBDHDD"))


varnames_old = c("id", "wt_z", "wt_x", "sex", "age", "race", "income", "edu", "paq", "fibe1", "sfat1", "chol1",
                 "bmi", "waist", "bpx", "diq", "glu", "wt_y", "y", "tg", "tc", "hdl") # 1999-2000, 2001-2002
varnames_new = c("id", "wt_z", "wt_x", "sex", "age", "race", "income", "edu", "paq", "fibe1", "sfat1", "chol1",
                 "fibe2", "sfat2", "chol2", "bmi", "waist", "bpx", "diq", "glu", "wt_y", "y", "tg", "tc", "hdl") # 2003-2004, 2015-2016, 2017-2020
varnames = c("id", "wt_z", "wt_x", "sex", "age", "race", "income", "edu", "paq", "fibe", "sfat", "chol",
             "bmi", "waist", "bpx", "diq", "glu", "wt_y", "y", "tg", "tc", "hdl")

varcodes = list(varcodes99, varcodes01, varcodes03, varcodes15, varcodes17)
varlist = list(varlist99, varlist01, varlist03, varlist15, varlist17)
nvec = vector(length = length(varcodes))
data = c()

for(i in 1:length(varcodes)){
  varcodes_i = varcodes[[i]]
  varlist_i = varlist[[i]]
  year_i = yearsvec[i]
  cat(years[year_i], "wave data extracting...\n")
  for(j in 1:length(varcodes_i)){
    data_i_new = read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s%s.XPT",
                                      years[year_i], varcodes_i[j], yearletters[year_i]))) %>% dplyr::select(c("SEQN", varlist_i[[j]]))
    if(j == 1) data_i = data_i_new
    if(j != 1) data_i = left_join(data_i, data_i_new, by = "SEQN")
  }
  if(i > 2) colnames(data_i) = varnames_new
  nvec[i] = nrow(data_i)
  data = rbind(data, data_i)

  if(i == 2){
    colnames(data) = varnames_old
    paqnum = (1:length(varnames_old))[varnames_old == "paq"]
    data = cbind(data[, 1:paqnum], data[, c("fibe1", "sfat1", "chol1")], data[, (paqnum+1):ncol(data)])
    colnames(data) = varnames_new
  }
}

n1 = sum(nvec[1:3])
n2 = sum(nvec[4:5])

data = data %>% # edu, fibe, sfat, chol
  mutate(edu = as.factor(case_when(
    edu <= 3 ~ 1,
    edu %in% c(4, 5) ~ 2)),
    fibe = (fibe1 + fibe2)/2,
    sfat = (sfat1 + sfat2)/2,
    chol = (chol1 + chol2)/2
  )

data = data %>% # diq : 2 groups
  mutate(diq = as.factor(case_when(
    diq %in% c(1, 3) ~ 1,
    diq == 2 ~ 2
  )))

data = data %>% # sex : 2 groups
  mutate(sex = as.factor(case_when(
    sex == 1 ~ 1,
    sex == 2 ~ 2)))

data = data %>% # paq : 2 groups
  mutate(paq = as.factor(case_when(
    paq == 1 ~ 1,
    paq == 2 ~ 2)))

data = data %>% # age : 4 groups
  mutate(age = as.factor(cut(age, breaks = c(19, 34, 49, 64, Inf), labels = 1:4)))

data = data %>% # race : 3 groups
  mutate(race = as.factor(case_when(
    race %in% c(1, 2) ~ 1, # Hispanic
    race == 3 ~ 2, # NH White
    race == 4 ~ 3))) # NH Black

data = data[, varnames]

data$wt_x = data$wt_x / 3
data$wt_y = data$wt_y / 3
data$wt_z = data$wt_z / 3

data <- data %>% mutate(ycat = as.factor(case_when(
  y < 100 ~ '0',
  y >= 100 & y < 130 ~ '1',
  y >= 130 & y < 160 ~ '2',
  y >= 160 & y < 190 ~ '3',
  y >= 190 ~ '4'
)))

data1 = data[1:n1, ]
data2 = data[(n1+1):(n1+n2), ]


#####################################################################################
#
#              Subgroup-specific overlap of continuous covariates
#
#####################################################################################

var_x <- c("tc", "bmi", "hdl")
var_x_labels <- c(tc = "Total cholesterol", bmi = "BMI", hdl = "HDL cholesterol")
group_variables <- c("sex", "age", "race")
analysis_covariates <- c(var_x, group_variables)

select_complete_analysis_covariates <- function(dataset) {
  complete_rows <- complete.cases(dataset[, analysis_covariates, drop = FALSE])
  finite_continuous_rows <- apply(
    dataset[, var_x, drop = FALSE],
    1,
    function(row) all(is.finite(row))
  )

  dataset[complete_rows & finite_continuous_rows, , drop = FALSE]
}

data1_complete_covariates <- select_complete_analysis_covariates(data1)
data2_complete_covariates <- select_complete_analysis_covariates(data2)

cat(
  "Complete analysis-covariate cases (response not required):\n",
  "t = 1: ", nrow(data1_complete_covariates), " / ", nrow(data1), "\n",
  "t = 2: ", nrow(data2_complete_covariates), " / ", nrow(data2), "\n",
  sep = ""
)

full_subgroup_grid <- expand_grid(
  race = levels(data$race),
  age = levels(data$age),
  sex = levels(data$sex)
) %>%
  mutate(
    subgroup_key = interaction(sex, age, race, drop = FALSE, sep = "."),
    subgroup = sprintf("Sex %s | Age %s | Race %s", sex, age, race)
  )

subgroup_levels <- full_subgroup_grid$subgroup
subgroup_lookup <- setNames(full_subgroup_grid$subgroup, full_subgroup_grid$subgroup_key)

make_overlap_data <- function(data1, data2, covariate) {
  bind_rows(
    data1 %>% mutate(dataset = "data1"),
    data2 %>% mutate(dataset = "data2")
  ) %>%
    transmute(
      dataset = factor(dataset, levels = c("data1", "data2")),
      sex = sex,
      age = age,
      race = race,
      value = .data[[covariate]],
      subgroup_key = as.character(interaction(sex, age, race, drop = FALSE, sep = "."))
    ) %>%
    filter(
      !is.na(sex), !is.na(age), !is.na(race),
      is.finite(value), subgroup_key %in% names(subgroup_lookup)
    ) %>%
    mutate(
      subgroup = factor(unname(subgroup_lookup[subgroup_key]), levels = subgroup_levels)
    )
}

validate_overlap_data <- function(plot_data, covariate) {
  cell_counts <- plot_data %>%
    count(subgroup, dataset, .drop = FALSE)

  invalid_cells <- cell_counts %>% filter(n < 2L)
  if (nrow(invalid_cells) > 0L) {
    stop(
      "Histogram requires at least two observations in every subgroup/time cell for ",
      covariate, ". Invalid cells: ",
      paste(
        paste(invalid_cells$subgroup, invalid_cells$dataset, sep = "/"),
        collapse = ", "
      )
    )
  }

  if (nlevels(plot_data$subgroup) != 24L) {
    stop("Expected 24 subgroups, but found ", nlevels(plot_data$subgroup), ".")
  }

  invisible(cell_counts)
}

make_mirrored_histogram_data <- function(plot_data, bins = 30L) {
  value_range <- range(plot_data$value, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) <= 0) {
    stop("Cannot construct histogram bins from the observed covariate range.")
  }

  breaks <- seq(value_range[1], value_range[2], length.out = bins + 1L)

  histogram_data <- plot_data %>%
    group_by(subgroup, dataset) %>%
    group_modify(~ {
      histogram <- hist(
        .x$value,
        breaks = breaks,
        plot = FALSE,
        include.lowest = TRUE,
        right = FALSE
      )

      data.frame(
        midpoint = histogram$mids,
        density = histogram$density
      )
    }) %>%
    ungroup() %>%
    mutate(
      signed_density = if_else(dataset == "t = 1", -density, density)
    )

  bin_width <- diff(breaks)[1]
  symmetric_limits <- histogram_data %>%
    group_by(subgroup) %>%
    summarise(limit = max(density) * 1.05, .groups = "drop") %>%
    crossing(sign = c(-1, 1)) %>%
    mutate(midpoint = mean(value_range), signed_density = sign * limit)

  list(
    histogram = histogram_data,
    symmetric_limits = symmetric_limits,
    bin_width = bin_width
  )
}

make_overlap_plot <- function(data1, data2, covariate, bins = 30L) {
  plot_data <- make_overlap_data(data1, data2, covariate)
  validate_overlap_data(plot_data, covariate)
  plot_data$dataset <- factor(
    plot_data$dataset,
    levels = c("data1", "data2"),
    labels = c("t = 1", "t = 2")
  )
  histogram_result <- make_mirrored_histogram_data(plot_data, bins = bins)

  ggplot(
    histogram_result$histogram,
    aes(x = midpoint, y = signed_density, fill = dataset)
  ) +
    geom_col(
      width = histogram_result$bin_width * 0.92,
      color = "white",
      linewidth = 0.12,
      alpha = 0.90
    ) +
    geom_blank(
      data = histogram_result$symmetric_limits,
      aes(x = midpoint, y = signed_density),
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey45",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    facet_wrap(
      ~ subgroup,
      ncol = 6,
      nrow = 4,
      scales = "free_y",
      drop = FALSE
    ) +
    scale_fill_manual(values = c("t = 1" = "#7fb6d8", "t = 2" = "#e98b8b")) +
    labs(
      title = NULL,
      x = covariate,
      y = "Density",
      fill = NULL
    ) +
    theme_test(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey94", color = "grey75"),
      strip.text = element_text(size = 8.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 11),
      axis.title.x = element_text(size = 16.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 16.5),
      legend.key.width = grid::unit(18, "pt"),
      legend.key.spacing.x = grid::unit(18, "pt"),
      plot.margin = margin(6, 6, 4, 6)
    )
}

overlap_plots <- setNames(lapply(var_x, function(covariate) {
  make_overlap_plot(
    data1_complete_covariates,
    data2_complete_covariates,
    covariate
  )
}), var_x)

output_pdfs <- setNames(
  file.path("../Figures", paste0("DA_WCP_NHANES_overlap_", var_x, ".pdf")),
  var_x
)

for (covariate in var_x) {
  grDevices::pdf(output_pdfs[[covariate]], width = 16, height = 10)
  print(overlap_plots[[covariate]])
  grDevices::dev.off()
}

cat(
  "NHANES subgroup mirrored-histogram figures saved to:\n",
  paste(output_pdfs, collapse = "\n"),
  "\n",
  sep = ""
)
