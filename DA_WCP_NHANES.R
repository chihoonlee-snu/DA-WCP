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

setwd("C:/Users/lchk0/OneDrive/바탕 화면/Research/DA-WCP/Code")
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
#                    Baseline characteristics of variables
#
#####################################################################################

# Population (Sum of sampling weights)

data1 %>%
  group_by(age) %>%
  summarise(weight_sum = sum(wt_z)) %>%
  mutate(prop = weight_sum / sum(weight_sum))

data2 %>%
  group_by(age) %>%
  summarise(weight_sum = sum(wt_z)) %>%
  mutate(prop = weight_sum / sum(weight_sum))

# sex : 0.488, 0.512 / 0.489, 0.511
# age : 0.217, 0.230, 0.150, 0.115, 0.288 / 0.206, 0.188, 0.195, 0.154, 0.256
# race : 0.148, 0.678, 0.120, 0.054 / 0.180, 0.600, 0.120, 0.101

# NHANES data

data1 %>%
  count(age) %>%
  mutate(prop = n / sum(n))

data2 %>%
  count(age) %>%
  mutate(prop = n / sum(n))

# sex : 0.488, 0.512 / 0.494, 0.506
# age : 0.133, 0.119, 0.101, 0.140, 0.507 / 0.141, 0.139, 0.158, 0.148, 0.414
# race : 0.326, 0.389, 0.243, 0.042 / 0.265, 0.327, 0.244, 0.165



### Continuous variables

var_conti = c("tc", "hdl", "bmi", "bpx")
varnames_y = c("y")

data1 %>%
  summarise(across(all_of(varnames_y),
                   list(mean  = ~ mean(.x, na.rm = TRUE),
                        wmean = ~ weighted.mean(.x, wt_y, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

data2 %>%
  summarise(across(all_of(varnames_y),
                   list(mean  = ~ mean(.x, na.rm = TRUE),
                        wmean = ~ weighted.mean(.x, wt_y, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))


#####################################################################################
#
#                        DA-WCP and GWCP to NHANES dataset
#
#####################################################################################

# Variables & Parameters
M = 100
beta = 0.6
var_x = c("tc", "bmi", "hdl")
var_y = c("y")
var_y_cat = c("ycat")
var_z = c()
var_g = c("sex", "age", "race")
var_h = c()
var_wt_x = "wt_x"
var_wt_y = "wt_y"

# Data preprocessing
d1_cov = data1[rowSums(is.na(data1[, c(var_x, var_z, var_g, var_h)])) == 0, 
               c(var_wt_x, var_wt_y, var_x, var_y, var_y_cat, var_z, var_g, var_h)]
d2_cov = data2[rowSums(is.na(data2[, c(var_x, var_z, var_g, var_h)])) == 0, 
               c(var_wt_x, var_wt_y, var_x, var_y, var_y_cat, var_z, var_g, var_h)]

lev_g <- levels(interaction(d1_cov[, var_g, drop = FALSE], drop = TRUE))
K_g <- length(lev_g)
d1_cov$g <- as.integer(factor(interaction(d1_cov[,  var_g, drop = FALSE], drop = TRUE), levels = lev_g))
d2_cov$g <- as.integer(factor(interaction(d2_cov[,  var_g, drop = FALSE], drop = TRUE), levels = lev_g))

ind1_y = as.numeric(!is.na(d1_tilde[[var_y]]))
ind2_y = as.numeric(!is.na(d2_tilde[[var_y]]))

d1_y = d1_cov[ind1_y == 1, ]
d2_y = d2_cov[ind2_y == 1, ]


### Marginal methods : eval.vec, length.vec
ar_m_da_wcp <- list()
ar_m_gwcp <- list()
cqr_m_da_wcp <- list()
cqr_m_gwcp <- list()
xgb_m_da_wcp <- list()
xgb_m_gwcp <- list()
rf_m_da_wcp <- list()
rf_m_gwcp <- list()
nn_m_da_wcp <- list()
nn_m_gwcp <- list()

rvec_marg <- vector("list", M)

ar_m_da_wcp$eval.vec <- vector(length = M)
ar_m_gwcp$eval.vec <- vector(length = M)
cqr_m_da_wcp$eval.vec <- vector(length = M)
cqr_m_gwcp$eval.vec <- vector(length = M)
xgb_m_da_wcp$eval.vec <- vector(length = M)
xgb_m_gwcp$eval.vec <- vector(length = M)
rf_m_da_wcp$eval.vec <- vector(length = M)
rf_m_gwcp$eval.vec <- vector(length = M)
nn_m_da_wcp$eval.vec <- vector(length = M)
nn_m_gwcp$eval.vec <- vector(length = M)

ar_m_da_wcp$length.vec <- vector(length = M)
ar_m_gwcp$length.vec <- vector(length = M)
cqr_m_da_wcp$length.vec <- vector(length = M)
cqr_m_gwcp$length.vec <- vector(length = M)

xgb_m_da_wcp$size.mat <- matrix(nrow = M, ncol = 5)
xgb_m_gwcp$size.mat <- matrix(nrow = M, ncol = 5)
rf_m_da_wcp$size.mat <- matrix(nrow = M, ncol = 5)
rf_m_gwcp$size.mat <- matrix(nrow = M, ncol = 5)
nn_m_da_wcp$size.mat <- matrix(nrow = M, ncol = 5)
nn_m_gwcp$size.mat <- matrix(nrow = M, ncol = 5)



### NHANES analysis

### Implementing marginal methods
for(m in 1:M){
  set.seed(123+m)
  cat("----------------------------------------------------------\n")
  cat("                ",  m, "th split progressing              \n", sep = "")
  cat("----------------------------------------------------------\n")
  
  ind_tra = replace(0*ind1_y, sample(which(ind1_y==1), sum(ind1_y==1)%/%2), 1)
  ind_cal = ind1_y - ind_tra
  ind_ttra = rep(1, nrow(d1_cov)) - ind_cal
  dtra = d1_cov[ind_tra == 1, ]
  dcal = d1_cov[ind_cal == 1, ]
  dttra = d1_cov[ind_ttra == 1, ]
  
  ind_test = replace(0*ind2_y, sample(which(ind2_y==1), ceiling(sum(ind2_y==1)*beta)), 1)
  ind_2 = rep(1, nrow(d2_cov)) - ind_test
  dtest = d2_cov[ind_test == 1, ]
  d2 = d2_cov[ind_2 == 1, ]
  
  # Continuous outcomes
  DA_WCP_ar_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "marg", 
                       score_type = "ar", gwcp = F, var_wt_x, var_wt_y, eval = T, rvec = T, rvec_cal = rvec_marg[[m]])
  # rvec_marg[[m]] = DA_WCP_ar_m$rvec_cal
  
  ar_m_da_wcp$eval.vec[m] = DA_WCP_ar_m$coverage
  ar_m_da_wcp$length.vec[m] = DA_WCP_ar_m$avglen
  
  DA_WCP_cqr_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "marg", 
                        score_type = "cqr", gwcp = F, var_wt_x, var_wt_y, eval = T, rvec = T, rvec_cal = rvec_marg[[m]])
  cqr_m_da_wcp$eval.vec[m] = DA_WCP_cqr_m$coverage
  cqr_m_da_wcp$length.vec[m] = DA_WCP_cqr_m$avglen
  
  GWCP_ar_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "marg", 
                     score_type = "ar", gwcp = T, var_wt_x, var_wt_y, eval = T)
  ar_m_gwcp$eval.vec[m] = GWCP_ar_m$coverage
  ar_m_gwcp$length.vec[m] = GWCP_ar_m$avglen
  
  GWCP_cqr_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "marg", 
                     score_type = "cqr", gwcp = T, var_wt_x, var_wt_y, eval = T)
  cqr_m_gwcp$eval.vec[m] = GWCP_cqr_m$coverage
  cqr_m_gwcp$length.vec[m] = GWCP_cqr_m$avglen
  
  # Categorical outcomes
  DA_WCP_xgb_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                     score_type = "xgb", gwcp = F, var_wt_x, var_wt_y, eval = T, rvec = T, rvec_cal = rvec_marg[[m]])
  xgb_m_da_wcp$eval.vec[m] = DA_WCP_xgb_m$coverage
  xgb_m_da_wcp$size.mat[m, ] = table(factor(rowSums(DA_WCP_xgb_m$result_mat), levels = 1:5))
  
  GWCP_xgb_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                        score_type = "xgb", gwcp = T, var_wt_x, var_wt_y, eval = T)
  xgb_m_gwcp$eval.vec[m] = GWCP_xgb_m$coverage
  xgb_m_gwcp$size.mat[m, ] = table(factor(rowSums(GWCP_xgb_m$result_mat), levels = 1:5))
  
  DA_WCP_rf_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                        score_type = "rf", gwcp = F, var_wt_x, var_wt_y, eval = T, rvec = T, rvec_cal = rvec_marg[[m]])
  rf_m_da_wcp$eval.vec[m] = DA_WCP_rf_m$coverage
  rf_m_da_wcp$size.mat[m, ] = table(factor(rowSums(DA_WCP_rf_m$result_mat), levels = 1:5))
  
  GWCP_rf_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                      score_type = "rf", gwcp = T, var_wt_x, var_wt_y, eval = T)
  rf_m_gwcp$eval.vec[m] = GWCP_rf_m$coverage
  rf_m_gwcp$size.mat[m, ] = table(factor(rowSums(GWCP_rf_m$result_mat), levels = 1:5))
  
  DA_WCP_nn_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                        score_type = "nn", gwcp = F, var_wt_x, var_wt_y, eval = T, rvec = T, rvec_cal = rvec_marg[[m]])
  nn_m_da_wcp$eval.vec[m] = DA_WCP_nn_m$coverage
  nn_m_da_wcp$size.mat[m, ] = table(factor(rowSums(DA_WCP_nn_m$result_mat), levels = 1:5))
  
  GWCP_nn_m = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "marg", 
                      score_type = "nn", gwcp = T, var_wt_x, var_wt_y, eval = T)
  nn_m_gwcp$eval.vec[m] = GWCP_nn_m$coverage
  nn_m_gwcp$size.mat[m, ] = table(factor(rowSums(GWCP_nn_m$result_mat), levels = 1:5))
  
  setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Code")
  save.image(file = "DA_WCP_NHANES_marg.RData")
}

# load("DA_WCP_NHANES_marg.RData")


### Visualization for continuous outcomes

eval.vec <- list(
  ar_m_da_wcp$eval.vec,
  ar_m_gwcp$eval.vec,
  cqr_m_da_wcp$eval.vec,
  cqr_m_gwcp$eval.vec
)

length.vec <- list(
  ar_m_da_wcp$length.vec,
  ar_m_gwcp$length.vec,
  cqr_m_da_wcp$length.vec,
  cqr_m_gwcp$length.vec
)

df_conti_m <- data.frame(
  coverage = do.call(base::c, eval.vec),
  length = do.call(base::c, length.vec),
  method = factor(rep(c("DA_WCP_AR_m", "GWCP_AR_m", "DA_WCP_CQR_m", "GWCP_CQR_m"), each = M))
)

setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Figures")

df_conti_m$method <- factor(df_conti_m$method, levels = c("GWCP_CQR_m", "DA_WCP_CQR_m", "GWCP_AR_m", "DA_WCP_AR_m"))

pdf(file = "LDL_continuous_coverage_m.pdf", width = 8, height = 6)
ggplot(df_conti_m, aes(x = coverage, y = method, fill = method)) +
  geom_boxplot(outlier.size = 1) +
  theme_test(base_size = 14) +
  geom_vline(xintercept = 0.8, color = "red", linewidth = 1) +
  xlim(c(0.55, 0.95)) +
  xlab("Coverage") +
  ylab(NULL) +
  ggtitle(NULL) +
  scale_fill_manual(
    values = c(
      "DA_WCP_AR_m"   = "#9ecae1",  
      "DA_WCP_CQR_m"  = "#3182bd",
      "GWCP_AR_m"   = "#fdae6b",  
      "GWCP_CQR_m"  = "#fd8d3c"   
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.position = "none")
dev.off()


pdf(file = "LDL_continuous_length_m.pdf", width = 8, height = 6)
ggplot(df_conti_m, aes(x = length, y = method, fill = method)) +
  geom_boxplot(outlier.size = 1) +
  theme_test(base_size = 14) +
  xlim(c(15, 30)) +
  xlab("Average interval length (mg/dL)") +
  ylab(NULL) +
  ggtitle(NULL) +
  scale_fill_manual(
    values = c(
      "DA_WCP_AR_m"   = "#9ecae1",  
      "GWCP_AR_m"   = "#fdae6b",  
      "DA_WCP_CQR_m"  = "#3182bd",
      "GWCP_CQR_m"  = "#fd8d3c"   
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.position = "none")
dev.off()


### Visualization for categorical outcomes

c(round(mean(xgb_m_da_wcp$eval.vec), 3), round(mean(rf_m_da_wcp$eval.vec), 3), round(mean(nn_m_da_wcp$eval.vec), 3),
  round(mean(xgb_m_gwcp$eval.vec), 3), round(mean(rf_m_gwcp$eval.vec), 3), round(mean(nn_m_gwcp$eval.vec), 3))

c(round(sd(xgb_m_da_wcp$eval.vec), 3), round(sd(rf_m_da_wcp$eval.vec), 3), round(sd(nn_m_da_wcp$eval.vec), 3),
  round(sd(xgb_m_gwcp$eval.vec), 3), round(sd(rf_m_gwcp$eval.vec), 3), round(sd(nn_m_gwcp$eval.vec), 3))


round(colMeans(sweep(xgb_m_da_wcp$size.mat, 1, rowSums(xgb_m_da_wcp$size.mat), "/")), 3)[1:3]
round(colMeans(sweep(rf_m_da_wcp$size.mat, 1, rowSums(rf_m_da_wcp$size.mat), "/")), 3)[1:3]
round(colMeans(sweep(nn_m_da_wcp$size.mat, 1, rowSums(nn_m_da_wcp$size.mat), "/")), 3)[1:3]
round(colMeans(sweep(xgb_m_gwcp$size.mat, 1, rowSums(xgb_m_gwcp$size.mat), "/")), 3)[1:3]
round(colMeans(sweep(rf_m_gwcp$size.mat, 1, rowSums(rf_m_gwcp$size.mat), "/")), 3)[1:3]
round(colMeans(sweep(nn_m_gwcp$size.mat, 1, rowSums(nn_m_gwcp$size.mat), "/")), 3)[1:3]

round(apply(sweep(xgb_m_da_wcp$size.mat, 1, rowSums(xgb_m_da_wcp$size.mat), "/"), 2, sd), 3)[1:3]
round(apply(sweep(rf_m_da_wcp$size.mat, 1, rowSums(rf_m_da_wcp$size.mat), "/"), 2, sd), 3)[1:3]
round(apply(sweep(nn_m_da_wcp$size.mat, 1, rowSums(nn_m_da_wcp$size.mat), "/"), 2, sd), 3)[1:3]
round(apply(sweep(xgb_m_gwcp$size.mat, 1, rowSums(xgb_m_gwcp$size.mat), "/"), 2, sd), 3)[1:3]
round(apply(sweep(rf_m_gwcp$size.mat, 1, rowSums(rf_m_gwcp$size.mat), "/"), 2, sd), 3)[1:3]
round(apply(sweep(nn_m_gwcp$size.mat, 1, rowSums(nn_m_gwcp$size.mat), "/"), 2, sd), 3)[1:3]


### Group-conditional methods
ar_g_da_wcp <- list()
ar_g_gwcp <- list()
cqr_g_da_wcp <- list()
cqr_g_gwcp <- list()
xgb_g_da_wcp <- list()
xgb_g_gwcp <- list()
rf_g_da_wcp <- list()
rf_g_gwcp <- list()
nn_g_da_wcp <- list()
nn_g_gwcp <- list()

rvec_group <- vector("list", M)

ar_g_da_wcp$eval.mat <- matrix(nrow = M, ncol = K_g)
ar_g_gwcp$eval.mat <- matrix(nrow = M, ncol = K_g)
cqr_g_da_wcp$eval.mat <- matrix(nrow = M, ncol = K_g)
cqr_g_gwcp$eval.mat <- matrix(nrow = M, ncol = K_g)
xgb_g_da_wcp$eval.mat <- matrix(nrow = M, ncol = K_g)
xgb_g_gwcp$eval.mat <- matrix(nrow = M, ncol = K_g)
rf_g_da_wcp$eval.mat <- matrix(nrow = M, ncol = K_g)
rf_g_gwcp$eval.mat <- matrix(nrow = M, ncol = K_g)
nn_g_da_wcp$eval.mat <- matrix(nrow = M, ncol = K_g)
nn_g_gwcp$eval.mat <- matrix(nrow = M, ncol = K_g)

ar_g_da_wcp$length.mat <- matrix(nrow = M, ncol = K_g)
ar_g_gwcp$length.mat <- matrix(nrow = M, ncol = K_g)
cqr_g_da_wcp$length.mat <- matrix(nrow = M, ncol = K_g)
cqr_g_gwcp$length.mat <- matrix(nrow = M, ncol = K_g)


### Implementing group-conditional methods

for(m in 1:M){
  set.seed(123+m)
  cat("----------------------------------------------------------\n")
  cat("                ",  m, "th split progressing              \n", sep = "")
  cat("----------------------------------------------------------\n")
  
  ind_tra = rep(0, nrow(d1_cov))
  ind_cal = rep(0, nrow(d1_cov))
  ind_ttra = rep(0, nrow(d1_cov))
  ind_test = rep(0, nrow(d2_cov))
  ind_2 = rep(0, nrow(d2_cov))
  
  for(k_g in 1:K_g){
    ind1_cov_g = as.numeric(d1_cov$g == k_g)
    ind1_y_g = as.numeric(d1_cov$g == k_g & ind1_y == 1)
    ind2_cov_g = as.numeric(d2_cov$g == k_g)
    ind2_y_g = as.numeric(d2_cov$g == k_g & ind2_y == 1)
    
    ind_tra_g = replace(0*ind1_y_g, sample(which(ind1_y_g==1), sum(ind1_y_g==1)%/%2), 1)
    ind_cal_g = ind1_y_g - ind_tra_g
    ind_ttra_g = ind1_cov_g - ind_cal_g
    
    ind_test_g = replace(0*ind2_y_g, sample(which(ind2_y_g==1), ceiling(sum(ind2_y_g==1)*beta)), 1)
    ind_2_g = ind2_cov_g - ind_test_g
    
    ind_tra = ind_tra + ind_tra_g
    ind_cal = ind_cal + ind_cal_g
    ind_ttra = ind_ttra + ind_ttra_g
    ind_test = ind_test + ind_test_g
    ind_2 = ind_2 + ind_2_g
  }
  
  dtra = d1_cov[ind_tra == 1, ]
  dcal = d1_cov[ind_cal == 1, ]
  dttra = d1_cov[ind_ttra == 1, ]
  dtest = d2_cov[ind_test == 1, ]
  d2 = d2_cov[ind_2 == 1, ]
  
  ### DA_WCP methods
  DA_WCP_ar_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "group", 
                       score_type = "ar", gwcp = F, var_wt_x, var_wt_y, rvec = T, rvec_cal = rvec_group[[m]])
  # rvec_group[[m]] = DA_WCP_ar_g$rvec_cal
  
  DA_WCP_cqr_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "group", 
                        score_type = "cqr", gwcp = F, var_wt_x, var_wt_y, rvec = T, rvec_cal = rvec_group[[m]])
  
  DA_WCP_xgb_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                        score_type = "xgb", gwcp = F, var_wt_x, var_wt_y, rvec = T, rvec_cal = rvec_group[[m]])
  
  DA_WCP_rf_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                        score_type = "rf", gwcp = F, var_wt_x, var_wt_y, rvec = T, rvec_cal = rvec_group[[m]])
  
  DA_WCP_nn_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                        score_type = "nn", gwcp = F, var_wt_x, var_wt_y, rvec = T, rvec_cal = rvec_group[[m]])
  
  ### GWCP methods
  GWCP_ar_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "group", 
                     score_type = "ar", gwcp = T, var_wt_x, var_wt_y)
  
  GWCP_cqr_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type = "group", 
                      score_type = "cqr", gwcp = T, var_wt_x, var_wt_y)
  
  GWCP_xgb_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                      score_type = "xgb", gwcp = T, var_wt_x, var_wt_y)
  
  GWCP_rf_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                      score_type = "rf", gwcp = T, var_wt_x, var_wt_y)
  
  GWCP_nn_g = DA_WCP(dtra, dcal, dttra, d2, dtest, var_x, var_y_cat, var_z, var_g, var_h, val_type = "group", 
                      score_type = "nn", gwcp = T, var_wt_x, var_wt_y)

  # Store results
  ar_g_da_wcp$eval.mat[m, ] = DA_WCP_ar_g$coverage
  ar_g_gwcp$eval.mat[m, ] = GWCP_ar_g$coverage
  cqr_g_da_wcp$eval.mat[m, ] = DA_WCP_cqr_g$coverage
  cqr_g_gwcp$eval.mat[m, ] = GWCP_cqr_g$coverage
  xgb_g_da_wcp$eval.mat[m, ] = DA_WCP_xgb_g$coverage
  xgb_g_gwcp$eval.mat[m, ] = GWCP_xgb_g$coverage
  rf_g_da_wcp$eval.mat[m, ] = DA_WCP_rf_g$coverage
  rf_g_gwcp$eval.mat[m, ] = GWCP_rf_g$coverage
  nn_g_da_wcp$eval.mat[m, ] = DA_WCP_nn_g$coverage
  nn_g_gwcp$eval.mat[m, ] = GWCP_nn_g$coverage
  
  ar_g_da_wcp$length.mat[m, ] = DA_WCP_ar_g$avglen
  ar_g_gwcp$length.mat[m, ] = GWCP_ar_g$avglen
  cqr_g_da_wcp$length.mat[m, ] = DA_WCP_cqr_g$avglen
  cqr_g_gwcp$length.mat[m, ] = GWCP_cqr_g$avglen
  
  setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Code")
  save.image(file = "DA_WCP_NHANES_group.RData")
}

load("DA_WCP_NHANES_group.RData")


### Visualization of group-conditional methods

eval.mat.conti <- list(
  ar_g_da_wcp$eval.mat, 
  ar_g_gwcp$eval.mat, 
  cqr_g_da_wcp$eval.mat, 
  cqr_g_gwcp$eval.mat
)

eval.mat.cat <- list(
  xgb_g_da_wcp$eval.mat, 
  xgb_g_gwcp$eval.mat, 
  rf_g_da_wcp$eval.mat, 
  rf_g_gwcp$eval.mat,
  nn_g_da_wcp$eval.mat, 
  nn_g_gwcp$eval.mat
)

length.mat.conti <- list(
  ar_g_da_wcp$length.mat, 
  ar_g_gwcp$length.mat, 
  cqr_g_da_wcp$length.mat, 
  cqr_g_gwcp$length.mat
)


df_conti_g <- data.frame(
  coverage = do.call(base::c, eval.mat.conti),
  length = do.call(base::c, length.mat.conti),
  method = factor(rep(c("DA_WCP_AR_g", "GWCP_AR_g", "DA_WCP_CQR_g", "GWCP_CQR_g"), each = M*K_g))
)

df_cat_g <- data.frame(
  coverage = do.call(base::c, eval.mat.cat),
  method = factor(rep(c("DA_WCP_XGB_g", "GWCP_XGB_g", "DA_WCP_RF_g", "GWCP_RF_g",
                        "DA_WCP_NN_g", "GWCP_NN_g"), each = M*K_g))
)


df_conti_g$subgroup <- factor(rep(seq_len(K_g), each = M, times = 4))
df_conti_g$method <- factor(df_conti_g$method, levels = c("DA_WCP_CQR_g", "GWCP_CQR_g", "DA_WCP_AR_g", "GWCP_AR_g"))

setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Figures")
pdf(file = "LDL_continuous_coverage_g.pdf", width = 13, height = 8)
ggplot(df_conti_g, aes(x = subgroup, y = coverage, fill = subgroup)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.8) +
  facet_wrap(~method, nrow = 2) +
  geom_hline(yintercept = 0.8,  color = "red", linewidth = 0.6) +
  ylim(0.3, 1.0) +
  labs(x = "Subroup", y = "Subgroup-wise coverage") +
  theme_test(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray95"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )
dev.off()


df_cat_g$subgroup <- factor(rep(seq_len(K_g), each = M, times = 6))
df_cat_g$method <- factor(df_cat_g$method, 
                          levels = c("DA_WCP_XGB_g", "GWCP_XGB_g", 
                                     "DA_WCP_RF_g", "GWCP_RF_g",
                                     "DA_WCP_NN_g", "GWCP_NN_g"))

pdf(file = "LDL_categorical_coverage_g.pdf", width = 13, height = 11)
ggplot(df_cat_g, aes(x = subgroup, y = coverage, fill = subgroup)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.8) +
  facet_wrap(~method, ncol = 2) +
  geom_hline(yintercept = 0.8,  color = "red", linewidth = 0.6) +
  ylim(0.4, 1.0) +
  labs(x = "Subroup", y = "Subgroup-wise coverage") +
  theme_test(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray95"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )
dev.off()






