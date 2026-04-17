# =============================================================================
# Script : Spatial cross-validation of the joint INLA model — block holdout by Site
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
# Install INLA :
# install.packages(
#   "INLA",
#   repos = c(getOption("repos"),
#             INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE
# )
library(dplyr)
library(tidyr)
library(sf)
library(INLA)
library(ggplot2)
library(fmesher)
library(pROC)
library(parallel)
library(doParallel)
library(foreach)
library(Matrix)


# -----------------------------------------------------------------------------
# 1. Importing data
# -----------------------------------------------------------------------------
substrate_df <- read.csv("01_Statistics/00_Data/02_Processed/CSV/substrate.csv") %>%
  select(Site, X, Y, Type, Essence, Diameter, Sap_max, P_Bux)

# -----------------------------------------------------------------------------
# 2. Full model based on all the data
# -----------------------------------------------------------------------------
out_inla <- "01_Statistics/00_Data/03_Results/05_joint_INLA"
dir.create(out_inla, recursive = TRUE, showWarnings = FALSE)
# Uncomment to recalculate; otherwise, load the RDS
# source("01_Statistics/01_R_Scripts/Functions/04_joint_INLA_fun.R")
# res <- fun_joint_inla(
#   data = substrate_df,
#   max.edge = c(5, 20), cutoff = 1, offset = c(10, 30),
#   prior.range1 = c(20, 0.9), prior.sigma1 = c(1, 0.5),
#   prior.range2 = c(50, 0.5), prior.sigma2 = c(1, 0.01)
# )

# saveRDS(res, file.path(out_inla, "res_joint_inla.rds"))
res <- readRDS(file.path(out_inla, "res_joint_inla.rds"))

# Verify that the site is a factor with the correct levels
res$data$Site <- factor(res$data$Site)
message("Sites : ", paste(levels(res$data$Site), collapse = ", "))
message("N observations : ", nrow(res$data))


# -----------------------------------------------------------------------------
# 3. Test on 1 replicate
# -----------------------------------------------------------------------------
source("01_Statistics/01_R_Scripts/Functions/05_joint_INLA_valid_fun.R")
cv_test <- run_parallel_cv(
  res = res,
  n_rep = 1,
  radius = 40,
  n_cores = 3,
  error_handling = "pass" # returns the error
)
print(cv_test)

# -----------------------------------------------------------------------------
# 4. Full run
# -----------------------------------------------------------------------------
# cv_results <- run_parallel_cv(
#   res = res,
#   n_rep = 100,
#   radius = 40,
#   n_cores = 9,
#   error_handling = "remove"
# )
# saveRDS(cv_results,"01_Statistics/00_Data/03_Results/06_joint_INLA_valid/cv_results_block40m_100rep.rds")

# -----------------------------------------------------------------------------
# 5. Visualization
# -----------------------------------------------------------------------------
# cv_results <- readRDS("01_Statistics/00_Data/03_Results/06_joint_INLA_valid/cv_results_block40m_100rep.rds")
#
# ggplot(cv_results, aes(x = auc)) +
#   geom_histogram(bins = 20, fill = "#1D9E75", color = "white", alpha = 0.8) +
#   geom_vline(xintercept = median(cv_results$auc, na.rm = TRUE),
#              linetype = "dashed", color = "gray30") +
#   labs(
#     title = "AUC distribution — Spatial cross validation (40 m radius)",
#     x = "AUC",
#     y = "Number of replicates"
#   ) +
#   theme_minimal()
