# =============================================================================
# Script : Joint INLA model — LGCP (deadwood) + Binomial (Buxbaumia viridis)
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
library(spdep)
library(fmesher)
library(pROC)
library(patchwork)
library(gtsummary)
library(gt)
setwd("C:/Users/pm83056/OneDrive - Office National des Forets/Bureau/Spatial_Buxbaumia_viridis/01_Statistics")

# -----------------------------------------------------------------------------
# 1. Importing data
# -----------------------------------------------------------------------------
substrate_df <- read.csv("00_Data/02_Processed/CSV/substrate.csv") %>%
  select(Site, X, Y, Type, Essence, Diameter, Sap_max, P_Bux)

# -----------------------------------------------------------------------------
# 2. Run joint INLA model
# -----------------------------------------------------------------------------
source("01_R_Scripts/Functions/04_joint_INLA_fun.R")

# For the final run, replace eb with ccd
res <- fun_joint_inla(
  data = substrate_df,
  max.edge = c(5, 20),
  cutoff = 1,
  offset = c(10, 30),
  
  # Priors w1 (deadwood spatial field)
  # Expected spatial autocorrelation range of deadwood distribution : ~10-20 m
  # within a forest plot. P(range < 20 m) = 0.9 is a strong prior that prevents
  # w1 from competing with w2. Without this constraint, w1 spontaneously converges
  # towards ~80 m — the same scale as w2 — creating competition between the two
  # fields that renders the model non-identifiable and artificially inflates all
  # fixed effects. The prior c(20, 0.9) cleanly separates the two processes by
  # assigning them distinct spatial scales, a necessary condition for joint model
  # identifiability.
  prior.range1 = c(20, 0.9),
  
  # In a forest, deadwood distribution is typically heterogeneous : some areas
  # concentrate windthrows and large deadwood, others have none. A field amplitude
  # σ ~ 1 on the log-intensity scale corresponds to an intensity ratio of
  # exp(1) ≈ 2.7 between dense and sparse areas, which is ecologically realistic.
  # This prior has little influence on Process 2 fixed effects (stable model).
  prior.sigma1 = c(1, 0.5),
  
  # Priors w2 (residual Buxbaumia field)
  # Prior centred on 50 m : residual autocorrelation of Buxbaumia after controlling
  # for covariates likely reflects spore dispersal and microclimatic effects at
  # plot scale. Epiphytic moss spores can disperse over several tens of metres.
  # This prior allows the data to freely estimate the actual range.
  prior.range2 = c(50, 0.5),
  
  # Weakly informative prior : σ can vary freely from near-zero to very large.
  # The amplitude of w2 is not constrained because it is a residual effect whose
  # magnitude is unknown a priori — a non-informative prior lets the data speak.
  prior.sigma2 = c(1, 0.01)
)

# PC (Penalized Complexity) priors on w1 were defined to constrain its range to
# the local scale of deadwood (~20 m), preventing competition with w2 and ensuring
# joint model identifiability. This constraint is justified by the nature of tree
# fall processes, which generate spatial structures at the scale of a few metres.
# The robustness of results was verified by a prior sensitivity analysis (Table X),
# showing that Process 2 fixed effects are stable across all tested configurations.

# -----------------------------------------------------------------------------
# 3. Export / Import
# -----------------------------------------------------------------------------
out_inla <- "00_Data/03_Results/05_joint_INLA"

saveRDS(res, file.path(out_inla, "res_joint_inla.rds"))

res <- readRDS(file.path(out_inla, "res_joint_inla.rds"))

# -----------------------------------------------------------------------------
# 4. Global summary
# -----------------------------------------------------------------------------
r <- res$inla
summary(r)

local({
  cat("\n--- Model metrics ---\n")
  cat("DIC  =", round(res$metrics$DIC,  2), "\n")
  cat("WAIC =", round(res$metrics$WAIC, 2), "\n")
  cat("CPO  =", round(res$metrics$CPO,  4), "\n")
  cat("\n--- Mesh ---\n")
  cat("Mesh nodes   :", res$mesh$n, "\n")
  cat("Observations :", nrow(res$data), "\n")
  cat("Ratio        :", round(res$mesh$n / nrow(res$data), 2), "\n")
})

# Predicted vs observed colonisation probability by site
pred_summary <- res$pred %>%
  group_by(Site) %>%
  summarise(
    P_predicted  = round(mean(p_buxbaumia), 3),
    P_observed   = round(mean(P_Bux), 3),
    N_substrates = n(),
    N_colonised  = sum(P_Bux)
  )

tbl_pred_summary <- pred_summary %>%
  gt() %>%
  tab_header(title    = "Joint INLA — predicted vs observed colonisation by site") %>%
  fmt_number(columns = c(P_predicted, P_observed), decimals = 3) %>%
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(pred_summary), 2))
  ) %>%
  tab_options(
    table.border.top.color            = "black",
    table.border.bottom.color         = "black",
    column_labels.border.bottom.color = "black",
    column_labels.font.weight         = "bold",
    table.font.size                   = 12
  )

tbl_pred_summary

# Save
#Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")
#gtsave(tbl_pred_summary, "pred_vs_observed_by_site.png", path = "02_Displays/Tables")


# --- Tables ---
# Fixed effects table
fixed_df <- res$inla$summary.fixed %>%
  tibble::rownames_to_column("Term") %>%
  rename(Mean = mean, SD = sd, Q025 = `0.025quant`,
         Median = `0.5quant`, Q975 = `0.975quant`,
         Mode = mode, KLD = kld)

tbl_joint_fixed <- fixed_df %>%
  gt() %>%
  tab_header(title = "Joint INLA model — fixed effects") %>%
  fmt_number(columns = c(Mean, SD, Q025, Median, Q975, Mode), decimals = 3) %>%
  fmt_scientific(columns = KLD, decimals = 3) %>%
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(fixed_df), 2))
  ) %>%
  tab_options(
    table.border.top.color            = "black",
    table.border.bottom.color         = "black",
    column_labels.border.bottom.color = "black",
    column_labels.font.weight         = "bold",
    table.font.size                   = 12
  )

tbl_joint_fixed

# Hyperparameters table (spatial fields + coupling)
hyper_df <- res$inla$summary.hyperpar %>%
  tibble::rownames_to_column("Parameter") %>%
  rename(Mean = mean, SD = sd, Q025 = `0.025quant`,
         Median = `0.5quant`, Q975 = `0.975quant`,
         Mode = mode) %>%
  # Add readable parameter group label
  mutate(Process = case_when(
    grepl("w1",   Parameter) & !grepl("copy", Parameter) ~ "Process 1 — Deadwood field (w1)",
    grepl("copy", Parameter) ~ "Coupling (β_s)",
    grepl("w2",   Parameter) ~ "Process 2 — Buxbaumia field (w2)"
  ))

tbl_joint_hyper <- hyper_df %>%
  select(Process, Parameter, Mean, SD, Q025, Median, Q975, Mode) %>%
  gt(groupname_col = "Process") %>%
  tab_header(title = "Joint INLA model — spatial field hyperparameters") %>%
  fmt_number(columns = c(Mean, SD, Q025, Median, Q975, Mode), decimals = 3) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(hyper_df), 2))
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    row_group.border.bottom.color = "black",
    column_labels.font.weight = "bold",
    table.font.size = 12
  )

tbl_joint_hyper

# Save tables
#Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")
#gtsave(tbl_joint_fixed, "joint_INLA_fixed.png",  path = "02_Displays/Tables")
#gtsave(tbl_joint_hyper, "joint_INLA_hyper.png",  path = "02_Displays/Tables")

# -----------------------------------------------------------------------------
# 5. Fixed effects - Process 2
# -----------------------------------------------------------------------------
fig_dir_joint <- "02_Displays/Figures/Joint_INLA"
save_fig <- FALSE # save figures in folder

fixed <- r$summary.fixed %>%
  mutate(param = rownames(.)) %>%
  filter(!param %in% c("intercept_p1", "intercept_p2")) %>%
  mutate(param = reorder(param, mean))

p_fixed <- ggplot(fixed, aes(x = mean, y = param)) +
  geom_point(size = 3, color = "#6161d9") +
  geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`),
                 height = 0.25, color = "#6161d9") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Joint INLA fixed effects — Process 2 (Buxbaumia viridis)",
       x = "Effect on logit(P(Buxbaumia))",
       y = "") +
  theme_minimal()

print(p_fixed)
if (save_fig) ggsave(file.path(fig_dir_joint, "fixed_effects_P2.png"),
                       plot = p_fixed, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 6. Spatial field hyperparameters
# -----------------------------------------------------------------------------
hw1 <- r$summary.hyperpar["Range for w1", ]
sw1 <- r$summary.hyperpar["Stdev for w1", ]
hw2 <- r$summary.hyperpar["Range for w2", ]
sw2 <- r$summary.hyperpar["Stdev for w2", ]

local({
  cat("\n--- Field w1 (deadwood distribution) ---\n")
  cat("Range :", round(hw1$mean, 1), "m",
      " [", round(hw1$`0.025quant`, 1), "—", round(hw1$`0.975quant`, 1), "m]\n")
  cat("Sigma :", round(sw1$mean, 3),
      " [", round(sw1$`0.025quant`, 3), "—", round(sw1$`0.975quant`, 3), "]\n")
  cat("\n--- Field w2 (residual Buxbaumia) ---\n")
  cat("Range :", round(hw2$mean, 1), "m",
      " [", round(hw2$`0.025quant`, 1), "—", round(hw2$`0.975quant`, 1), "m]\n")
  cat("Sigma :", round(sw2$mean, 3),
      " [", round(sw2$`0.025quant`, 3), "—", round(sw2$`0.975quant`, 3), "]\n")
})

# -----------------------------------------------------------------------------
# 7. β_s - coupling coefficient
# -----------------------------------------------------------------------------
beta_s_marg <- r$marginals.hyperpar[[grep("Beta for w1", names(r$marginals.hyperpar))]]
beta_s_df   <- as.data.frame(beta_s_marg)
colnames(beta_s_df) <- c("x", "density")
p_positive  <- 1 - inla.pmarginal(0, beta_s_marg)

local({
  cat("\n--- β_s ---\n")
  cat("Mean   :", round(r$summary.hyperpar["Beta for w1.copy", "mean"],       4), "\n")
  cat("95% CI : [",
      round(r$summary.hyperpar["Beta for w1.copy", "0.025quant"], 4), ",",
      round(r$summary.hyperpar["Beta for w1.copy", "0.975quant"], 4), "]\n")
  cat("P(β_s > 0) :", round(p_positive, 3), "\n")
})

p_beta <- ggplot(beta_s_df, aes(x, density)) +
  geom_line(linewidth = 1, color = "#6161d9") +
  geom_area(data = subset(beta_s_df, x > 0), fill = "#6161d9", alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  annotate("text",
           x = max(beta_s_df$x) * 0.5,
           y = max(beta_s_df$density) * 0.85,
           label = paste0("P(β_s > 0) = ", round(p_positive, 3)),
           size = 4, color = "gray30") +
  labs(title = "Posterior distribution of β_s",
       subtitle = "Coupling : deadwood density → Buxbaumia viridis",
       x = "β_s", y = "Posterior density") +
  theme_minimal()

print(p_beta)
if (save_fig) ggsave(file.path(fig_dir_joint, "beta_s_posterior.png"),
                       plot = p_beta, width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. Predictions - calibration and discrimination
# -----------------------------------------------------------------------------

# Distribution of predicted probabilities by site
pred_distrib <- ggplot(res$pred, aes(x = p_buxbaumia, fill = factor(P_Bux))) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 40) +
  facet_wrap(~ Site) +
  scale_fill_manual(values = c("0" = "#888780", "1" = "#6161d9"),
                    labels = c("Absent", "Present")) +
  labs(title = "Distribution of predicted probabilities by site",
       x = "Predicted P(Buxbaumia)",
       y = "Number of substrates",
       fill = "Observation") +
  theme_minimal()
print(pred_distrib)

if (save_fig) ggsave(file.path(fig_dir_joint, "predictions_distribution.png"),
                     plot = pred_distrib, width = 10, height = 7, dpi = 300)

# ROC curve and AUC
roc_obj <- pROC::roc(as.numeric(res$pred$P_Bux), res$pred$p_buxbaumia, quiet = TRUE)
auc_val <- round(pROC::auc(roc_obj), 3)
cat(paste0("\nAUC = ", auc_val, "\n"))

# Convert to dataframe for ggplot
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#6161d9", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.75, y = 0.1,
           label = paste0("AUC = ", auc_val),
           size = 4.5, color = "gray30") +
  labs(title = "ROC curve — joint INLA model",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal()

print(p_roc)
if (save_fig) ggsave(file.path(fig_dir_joint, "ROC_curve.png"),
                       plot = p_roc, width = 7, height = 7, dpi = 300)
# -----------------------------------------------------------------------------
# 9. CPO / PIT
# -----------------------------------------------------------------------------
cpo_vals <- r$cpo$cpo
pit_vals <- r$cpo$pit

# Filter on binomial observations only
# The first n_nodes * n_sites values correspond to LGCP mesh nodes
n_lgcp <- res$mesh$n * nlevels(res$data$Site)
idx_binom <- (n_lgcp + 1):length(cpo_vals)

cpo_binom <- cpo_vals[idx_binom]
pit_binom <- pit_vals[idx_binom]

cpo_clean <- cpo_binom[!is.na(cpo_binom) & cpo_binom > 1e-6]
cat(paste0("Valid CPO values : ", length(cpo_clean), " / ", length(cpo_binom), "\n",
           "CPO log-score   : ", round(-mean(log(cpo_clean), na.rm = TRUE), 4), "\n"))
# CPO ~ 0.3 indicates good predictive performance

# PIT histogram : should be uniform if the model is well calibrated
pit_binom <- res$inla$cpo$pit[idx_binom]
pit_df <- data.frame(pit = pit_binom[!is.na(pit_binom)])
expected_height <- nrow(pit_df) / 20

p_pit <- ggplot(pit_df, aes(x = pit)) +
  geom_histogram(breaks = seq(0, 1, length.out = 21),
                 fill = "#9595dc", color = "white") +
  geom_hline(yintercept = expected_height,
             linetype = "dashed", color = "#6161d9") +
  labs(title = "PIT — calibration (Buxbaumia observations only)",
       x = "PIT", y = "Count") +
  theme_minimal()

# Note : concentration near 1 indicates the model underestimates presences.
# This pattern is expected when modelling rare species occurrences —
# see literature on the limits of SDMs for rare species.

print(p_pit)
if (save_fig) ggsave(file.path(fig_dir_joint, "PIT_calibration.png"),
                     plot = p_pit, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 10. Spatial field and prediction maps by site
# -----------------------------------------------------------------------------
w1_mean <- r$summary.random$w1$mean
w2_mean <- r$summary.random$w2$mean
n_nodes <- res$mesh$n

if ("package:patchwork" %in% search()) detach("package:patchwork", unload = TRUE)

for (k in 1:nlevels(res$data$Site)) {
  
  site_k    <- levels(res$data$Site)[k]
  df_k      <- res$data %>% filter(site_id == k)
  df_pred_k <- res$pred %>% filter(site_id == k)
  idx_k     <- ((k - 1) * n_nodes + 1):(k * n_nodes)
  
  # Prediction grid
  x_seq <- seq(min(df_k$X) - 5, max(df_k$X) + 5, length.out = 100)
  y_seq <- seq(min(df_k$Y) - 5, max(df_k$Y) + 5, length.out = 100)
  grid  <- expand.grid(X = x_seq, Y = y_seq)
  grid$X_loc <- grid$X - mean(df_k$X)
  grid$Y_loc <- grid$Y - mean(df_k$Y)
  
  A_grid  <- inla.spde.make.A(res$mesh, loc = as.matrix(grid[, c("X_loc", "Y_loc")]))
  grid$w1 <- as.vector(A_grid %*% w1_mean[idx_k])
  grid$w2 <- as.vector(A_grid %*% w2_mean[idx_k])
  
  # w1 map : deadwood spatial field
  p_w1 <- ggplot(grid, aes(X, Y, fill = w1)) +
    geom_raster() + coord_equal() +
    scale_fill_viridis_c(name = "w1") +
    geom_point(data = df_k, aes(X, Y, shape = factor(P_Bux)),
               inherit.aes = FALSE, color = "darkred", size = 1.5) +
    scale_shape_manual(values = c("0" = 1, "1" = 16),
                       labels = c("Absent", "Present"),
                       name   = "Buxbaumia") +
    labs(title = paste("w1 — deadwood density,", site_k)) +
    theme_minimal()
  
  # w2 map : residual Buxbaumia spatial field
  p_w2 <- ggplot(grid, aes(X, Y, fill = w2)) +
    geom_raster() + coord_equal() +
    scale_fill_viridis_c(option = "magma", name = "w2") +
    geom_point(data = df_k, aes(X, Y, shape = factor(P_Bux)),
               inherit.aes = FALSE, color = "white", size = 1.5) +
    scale_shape_manual(values = c("0" = 1, "1" = 16),
                       labels = c("Absent", "Present"),
                       name   = "Buxbaumia") +
    labs(title = paste("w2 — residual Buxbaumia spatial field,", site_k)) +
    theme_minimal()
  
  # P(Buxbaumia) prediction map
  p_pred <- ggplot(df_pred_k, aes(X, Y, color = p_buxbaumia, shape = factor(P_Bux))) +
    geom_point(size = 2.5) + coord_equal() +
    scale_color_viridis_c(name = "P(BUX)") +
    scale_shape_manual(values = c("0" = 1, "1" = 16),
                       labels = c("Absent", "Present"),
                       name   = "Observation") +
    labs(title = paste("Predicted P(Buxbaumia),", site_k)) +
    theme_minimal()
  
  print(p_w1)
  print(p_w2)
  print(p_pred)
  
  if (save_fig) {
    ggsave(file.path(fig_dir_joint, paste0("w1_", site_k, ".png")),
           plot = p_w1,  width = 7, height = 6, dpi = 300)
    ggsave(file.path(fig_dir_joint, paste0("w2_", site_k, ".png")),
           plot = p_w2,  width = 7, height = 6, dpi = 300)
    ggsave(file.path(fig_dir_joint, paste0("pred_", site_k, ".png")),
           plot = p_pred, width = 7, height = 6, dpi = 300)
  }
}

# -----------------------------------------------------------------------------
# 11. Global prediction map (all sites)
# -----------------------------------------------------------------------------

plots_sites <- lapply(levels(res$pred$Site), function(s) {
  df_s <- res$pred %>% filter(Site == s)
  ggplot(df_s, aes(X, Y, color = p_buxbaumia, shape = factor(P_Bux))) +
    geom_point(size = 2) +
    scale_color_viridis_c(name   = "P(BUX)",
                          limits = range(df_s$p_buxbaumia)) +
    scale_shape_manual(values = c("0" = 1,    "1" = 16),
                       labels = c("Absent", "Present"),
                       name   = "Observation") +
    coord_equal() +
    labs(title = s, x = "X", y = "Y") +
    theme_minimal()
})

library(patchwork)
p_all_sites <- wrap_plots(plots_sites, ncol = 2) +
  plot_annotation(title = "Predicted presence probabilities - all sites")

print(p_all_sites)

if (save_fig) ggsave(file.path(fig_dir_joint, "pred_all_sites.png"),
                     plot = p_all_sites, width = 12, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 12. Mesh visualization
# -----------------------------------------------------------------------------
# Mesh plot 
res$mesh_plot()
if (save_fig) {
  png(file.path(fig_dir_joint, "mesh.png"), width = 1800, height = 1600, res = 150)
  res$mesh_plot()
  dev.off()
}

# Voronoi plot — returns a ggplot object, use ggsave()
p_voronoi <- res$voronoi_plot()
print(p_voronoi)
if (save_fig) ggsave(file.path(fig_dir_joint, "voronoi.png"),
                     plot = p_voronoi, width = 15, height = 12, dpi = 300)

# -----------------------------------------------------------------------------
# 13. Absence of coupling signal - verification
# -----------------------------------------------------------------------------
radius <- 10
coords_mat <- as.matrix(res$data[, c("X_loc", "Y_loc")])
density_bm <- sapply(seq_len(nrow(res$data)), function(i) {
  dists <- sqrt(rowSums((coords_mat - coords_mat[i, ])^2))
  sum(dists > 0 & dists <= radius)
})
res$data$density_bm <- density_bm

# Spearman correlation by site
spearman_site <- res$data %>%
  group_by(Site) %>%
  summarise(
    rho   = cor(density_bm, P_Bux, method = "spearman"),
    p_val = cor.test(density_bm, P_Bux, method = "spearman")$p.value,
    n     = n()
  )

# Global Spearman test
spearman_global <- cor.test(res$data$density_bm, res$data$P_Bux, method = "spearman")

# Add global row
spearman_df <- bind_rows(
  spearman_site,
  data.frame(Site  = "All sites",
             rho   = spearman_global$estimate,
             p_val = spearman_global$p.value,
             n     = nrow(res$data))
)

tbl_spearman <- spearman_df %>%
  rename(Site = Site, Rho = rho, P_value = p_val, N = n) %>%
  gt() %>%
  tab_header(
    title    = "Spearman correlation : local deadwood density vs. Buxbaumia presence",
    subtitle = paste0("Local density = number of substrates within ", radius, " m radius")
  ) %>%
  fmt_number(columns = Rho, decimals = 3) %>%
  fmt_scientific(columns = P_value, decimals = 3) %>%
  tab_style(
    style = cell_fill(color = "#f5f5f5"),
    locations = cells_body(rows = seq(1, nrow(spearman_df), 2))
  ) %>%
  # Visual separator before global row
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(2)),
    locations = cells_body(rows = Site == "All sites")
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.font.weight = "bold",
    table.font.size = 12
  )

tbl_spearman

# Save table
#Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")
#gtsave(tbl_spearman, "spearman_coupling.png", path = "02_Displays/Tables")



