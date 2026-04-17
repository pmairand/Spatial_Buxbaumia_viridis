# =============================================================================
# Script : Spatial point pattern analysis — PCF, MCF and PPM (site-by-site)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(GGally)
library(spatstat)
library(spatstat.explore)
library(spatstat.geom)
library(patchwork)
library(ggplot2)
library(gtsummary)
library(gt)


# -----------------------------------------------------------------------------
# 1. Importing data
# -----------------------------------------------------------------------------
substrate_df <- read.csv("01_Statistics/00_Data/02_Processed/CSV/substrate.csv") %>%
  select(Site, X, Y, Z, Type, Essence, Position, Contact, Diameter, Sap_max, P_Bux)

# -----------------------------------------------------------------------------
# 2. Correlation plots
# -----------------------------------------------------------------------------
ggpairs(
  substrate_df %>% select(Z, Type, Essence, Position, Contact, Diameter, Sap_max, P_Bux),
  upper = list(continuous = wrap("cor", size = 4, method = "spearman")),
  lower = list(continuous = wrap("points", alpha = 0.4, size = 0.8)),
  diag  = list(continuous = "densityDiag")
) + theme_bw()

# Save plot to figures folder
# output_dir <- "01_Statistics/02_Displays/Figures"
# ggsave(
#  filename = file.path(output_dir, "substrate_var_correlations.png"),
#  width = 18, height = 9, dpi = 300)

# -----------------------------------------------------------------------------
# 3. Separation of spatial coordinates by site
# -----------------------------------------------------------------------------
substrate_pts <- substrate_df %>%
  select(Site, X, Y, P_Bux)

sites <- unique(substrate_pts$Site)

for (s in sites) {
  assign(
    paste0("points_", s),
    substrate_pts %>%
      filter(Site == s) %>%
      select(-Site)
  )
}

# -----------------------------------------------------------------------------
# 4. Point pattern objects (ppp) - Ripley-Rasson window estimation
# -----------------------------------------------------------------------------
for (s in sites) {
  pts_df <- get(paste0("points_", s))

  # Window estimation (Ripley-Rasson + dilation to avoid edge effects)
  X_tmp <- ppp(
    x = pts_df$X, y = pts_df$Y,
    window = owin(range(pts_df$X), range(pts_df$Y))
  )
  W_poly <- dilation(ripras(X_tmp), r = 1)
  assign(paste0("W_", s, "_poly"), W_poly)

  # ppp object with presence/absence marks
  ppp_obj <- ppp(x = pts_df$X, y = pts_df$Y, window = W_poly)
  marks(ppp_obj) <- factor(pts_df$P_Bux, levels = c(0, 1), labels = c("abs", "pres"))
  assign(paste0(s, "_ppp"), ppp_obj)

  # Summary
  cat("\n========================================\n")
  cat("Summary of", paste0(s, "_ppp"), "\n")
  cat("========================================\n")
  print(summary(ppp_obj))

  # Convert ppp to dataframe for ggplot
  ppp_df <- as.data.frame(ppp_obj)

  # Extract window boundary for plotting
  win_df <- as.data.frame(W_poly$bdry[[1]])

  p <- ggplot(ppp_df, aes(x = x, y = y, color = marks, shape = marks)) +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("abs" = "steelblue", "pres" = "firebrick")) +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste(s, "- Ripley-Rasson Window"),
      x = "X", y = "Y",
      color = "Status", shape = "Status"
    )
  print(p)

  # Save plots
  # output_dir <- "01_Statistics/02_Displays/Figures"
  # ggsave(
  #  filename = file.path(output_dir, paste0("Rip-Ras_ppp_", s, ".png")),
  #  plot = p,
  #  width = 8, height = 8, dpi = 300)
}

# M2 : 2.2 ha | 7%  presence
# M3 : 3.4 ha | 22% presence
# M4 : 2.8 ha | 8%  presence
# M5 : 1.6 ha | 6%  presence

# -----------------------------------------------------------------------------
# 5. Site-by-site PCF & MCF analysis
# -----------------------------------------------------------------------------
source("01_Statistics/01_R_Scripts/Functions/01_PCF_fun.R")
source("01_Statistics/01_R_Scripts/Functions/02_MCF_fun.R")

names_ppp <- c("M2_ppp", "M3_ppp", "M4_ppp", "M5_ppp")

par(mfrow = c(1, 2))
for (nm in names_ppp) {
  s <- sub("_ppp", "", nm)
  message("\nPCF for site : ", s)
  PCF_function(s = s, recalculate = FALSE, save_plot = TRUE)
  message("\nMCF for site : ", s)
  MCF_function(s = s, recalculate = FALSE, save_plot = TRUE)
}

# -----------------------------------------------------------------------------
# 6. Simple PPM - site-by-site
# -----------------------------------------------------------------------------
# Effects plots function
make_effect_plot <- function(ef, xvar, xlab, rug_data, mean_val, sd_val) {
  ef_df <- as.data.frame(ef)
  ef_df$x_orig <- ef_df[[xvar]] * sd_val + mean_val

  ggplot(ef_df, aes(x = x_orig)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = lambda), color = "steelblue", linewidth = 0.8) +
    geom_rug(
      data = data.frame(x = rug_data), aes(x = x),
      inherit.aes = FALSE, alpha = 0.4
    ) +
    theme_bw() +
    labs(x = xlab, y = "Intensity", title = paste("Effect of", xlab))
}
# Note: rug marks may extend beyond the curve range. This is expected: effectfun()
# evaluates the model over the covariate image extent, which is slightly narrower
# than the raw data range due to kernel smoothing contraction at boundaries (Smooth.ppp, sigma = 15).

# PPM loop
results_ppm <- data.frame()

output_dir <- "01_Statistics/02_Displays/Figures/PPM"
save_ppm_plots <- FALSE # Save plots in Figures folder

for (s in sites) {
  # --- 6.1 Data preparation --------------------------------------------------
  df_site <- substrate_df %>%
    filter(Site == s) %>%
    select(X, Y, P_Bux, Z, Diameter, Sap_max) %>%
    drop_na()

  # --- 6.2 Point pattern objects ---------------------------------------------
  W <- get(paste0("W_", s, "_poly"))
  all_pts <- ppp(df_site$X, df_site$Y, window = W)
  df_pres <- df_site %>% filter(P_Bux == 1)
  pres_pts <- ppp(df_pres$X, df_pres$Y, window = W)
  win_df <- as.data.frame(W$bdry[[1]])

  # Plot colonised substrates
  p_pres <- ggplot() +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(
      data = df_pres, aes(x = X, y = Y),
      color = "firebrick", size = 1.5, alpha = 0.7
    ) +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - Colonised substrates (P_Bux = 1)"),
      x = "X", y = "Y"
    )

  # --- 6.3 Spatial covariates ------------------------------------------------

  # D10(u) : store mean/sd before scaling
  dens10_im <- density.ppp(all_pts, sigma = 10, edge = TRUE, weights = df_site$Diameter)
  vals <- as.vector(dens10_im)
  mean_d10 <- mean(vals, na.rm = TRUE)
  sd_d10 <- sd(vals, na.rm = TRUE)
  dens10_im_sc <- eval.im((dens10_im - mean_d10) / sd_d10)
  df_site$dens10 <- spatstat.geom::lookup.im(dens10_im_sc, x = df_site$X, y = df_site$Y)
  dens_df <- as.data.frame(dens10_im_sc)

  p_dens <- ggplot(dens_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(
      data = df_site, aes(x = X, y = Y),
      inherit.aes = FALSE, size = 0.4, alpha = 0.5
    ) +
    scale_fill_viridis_c() +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - Local deadwood density (sigma = 10 m)"),
      x = "X", y = "Y", fill = "Scaled\ndensity"
    )

  # A(u) : store mean/sd before scaling
  X_Z <- all_pts
  marks(X_Z) <- df_site$Z
  alt_im <- Smooth.ppp(X_Z, sigma = 15)
  vals_alt <- as.vector(alt_im)
  mean_alt <- mean(vals_alt, na.rm = TRUE)
  sd_alt <- sd(vals_alt, na.rm = TRUE)
  alt_im_sc <- eval.im((alt_im - mean_alt) / sd_alt)
  alt_df <- as.data.frame(alt_im_sc)

  p_alt <- ggplot(alt_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    scale_fill_viridis_c(option = "magma") +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - Smoothed altitude"),
      x = "X", y = "Y", fill = "Scaled\naltitude"
    )

  # S(u) : store mean/sd before scaling
  X_S <- all_pts
  marks(X_S) <- df_site$Sap_max
  sap_smooth <- Smooth.ppp(X_S, sigma = 3)
  dens_tiny <- density.ppp(all_pts, sigma = 0.5, edge = TRUE)
  sap_im <- eval.im(ifelse(dens_tiny > 1e-6, sap_smooth, 0))
  vals_sap <- as.vector(sap_im)
  mean_sap <- mean(vals_sap, na.rm = TRUE)
  sd_sap <- sd(vals_sap, na.rm = TRUE)
  sap_im_sc <- eval.im((sap_im - mean_sap) / sd_sap)
  sap_df <- as.data.frame(sap_im_sc)

  p_sap <- ggplot(sap_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    scale_fill_viridis_c(option = "cividis") +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - Smoothed saproxylic index"),
      x = "X", y = "Y", fill = "Scaled\nsap. index"
    )

  # Save covariates panel
  p_covs <- p_pres + p_dens + p_alt + p_sap +
    plot_annotation(title = paste("Site", s, "- Covariates"))
  if (save_ppm_plots) {
    ggsave(file.path(output_dir, paste0(s, "_covariates.png")),
      plot = p_covs, width = 16, height = 8, dpi = 300
    )
  }


  # --- 6.4 Point Process Model -----------------------------------------------
  covs <- list(dens10 = dens10_im_sc, alt = alt_im_sc, sap = sap_im_sc)
  fit_ppm <- ppm(pres_pts ~ dens10 + alt + sap, covariates = covs)

  # Coefficients
  coef_ppm <- coef(fit_ppm)
  se_ppm <- sqrt(diag(vcov(fit_ppm)))
  results_ppm <- bind_rows(
    results_ppm,
    data.frame(
      Site      = s,
      Term      = names(coef_ppm),
      Estimate  = round(coef_ppm, 4),
      Std_Error = round(se_ppm, 4),
      Z_value   = round(coef_ppm / se_ppm, 3),
      P_value   = round(2 * pnorm(-abs(coef_ppm / se_ppm)), 4)
    )
  )

  # Residuals
  res <- residuals(fit_ppm)

  res_df <- data.frame(
    x     = res$loc$x,
    y     = res$loc$y,
    value = res$val
  )

  p_resid <- ggplot(res_df, aes(x = x, y = y, color = value)) +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_gradient2(low = "#5D0C96", mid = "grey90", high = "#45960C", midpoint = 0) +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - PPM residuals"),
      x = "X", y = "Y", color = "Residuals"
    )

  # Predicted intensity
  pred_df <- as.data.frame(predict(fit_ppm, type = "trend"))
  p_pred <- ggplot(pred_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_polygon(
      data = win_df, aes(x = x, y = y),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(
      data = df_pres, aes(x = X, y = Y),
      inherit.aes = FALSE, color = "white", size = 0.8, alpha = 0.7
    ) +
    scale_fill_gradientn(colors = terrain.colors(100)) +
    coord_equal() +
    theme_bw() +
    labs(
      title = paste0(s, " - PPM predicted intensity"),
      x = "X", y = "Y", fill = "Intensity"
    )

  # Save residuals + predicted panel
  p_model <- p_resid + p_pred +
    plot_annotation(title = paste("Site", s, "- PPM diagnostics"))
  if (save_ppm_plots) {
    ggsave(file.path(output_dir, paste0(s, "_PPM_diagnostics.png")),
      plot = p_model, width = 14, height = 6, dpi = 300
    )
  }

  # Marginal effects
  ef1 <- effectfun(fit_ppm, "sap", alt = 0, dens10 = 0, se.fit = TRUE)
  ef2 <- effectfun(fit_ppm, "alt", sap = 0, dens10 = 0, se.fit = TRUE)
  ef3 <- effectfun(fit_ppm, "dens10", sap = 0, alt = 0, se.fit = TRUE)

  p_ef1 <- make_effect_plot(ef1, "sap", "Saproxylic index", df_site$Sap_max, mean_sap, sd_sap)
  p_ef2 <- make_effect_plot(ef2, "alt", "Altitude", df_site$Z, mean_alt, sd_alt)
  p_ef3 <- make_effect_plot(ef3, "dens10", "Deadwood density", df_site$dens10, mean_d10, sd_d10)

  p_effects <- p_ef1 + p_ef2 + p_ef3 +
    plot_annotation(title = paste("Site", s, "- Marginal effects"))

  if (save_ppm_plots) {
    ggsave(file.path(output_dir, paste0(s, "_marginal_effects.png")),
      plot = p_effects, width = 14, height = 5, dpi = 300
    )
  }

  print(p_covs)
  print(p_model)
  print(p_effects)
}

# -----------------------------------------------------------------------------
# 7. Results summary
# -----------------------------------------------------------------------------
make_ppm_table <- function(data, sites, title) {
  df <- data %>% filter(Site %in% sites)
  df %>%
    gt(groupname_col = "Site") %>%
    tab_header(title = title) %>%
    fmt_number(columns = c(Estimate, Std_Error, Z_value), decimals = 3) %>%
    fmt_number(columns = P_value, decimals = 4) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    tab_style(
      style = cell_fill(color = "#f5f5f5"),
      locations = cells_body(rows = seq(1, nrow(df), 2))
    ) %>%
    tab_options(
      table.border.top.color            = "black",
      table.border.bottom.color         = "black",
      column_labels.border.bottom.color = "black",
      row_group.border.bottom.color     = "black",
      column_labels.font.weight         = "bold",
      table.font.size                   = 12
    )
}

tbl_M2M3 <- make_ppm_table(results_ppm, c("M2", "M3"), "PPM coefficients — M2 & M3")
tbl_M4M5 <- make_ppm_table(results_ppm, c("M4", "M5"), "PPM coefficients — M4 & M5")

tbl_M2M3
tbl_M4M5

# Save tables
# Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")
# gtsave(tbl_M2M3, "results_ppm_M2M3.png", path = "01_Statistics/02_Displays/Tables")
# gtsave(tbl_M4M5, "results_ppm_M4M5.png", path = "01_Statistics/02_Displays/Tables")
