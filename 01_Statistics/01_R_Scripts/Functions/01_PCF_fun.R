# -----------------------------------------------------------------------------
# PCF_function : Pair Correlation Function under CSR
# -----------------------------------------------------------------------------
# Computes or imports a PCF envelope for a given site under Complete Spatial
# Randomness (CSR) hypothesis, then plots the result.
#
# Arguments:
#   s : site name (e.g. "M2")
#   recalculate : if TRUE, recomputes and saves the envelope ;
#                 if FALSE, imports the previously saved .rds file
#   correction : edge correction method passed to envelope()
#   nsim : number of simulations for the envelope
#   nrank : rank of the envelope, i.e. number of simulations ignored on
#           each side to build the CI
#   save_plot : saving plots in Figures folder
#
# Returns a list with
#   $PCF : the envelope object
#   $file : path to the .rds file
#   $ppp_name : name of the ppp object used
# -----------------------------------------------------------------------------

PCF_function <- function(s,
                         recalculate = FALSE,
                         correction = "best",
                         nsim = 199,
                         nrank = 5,
                         save_plot = FALSE) {
  # Get ppp object
  ppp_name <- paste0(s, "_ppp")
  X <- get(ppp_name, envir = .GlobalEnv)

  # File path
  fp_PCF <- file.path(paste0("01_Statistics/00_Data/03_Results/01_PCF_CSR/", ppp_name, ".rds"))

  # Recalculate or import
  if (recalculate) {
    message("Calculating PCF CSR for ", ppp_name)
    PCF <- envelope(X,
      fun = pcf, correction = correction,
      nsim = nsim, nrank = nrank
    )
    saveRDS(PCF, fp_PCF)
  } else {
    message("Importing PCF CSR from :\n", fp_PCF)
    PCF <- readRDS(fp_PCF)
  }

  # Convert envelope to dataframe for ggplot
  pcf_df <- as.data.frame(PCF)

  p <- ggplot(pcf_df, aes(x = r)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.6) +
    geom_line(aes(y = obs), color = "darkorange", linewidth = 0.7) +
    geom_line(aes(y = theo), color = "red", linetype = "dashed", linewidth = 0.6) +
    theme_bw() +
    labs(
      title = paste0("PCF CSR - ", s),
      x = "r",
      y = "g(r)"
    )

  print(p)

  # Save plot
  if (save_plot) {
    output_dir <- "01_Statistics/02_Displays/Figures/PCF_MCF"
    ggsave(
      filename = file.path(output_dir, paste0("PCF_CSR_", s, ".png")),
      plot = p,
      width = 8, height = 5, dpi = 300
    )
  }

  invisible(list(PCF = PCF, file = fp_PCF, ppp_name = ppp_name))
}
