# =============================================================================
# Validation en deux étapes — Vérification du β_s du modèle joint
#
# Logique :
#   Étape 1 : LGCP seul → estimation de w1
#   Étape 2 : Tirage de réalisations de la post. de w1 → covariable dans
#             un modèle binomial classique → estimation de β_s séquentiel
#
# Si β_s (joint) ≈ β_s (séquentiel) : le couplage est cohérent
# Si divergence importante : problème d'identifiabilité ou de propagation
#                            d'incertitude dans le modèle joint
# =============================================================================

# On suppose que res (output de fun_joint_inla) et substrate_df sont en mémoire.
# Ajuster les chemins si nécessaire.

library(INLA)
library(dplyr)
library(fmesher)
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. Récupération des objets du modèle joint
# -----------------------------------------------------------------------------
out_inla <- "01_Statistics/00_Data/03_Results/05_joint_INLA"
res <- readRDS(file.path(out_inla, "res_joint_inla.rds"))
# Réutilise le mesh, les SPDE, les données préparées de res
df       <- res$data          # données nettoyées (sans Indéterminé)
mesh     <- res$mesh
r        <- res$inla

n_sites  <- nlevels(df$Site)
n_nodes  <- mesh$n
coords_loc <- as.matrix(df[, c("X_loc", "Y_loc")])

# Reconstruire spde1 avec les mêmes hyperpriors que dans fun_joint_inla
spde1_2s <- inla.spde2.pcmatern(
  mesh        = mesh,
  alpha       = 2,
  prior.range = c(20, 0.9),
  prior.sigma = c(1, 0.5),
  constr      = TRUE
)

spde2_2s <- inla.spde2.pcmatern(
  mesh        = mesh,
  alpha       = 2,
  prior.range = c(50, 0.5),
  prior.sigma = c(1, 0.01)
)

# =============================================================================
# ÉTAPE 1 — LGCP seul : estimation de w1
# =============================================================================
message("=== ÉTAPE 1 : LGCP seul ===")

# Voronoi areas (identiques à celles du modèle joint)
voronoi_areas <- res$voronoi_areas

# Index w1 avec réplicats par site
idx_w1_s1 <- inla.spde.make.index("w1", n.spde = spde1_2s$n.spde, n.repl = n_sites)

# A matrix bloc-diagonale (identité par site, comme dans fun_joint_inla)
A_lgcp_s1 <- Matrix::bdiag(lapply(seq_len(n_sites), function(k) Matrix::Diagonal(n = n_nodes)))

n_nodes_total <- n_nodes * n_sites
E_lgcp_s1     <- rep(voronoi_areas, n_sites)

stack_lgcp_s1 <- inla.stack(
  data    = list(y_lgcp = rep(0L, n_nodes_total), E_vec = E_lgcp_s1),
  A       = list(1, A_lgcp_s1),
  effects = list(
    data.frame(intercept_p1 = rep(1, n_nodes_total)),
    w1 = idx_w1_s1
  ),
  tag = "lgcp_s1"
)

formula_lgcp <- y_lgcp ~ -1 +
  intercept_p1 +
  f(w1, model = spde1_2s, replicate = w1.repl)

dat_lgcp <- inla.stack.data(stack_lgcp_s1)
A_lgcp   <- inla.stack.A(stack_lgcp_s1)

res_lgcp <- inla(
  formula_lgcp,
  family  = "poisson",
  data    = dat_lgcp,
  E       = dat_lgcp$E_vec,
  control.predictor = list(A = A_lgcp, compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.inla      = list(int.strategy = "eb"),
  verbose = FALSE
)

message(paste0("LGCP DIC = ", round(res_lgcp$dic$dic, 2)))
message(paste0("Range w1 (séquentiel) = ",
               round(res_lgcp$summary.hyperpar["Range for w1", "mean"], 2), " m"))
message(paste0("Sigma w1 (séquentiel) = ",
               round(res_lgcp$summary.hyperpar["Stdev for w1", "mean"], 3)))

# =============================================================================
# ÉTAPE 2 — Tirage de réalisations de w1 et modèle binomial
# =============================================================================
message("=== ÉTAPE 2 : Tirage de réalisations de w1 ===")

# ---- 2a. Tirer N réalisations de la posterieure de w1 ----------------------
# inla.posterior.sample() tire des vecteurs complets (latent field + hyperpar)
# depuis la posterieure jointe approximée par INLA.
set.seed(42)
N_samples <- 100  # 100 suffit pour estimer la variabilité ; augmenter si besoin

samples <- inla.posterior.sample(N_samples, res_lgcp)

# Identifier les indices de w1 dans le vecteur latent
# Le vecteur latent INLA est ordonné : (effets fixes) | (random effects)
# On cherche les lignes dont le nom commence par "w1:"
latent_names <- rownames(samples[[1]]$latent)
idx_w1_latent <- grep("^w1:", latent_names)
message(paste0("Indices w1 dans le vecteur latent : ", length(idx_w1_latent),
               " (attendu : ", n_nodes * n_sites, ")"))

# ---- 2b. Interpoler w1 aux positions des observations ----------------------
# Pour chaque site k, extraire les n_nodes valeurs de w1 correspondantes,
# puis projeter vers les observations via la matrice A.

A_obs_list <- lapply(seq_len(n_sites), function(k) {
  df_k   <- df %>% filter(site_id == k)
  coords_k <- as.matrix(df_k[, c("X_loc", "Y_loc")])
  inla.spde.make.A(mesh = mesh, loc = coords_k)
})

# Fonction : pour un tirage donné, retourner w1 interpolé aux observations
get_w1_obs <- function(samp) {
  w1_full <- samp$latent[idx_w1_latent, 1]  # vecteur de longueur n_nodes * n_sites
  
  # w1 est stocké par blocs de n_nodes, ordonnés par site (site_id 1..n_sites)
  w1_obs_all <- numeric(nrow(df))
  for (k in seq_len(n_sites)) {
    idx_block <- ((k - 1) * n_nodes + 1):(k * n_nodes)
    w1_k      <- w1_full[idx_block]
    df_k_rows <- which(df$site_id == k)
    w1_obs_all[df_k_rows] <- as.vector(A_obs_list[[k]] %*% w1_k)
  }
  w1_obs_all
}

# Matrice des w1 interpolés : nrow(df) × N_samples
W1_mat <- sapply(samples, get_w1_obs)
# Moyenne des réalisations → covariable ponctuelle pour le modèle
w1_mean_sampled <- rowMeans(W1_mat)

message(paste0("w1 interpolé : mean = ", round(mean(w1_mean_sampled), 4),
               ", sd = ", round(sd(w1_mean_sampled), 4)))

# ---- 2c. Modèle binomial avec w1 comme covariable fixe ---------------------
message("=== Modèle binomial séquentiel ===")

# Encodage catégoriel identique à fun_joint_inla
df_binom_s2 <- df %>%
  mutate(
    Sap_max2        = as.integer(Sap_max == "2"),
    Sap_max3        = as.integer(Sap_max == "3"),
    Sap_max4        = as.integer(Sap_max == "4"),
    Sap_max5        = as.integer(Sap_max == "5"),
    Essence_Conifer = as.integer(Essence == "Résineux"),
    Type_Stump      = as.integer(Type    == "Souche"),
    intercept_p2    = 1,
    w1_cov          = w1_mean_sampled   # w1 moyen comme covariable fixe
  )

# Matrice A pour le champ résiduel w2
A_obs_w2_s2 <- inla.spde.make.A(mesh = mesh, loc = coords_loc, repl = df$site_id)
idx_w2_s2   <- inla.spde.make.index("w2", n.spde = spde2_2s$n.spde, n.repl = n_sites)

effects_binom_s2 <- df_binom_s2 %>%
  select(intercept_p2, Diam_sc,
         Sap_max2, Sap_max3, Sap_max4, Sap_max5,
         Essence_Conifer, Type_Stump,
         w1_cov)         # <-- w1 comme effet fixe ordinaire

stack_binom_s2 <- inla.stack(
  data    = list(y_binom = df_binom_s2$P_Bux),
  A       = list(1, A_obs_w2_s2),
  effects = list(effects_binom_s2, w2 = idx_w2_s2),
  tag     = "binom_s2"
)

formula_binom_s2 <- y_binom ~ -1 +
  intercept_p2 +
  Diam_sc +
  Sap_max2 + Sap_max3 + Sap_max4 + Sap_max5 +
  Essence_Conifer +
  Type_Stump +
  w1_cov +                                         # β_s séquentiel
  f(w2, model = spde2_2s, replicate = w2.repl)     # champ résiduel

dat_binom_s2 <- inla.stack.data(stack_binom_s2)
A_binom_s2   <- inla.stack.A(stack_binom_s2)

res_binom_s2 <- inla(
  formula_binom_s2,
  family  = "binomial",
  data    = dat_binom_s2,
  Ntrials = rep(1L, nrow(df_binom_s2)),
  control.predictor = list(A = A_binom_s2, compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla      = list(int.strategy = "eb"),
  verbose = FALSE
)

message("Modèle binomial séquentiel — terminé.")

# =============================================================================
# COMPARAISON : β_s joint vs séquentiel
# =============================================================================

# β_s du modèle joint (hyperparamètre "Beta for w1.copy")
beta_s_joint <- r$summary.hyperpar[
  grep("Beta for w1", rownames(r$summary.hyperpar)), ]

# β_s séquentiel (effet fixe w1_cov dans le modèle binomial)
beta_s_seq <- res_binom_s2$summary.fixed["w1_cov", ]

cat("\n")
cat("=================================================================\n")
cat("COMPARAISON β_s : modèle joint vs inférence séquentielle\n")
cat("=================================================================\n")
cat(sprintf("  Joint      : mean = %6.4f | 95%% CI [%6.4f, %6.4f]\n",
            beta_s_joint$mean, beta_s_joint$`0.025quant`, beta_s_joint$`0.975quant`))
cat(sprintf("  Séquentiel : mean = %6.4f | 95%% CI [%6.4f, %6.4f]\n",
            beta_s_seq$mean, beta_s_seq$`0.025quant`, beta_s_seq$`0.975quant`))
cat("=================================================================\n\n")

# Comparaison des autres effets fixes
cat("--- Effets fixes : joint vs séquentiel ---\n")
fixed_joint <- r$summary.fixed
fixed_seq   <- res_binom_s2$summary.fixed

# Termes communs (hors intercepts et w1_cov vs β_s)
common_terms <- intersect(
  rownames(fixed_joint)[!rownames(fixed_joint) %in% c("intercept_p1", "intercept_p2")],
  rownames(fixed_seq)[rownames(fixed_seq) != "intercept_p2"]
)

compare_df <- data.frame(
  Term          = common_terms,
  Joint_mean    = round(fixed_joint[common_terms, "mean"], 4),
  Joint_q025    = round(fixed_joint[common_terms, "0.025quant"], 4),
  Joint_q975    = round(fixed_joint[common_terms, "0.975quant"], 4),
  Seq_mean      = round(fixed_seq[common_terms, "mean"], 4),
  Seq_q025      = round(fixed_seq[common_terms, "0.025quant"], 4),
  Seq_q975      = round(fixed_seq[common_terms, "0.975quant"], 4),
  Delta_mean    = round(fixed_seq[common_terms, "mean"] -
                          fixed_joint[common_terms, "mean"], 4)
)

print(compare_df)

# =============================================================================
# VISUALISATION : comparaison des posterieures
# =============================================================================

# ---- Figure 1 : β_s joint vs séquentiel ------------------------------------
beta_s_marg_joint <- r$marginals.hyperpar[[
  grep("Beta for w1", names(r$marginals.hyperpar))]]
beta_s_df_joint <- as.data.frame(beta_s_marg_joint)
colnames(beta_s_df_joint) <- c("x", "density")
beta_s_df_joint$model <- "Joint"

beta_s_marg_seq <- res_binom_s2$marginals.fixed$w1_cov
beta_s_df_seq   <- as.data.frame(beta_s_marg_seq)
colnames(beta_s_df_seq) <- c("x", "density")
beta_s_df_seq$model <- "Séquentiel"

beta_compare_df <- bind_rows(beta_s_df_joint, beta_s_df_seq)

p_beta_compare <- ggplot(beta_compare_df, aes(x, density, color = model, linetype = model)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Joint" = "#6161d9", "Séquentiel" = "#d96161")) +
  labs(
    title    = "Postérieure de β_s : modèle joint vs inférence séquentielle",
    subtitle = "Convergence des deux distributions → couplage cohérent",
    x        = "β_s",
    y        = "Densité postérieure",
    color    = "Modèle", linetype = "Modèle"
  ) +
  theme_minimal()

print(p_beta_compare)

# ---- Figure 2 : forêt — effets fixes joint vs séquentiel -------------------
plot_terms <- common_terms  # termes partagés entre les deux modèles

fe_joint <- fixed_joint[plot_terms, ] %>%
  tibble::rownames_to_column("Term") %>%
  mutate(model = "Joint")

fe_seq <- fixed_seq[plot_terms, ] %>%
  tibble::rownames_to_column("Term") %>%
  mutate(model = "Séquentiel")

fe_all <- bind_rows(fe_joint, fe_seq) %>%
  mutate(Term = reorder(Term, mean))

p_forest <- ggplot(fe_all, aes(x = mean, y = Term, color = model, shape = model)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbarh(
    aes(xmin = `0.025quant`, xmax = `0.975quant`),
    height = 0.25, position = position_dodge(0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Joint" = "#6161d9", "Séquentiel" = "#d96161")) +
  labs(
    title  = "Effets fixes : modèle joint vs inférence séquentielle",
    x      = "Effet sur logit(P(Buxbaumia))",
    y      = "",
    color  = "Modèle", shape = "Modèle"
  ) +
  theme_minimal()

print(p_forest)

# =============================================================================
# BONUS : variabilité du β_s selon les réalisations de w1
# =============================================================================
# Au lieu de la moyenne des tirages, on peut estimer β_s pour chaque réalisation
# de w1 séparément (sans champ résiduel w2, pour la rapidité), et visualiser
# la distribution de β_s induite par l'incertitude sur w1.
#
# Cela mesure si la propagation d'incertitude du modèle joint est bien capturée.

message("=== Variabilité de β_s selon les réalisations de w1 ===")
message("(Modèle simplifié sans w2 pour la rapidité)")

# Sous-ensemble de tirages pour la rapidité
results_per_sample <- lapply(seq_len(N_sub), function(j) {
  w1_j <- W1_mat[, j]
  
  df_j <- df %>%
    mutate(
      Sap_max2        = as.integer(Sap_max == "2"),
      Sap_max3        = as.integer(Sap_max == "3"),
      Sap_max4        = as.integer(Sap_max == "4"),
      Sap_max5        = as.integer(Sap_max == "5"),
      Essence_Conifer = as.integer(Essence == "Résineux"),
      Type_Stump      = as.integer(Type    == "Souche"),
      intercept_p2    = 1,
      w1_cov          = w1_j
    )
  
  res_j <- tryCatch(
    inla(
      P_Bux ~ -1 + intercept_p2 + Diam_sc +
        Sap_max2 + Sap_max3 + Sap_max4 + Sap_max5 +
        Essence_Conifer + Type_Stump + w1_cov,
      family  = "binomial",
      data    = df_j,
      Ntrials = rep(1L, nrow(df_j)),
      control.inla    = list(int.strategy = "eb"),
      control.compute = list(dic = FALSE),
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(res_j)) {
    fx <- res_j$summary.fixed["w1_cov", ]
    # P(β_s > 0) depuis la marginale
    marg <- res_j$marginals.fixed$w1_cov
    p_pos <- 1 - inla.pmarginal(0, marg)
    data.frame(
      tirage  = j,
      mean    = fx$mean,
      q025    = fx$`0.025quant`,
      q975    = fx$`0.975quant`,
      p_pos   = p_pos,                        # P(β_s > 0)
      sign    = fx$`0.025quant` > 0 |         # IC exclut zéro ?
        fx$`0.975quant` < 0
    )
  } else {
    data.frame(tirage = j, mean = NA, q025 = NA, 
               q975 = NA, p_pos = NA, sign = NA)
  }
})

results_df <- do.call(rbind, results_per_sample)
print(results_df)

# Résumé
cat(sprintf("\nProportion de tirages où IC95 exclut zéro : %d / %d\n",
            sum(results_df$sign, na.rm = TRUE), N_sub))
cat(sprintf("P(β_s > 0) médiane sur les tirages : %.3f\n",
            median(results_df$p_pos, na.rm = TRUE)))

p_beta_dist <- ggplot(data.frame(beta_s = beta_s_per_sample[!is.na(beta_s_per_sample)]),
                      aes(x = beta_s)) +
  geom_histogram(fill = "#d96161", color = "white", bins = 10, alpha = 0.8) +
  geom_vline(
    xintercept = beta_s_joint$mean,
    color = "#6161d9", linewidth = 1.2, linetype = "dashed"
  ) +
  annotate("text",
           x = beta_s_joint$mean, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("β_s joint = %.4f", beta_s_joint$mean),
           color = "#6161d9", size = 3.5) +
  labs(
    title    = sprintf("Distribution de β_s sur %d réalisations de w1 (sans w2)", N_sub),
    subtitle = "Chaque barre = un tirage de la postérieure de w1",
    x        = "β_s (par tirage)",
    y        = "Nombre de tirages"
  ) +
  theme_minimal()

print(p_beta_dist)

message("=== Validation deux étapes terminée ===")
