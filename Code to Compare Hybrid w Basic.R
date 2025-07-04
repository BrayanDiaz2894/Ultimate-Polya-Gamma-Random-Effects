#------------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
# Only alpha is sampled; beta and V_alpha are fixed to the true population values.
#Version with two augmentations latent and pg.
#------------------------------------------------------------------

rm(list = ls())  # Remove all variables
gc()             # Free memory

start_time <- Sys.time()

# I. Load required libraries
library(BayesLogit)  # For Polya-Gamma sampling (rpg)
library(MASS)        # For mvrnorm
library(ggplot2)
library(MCMCpack)    # For riwish() (not used now)
library(coda)
library(truncnorm)   # For truncated normal sampling
library(coda)
library(ggplot2)
library(tidyr)   
library(tidyverse)

#Pick a random individual
pick_random_real <- function(n = 1, min = 0, max = 100) {
  runif(n, min = min, max = max)
}
set.seed(123)
rnumb <- ceiling(pick_random_real())  # un número real


# II. Set seed for reproducibility
set.seed(123)

# III. Set parameters
n <- 100       # Number of individuals
t_obs <- 30    # Observations per individual
p <- 3         # Number of predictors (fixed effects)
q <- 2         # Number of individual-level covariates (random effects)
reg <- 0       # Regularization parameter (if needed)

# IV. Simulate data
# i. Basic data:
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))  # Predictor array (n x t_obs x p)
X[, , 1] <- 1   # Force intercept = 1
Z <- matrix(rnorm(n * q), nrow = n, ncol = q)          # Individual-level covariates (n x q)

# ii. Additional elements:
V_alpha_true <- matrix(c(1, 0.5,
                         0.5, 3), nrow = q, ncol = q)    # True covariance for random effects
alpha_true <- mvrnorm(n, mu = rep(0, q), Sigma = V_alpha_true)  # True random effects (n x q)
beta_true <- c(-1, 2, -1.5)         

load("Output/V_sampleslatentvar.RData")
load("Output/V_samples.RData") 
load("Output/V_samplesupg.RData") 
#load("V_sampleslatentvarrealvlpha.RData") 
#load("V_sampleslatentvarrealalpha.RData") 
#load("V_samplesnlatentvarrealalpha.RData") 
#load("V_samplesauglatentvarrealvlpha.RData") 

#load("V_samplesnlatentvarrealalphaeqsize.RData") 
#load("V_sampleslatentvarrealalphaeqsize.RData") 



# Load required package

# Define the density_overlay function with a maximum of 4 input vectors
density_overlay <- function(..., 
                            vertical_line = NULL, 
                            vline_color = "black", 
                            vline_linetype = "dashed", 
                            plot_title = "Overlay Density Plot", 
                            xlab = "Value", 
                            ylab = "Density",
                            labels = NULL) {
  # Capture the list of vectors provided via ellipsis
  vec_list <- list(...)
  
  # Check that at least two vectors were provided
  if (length(vec_list) < 2) {
    stop("Please provide at least two vectors for comparison.")
  }
  
  # Limit the number of vectors to 4.
  if (length(vec_list) > 4) {
    warning("More than 4 vectors provided. Only the first 4 will be used.")
    vec_list <- vec_list[1:4]
  }
  
  # Determine names/labels for each vector
  if (!is.null(labels)) {
    if (length(labels) != length(vec_list)) {
      warning("Length of labels is not equal to the number of vectors provided. Using default names.")
      vec_names <- names(vec_list)
      if (is.null(vec_names)) {
        vec_names <- paste("Vector", seq_along(vec_list))
      }
    } else {
      vec_names <- labels
    }
  } else {
    vec_names <- names(vec_list)
    if (is.null(vec_names)) {
      vec_names <- paste("Vector", seq_along(vec_list))
    }
  }
  
  # Combine the vectors into one data frame.
  df <- data.frame(value = numeric(0), vector = character(0), stringsAsFactors = FALSE)
  for (i in seq_along(vec_list)) {
    tmp <- data.frame(value = vec_list[[i]], 
                      vector = rep(vec_names[i], length(vec_list[[i]])),
                      stringsAsFactors = FALSE)
    df <- rbind(df, tmp)
  }
  
  # Create the density plot with overlaid densities
  p <- ggplot(df, aes(x = value, color = vector)) +
    geom_density(size = 1) +
    labs(title = plot_title, x = xlab, y = ylab, color = "Vector") +
    theme_minimal()
  
  # Add vertical line(s) if the 'vertical_line' parameter is provided
  if (!is.null(vertical_line)) {
    for (v in vertical_line) {
      p <- p + geom_vline(xintercept = v, color = vline_color, 
                          linetype = vline_linetype, size = 1)
    }
  }
  
  return(p)
}


#### Some basic Statistics. 

#SE and posterior mean 

# Función para resumir una cadena MCMC
# Requiere el paquete coda
# Función para resumir una cadena MCMC e incluir prueba de Geweke
# Requiere el paquete coda
summary_mcmc <- function(draws) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Please install the 'coda' package to use this function.")
  }
  
  # Convertir a mcmc y extraer vector numérico (si fuera multivariante, toma la primera columna)
  draws_mcmc <- coda::as.mcmc(draws)
  draws_vec  <- as.numeric(draws_mcmc)
  
  # Estadísticos básicos
  m   <- mean(draws_vec)
  mn  <- min(draws_vec)
  mx  <- max(draws_vec)
  
  # Effective sample size
  n_eff <- coda::effectiveSize(draws_mcmc)
  # Desviación estándar de la cadena
  sd_chain <- sd(draws_vec)
  # Monte Carlo standard error
  mcse <- sd_chain / sqrt(n_eff)
  
  # Geweke diagnostic (z-score)
  gw <- coda::geweke.diag(draws_mcmc)
  z_geweke <- as.numeric(gw$z)[1]
  
  # Devolver lista con resultados
  list(
    mean        = m,
    minimum     = mn,
    maximum     = mx,
    mcse        = mcse,
    geweke_z    = z_geweke
  )
}
##############
#############

summary_mcmc(beta_1_samplesnaug)
summary_mcmc(beta_2_samplesnaug)
summary_mcmc(beta_3_samplesnaug)
summary_mcmc(V1_naug)
summary_mcmc(V2_naug)
summary_mcmc(V21_naug)
summary_mcmc(alpha_samples_naug[, rnumb, 1])
summary_mcmc(alpha_samples_naug[, rnumb, 2])

##############
#############

summary_mcmc(beta_1_samplesaug)
summary_mcmc(beta_2_samplesaug)
summary_mcmc(beta_3_samplesaug)
summary_mcmc(V1_aug)
summary_mcmc(V2_aug)
summary_mcmc(V21_aug)
summary_mcmc(alpha_samples_aug[, rnumb, 1])
summary_mcmc(alpha_samples_aug[, rnumb, 2])


#### Plots for beta. 

#1. Compare Behaviors. 

density_plot <- density_overlay(vec1 = beta_1_samplesnaug, 
                                vec2 = beta_1_samplesaug,
                                vec3 = beta_1_samplesupg,
                                vertical_line = beta_true[1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent", "UPG"))

print(density_plot)

density_plot <- density_overlay(vec1 = beta_2_samplesnaug, 
                                vec2 = beta_2_samplesaug, 
                                vertical_line = beta_true[2],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

density_plot <- density_overlay(vec1 = beta_3_samplesnaug, 
                                vec2 = beta_3_samplesaug, 
                                vertical_line = beta_true[3],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

#2. Get the Computational Efficiency. 

#ESS. 
beta_1_samplesnaug <- as.mcmc(beta_1_samplesnaug)
essb1naug <- effectiveSize(beta_1_samplesnaug)/as.numeric(timesamples_naug, units = "secs")
beta_1_samplesaug <- as.mcmc(beta_1_samplesaug)
essb1aug <- effectiveSize(beta_1_samplesaug)/as.numeric(timesamples_aug, units = "secs")

beta_2_samplesnaug <- as.mcmc(beta_2_samplesnaug)
essb2naug <- effectiveSize(beta_2_samplesnaug)/as.numeric(timesamples_naug, units = "secs")
beta_2_samplesaug <- as.mcmc(beta_2_samplesaug)
essb2aug <- effectiveSize(beta_2_samplesaug)/as.numeric(timesamples_aug, units = "secs")

beta_3_samplesnaug <- as.mcmc(beta_3_samplesnaug)
essb3naug <- effectiveSize(beta_3_samplesnaug)/as.numeric(timesamples_naug, units = "secs")
beta_3_samplesaug <- as.mcmc(beta_3_samplesaug)
essb3aug <- effectiveSize(beta_3_samplesaug)/as.numeric(timesamples_aug, units = "secs")

# build a data frame with correct factor levels
ess_df <- data.frame(
  beta   = rep(paste0("Beta", 1:3), each = 2),
  method = factor(
    rep(c("Polya Gamma Augmentation", "Polya Gamma + Latent Variable Augmentation"), times = 3),
    levels = c("Polya Gamma Augmentation", "Polya Gamma + Latent Variable Augmentation")  # legend order
  ),
  ess = c(
    essb1naug, essb1aug,
    essb2naug, essb2aug,
    essb3naug, essb3aug
  )
)

# grouped bar chart with legend control
ggplot(ess_df, aes(x = beta, y = ess, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = c(
      "Polya Gamma Augmentation" = "blue",
      "Polya Gamma + Latent Variable Augmentation" = "red"
    )
  ) +
  scale_x_discrete(labels = c(
    Beta1 = expression(beta[1]),
    Beta2 = expression(beta[2]),
    Beta3 = expression(beta[3])
  )) +
  labs(
    x    = NULL,
    y    = "ESS per second",
    fill = "Sampler"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


ggsave(filename = "Output/Results/essbeta.pdf", width = 8, height = 4, dpi = 300)

#3. Tracepots.

plot_beta_trace <- function(idx, k = NULL) {
  # build the object names
  name_naug <- paste0("beta_", idx, "_samplesnaug")
  name_aug  <- paste0("beta_", idx, "_samplesaug")
  
  # grab them from the environment
  b_naug <- get(name_naug)
  b_aug  <- get(name_aug)
  
  # optionally restrict to first k draws
  if (!is.null(k)) {
    b_naug <- head(b_naug, k)
    b_aug  <- head(b_aug,  k)
  }
  
  # build a data.frame
  df <- data.frame(
    iteration = c(seq_along(b_naug), seq_along(b_aug)),
    value     = c(b_naug, b_aug),
    sampler   = factor(
      c(
        rep("Polya Gamma Augmentation", length(b_naug)),
        rep("Polya Gamma + Latent Variable", length(b_aug))
      ),
      levels = c("Polya Gamma Augmentation", "Polya Gamma + Latent Variable")
    )
  )
  
  # plot with horizontal line for true value
  ggplot(df, aes(x = iteration, y = value, color = sampler)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = beta_true[idx], 
               linetype = "dashed", 
               color = "black", 
               size = 0.9) +
    scale_color_manual(
      values = c(
        "Polya Gamma Augmentation"         = "blue", 
        "Polya Gamma + Latent Variable"    = "red"
      )
    ) +
    labs(
      x     = "Iteration",
      y     = bquote(beta[.(idx)]),
      color = "Sampler"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


plot_beta_trace(1, k = 1000) 
ggsave(filename = "Output/Results/tracebeta1.pdf", width = 8, height = 4, dpi = 300)
plot_beta_trace(2, k = 1000) 
ggsave(filename = "Output/Results/tracebeta2.pdf", width = 8, height = 4, dpi = 300)
plot_beta_trace(3, k = 1000) 
ggsave(filename = "Output/Results/tracebeta3.pdf", width = 8, height = 4, dpi = 300)

#### Plots for V_alpha

#1. Compare Behaviors.
density_plot <- density_overlay(vec1 = V1_naug,
                                vec2 = V1_aug,
                                #vec3 = V1_naugreal,
                                #vec3 = V1_augreal,
                                vertical_line = V_alpha_true[1,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

density_plot <- density_overlay(vec1 = V2_naug,
                                vec2 = V2_aug,
                                #vec3 = V1_naugreal,
                                vertical_line = V_alpha_true[2,2],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

density_plot <- density_overlay(vec1 = V21_naug,
                                vec2 = V21_aug,
                                #vec3 = V1_naugreal,
                                vertical_line = V_alpha_true[2,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

#2. Get the Computational Efficiency.

#ESS.
V1_naug <- as.mcmc(V1_naug)
essV1_naug <- effectiveSize(V1_naug)/as.numeric(timesamples_naug, units = "secs")
V1_aug <- as.mcmc(V1_aug)
essV1_aug <- effectiveSize(V1_aug)/as.numeric(timesamples_aug, units = "secs")

V2_naug <- as.mcmc(V2_naug)
essV2_naug <- effectiveSize(V2_naug)/as.numeric(timesamples_naug, units = "secs")
V2_aug <- as.mcmc(V2_aug)
essV2_aug <- effectiveSize(V2_aug)/as.numeric(timesamples_aug, units = "secs")

V21_naug <- as.mcmc(V21_naug)
essV21_naug <- effectiveSize(V21_naug)/as.numeric(timesamples_naug, units = "secs")
V21_aug <- as.mcmc(V21_aug)
essV21_aug <- effectiveSize(V21_aug)/as.numeric(timesamples_aug, units = "secs")

# 1) Build the data frame
essV_df <- data.frame(
  param  = rep(c("V11", "V22", "V21"), each = 2),
  method = rep(
    c("Polya Gamma Augmentation",
      "Polya Gamma + Latent Variable Augmentation"),
    times = 3
  ),
  ess    = c(
    essV1_naug,  essV1_aug,
    essV2_naug,  essV2_aug,
    essV21_naug, essV21_aug
  )
)

# 2) Plot grouped bar chart
# Asegurar niveles en el orden deseado para la leyenda
essV_df$method <- factor(
  essV_df$method,
  levels = c("Polya Gamma Augmentation", "Polya Gamma + Latent Variable Augmentation")
)

# Gráfico
ggplot(essV_df, aes(x = param, y = ess, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = c(
      "Polya Gamma Augmentation" = "blue",
      "Polya Gamma + Latent Variable Augmentation" = "red"
    )
  ) +
  scale_x_discrete(labels = c(
    V11 = expression(V[alpha][1]),
    V22 = expression(V[alpha][2]),
    V21 = expression(V[alpha][2][1])
  )) +
  labs(
    x    = NULL,
    y    = "ESS per second",
    fill = "Sampler"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = "Output/Results/essvalpha.pdf", width = 8, height = 4, dpi = 300)

#3) traceplots

plot_V_trace <- function(idx, k = NULL) {
  # map idx→object names
  naug_names <- c("V1_naug", "V2_naug", "V21_naug")
  aug_names  <- c("V1_aug",  "V2_aug",  "V21_aug")
  
  # grab the correct vectors
  v_naug <- as.numeric(get(naug_names[idx]))
  v_aug  <- as.numeric(get(aug_names[idx]))
  
  # optional head()
  if (!is.null(k)) {
    v_naug <- head(v_naug, k)
    v_aug  <- head(v_aug,  k)
  }
  
  # assemble data
  df <- data.frame(
    iteration = c(seq_along(v_naug), seq_along(v_aug)),
    value     = c(v_naug, v_aug),
    sampler   = factor(
      c(rep("Polya Gamma Augmentation", length(v_naug)),
        rep("Polya Gamma + Latent Variable", length(v_aug))),
      levels = c("Polya Gamma Augmentation","Polya Gamma + Latent Variable")
    )
  )
  
  # nice y-labels for each entry of V_alpha
  ylabels <- list(
    expression(V[alpha][1]),       # idx=1
    expression(V[alpha][2]),       # idx=2
    expression(V[alpha][2][1])     # idx=3 (off-diagonal)
  )
  
  # true population value from V_alpha_true
  pop_val <- switch(idx,
                    V_alpha_true[1, 1],  # idx = 1
                    V_alpha_true[2, 2],  # idx = 2
                    V_alpha_true[2, 1]   # idx = 3
  )
  
  ggplot(df, aes(x = iteration, y = value, color = sampler)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = pop_val,
               linetype = "dashed",
               color = "black",
               size = 0.9) +
    scale_color_manual(
      values = c(
        "Polya Gamma Augmentation"      = "blue",
        "Polya Gamma + Latent Variable" = "red"
      )
    ) +
    labs(
      x     = "Iteration",
      y     = ylabels[[idx]],
      color = "Sampler"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


# Examples:
plot_V_trace(1, k = 1000)         # V[alpha][1,1]
ggsave(filename = "Output/Results/tracevalpha1.pdf", width = 8, height = 4, dpi = 300)
plot_V_trace(2, k = 1000) # first 200 draws for V[alpha][2,2]
ggsave(filename = "Output/Results/tracevalpha2.pdf", width = 8, height = 4, dpi = 300)
plot_V_trace(3, k = 1000)      # the 'third' V_alpha entry (correlation)
ggsave(filename = "Output/Results/tracevalpha21.pdf", width = 8, height = 4, dpi = 300)



#### Plots for alpha

#1. Compare Behaviors.

density_plot <- density_overlay(vec1 = a1_nag[,1],
                                vec2 = a1_ag[,1],
                                vertical_line = alpha_true[1,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

density_plot <- density_overlay(vec1 = a2_nag[,2],
                                vec2 = a2_ag[,2],
                                vertical_line = alpha_true[2,2],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent"))

print(density_plot)

#2. Get the Computational Efficiency.

#ESS.
test <- as.matrix(alpha_samples_aug[,1,1])
dim(test)

#Extract the ESS for all alphas.
n <- dim(alpha_samples_aug)[2]          # nº de individuos

## --- 2 · data.frame vacío ----------------------------------------------------
df_ess <- data.frame(
  ind            = integer(n),
  ess_alpha1naug = numeric(n),
  ess_alpha1aug  = numeric(n),
  ess_alpha2naug = numeric(n),
  ess_alpha2aug  = numeric(n),
  meanalpha1_max = numeric(n),
  meanalpha2_max = numeric(n)
)

## --- 3 · bucle ----------------------------------------------------------------
for (j in seq_len(n)) {
  
  ## ··· cadenas del individuo j ···
  chain1aug  <- alpha_samples_aug[ , j, 1]
  chain2aug  <- alpha_samples_aug[ , j, 2]
  chain1naug <- alpha_samples_naug[, j, 1]
  chain2naug <- alpha_samples_naug[, j, 2]
  
  ## ··· ESS por segundo ···
  ess1aug <- effectiveSize(as.mcmc(chain1aug))  /
    as.numeric(timesamples_aug,  units = "secs")
  ess2aug <- effectiveSize(as.mcmc(chain2aug))  /
    as.numeric(timesamples_aug,  units = "secs")
  ess1nau <- effectiveSize(as.mcmc(chain1naug)) /
    as.numeric(timesamples_naug, units = "secs")
  ess2nau <- effectiveSize(as.mcmc(chain2naug)) /
    as.numeric(timesamples_naug, units = "secs")
  
  ## ··· medias ···
  m1aug <- mean(chain1aug);  m2aug <- mean(chain2aug)
  m1nau <- mean(chain1naug); m2nau <- mean(chain2naug)
  
  ## ··· rellenar fila j ···
  df_ess[j, ] <- list(
    ind            = j,
    ess_alpha1naug = ess1nau,
    ess_alpha1aug  = ess1aug,
    ess_alpha2naug = ess2nau,
    ess_alpha2aug  = ess2aug,
    meanalpha1_max = max(m1aug, m1nau),
    meanalpha2_max = max(m2aug, m2nau)
  )
}

plot_ess <- function(df, alpha = 1, n = nrow(df)) {
  stopifnot(alpha %in% c(1, 2))
  stopifnot(n >= 1, n <= nrow(df))
  
  # Pick relevant columns ------------------------------------------------------
  col_naug <- paste0("ess_alpha", alpha, "naug")
  col_aug  <- paste0("ess_alpha", alpha, "aug")
  col_mean <- paste0("meanalpha", alpha, "_max")
  
  # Sub-data frame ordered by the chosen posterior mean ------------------------
  df_sub <- df[order(df[[col_mean]]), ][1:n, ]
  df_sub$ind <- seq_len(n)  # reset x-axis index after ordering
  
  # ---------------------------------------------------------------------------
  ggplot(df_sub, aes(x = ind)) +
    geom_col(aes(y = .data[[col_naug]],
                 fill = "Polya Gamma Augmentation"),
             alpha = 0.8) +
    geom_col(aes(y = .data[[col_aug]],
                 fill = "PG + Latent Utility Augmentation"),
             alpha = 0.5) +
    scale_fill_manual(
      name   = "Sampler",
      values = c("Polya Gamma Augmentation"           = "blue",
                 "PG + Latent Utility Augmentation" = "red")
    ) +
    scale_x_continuous(breaks = seq(0, n, by = max(1, floor(n/10)))) +
    labs(x = "Index ordered by posterior mean",
         y = "ESS / s") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_ess_diff <- function(df, alpha = 1, n = nrow(df)) {
  stopifnot(alpha %in% c(1, 2))
  stopifnot(n >= 1, n <= nrow(df))
  
  ## --- columnas que dependen del componente α -------------------------------
  col_naug <- paste0("ess_alpha", alpha, "naug")
  col_aug  <- paste0("ess_alpha", alpha, "aug")
  col_mean <- paste0("meanalpha", alpha, "_max")
  
  ## --- sub-data.frame ordenado por la media posterior -----------------------
  df_sub <- df[order(df[[col_mean]]), ][1:n, ]
  df_sub$ind  <- seq_len(n)                         # re-indexar eje x
  df_sub$diff <- abs(df_sub[[col_aug]] - df_sub[[col_naug]])  # Δ ESS/s
  
  ## --- gráfico --------------------------------------------------------------
  ggplot(df_sub, aes(x = ind, y = diff,
                     fill = diff > 0)) +             # TRUE ⇒ aug mejor
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(
      name   = "",
      values = c("TRUE"  = "lightcoral",   # PG + utility > PG
                 "FALSE" = "blue"), # PG ≥ PG + utility
      labels = c("TRUE"  = "",
                 "FALSE" = "")
    ) +
    scale_x_continuous(breaks = seq(0, n, by = max(1, floor(n/10)))) +
    labs(x = "Index ordered by posterior mean",
         y = "Absolute ValueΔ ESS per second") +
    theme_minimal() +
    theme(legend.position = "bottom")
}


plot_ess(df_ess, alpha = 1, n = 100)
ggsave("Output/Results/ess_alpha1.pdf", width = 8, height = 4, dpi = 300)

plot_ess(df_ess, alpha = 2, n = 100)
ggsave("Output/Results/ess_alpha2.pdf", width = 8, height = 4, dpi = 300)

plot_ess_diff(df_ess, alpha = 1, n = 100)
ggsave("Output/Results/ess_alpha1diff.pdf", width = 8, height = 4, dpi = 300)
plot_ess_diff(df_ess, alpha = 2, n = 100)
ggsave("Output/Results/ess_alpha2diff.pdf", width = 8, height = 4, dpi = 300)


#3) Traceplots.


plot_alpha_trace <- function(rnumb, chain_num = 1, n_iter = NULL) {
  # Validación del argumento
  if (!(chain_num %in% c(1, 2))) {
    stop("El argumento 'chain_num' debe ser 1 o 2.")
  }
  
  # Seleccionar las cadenas correctas
  chain_aug  <- alpha_samples_aug[, rnumb, chain_num]
  chain_naug <- alpha_samples_naug[, rnumb, chain_num]
  
  # Limitar número de iteraciones si se especifica
  if (!is.null(n_iter)) {
    chain_aug  <- head(chain_aug, n_iter)
    chain_naug <- head(chain_naug, n_iter)
  }
  
  # Valor verdadero de α
  alpha_true_val <- alpha_true[rnumb, chain_num]
  
  # Crear dataframe
  df_trace <- tibble(
    Iteration     = 1:length(chain_aug),
    `Polya Gamma Augmentation`         = chain_naug,
    `Polya Gamma + Latent Variable`    = chain_aug
  )
  
  # Formato largo
  df_long <- df_trace %>%
    pivot_longer(cols = -Iteration, names_to = "Sampler", values_to = "Value") %>%
    mutate(Sampler = factor(Sampler, levels = c(
      "Polya Gamma Augmentation", "Polya Gamma + Latent Variable"
    )))
  
  # Graficar
  ggplot(df_long, aes(x = Iteration, y = Value, color = Sampler)) +
    geom_line(alpha = 0.8) +
    geom_hline(yintercept = alpha_true_val, 
               linetype = "dashed", color = "black", size = 0.9) +
    scale_color_manual(
      values = c(
        "Polya Gamma Augmentation"      = "blue",
        "Polya Gamma + Latent Variable" = "red"
      )
    ) +
    labs(
      title = paste(""),
      y     = bquote(alpha[.(chain_num)]),
      x     = "Iteration",
      color = "Sampler"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_alpha_trace(rnumb,chain_num = 1, n_iter = 1000)
ggsave(filename = "Output/Results/alpharandomtrace.pdf", width = 8, height = 4, dpi = 300)

