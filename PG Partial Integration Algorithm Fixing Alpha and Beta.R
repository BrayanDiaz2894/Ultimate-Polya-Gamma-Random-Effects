#-------------------------------------------------------------
# Gibbs Sampler – betas y alphas fijos en sus valores “true”
#-------------------------------------------------------------

library(BayesLogit)   # rpg()
library(MASS)         # mvrnorm()
library(MCMCpack)     # riwish()
library(truncnorm)    # rtruncnorm()

set.seed(123)

## --- 1. Dimensiones y datos simulados -----------------------
n      <- 500
t_obs  <- 30
p      <- 3
q      <- 2

# Diseño X
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))
X[ , , 1] <- 1

# Covariables Z
Z <- matrix(rnorm(n * q), n, q)

# Valores poblacionales “verdaderos”
V_alpha_true <- matrix(c(1, .5,
                         .5, 3), q, q)
alpha_true <- mvrnorm(n, rep(0, q), V_alpha_true)
beta_true  <- c(-1, 2, -1.5)

# η y  Y
eta <- matrix(0, n, t_obs)
for (i in 1:n)
  eta[i, ] <- X[i, , ] %*% beta_true + as.numeric(Z[i, ] %*% alpha_true[i, ])

prob <- 1 / (1 + exp(-eta))
Y    <- matrix(rbinom(n * t_obs, 1, prob), n, t_obs)

## --- 2. Priors ---------------------------------------------
mu0        <- rep(0, p)          # ya no se usan (β fijo) pero los dejamos
Sigma0     <- diag(p) * 10
Sigma0_inv <- solve(Sigma0)
G0         <- 10                 # prior de γ
nu0        <- q + 2
Lambda0    <- diag(q)

## --- 3. Valores iniciales (β y α fijos) ---------------------
Beta    <- beta_true             # <<– FIJO
Alpha   <- alpha_true            # <<– FIJO
V_alpha <- diag(q)
gamma   <- 0

##  Matriz X “apilada” (útil para ω, γ)
N_total <- n * t_obs
X_all   <- matrix(NA, N_total, p)
for (i in 1:n)
  X_all[((i-1)*t_obs + 1):(i*t_obs), ] <- X[i, , ]

## --- 4. Almacenamiento --------------------------------------
iters   <- 2000
burnin  <- floor(.1 * iters)

Beta_save    <- matrix(rep(beta_true, each = iters), iters, p)  # copia fija
V_alpha_save <- vector("list", iters)
gamma_save   <- numeric(iters)

## --- 5. Gibbs loop ------------------------------------------
for (it in 1:iters) {
  print(it)
  ## 5.1  Latent utilities h_ij
  eta_cur <- matrix(0, n, t_obs)
  for (i in 1:n)
    eta_cur[i, ] <- X[i, , ] %*% Beta + as.numeric(Z[i, ] %*% Alpha[i, ])
  
  lambda_mat <- exp(eta_cur)
  pi_mat     <- lambda_mat / (1 + lambda_mat)
  
  U_mat <- matrix(runif(N_total, .Machine$double.eps,
                        1 - .Machine$double.eps), n, t_obs)
  V_mat <- Y + U_mat * (1 - Y - pi_mat)
  eps_mat <- log( V_mat / (1 - V_mat) )
  h_mat   <- eta_cur + eps_mat
  
  ## 5.2  Polya–Gamma ω_ij
  diff_vec   <- as.vector(abs(h_mat - eta_cur))
  omega_vec  <- rpg(N_total, 2, diff_vec)
  omega_mat  <- matrix(omega_vec, n, t_obs)
  
  ## 5.3  Paso de γ  (iMDA)
  gamma_tilde <- rnorm(1, 0, sqrt(G0))
  h_tilde     <- h_mat + gamma_tilde
  L_bound <- max(h_tilde[Y == 0]); if (is.infinite(L_bound)) L_bound <- -Inf
  U_bound <- min(h_tilde[Y == 1]); if (is.infinite(U_bound)) U_bound <-  Inf
  
  # términos A y B
  mu_beta1_star <- sum(omega_vec * h_tilde) + as.numeric(Sigma0_inv %*% mu0)
  mu_beta2_star <- crossprod(X_all, omega_vec)
  Sigma_beta_star <- crossprod(X_all * omega_vec, X_all) + Sigma0_inv
  
  A_num <- sum(omega_vec * h_tilde) +
    t(mu_beta2_star) %*% solve(Sigma_beta_star) %*% mu_beta1_star
  B_den <- sum(omega_vec) -
    t(mu_beta2_star) %*% solve(Sigma_beta_star) %*% mu_beta2_star +
    1 / G0
  
  mu_gamma  <- as.numeric(A_num / B_den)
  var_gamma <- 1 / as.numeric(B_den)
  
  gamma <- rtruncnorm(1, a = L_bound, b = U_bound,
                      mean = mu_gamma, sd = sqrt(var_gamma))
  
  h_mat <- h_tilde - gamma   # recentrado
  
  ## 5.4  V_alpha | Alpha  (Alpha fijo)
  S_post   <- Lambda0 + t(Alpha) %*% Alpha
  V_alpha  <- riwish(nu0 + n, S_post)
  
  ## 5.5  Guardar
  V_alpha_save[[it]] <- V_alpha
  gamma_save[it]     <- gamma
}

## --- 6. Posterior después de burn-in ------------------------
V_alpha_samps <- V_alpha_save[(burnin+1):iters]
gamma_samps   <- gamma_save[(burnin+1):iters]

cat("Posterior mean Beta  (fijo):\n"); print(beta_true)
cat("Posterior mean V_alpha:\n")
print( Reduce("+", V_alpha_samps) / length(V_alpha_samps) )

##IT should converges to
denom <- (q + 2 + n) - q - 1        # 101
S_post <- Lambda0 + t(alpha_true) %*% alpha_true
post_mean_theory <- S_post / denom
post_mean_theory

