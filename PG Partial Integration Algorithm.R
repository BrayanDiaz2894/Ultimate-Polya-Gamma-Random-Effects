#-------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
#-------------------------------------------------------------

# I. Load required libraries
library(BayesLogit)   # For Polya-Gamma sampling (rpg)
library(MASS)         # For multivariate normal (mvrnorm)
library(MCMCpack)     # For inverse-Wishart sampling (riwish)
library(truncnorm)    # For truncated normal sampling
library(coda)         # For MCMC diagnostics (if needed)
library(ggplot2)      # For plotting (if needed)

# II. Set seed for reproducibility
set.seed(123)

# III. Set simulation parameters
n      <- 100  # Number of individuals
t_obs  <- 30    # Observations per individual
p      <- 3     # Number of fixed-effect predictors (including intercept)
q      <- 2     # Number of random-effect covariates (random effects dimension)

# IV. Simulate data
# 1. Simulate predictor matrix X (n x t_obs x p) and ensure first column is intercept = 1
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))
X[,,1] <- 1  # first predictor is intercept (set to 1 for all observations)

# 2. Simulate individual-level covariates Z for random effects (n x q)
Z <- matrix(rnorm(n * q), nrow = n, ncol = q)

# 3. Set "true" covariance for random effects and simulate true random effects alpha_true

V_alpha_true <- matrix(c(1, 0.5,   # Covariance matrix (q x q) for random effects
                         0.5, 3), nrow = q, ncol = q)
alpha_true <- mvrnorm(n = n, mu = rep(0, q), Sigma = V_alpha_true)  # random effects for each individual (n x q)

# 4. Set "true" fixed effects beta_true
beta_true <- c(-1, 2, -1.5)  # length p = 3 (intercept, X2, X3)

# 5. Generate latent linear predictor eta and binary response Y
eta <- matrix(0, nrow = n, ncol = t_obs)  # linear predictor for each observation
for (i in 1:n) {
  # Calculate eta[i,] = X[i,,] %*% beta_true + Z[i,] %*% alpha_true[i,]
  # (Note: Z[i,] %*% alpha_true[i,] gives a scalar, apply to all obs of individual i)
  eta[i, ] <- X[i,, ] %*% beta_true + as.numeric(Z[i, ] %*% alpha_true[i, ])
}
# Generate binary outcomes Y from logistic model: P(Y=1) = logistic(eta)
prob <- 1 / (1 + exp(-eta))
Y <- matrix(rbinom(n * t_obs, size = 1, prob = prob), nrow = n, ncol = t_obs)

# V. Set up priors and initial values for the Gibbs sampler
# 1. Priors for fixed effects beta ~ N(mu0, Sigma0)
mu0     <- rep(0, p)
Sigma0  <- diag(p)          # relatively diffuse prior (variance 10 on diagonal)
Sigma0_inv <- solve(Sigma0)

# 2. Prior for gamma (global location parameter) ~ N(0, G0)
G0    <- 10                      # variance of gamma's prior (could be tuned as needed)

# 3. Prior for random effects covariance V_alpha ~ Inverse-Wishart(nu0, Lambda0)
nu0     <- q + 3                 # degrees of freedom for IW prior (q+2 ensures finite mean)
Lambda0 <- diag(2, q)               # scale matrix for IW prior (identity for simplicity)

# 4. Initial values for parameters
Beta    <- rep(0, p)             # start fixed effects at 0
Alpha   <- matrix(0, n, q)       # start random effects at 0 for all individuals
V_alpha <- diag(q)               # start covariance as identity matrix
gamma   <- 0                     # start gamma (location expansion parameter) at 0

# 5. Pre-compute X_all (flattened design matrix) for efficiency in Beta updates
N_total <- n * t_obs
X_all   <- matrix(NA, nrow = N_total, ncol = p)
for (i in 1:n) {
  # Fill rows for individual i
  idx_start <- (i-1)*t_obs + 1
  idx_end   <- i * t_obs
  X_all[idx_start:idx_end, ] <- X[i,, ]
}

# VI. Prepare storage for MCMC samples
iterations <- 10000
burnin     <- floor(iterations * 0.1)    # 10% burn-in
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)
gamma_save   <- numeric(iterations)                       # store gamma draws

# VII. Gibbs Sampler loop


matriz <- matrix(data = 1:(3 * iterations), nrow = iterations, ncol = 3)


#iter <- 2
for (iter in 1:iterations) {
  
  print(iter)
  ## ---------------------------------------------------------------
  ## STEP 1 : draw latent utilities h_{ij}
  ## ---------------------------------------------------------------
  # 1. linear predictor η_ij
  ## ----------  STEP (Z)  : sample latent utilities  z_ij  -----------------
  ## Linear predictor η_ij  (same as before)
  eta_current <- matrix(0, n, t_obs)
  for (i in 1:n)
    eta_current[i, ] <- X[i, , ] %*% Beta + as.numeric(Z[i, ] %*% Alpha[i, ])
  
  ## λ_ij   and   π_ij  = F(η_ij)
  lambda_mat <- exp(eta_current)                     # λ_ij = exp(η_ij)
  pi_mat     <- lambda_mat / (1 + lambda_mat)        # π_ij = λ /(1+λ) = plogis(η)
  
  ## Draw U_ij  ∼  U(0,1)   (avoid the exact endpoints)
  epsU  <- .Machine$double.eps
  U_mat <- matrix(runif(n * t_obs, epsU, 1 - epsU), n, t_obs)
  
  ## Construct V_ij  as in the paper:
  ##   V_ij = y_ij + U_ij * ( 1 - y_ij - π_ij )
  V_mat <- Y + U_mat * (1 - Y - pi_mat)
  
  ## In words:
  ##   y = 0  →  V ∼ U(0 , 1-π)
  ##   y = 1  →  V ∼ U(1-π , 1)
  
  ## Logistic inverse-CDF  ε_ij  and latent utility  h_ij
  eps_mat <- log( V_mat / (1 - V_mat) )          # F^{-1}_ε (V)
  h_mat   <- eta_current + eps_mat               # h_ij = η_ij + ε_ij
  # ------------------------------------------------------------------------
  
  ## (Optional sanity check)
  all( h_mat[Y == 1]  > 0 - 1e-12 )   # should be TRUE
  all( h_mat[Y == 0] <= 0 + 1e-12 )   # should be TRUE
  
  
  
  # 2. **Sample Polya-Gamma ω_{ij} | h, Beta, Alpha** for each observation
  #    ω_{ij} ~ PG(2, |h_{ij} - X_{ij}^T Beta - Z_{ij}^T Alpha_i|).
  #    Note: X_{ij}^T Beta + Z_{ij}^T Alpha_i = eta_current, so use |h - eta|.
  diff_mat <- abs(h_mat - eta_current)
  # Flatten diff_mat to vector and draw PG in one call for efficiency
  diff_vec <- as.vector(diff_mat)
  omega_vec <- rpg(length(diff_vec), 2, diff_vec)   #  h = 2 is correct here
  omega_mat <- matrix(omega_vec, nrow = n, ncol = t_obs)
  
  # 3. **Sample auxiliary location parameter gamma_tilde ~ N(0, G0)**
  gamma_tilde <- rnorm(1, mean = 0, sd = sqrt(G0))
  
  # 4. **Propose shifted utilities:** z_tilde (or h_tilde) = h + gamma_tilde
  h_tilde <- h_mat + gamma_tilde
  
  #    Compute truncation bounds for gamma_new:
  #    L = max_{i,j: Y_{ij}=0} h_tilde_{ij},  U = min_{i,j: Y_{ij}=1} h_tilde_{ij}.
  L_bound <- max(h_tilde[Y == 0])   # if no Y==0, this will be -Inf
  U_bound <- min(h_tilde[Y == 1])   # if no Y==1, this will be Inf
  if (is.infinite(L_bound)) L_bound <- -Inf
  if (is.infinite(U_bound)) U_bound <- Inf
  
  # 5. **Sample gamma_new | ω, h_tilde, Y** from N(g_N, G_N) truncated to [L, U].
  
  # mu_beta1_star = sum_ij omega_ij * h_tilde_ij + Sigma0_inv %*% mu0
  omega_vec <- as.vector(omega_mat)              # vectorizar omega
  mu_beta1_star <- sum(omega_vec * h_tilde) + solve(Sigma0) %*% mu0
  # mu_beta2_star = sum_ij omega_ij * X_ij^T  (X_all' * omega)
  mu_beta2_star <- t(X_all) %*% omega_vec        # X_all es (n*t_obs x p)
  # Sigma_beta_star = sum_ij omega_ij * X_ij X_ij^T + Sigma0_inv
  # Equivalent to: X_all' diag(omega) X_all + Sigma0_inv
  Sigma_beta_star <- t(X_all) %*% (X_all * omega_vec) + solve(Sigma0)
  
  # Linear Term
  A <- sum(omega_vec * h_tilde) +
    t(mu_beta2_star) %*% solve(Sigma_beta_star) %*% mu_beta1_star
  
  # Squared Term
  B <- sum(omega_vec) -
    t(mu_beta2_star) %*% solve(Sigma_beta_star) %*% mu_beta2_star +
    1 / G0
  
  matriz[iter,1] <- sum(omega_vec)
  matriz[iter,2] <-  t(mu_beta2_star) %*% solve(Sigma_beta_star) %*% mu_beta2_star
  matriz[iter,3] <- 1 / G0
  
  # 3. Media y varianza
  mu_gamma     <- as.numeric(A / B)
  var_gamma    <- 1 / as.numeric(B)
  
  gamma_new <- rtruncnorm(1, a = L_bound, b = U_bound, mean = mu_gamma, sd = sqrt(var_gamma))
  
  # 6. **Shift utilities:** define new h = h_tilde - gamma_new (to center latent utilities back around 0 threshold)
  h_mat <- h_tilde - gamma_new
  #    Update gamma to the newly sampled value (store for output)
  gamma <- gamma_new
  
  # 7. **Sample Beta | h, ω, Alpha** from its full conditional (Gaussian).
  ## ---------------------------------------------------------------
  ## Step to sample β     (clear version, non-vectorized)
  ## ---------------------------------------------------------------
  #
  #   We want:    Q_beta  =  Σ_{ij} ω_ij X_ij X_ijᵀ  +  Σ0^{-1}
  #               s       =  Σ_{ij} ω_ij X_ij (h_ij - Z_ijᵀ α_i)  +  Σ0^{-1} μ0
  #               m_beta  =  Q_beta^{-1} s
  #               β ~ N( m_beta ,  Q_beta^{-1} )
  #
  
  # 0. Initialize accumulators
  Q_beta <- matrix(0, nrow = p, ncol = p)      # p × p
  s_vec  <- rep(0, p)                          # p × 1
  
  # 1. Loop over each individual i and observation j
  obs_index <- 0           # to keep track of linear position in X_all and omega_vec
  for (i in 1:n) {
    
    # — random effect for individual i (scalar): Z_i %*% Alpha_i
    zalpha_i <- sum( Z[i, ] * Alpha[i, ] )
    
    for (j in 1:t_obs) {
      obs_index <- obs_index + 1               # 1 … N_total
      
      # Observation (i,j) data
      x_ij     <- X[i, j, ]                    # p-dimensional vector
      omega_ij <- omega_vec[ obs_index ]       # scalar
      h_ij     <- h_mat   [ i, j ]             # scalar
      
      # Residual    r_ij = h_ij - Zα_i
      r_ij <- h_ij - zalpha_i
      
      # Update Σ ω X Xᵀ
      Q_beta <- Q_beta + omega_ij * ( x_ij %*% t(x_ij) )
      
      # Update Σ ω X r
      s_vec  <- s_vec  + omega_ij * x_ij * r_ij
    }
  }
  
  # 2. Add the prior contribution (Σ0^{-1} and Σ0^{-1} μ0)
  Q_beta <- Q_beta + Sigma0_inv
  s_vec  <- s_vec  + as.numeric( Sigma0_inv %*% mu0 )
  
  # 3. Posterior mean and covariance
  Q_inv   <- solve(Q_beta)             # (p × p)
  m_beta  <- Q_inv %*% s_vec           # (p × 1)
  
  # 4. Sample  β  ∼  N(m_beta ,  Q_beta^{-1})
  Beta <- as.numeric( mvrnorm(1, mu = m_beta, Sigma = Q_inv) )
  
  
  # 8. **Sample Alpha_i | Beta, h, ω for each i** from its Gaussian full conditional.
  #    For each individual i:
  for (i in 1:n) {
    # Extract vectors for this individual's observations
    omega_i <- omega_mat[i, ]            # length t_obs
    h_i     <- h_mat[i, ]               # length t_obs
    X_i     <- X[i,, ]                  # (t_obs x p) matrix
    # Compute precision: Q_alpha_i = V_alpha^{-1} + ∑_j ω_ij Z_ij Z_ij^T.
    # (If Z is constant per individual, Z_ij = Z[i,] for all j.)
    Z_i <- matrix(Z[i, ], nrow = q, ncol = 1)    # (q x 1) vector for this individual
    # Sum over j: since Z_i is constant, ∑_j ω_ij * Z_i Z_i^T = (∑_j ω_ij) * (Z_i %*% t(Z_i))
    sum_omega_i <- sum(omega_i)
    Q_alpha_i <- solve(V_alpha) + sum_omega_i * (Z_i %*% t(Z_i))
    # Compute mean: m_alpha_i = Q_alpha_i^{-1} [ ∑_j ω_ij * Z_ij * (h_ij - X_ij^T Beta) ]
    # (Note: h_ij - X_ij^T Beta is the part of latent utility not explained by fixed effects.)
    resid_i <- h_i - (X_i %*% Beta)     # vector length t_obs
    # Since Z_i is constant, ∑_j ω_ij * (h_ij - X_ij β) * Z_i = (Z_i) * ∑_j ω_ij * resid_i
    S_i <- Z_i * sum(omega_i * resid_i)   # (q x 1) vector
    m_alpha_i <- solve(Q_alpha_i, S_i)
    # Draw Alpha_i ~ MVN(m_alpha_i, Q_alpha_i^{-1})
    cov_alpha_i <- solve(Q_alpha_i)
    Alpha[i, ] <- as.numeric(mvrnorm(1, mu = m_alpha_i, Sigma = cov_alpha_i))
  }
  
  # 9. **Sample V_alpha | Alpha** from Inverse-Wishart(ν0 + n, Λ0 + ∑_{i=1}^n α_i α_i^T).
  S_post <- Lambda0 + t(Alpha) %*% Alpha   # q x q scale matrix (sum of outer products of Alpha_i plus Lambda0)
  V_alpha <- riwish(nu0 + n, S_post)
  
  # -- Save current draws
  Beta_save[iter, ]    <- Beta
  Alpha_save[iter, , ] <- Alpha
  V_alpha_save[[iter]] <- V_alpha
  gamma_save[iter]     <- as.numeric(gamma)
}

# VIII. Discard burn-in samples
Beta_samples    <- Beta_save[(burnin+1):iterations, , drop=FALSE]
Alpha_samples   <- Alpha_save[(burnin+1):iterations, , , drop=FALSE]
V_alpha_samples <- V_alpha_save[(burnin+1):iterations]
gamma_samples   <- gamma_save[(burnin+1):iterations]

# The variables Beta_samples, Alpha_samples, V_alpha_samples, and gamma_samples 
# now contain the posterior samples after burn-in.
# For example, we can examine the posterior means:
colMeans(Beta_samples)        # posterior mean of beta
apply(Alpha_samples, 3, mean) # posterior mean of alpha (averaged over individuals, for each component)
Reduce("+", V_alpha_samples) / length(V_alpha_samples)  # posterior mean of V_alpha

# (Optional) Use coda to summarize or plot the MCMC:
beta_mcmc <- mcmc(Beta_samples)
summary(beta_mcmc)
# Example trace plot for beta:
plot(beta_mcmc) 


## ------------------------------------------------------------------
## 1.  Extraer columnas de Beta_samples  →  vectores individuales
## ------------------------------------------------------------------
beta_1_samples <- Beta_samples[ , 1]                  # (iter × 1)
beta_2_samples <- Beta_samples[ , 2]
beta_3_samples <- Beta_samples[ , 3]

## Etiquetas “aug” para ser coherentes
beta_1_samplespupg <- beta_1_samples
beta_2_samplespupg <- beta_2_samples
beta_3_samplespupg <- beta_3_samples

## ------------------------------------------------------------------
## 2.  Convertir lista de V_alpha_samples en array 3-D y sacar elementos
##      (q = 2 ⇒ cada matriz es 2×2, así que simplificamos a 3-D)
## ------------------------------------------------------------------
V_arr <- simplify2array(V_alpha_samples)   # dims: 2 × 2 × draws

V1_pupg  <- V_arr[1, 1, ]   # trayectoria de V[1,1]
V2_pupg  <- V_arr[2, 2, ]   # trayectoria de V[2,2]
V21_pupg <- V_arr[2, 1, ]   # trayectoria de V[2,1]  (equivalente a V[1,2])

## ------------------------------------------------------------------
## 3.  Sacar ejemplos de alpha:
##     a1_ag = componente 1 de α para el INDIVIDUO 1  a lo largo de draws
##     a2_ag = componente 2 de α para el INDIVIDUO 1  a lo largo de draws
##     (Ajusta el índice “1” si te interesa otro individuo)
## ------------------------------------------------------------------
a1_pupg <- Alpha_samples[ , 1, 1]   # (iter × 1)
a2_pupg <- Alpha_samples[ , 1, 2]

## Mantén el array completo por si lo necesitas
alpha_samples_pupg <- Alpha_samples

## ------------------------------------------------------------------
## 4.  Si quieres guardar también γ (opcional: lo llamo timesamples_aug
##     para respetar tu variable original; renómbralo si prefieres)
## ------------------------------------------------------------------
timesamples_pupg <- gamma_samples   # vector (iter × 1)

## ------------------------------------------------------------------
## 5.  Guardar TODO en un único .RData, idéntico a tu otro script
## ------------------------------------------------------------------
if (!dir.exists("Output")) dir.create("Output")

save(a1_pupg, a2_pupg,
     V21_pupg, V2_pupg, V1_pupg,
     beta_1_samplespupg, beta_2_samplespupg, beta_3_samplespupg,
     alpha_samples_pupg, timesamples_pupg,
     file = "Output/V_samplespupg.RData")

