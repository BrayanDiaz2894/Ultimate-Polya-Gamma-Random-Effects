###############################################################################
##  draw_gamma()  –  single draw of γ integrating out β and α
###############################################################################
draw_gamma <- function(omega_mat, h_mat, X, Z,
                       Sigma0_inv, mu0, V_alpha_inv, G0)
{
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required")
  
  ## ---------- basic dimensions ---------------------------------------------
  n     <- dim(X)[1]
  t_obs <- dim(X)[2]
  p     <- dim(X)[3]
  q     <- ncol(Z)
  Ntot  <- n * t_obs
  
  omega_vec <- as.vector(omega_mat)         # length Ntot
  h_vec     <- as.vector(h_mat)
  
  ## ---------- 1. build Σ_{α_i}^{-1},  d_i = Σ ω Z ,  c_i = Σ ω Z h ----------
  Sigma_blocks <- vector("list", n)         # each (q×q)
  d_mat <- matrix(0, n, q)                  # rows = d_i^T  (∑ ω Z)
  c_mat <- matrix(0, n, q)                  # rows = c_i^T  (∑ ω Z h)
  
  for (i in seq_len(n)) {
    Z_i <- matrix(Z[i, ], q, 1)             # q × 1   (constant over j)
    w_i <- omega_mat[i, ]                   # length t_obs
    h_i <- h_mat[i, ]
    
    sum_w_i <- sum(w_i)
    Sigma_ai <- sum_w_i * (Z_i %*% t(Z_i)) + V_alpha_inv
    Sigma_blocks[[i]] <- solve(Sigma_ai)    # (Σ_{α_i})⁻¹
    
    d_mat[i, ] <- sum_w_i       * t(Z_i)    # d_i = Σ ω Z
    c_mat[i, ] <- sum(w_i*h_i)  * t(Z_i)    # c_i = Σ ω Z h
  }
  
  ## stacked vectors and block-diag Pα⁻¹
  P_inv <- as.matrix(Matrix::bdiag(Sigma_blocks))          # (nq × nq)
  d_vec <- as.vector(t(d_mat))                             # (nq) × 1
  c_vec <- as.vector(t(c_mat))                             # (nq) × 1
  
  ## ---------- 2. flatten X once --------------------------------------------
  X_all <- matrix(NA_real_, Ntot, p)
  row_pt <- 0L
  for (i in seq_len(n)) {
    rng            <- row_pt + seq_len(t_obs)
    row_pt         <- row_pt + t_obs
    X_all[rng,  ]  <- X[i, , ]
  }
  
  ## ---------- 3. matrices to integrate α -----------------------------------
  eta_mat <- matrix(0, p, q)
  P_small <- matrix(0, q, q)
  
  for (i in seq_len(n)) {
    Z_i <- matrix(Z[i, ], q, 1)
    for (j in seq_len(t_obs)) {
      idx   <- (i-1L)*t_obs + j
      X_ij  <- matrix(X_all[idx, ], ncol = 1)
      w_ij  <- omega_vec[idx]
      
      eta_mat <- eta_mat + w_ij * (X_ij %*% t(Z_i))
      P_small <- P_small + w_ij * (Z_i  %*% t(Z_i))
    }
  }
  P_small <- P_small + V_alpha_inv
  P_small_inv <- solve(P_small)
  eta_t <- t(eta_mat)
  
  ## ---------- 4. integrate β  ----------------------------------------------
  Sigma_beta_star <- crossprod(X_all * sqrt(omega_vec)) +
    Sigma0_inv -
    (eta_mat %*% P_small_inv %*% eta_t)
  
  # helpers based on previous sums (length q)
  omega_d      <- colSums(d_mat)               # Σ ω Z          (identical to rowSums)
  omega_c      <- colSums(c_mat)               # Σ ω Z h
  
  # helpers in X (length p)
  omega_h_X <- crossprod(X_all, omega_vec * h_vec)  # p × 1
  omega_X   <- crossprod(X_all, omega_vec)          # p × 1
  
  mu_beta1  <- as.vector( omega_h_X -
                            eta_mat %*% P_small_inv %*% omega_c +
                            Sigma0_inv %*% mu0 )
  mu_beta2  <- as.vector( omega_X -
                            eta_mat %*% P_small_inv %*% omega_d )
  
  ## ---------- 5. pre-compute quadratic pieces ------------------------------
  dPd  <- as.numeric( t(d_vec) %*% P_inv %*% d_vec )   # dᵀ P⁻¹ d
  dPc  <- as.numeric( t(d_vec) %*% P_inv %*% c_vec )   # dᵀ P⁻¹ c
  cPc  <- as.numeric( t(c_vec) %*% P_inv %*% c_vec )   # cᵀ P⁻¹ c
  
  Sigma_beta_inv <- solve(Sigma_beta_star)
  mu_beta1_vec   <- matrix(mu_beta1, ncol = 1)
  mu_beta2_vec   <- matrix(mu_beta2, ncol = 1)
  
  bPb <- as.numeric( t(mu_beta2_vec) %*% Sigma_beta_inv %*% mu_beta2_vec )
  bPa <- as.numeric( t(mu_beta2_vec) %*% Sigma_beta_inv %*% mu_beta1_vec )
  
  ## ---------- 6. coefficients of the quadratic in γ ------------------------
  S0 <- sum(omega_vec)                           # Σ ω
  S1 <- sum(omega_vec * h_vec)                   # Σ ω h
  
  # Precision (coefficient of γ²) and mean
  prec_gamma <- S0 + dPd + 1/G0                 # a
  mean_num   <- S1 + dPc - bPa                  # b
  #  (signs match the algebra: see explanation en texto)
  
  mu_gamma     <- mean_num / prec_gamma
  sigma2_gamma <- 1 / prec_gamma
  
  ## ---------- 7. return a draw ---------------------------------------------
  rnorm(1, mu_gamma, sqrt(sigma2_gamma))
}
