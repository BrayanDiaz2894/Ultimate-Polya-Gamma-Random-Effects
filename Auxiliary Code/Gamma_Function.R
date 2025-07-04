###############################################################################
##  draw_gamma()  –  single draw of γ integrating out β and α
##  Versión con  η  apilada por individuo:
##      η_stack  =  rbind_i  Σ_j ω_ij Z_ij X_ijᵀ   ∈ ℝ^{(nq)×p}
###############################################################################
draw_gamma <- function(omega_mat, h_mat, X, Z,
                       Sigma0_inv, mu0, V_alpha_inv, G0)
{
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required")
  
  ## ---------- 0. dimensiones básicas ---------------------------------------
  n     <- dim(X)[1]
  t_obs <- dim(X)[2]
  p     <- dim(X)[3]
  q     <- ncol(Z)
  Ntot  <- n * t_obs
  
  omega_vec <- as.vector(omega_mat)          # longitud Ntot
  h_vec     <- as.vector(h_mat)
  
  ## ---------- 1. bloques Σ_{α_i}^{-1},  d_i,  c_i ---------------------------
  Sigma_blocks <- vector("list", n)          # cada uno (q×q)
  d_mat <- matrix(0, n, q)                   # filas = d_iᵀ = Σ_j ω_ij Z_ij
  c_mat <- matrix(0, n, q)                   # filas = c_iᵀ = Σ_j ω_ij Z_ij h_ij
  
  for (i in seq_len(n)) {
    Z_i <- matrix(Z[i, ], q, 1)              # q × 1  (constante en j)
    w_i <- omega_mat[i, ]                    # longitud t_obs
    h_i <- h_mat[i, ]
    
    sum_w_i <- sum(w_i)
    Sigma_ai      <- sum_w_i * (Z_i %*% t(Z_i)) + V_alpha_inv
    Sigma_blocks[[i]] <- solve(Sigma_ai)     # (Σ_{α_i})⁻¹
    
    d_mat[i, ] <- sum_w_i      * t(Z_i)      # d_i
    c_mat[i, ] <- sum(w_i*h_i) * t(Z_i)      # c_i
  }
  
  ## vectores apilados y bloque diagonal P⁻¹
  P_inv <- as.matrix(Matrix::bdiag(Sigma_blocks))        # (nq × nq)
  d_vec <- as.vector(t(d_mat))                           # (nq) × 1
  c_vec <- as.vector(t(c_mat))                           # (nq) × 1
  
  ## ---------- 2. aplanar X una sola vez ------------------------------------
  X_all <- matrix(NA_real_, Ntot, p)
  row_pt <- 0L
  for (i in seq_len(n)) {
    rng            <- row_pt + seq_len(t_obs)
    row_pt         <- row_pt + t_obs
    X_all[rng,  ]  <- X[i, , ]
  }
  
  ## ---------- 3. construir η apilada por individuo -------------------------
  eta_stack <- matrix(0, n * q, p)                       # (nq × p)
  
  for (i in seq_len(n)) {
    Z_i    <- matrix(Z[i, ], q, 1)                       # q × 1
    eta_i  <- matrix(0, q, p)                            # bloque q × p
    
    for (j in seq_len(t_obs)) {
      idx   <- (i - 1L) * t_obs + j
      X_ij  <- matrix(X_all[idx, ], ncol = 1)            # p × 1
      w_ij  <- omega_vec[idx]
      
      eta_i <- eta_i + w_ij * (Z_i %*% t(X_ij))          # q × p
    }
    
    rows_i <- ((i - 1) * q + 1):(i * q)                  # filas del bloque
    eta_stack[rows_i, ] <- eta_i                         # colocar bloque
  }
  
  ## ---------- 4. integrar β -------------------------------------------------
  ##   Σ_β⁻¹ = Xᵀ Ω X  + Σ₀⁻¹  −  ηᵀ P⁻¹ η
  Sigma_beta_star <-
    crossprod(X_all * sqrt(omega_vec)) + #This is because \text{crossprod}(X_{\text{all}} * \sqrt{\omega}) = X^\top \, \mathrm{diag}(\omega) \, X = \sum_{i=1}^{N} \omega_i \, X_i X_i^\top
    Sigma0_inv -
    (t(eta_stack) %*% P_inv %*% eta_stack)
  
  ## auxiliares en X (p × 1) y en Z (nq × 1)
  omega_h_X <- crossprod(X_all, omega_vec * h_vec)       # p × 1
  omega_X   <- crossprod(X_all, omega_vec)               # p × 1
  
  ## medias condicionales de β
  mu_beta1 <- as.vector(
    omega_h_X -
      t(eta_stack) %*% P_inv %*% c_vec +
      Sigma0_inv %*% mu0)
  
  mu_beta2 <- as.vector(
    omega_X -
      t(eta_stack) %*% P_inv %*% d_vec)
  
  ## ---------- 5. piezas cuadráticas previas --------------------------------
  dPd <- as.numeric( t(d_vec) %*% P_inv %*% d_vec )      # dᵀ P⁻¹ d
  dPc <- as.numeric( t(d_vec) %*% P_inv %*% c_vec )      # dᵀ P⁻¹ c
  
  Sigma_beta_inv <- solve(Sigma_beta_star)
  mu_beta1_vec   <- matrix(mu_beta1, ncol = 1)
  mu_beta2_vec   <- matrix(mu_beta2, ncol = 1)
  
  bPb <- as.numeric( t(mu_beta2_vec) %*% Sigma_beta_inv %*% mu_beta2_vec )
  bPa <- as.numeric( t(mu_beta2_vec) %*% Sigma_beta_inv %*% mu_beta1_vec )
  
  ## ---------- 6. coeficientes del cuadrático en γ --------------------------
  S0 <- sum(omega_vec)                              # Σ ω
  S1 <- sum(omega_vec * h_vec)                      # Σ ω h
  
  prec_gamma <- S0 + dPd + 1 / G0                  # coef. γ²
  mean_num   <- S1 + dPc - bPa                     # coef. γ
  
  mu_gamma     <- mean_num / prec_gamma
  sigma2_gamma <- 1 / prec_gamma
  
  ## ---------- 7. devolver un draw ------------------------------------------
  rnorm(1, mu_gamma, sqrt(sigma2_gamma))
}
