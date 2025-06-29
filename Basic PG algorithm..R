#Include Random Effects
#Includes Explicit Sampling for the V_\alpha parameter.
#More efficient. 
#LAST VERSION. This version introduces more vectorization.

rm(list = ls())  # Removes all variables from memory
gc()             # Frees up unused memory



#0 Background. 

  #I. Load necessary libraries
  library(BayesLogit)
  library(MASS)
  library(ggplot2)
  library(MCMCpack)  # for riwish()
library(coda)
  library(ggplot2)
  #II. Set seed for reproducibility
  set.seed(123)
  
  #III. set parameters
  n <- 100  # Number of individuals
  t <- 30    # Number of observations per individual
  p <- 3    # Number of predictors
  q <- 2    # Number of individual-level covariates (random effects)
  reg <- 0 ##Regularization parameter in case that the inversechol of B_i is not invertible.
  
  #*IV. simulate the data.
  
    #i. basic data
    X <- array(rnorm(n * t * p), dim = c(n, t, p))  # Predictor matrix
    X[, , 1] <- 1
    Z <- matrix(rnorm(n * q), nrow = n, ncol = q)      # Individual-level covariates (time-invariant)
    
    # ii. additional elements
    V_alpha_true <- matrix(c(1, 0.5, 0.5, 3), nrow = q, ncol = q)  # True variance-covariance matrix #Needs to be symmetric.
    alpha_true <- mvrnorm(n, mu = rep(0, q), Sigma = V_alpha_true)  # Generate random effects
    beta_true <- c(-1, 2, -1.5)  # True beta coefficients
  
    #iii. Generate response variable Y
    eta <- array(0, dim = c(n, t))
    for (i in 1:n) {
      eta[i, ] <-  t(X[i, , ] %*% beta_true) + as.vector(Z[i, ] %*% alpha_true[i, ])
    }
    prob <- 1 / (1 + exp(-eta))
    Y <- matrix(rbinom(n * t, size = 1, prob = prob), nrow = n, ncol = t)
    
    #Lets define a funcion which invert matrixes based on cholesky decomp
    inversechol <- function(matrix,reg=0) { #matrix and regularization if any
      matrix_reg <- matrix + reg * diag(nrow(matrix))
      L <- chol(matrix_reg)  # L such that t(L) %*% L = Bi_reg
      inverted <- chol2inv(L)
      return(inverted)
    }

  #V. Gibbs Sampler settings
  num_iter <- 10000  # Number of MCMC iterations
  burn_in <- num_iter*0.1     # Burn-in period

  #VI. Store Parameters
  beta <- matrix(1, nrow = num_iter, ncol = p)  # Store sampled betas
  alpha <- array(1, dim = c(num_iter, n, q))    # Store sampled alphas
  omega_store <- array(1, dim = c(num_iter, n, t))  # Store omega samples
  V_alpha_store <- array(1, dim = c(num_iter, q, q))  # Store sampled V_alpha matrices
  
  #Start Parameters
  beta_current <- rep(0, p)  # Initial values for beta in zero such that it explores both positive and negative domain. 
  alpha_current <- mvrnorm(n, mu = rep(0, q), Sigma = diag(1, q)) #start close to zero, which is kind of informative, and diag(1) vcov
  V_alpha_current <- diag(1, q) #Diagonal matrix.
  V_alpha_inv <- inversechol(V_alpha_current)
   

  # Prior specification for beta #Normal prior for beta.
  mu_0 <- rep(0, p)
  Sigma_0_inv <- diag(1, p)

  # Prior for V_alpha: Inverse-Wishart parameters
  v0 <- 5 ##it needs to be larger than q-1
  Lambda0 <- diag(2, q) #We know that E[V_\alpha] = \Lambda_0 / (v_0 - q - 1) therefore with q=1  and v_0=5
  #We have that E[V_\alpha] = Diag(1,1)
  
#Stuff that can be computed before hand 
  
  # Precompute Z_i Z_i^T for each subject i
  ZZ_precomp <- array(NA, dim = c(n, q, q))
  for (i in 1:n) {
    ZZ_precomp[i,,] <- tcrossprod(Z[i, ])  # outer product: Z_i Z_i^T
  }
  
  X_flat <- matrix(X, nrow = n * t, ncol = p) #Used to vectorize beta. 
  X_flatt_cross <- t(X_flat)%*%X_flat
  row_id    <- rep(1:n, times = t)     
  #X is an array nxtxp. The code before basically what is doing is reducing one dimension and having on the rows, 
  #one row for each individual time. 
  
#Lets define a funcion which invert matrixes based on cholesky decomp
  
  
#2.  Gibbs Sampling
 start_time <- Sys.time()
for (iter in 1:num_iter) {
  print(iter*100/num_iter)
  
  omega_prev <- if (iter == 1) matrix(1, nrow = n, ncol = t) else omega_store[iter - 1, , ]
  
  # Step 1: Sample alpha_i (random effects)
  for (i in 1:n) {
    
    # 1. Posterior precision
    B_i <- V_alpha_inv + sum(omega_prev[i, ]) * ZZ_precomp[i,,]
    
    # 2. Posterior mean component
    temp <- Y[i, ] - 0.5 - omega_prev[i, ] * as.vector(X[i, , ] %*% beta_current)
    A_i <- sum(temp) * Z[i, ]
    
    # 3. Posterior covariance and mean
    cov_alpha <- inversechol(B_i, 1e-6)  # Σ = B⁻¹
    mu_i <- cov_alpha %*% A_i           # μ = Σ * A_i
    
    # 4. Cholesky-based sampling
    R <- chol(cov_alpha)                # chol() returns R such that Σ = R'R
    z <- rnorm(q)                       # standard normal vector
    alpha_current[i, ] <- as.numeric(mu_i + t(R) %*% z)  # αᵢ = μ + Rᵗ z
  }
  
  alpha[iter, , ] <- alpha_current
  

  # Step 2: Sample beta (fixed effects)

  # Vectorized computation for Sigma_beta_inv:
  omega_vec <- as.vector(omega_prev)
  # This first trick comes from this identity: 
  # \left(\sqrt{\boldsymbol{\omega}} \cdot \boldsymbol{X}\right)^\top\left(\sqrt{\boldsymbol{\omega}} \cdot \boldsymbol{X}\right)\;=\;\sum_{i,t} \omega_{it} \, \boldsymbol{X}_{it} \boldsymbol{X}_{it}^\top
  Sigma_beta_inv_data <- t(X_flat * sqrt(omega_vec)) %*% (X_flat * sqrt(omega_vec))
  Sigma_beta_inv <- Sigma_0_inv + Sigma_beta_inv_data
  
  # For mu_beta, first compute the contribution from the data.
  #We aim to \mu_\beta = \Sigma_\beta \left( \sum_{i=1}^n \sum_{j=1}^t \left( y_{ij} - \frac{1}{2} - \omega_{ij} \cdot (Z_i^\top \alpha_i) \right) X_{ij} + \Sigma_0^{-1} \mu_0 \right)

  # For each individual, compute the scalar sum_{k} (Z[i,] * alpha_current[i,]) #Note that here only the dimmension j is being summed,
  # So in the double summation we can take it out. 
  dot_Z_alpha <- rowSums(Z * alpha_current)  # length n
  
  # Replicate dot_Z_alpha across t columns to match Y and omega_prev.
  # Each row i is repeated t times.
  dot_Z_alpha_mat <- matrix(rep(dot_Z_alpha, each = t), nrow = n, ncol = t, byrow = TRUE)
  #This cretes a matrix where dot_Z_apha is the same for all individuals across t times. We will use it in a moment. 
  
  # Compute the t-by-n matrix of scalars for each (i,j):
  #Here y, omega_prev and dot_Z_alpha_mat are nxt, but dot_Z_alpha_mat has the same values across all columns (time invariant)
  #Then the next expression is 0.5*(2*Y_ij-1)*omega_i,j*Z_ij this is, we have a matrix where the rows is the individual and columns time, for the expression mentioned before.
  temp_mat <- 0.5 * (2 * Y - 1) - omega_prev * dot_Z_alpha_mat
  
  # Flatten temp_mat to a vector of length n*t:
  #We vectorize it, so we dont have 5 collumns any more.
  temp_vec <- as.vector(temp_mat)
  
  # Data contribution to mu_beta: sum_{i,j} (temp_scalar[i,j] * X[i,j,])
  mu_beta_data <- t(X_flat) %*% temp_vec
  
  # Finally, combine with the prior contribution:
  mu_beta <- inversechol(Sigma_beta_inv) %*% (mu_beta_data + Sigma_0_inv %*% mu_0)
  Sigma_beta <- inversechol(Sigma_beta_inv, 1e-6)  # posterior covariance
  R <- chol(Sigma_beta)                            # upper-triangular: Sigma_beta = R'R
  z <- rnorm(p)                                    # standard normal vector
  beta_current <- as.numeric(mu_beta + t(R) %*% z) # correct draw
  beta[iter, ] <- beta_current
  
  # Step 3: Sample omega_ij
  Zalpha_vec <- dot_Z_alpha[row_id]             # ✓ FIX  use 'times = t', not 'each = t'
  XB_vec <- X_flat %*% beta_current
  psi_vec      <- XB_vec + Zalpha_vec
  omega_vec    <- rpg(n*t, h = 1, z = psi_vec)     # BayesLogit::rpg is vectorised
  omega_matrix <- matrix(omega_vec, nrow = n, ncol = t, byrow = FALSE)
  omega_store[iter, , ] <- omega_matrix
  omega_prev            <- omega_matrix
  
  # Step 4: Sample V_alpha
  V_alpha_current <- riwish(v0 + n, Lambda0 + t(alpha_current) %*% alpha_current)
  V_alpha_store[iter, , ] <- V_alpha_current
  V_alpha_inv <- inversechol(V_alpha_current)
}
  end_time <- Sys.time()
  print(end_time - start_time)
  timesamples_naug <- end_time - start_time 
  
  
#3. Plot alpha
   # Extract alpha samples for the first individual
   alpha_first <- alpha[(burn_in + 1):num_iter, 3, 1]
   alpha_first_true <- alpha_true[3,1]
   
   # Plot trace plot comparing MCMC draws to population parameter
   alpha_trace <- data.frame(iter = (burn_in + 1):num_iter, alpha_first = alpha_first)
    
   ggplot(alpha_trace, aes(x = iter, y = alpha_first)) +
     geom_line(color = "blue") +
     geom_hline(yintercept = alpha_first_true, linetype = "dashed", color = "red") +
     labs(title = "Trace Plot of First Alpha Element for First Individual",
         x = "Iteration",
          y = "Alpha First Element") +
     theme_minimal()
    
#4. Sample V_\alpha

   # Extract samples for each element of V_alpha
   V_11_samples <- V_alpha_store[(burn_in + 1):num_iter, 1, 1]
   V_12_samples <- V_alpha_store[(burn_in + 1):num_iter, 1, 2]
   V_21_samples <- V_alpha_store[(burn_in + 1):num_iter, 2, 1]
   V_22_samples <- V_alpha_store[(burn_in + 1):num_iter, 2, 2]

   # Create iteration sequence
   iterations <- (burn_in + 1):num_iter

   # Convert to a data frame for ggplot
   V_alpha_df <- data.frame(
     iteration = rep(iterations, 4),
     value = c(V_11_samples, V_12_samples, V_21_samples, V_22_samples),
     element = rep(c("V[1,1]", "V[1,2]", "V[2,1]", "V[2,2]"), each = length(iterations))
   )

   # Create a mapping for true values
   true_values <- data.frame(
     element = c("V[1,1]", "V[1,2]", "V[2,1]", "V[2,2]"),
     true_value = c(V_alpha_true[1,1], V_alpha_true[1,2], V_alpha_true[2,1], V_alpha_true[2,2])
   )

   # Generate trace plots for each element of V_alpha
   ggplot(V_alpha_df, aes(x = iteration, y = value, color = element)) +
     geom_line(alpha = 0.6) +  # MCMC samples
     facet_wrap(~element, scales = "free_y") +  # Separate plots for each element
     geom_hline(data = true_values, aes(yintercept = true_value, color = element),
                linetype = "dashed", size = 1) +  # True values in matching colors
     scale_color_manual(values = c("blue", "red", "green", "purple")) +  # Distinct colors for elements
     labs(title = "Trace Plots for Each Element of V_alpha",
          x = "Iteration",
          y = "Value") +
     theme_minimal() +
     theme(legend.position = "none")  # Remove redundant legend

#5. Plot Beta


   # Extract beta samples for each parameter
   beta_1_samples <- beta[(burn_in + 1):num_iter, 1]
   beta_2_samples <- beta[(burn_in + 1):num_iter, 2]
   beta_3_samples <- beta[(burn_in + 1):num_iter, 3]

   # Create iteration sequence
   iterations <- (burn_in + 1):num_iter

   # Convert to a data frame for ggplot
   beta_df <- data.frame(
     iteration = rep(iterations, 3),
     value = c(beta_1_samples, beta_2_samples, beta_3_samples),
     element = rep(c("Beta[1]", "Beta[2]", "Beta[3]"), each = length(iterations))
   )

   # Create a mapping for true values
   true_values <- data.frame(
     element = c("Beta[1]", "Beta[2]", "Beta[3]"),
     true_value = beta_true
   )

   # Generate trace plots for each element of beta
   ggplot(beta_df, aes(x = iteration, y = value, color = element)) +
     geom_line(alpha = 0.6) +  # MCMC samples
     facet_wrap(~element, scales = "free_y") +  # Separate plots for each element
     geom_hline(data = true_values, aes(yintercept = true_value, color = element),
                linetype = "dashed", size = 1) +  # True values in matching colors
     scale_color_manual(values = c("blue", "red", "green")) +  # Distinct colors for elements
     labs(title = "Trace Plots for Each Element of Beta",
          x = "Iteration",
          y = "Value") +
     theme_minimal() +
     theme(legend.position = "none")  # Remove redundant legend
   

   a1_nag <- alpha[(burn_in + 1):num_iter, 1, ]
   a2_nag <- alpha[(burn_in + 1):num_iter, 2, ]
   beta_1_samplesnaug <- beta_1_samples
   beta_2_samplesnaug <- beta_2_samples
   beta_3_samplesnaug <- beta_3_samples
   V1_naug <- V_11_samples
   V2_naug <- V_22_samples
   V21_naug <- V_21_samples
   alpha_samples_naug <- alpha
   save(a1_nag, a2_nag, V21_naug, V2_naug, V1_naug, beta_1_samplesnaug, beta_2_samplesnaug, beta_3_samplesnaug,timesamples_naug, alpha_samples_naug, file = "Output/V_samples.RData")
 