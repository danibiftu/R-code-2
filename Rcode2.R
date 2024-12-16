###### R code for simulation of spatial latent trait model with Bayesian rank likelihood  #########################################################
# Load necessary libraries
library("LaplacesDemon")
library("geoR")
library("mvtnorm")
library("MCMCpack")
library("msm")
# Set seed for reproducibility
set.seed(12345)
# Define constants
n <- 400  # Sample size
groups <- 11
beta <- matrix(c(-1, 0.8, 0.6, 0.2), ncol = 1)
sigma2_random_effects <- 0.83
tau2_spatial_effects <- 0.67
lamda2 <- 0.5
psi2 <- 0.4
omega <- 5
eta <- 5
n_r <- 3
r_spatial <- c(0.2, 0.4, 0.1)
mean_spatial <- sum(r_spatial) / n_r
# Create design matrices
x_matt <- list()
for (i in 1:groups) {
  x_matt[[i]] <- matrix(NA, nrow = n, ncol = length(beta))
  x_matt[[i]][, 1] <- rnorm(n, 15, 0.1)
  x_matt[[i]][, 2] <- rpois(n, 2)
  x_matt[[i]][, 3] <- rnorm(n, 23, 0.2)
  x_matt[[i]][, 4] <- rpois(n, 2)
}
# Generate epsilon
epsilon <- rnorm(n, 0, 1)
# Initialize matrices and parameters
latent <- matrix(NA, n, groups)
random_effects <- rnorm(groups, 0, sigma2_random_effects)
spatial_effects <- rnorm(groups, mean_spatial, tau2_spatial_effects / n_r)
# Calculate latent variables
for (j in 1:groups) {
  latent[, j] <- x_matt[[j]] %*% beta + random_effects[j] + spatial_effects[j] + epsilon
}
# Convert latent variables to categorical data
latent_cat <- ceiling(latent)
for (i in 1:nrow(latent_cat)) {
  for (j in 1:ncol(latent_cat)) {
    if (latent_cat[i, j] >= -3 & latent_cat[i, j] < 0) {
      latent_cat[i, j] <- (-3 + 0) / 2
    } else if (latent_cat[i, j] >= 0 & latent_cat[i, j] < 3) {
      latent_cat[i, j] <- (0 + 3) / 2
    } else if (latent_cat[i, j] >= 3 & latent_cat[i, j] < 6) {
      latent_cat[i, j] <- (3 + 6) / 2
    } else if (latent_cat[i, j] >= 6 & latent_cat[i, j] < 9) {
      latent_cat[i, j] <- (6 + 9) / 2
    }
  }
}
# Plot latent categorical variables
par(mfrow = c(3, 1))
plot(latent_cat[, 1], type = "l", main = "Latent Category 1")
plot(latent_cat[, 5], type = "l", main = "Latent Category 5")
plot(latent_cat[, 11], type = "l", main = "Latent Category 11")
# Prepare data for further analysis
Y <- latent_cat
n_obs <- nrow(Y) - 1
Y <- Y[1:n_obs, ]
Y_cat <- ceiling(Y)
# Scale covariates
scaled_x_matt <- list()
for (i in 1:groups) {
  scaled_x_matt[[i]] <- scale(x_matt[[i]], center = TRUE, scale = TRUE)
}
# Combine design matrices
X <- array(cbind(scaled_x_matt[[1]][1:n_obs,], scaled_x_matt[[2]][1:n_obs,], scaled_x_matt[[3]][1:n_obs,], scaled_x_matt[[4]][1:n_obs,],
                 scaled_x_matt[[5]][1:n_obs,], scaled_x_matt[[6]][1:n_obs,], scaled_x_matt[[7]][1:n_obs,], scaled_x_matt[[8]][1:n_obs,],
                 scaled_x_matt[[9]][1:n_obs,], scaled_x_matt[[10]][1:n_obs,], scaled_x_matt[[11]][1:n_obs,]), c(n_obs, ncol(scaled_x_matt[[1]]), groups))
X_init <- array(cbind(scaled_x_matt[[1]], scaled_x_matt[[2]], scaled_x_matt[[3]], scaled_x_matt[[4]], scaled_x_matt[[5]],
                      scaled_x_matt[[6]], scaled_x_matt[[7]], scaled_x_matt[[8]], scaled_x_matt[[9]], scaled_x_matt[[10]], scaled_x_matt[[11]]), 
                c((n_obs + 1), ncol(scaled_x_matt[[1]]), groups))
# Prior hyperparameters
lamda2_prior <- 5
psi2_prior <- 10
omega_upsilon_prior <- 10
eta_gamma_prior <- 15
# Initial values for MCMC
BETA <- matrix(0, ncol(X[, , 1]), 1)
z <- matrix(0.5, nrow = n_obs, ncol = groups)
a <- b <- matrix(NA, nrow = n_obs, ncol = groups)
# MCMC settings
n_iter <- 50000
BETA_POST_3<- gamma_POST_3 <- upsilon_POST_3 <- tausq_gamma_POST_3 <- sigmasq_upsilon_POST_3<- NULL
upsilon_posterior <- rep(0.1, groups)
gamma_posterior <- rep(0.1, groups)
sigmasq_upsilon_post <- NULL
tausq_gamma_post <- NULL
z_all_3 <- NULL
# Gibbs Sampler
for (i in 1:n_iter) {
  for (j in 1:groups) {
    # Update degrees of freedom and scale parameters for random and spatial effects
    omega_upsilon_posterior <- omega_upsilon_prior + 1
    lamda2_posterior <- lamda2_prior + (upsilon_posterior[j]^2) / omega_upsilon_posterior 
    eta_gamma_posterior <- eta_gamma_prior + 1
    psi2_posterior <- psi2_prior + (n_r * (gamma_posterior[j] - mean_spatial)^2) / eta_gamma_posterior
    
    # Sample variance of random and spatial effects
    sigmasq_upsilon_post[j] <- rinvchisq(1, omega_upsilon_posterior, lamda2_posterior)
    tausq_gamma_post[j] <- rinvchisq(1, eta_gamma_posterior, psi2_posterior)
  }
  # Sample from random and spatial effects
  mean_upsilon_posterior <- array(NA, c(n_obs, 1, groups))
  mean_gamma_posterior <- array(NA, c(n_obs, 1, groups)) 
  for (t in 1:n_obs) {
    for (j in 1:groups) {
      mean_upsilon_posterior[t,,j] <- (z[t, j] - X[t, , j] %*% BETA)
      mean_gamma_posterior[t,,j] <- (z[t, j] - X[t, , j] %*% BETA)
    }
  }
  for (j in 1:groups) {
    var_upsilon_posterior <- 1 / (n_obs + (1 / sigmasq_upsilon_post[j]))
    var_gamma_posterior <- 1 / (n_obs + (n_obs / tausq_gamma_post[j])) 
    mean_upsilon_posterior_sum <- sum(mean_upsilon_posterior[,,j])
    mean_gamma_posterior_sum <- sum(mean_gamma_posterior[,,j]) 
    mean_upsilon_posterior_sum <- (mean_upsilon_posterior_sum - n_obs * gamma_posterior[j]) / (n_obs + 1 / sigmasq_upsilon_post[j])
    upsilon_posterior[j] <- rnorm(1, mean_upsilon_posterior_sum, sqrt(var_upsilon_posterior))
    mean_gamma_posterior_sum <- (r_spatial_sum + tausq_gamma_post[j] * mean_gamma_posterior_sum - n_obs * upsilon_posterior[j]) / (n_obs * tausq_gamma_post[j] + n_obs)
    gamma_posterior[j] <- rnorm(1, mean_gamma_posterior_sum, sqrt(var_gamma_posterior))
  }
  # Update latent variable z
  upsilon_rep <- array(rep(upsilon_posterior, each = n_obs), c(n_obs, 1, groups))
  gamma_rep <- array(rep(gamma_posterior, each = n_obs), c(n_obs, 1, groups))
  for (t in 1:n_obs) {
    for (j in 1:groups) {
      a[t, j] <- max(z[t, Y_cat[t,] < Y_cat[t, j]])
      b[t, j] <- min(z[t, Y_cat[t, j] < Y_cat[t,]])
      z[t, j] <- rtnorm(1, X[t, , j] %*% BETA + upsilon_rep[t, 1, j] + gamma_rep[t, 1, j], 1, a[t, j], b[t, j])
    }
  }
  # Update BETA
  muj1 <- array(NA, c(ncol(X[, , 1]), ncol(X[, , 1]), groups))
  muj2 <- array(NA, c(ncol(X[, , 1]), 1, groups))
  for (j in 1:groups) {
    muj1[,,j] <- t(X[,,j]) %*% X[,,j]
    muj2[,,j] <- t(X[,,j]) %*% (z[,j] - upsilon_rep[,1,j] - gamma_rep[,1,j])
  }
  muj1_sum <- apply(muj1, c(1, 2), sum)
  muj2_sum <- apply(muj2, c(1, 2), sum)
  muj_final <- solve(muj1_sum) %*% muj2_sum
  varj_final <- solve(muj1_sum)
  BETA <- mvrnorm(1, muj_final, varj_final)
  # Store results
  BETA_POST_3<- rbind(BETA_POST_3, BETA)
  upsilon_POST_3 <- rbind(upsilon_POST_3, upsilon_posterior)
  gamma_POST_3 <- rbind(gamma_POST_5, gamma_posterior)
  sigmasq_upsilon_POST_3 <- rbind(sigmasq_upsilon_POST_3, sigmasq_upsilon_post)
  tausq_gamma_POST_3 <- rbind(tausq_gamma_POST_3, tausq_gamma_post)
  z_all_3 <- rbind(z_all_3, z) 
  print(paste("Iteration", i))
}
# Load necessary libraries
library("coda")
library("ggplot2")
library("bayesplot")
# Define the chains for BETA parameters
chain_1 <- BETA_POST_1
chain_2 <- BETA_POST_2
chain_3 <- BETA_POST_3
# Combine specific chains for each BETA parameter
beta1 <- cbind(chain_1[1:50000, 1], chain_2[1:50000, 1], chain_3[1:50000, 1])
beta2 <- cbind(chain_1[1:50000, 2], chain_2[1:50000, 2], chain_3[1:50000, 2])
beta3 <- cbind(chain_1[1:50000, 3], chain_2[1:50000, 3], chain_3[1:50000, 3])
beta4 <- cbind(chain_1[1:50000, 4], chain_2[1:50000, 4], chain_3[1:50000, 4])
# Combine all BETA parameters into a single array
chains_BETA <- array(cbind(beta1, beta2, beta3, beta4), c(50000, 3, 4))
# Set new column and parameter names
dimnames(chains_BETA)[[2]] <- c("chain1", "chain2", "chain3")
dimnames(chains_BETA)[[3]] <- c("\u03B21", "\u03B22", "\u03B23", "\u03B24")
# Plot MCMC sensitivity traces for BETA parameters
color_scheme_set("viridis")
mcmc_trace(chains_BETA, pars = c("\u03B21", "\u03B22", "\u03B23", "\u03B24"), 
           facet_args = list(ncol = 1, strip.position = "left"))
# Define the chains for upsilon parameters
chain_11 <- upsilon_POST_1
chain_22<- upsilon_POST_2
chain_33<- upsilon_POST_3
# Combine specific chains for each upsilon parameter
upsilon1 <- cbind(chain_11[1:50000, 1], chain_22[1:50000, 1], chain_33[1:50000, 1])
upsilon2 <- cbind(chain_11[1:50000, 2], chain_22[1:50000, 2], chain_33[1:50000, 2])
upsilon3 <- cbind(chain_11[1:50000, 3], chain_22[1:50000, 3], chain_33[1:50000, 3])
upsilon4 <- cbind(chain_11[1:50000, 4], chain_22[1:50000, 4], chain_33[1:50000, 4])
upsilon5 <- cbind(chain_11[1:50000, 5], chain_22[1:50000, 5], chain_33[1:50000, 5])
upsilon6 <- cbind(chain_11[1:50000, 6], chain_22[1:50000, 6], chain_33[1:50000, 6])
upsilon7 <- cbind(chain_11[1:50000, 7], chain_22[1:50000, 7], chain_33[1:50000, 7])
upsilon8 <- cbind(chain_11[1:50000, 8], chain_22[1:50000, 8], chain_33[1:50000, 8])
upsilon9 <- cbind(chain_11[1:50000, 9], chain_22[1:50000, 9], chain_33[1:50000, 9])
upsilon10 <- cbind(chain_11[1:50000, 10], chain_22[1:50000, 10], chain_33[1:50000, 10])
upsilon11 <- cbind(chain_11[1:50000, 11], chain_22[1:50000, 11], chain_33[1:50000, 11])
# Combine all upsilon parameters into a single array
chains_upsilon <- array(cbind(upsilon1, upsilon2, upsilon3, upsilon4, upsilon5, upsilon6, 
                              upsilon7, upsilon8, upsilon9, upsilon10, upsilon11), c(50000, 3, 11))
# Set new column and parameter names
dimnames(chains_upsilon)[[2]] <- c("chain1", "chain2", "chain3")
dimnames(chains_upsilon)[[3]] <- c("\u03C51", "\u03C52", "\u03C53", "\u03C54", "\u03C55", 
                                   "\u03C56", "\u03C57", "\u03C58", "\u03C59", "\u03C510", "\u03C511")
# Plot MCMC sensitivity traces for upsilon parameters
color_scheme_set("viridis")
mcmc_trace(chains_upsilon, pars = c("\u03C51", "\u03C52", "\u03C53", "\u03C54", "\u03C55", 
                                    "\u03C56", "\u03C57", "\u03C58", "\u03C59", "\u03C510", "\u03C511"), 
           facet_args = list(ncol = 2, strip.position = "left"))
# Define the chains for gamma parameters
chain_111<- gamma_POST_1
chain_222 <- gamma_POST_2
chain_333<- gamma_POST_3
# Combine specific chains for each gamma parameter
gamma1 <- cbind(chain_111[1:50000, 1], chain_222[1:50000, 1], chain_333[1:50000, 1])
gamma2 <- cbind(chain_111[1:50000, 2], chain_222[1:50000, 2], chain_333[1:50000, 2])
gamma3 <- cbind(chain_111[1:50000, 3], chain_222[1:50000, 3], chain_333[1:50000, 3])
gamma4 <- cbind(chain_111[1:50000, 4], chain_222[1:50000, 4], chain_333[1:50000, 4])
gamma5 <- cbind(chain_111[1:50000, 5], chain_222[1:50000, 5], chain_333[1:50000, 5])
gamma6 <- cbind(chain_111[1:50000, 6], chain_222[1:50000, 6], chain_333[1:50000, 6])
gamma7 <- cbind(chain_111[1:50000, 7], chain_222[1:50000, 7], chain_333[1:50000, 7])
gamma8 <- cbind(chain_111[1:50000, 8], chain_222[1:50000, 8], chain_333[1:50000, 8])
gamma9 <- cbind(chain_111[1:50000, 9], chain_222[1:50000, 9], chain_333[1:50000, 9])
gamma10 <- cbind(chain_111[1:50000, 10], chain_222[1:50000, 10], chain_333[1:50000, 10])
gamma11 <- cbind(chain_111[1:50000, 11], chain_222[1:50000, 11], chain_333[1:50000, 11])
# Combine all gamma parameters into a single array
chains_gamma <- array(cbind(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, 
                            gamma7, gamma8, gamma9, gamma10, gamma11), c(50000, 3, 11))
# Set new column and parameter names
dimnames(chains_gamma)[[2]] <- c("chain1", "chain2", "chain3")
dimnames(chains_gamma)[[3]] <- c("\u03B31", "\u03B32", "\u03B33", "\u03B34", "\u03B35", 
                                 "\u03B36", "\u03B37", "\u03B38", "\u03B39", "\u03B310", "\u03B311")
# Plot MCMC sensitivity traces for gamma parameters
color_scheme_set("viridis")
mcmc_trace(chains_gamma, pars = c("\u03B31", "\u03B32", "\u03B33", "\u03B34", "\u03B35", 
                                  "\u03B36", "\u03B37", "\u03B38", "\u03B39", "\u03B310", "\u03B311"), 
           facet_args = list(ncol = 2, strip.position = "left"))
# Define true parameter values
beta <- matrix(BETA, 4, 1)
upsilon <- matrix(upsilon_posterior, m, 1)
gamma <- matrix(gamma_posterior, m, 1)
# Compute RMSE & MAE for beta chain 2
posterior_means_beta_1 <- matrix(colMeans(BETA_POST_1), 4, 1)
mse_beta_1 <- mean((posterior_means_beta_1 - beta)^2)
RMSE_beta_1 <- sqrt(mse_beta_1)
MAE_beta_1 <- mean(abs(posterior_means_beta_1 - beta))
# Compute RMSE & MAE for upsilon chain 1
posterior_means_upsilon_1 <- matrix(colMeans(upsilon_POST_1), m, 1)
mse_upsilon_1 <- mean((posterior_means_upsilon_1 - upsilon)^2)
RMSE_upsilon_1 <- sqrt(mse_upsilon_1)
MAE_upsilon_1 <- mean(abs(posterior_means_upsilon_1 - upsilon))
# Compute RMSE & MAE for gamma chain 1
posterior_means_gamma_ <- matrix(colMeans(gamma_POST_1), m, 1)
mse_gamma_1 <- mean((posterior_means_gamma_1 - gamma)^2)
RMSE_gamma_1 <- sqrt(mse_gamma_1)
MAE_gamma_1 <- mean(abs(posterior_means_gamma_1 - gamma))
# Compute RMSE & MAE for beta chain 2
posterior_means_beta_2 <- matrix(colMeans(BETA_POST_2), 4, 1)
mse_beta_2 <- mean((posterior_means_beta_2 - beta)^2)
RMSE_beta_2 <- sqrt(mse_beta_2)
MAE_beta_2<- mean(abs(posterior_means_beta_2 - beta)))
# Compute RMSE & MAE for upsilon chain 2
posterior_means_upsilon_2 <- matrix(colMeans(upsilon_POST_2), m, 1)
mse_upsilon_2<- mean((posterior_means_upsilon_2 - upsilon)^2)
RMSE_upsilon_2<- sqrt(mse_upsilon_2)
MAE_upsilon_2<- mean(abs(posterior_means_upsilon_2 - upsilon))
# Compute RMSE & MAE for gamma chain 2
posterior_means_gamma_2 <- matrix(colMeans(gamma_POST_2), m, 1)
mse_gamma_2 <- mean((posterior_means_gamma_2 - gamma)^2)
RMSE_gamma_2 <- sqrt(mse_gamma_2)
MAE_gamma_2<- mean(abs(posterior_means_gamma_2 - gamma))
# Compute RMSE & MAE for beta chain 3
posterior_means_beta_3 <- matrix(colMeans(BETA_POST_3), 4, 1)
mse_beta_3 <- mean((posterior_means_beta_3 - beta)^2)
RMSE_beta_3 <- sqrt(mse_beta_3)
MAE_beta_3 <- mean(abs(posterior_means_beta_3 - beta)))
# Compute RMSE & MAE for upsilon chain 3
posterior_means_upsilon_3 <- matrix(colMeans(upsilon_POST_3), m, 1)
mse_upsilon_3 <- mean((posterior_means_upsilon_3 - upsilon)^2)
RMSE_upsilon_3 <- sqrt(mse_upsilon_3)
MAE_upsilon_3 <- mean(abs(posterior_means_upsilon_3 - upsilon)))
# Compute RMSE & MAE for gamma chain 3
posterior_means_gamma_3 <- matrix(colMeans(gamma_POST_3), m, 1)
mse_gamma_3<- mean((posterior_means_gamma_3 - gamma)^2)
RMSE_gamma_3 <- sqrt(mse_gamma_3)
MAE_gamma_3 <- mean(abs(posterior_means_gamma_3- gamma)))
# PSRF/Gelman-Rubin convergence test
library(stableGR)
library(coda)
library(bayesplot)
psrf_1_beta_test <- stable.GR(BETA_POST_1)
psrf_2_beta_test <- stable.GR(BETA_POST_2)
psrf_3_beta_test <- stable.GR(BETA_POST_3)
psrf_1_upsilon_test <- stable.GR(upsilon_POST_1)
psrf_2_upsilon_test <- stable.GR(upsilon_POST_2)
psrf_3_upsilon_test <- stable.GR(upsilon_POST_3)
psrf_1_gamma_test <- stable.GR(gamma_POST_1)
psrf_2_gamma_test <- stable.GR(gamma_POST_2)
psrf_3_gamma_test <- stable.GR(gamma_POST_3)
# Checking the convergence with burn
Beta_gibbis_thinned <- BETA_POST_5[seq(1, nrow(BETA_POST_5[10001:50000, ]), 10), ]
gamma_gibbis_thinned <- gamma_POST_5[seq(1, nrow(gamma_POST_5[10001:50000, ]), 10), ]
upsilon_gibbis_thinned <- upsilon_POST_5[seq(1, nrow(upsilon_POST_5[10001:50000, ]), 10), ]
# Trace plot and density plot for Beta
par(mfrow = c(4, 1), mar = c(3, 3, 2, 1))
for (i in 1:4) {
  traceplot(as.mcmc(Beta_gibbis_thinned[, i]), main = substitute(paste(beta[i]), list(i = i)))
}
par(mfrow = c(2, 2), mar = c(7, 3, 3, 1))
for (i in 1:4) {
  densplot(as.mcmc(Beta_gibbis_thinned[, i]), main = substitute(paste(beta[i]), list(i = i)))
}
# Autocorrelation plot for Beta
color_scheme_set("mix-blue-red")
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
for (i in 1:4) {
  acf(as.mcmc(Beta_gibbis_thinned[, i]), ylab = "Autocorrelation", ci = FALSE, lag.max = 4000, main = substitute(paste(beta[i]), list(i = i)))
}
# Posterior means of Beta
cat("Posterior means of Beta:\n")
for (i in 1:4) {
  cat(paste("Beta[", i, "]: ", mean(Beta_gibbis_thinned[, i]), "\n", sep = ""))
}
# Posterior standard deviations of Beta
cat("\nPosterior standard deviations of Beta:\n")
for (i in 1:4) {
  cat(paste("Beta[", i, "]: ", sd(Beta_gibbis_thinned[, i]), "\n", sep = ""))
}
# For upsilon
par(mfrow = c(6, 1), mar = c(3, 3, 2, 2))
for (i in 1:6) {
  traceplot(as.mcmc(upsilon_gibbis_thinned[, i]), main = substitute(paste(upsilon[i]), list(i = i)))
}
par(mfrow = c(5, 1), mar = c(3, 3, 2, 2))
for (i in 7:11) {
  traceplot(as.mcmc(upsilon_gibbis_thinned[, i]), main = substitute(paste(upsilon[i]), list(i = i)))
}
par(mfrow = c(3, 2), mar = c(7, 3, 3, 1))
for (i in 1:6) {
  densplot(as.mcmc(upsilon_gibbis_thinned[, i]), main = substitute(paste(upsilon[i]), list(i = i)))
}
for (i in 7:11) {
  densplot(as.mcmc(upsilon_gibbis_thinned[, i]), main = substitute(paste(upsilon[i]), list(i = i)))
}
par(mfrow = c(3, 2), mar = c(7, 4, 3, 1))
for (i in 1:6) {
  acf(as.mcmc(upsilon_gibbis_thinned[, i]), ylab = "Autocorrelation", ci = F, lag.max = 4000, main = substitute(paste(upsilon[i]), list(i = i)))
}
for (i in 7:11) {
  acf(as.mcmc(upsilon_gibbis_thinned[, i]), ylab = "Autocorrelation", ci = F, lag.max = 4000, main = substitute(paste(upsilon[i]), list(i = i)))
}
# Posterior means of upsilon
cat("Posterior means of upsilon:\n")
for (i in 1:11) {
  cat(paste("upsilon[", i, "]: ", mean(upsilon_gibbis_thinned[, i]), "\n", sep = ""))
}
# Posterior standard deviations of upsilon
cat("\nPosterior standard deviations of upsilon:\n")
for (i in 1:11) {
  cat(paste("upsilon[", i, "]: ", sd(upsilon_gibbis_thinned[, i]), "\n", sep = ""))
}
# For gamma
par(mfrow = c(6, 1), mar = c(3, 3, 2, 2))
for (i in 1:6) {
  traceplot(as.mcmc(gamma_gibbis_thinned[, i]), main = substitute(paste(gamma[i]), list(i = i)))
}
par(mfrow = c(5, 1), mar = c(3, 3, 2, 2))
for (i in 7:11) {
  traceplot(as.mcmc(gamma_gibbis_thinned[, i]), main = substitute(paste(gamma[i]), list(i = i)))
}

par(mfrow = c(3, 2), mar = c(7, 3, 3, 1))
for (i in 1:6) {
  densplot(as.mcmc(gamma_gibbis_thinned[, i]), main = substitute(paste(gamma[i]), list(i = i)))
}
for (i in 7:11) {
  densplot(as.mcmc(gamma_gibbis_thinned[, i]), main = substitute(paste(gamma[i]), list(i = i)))
}

par(mfrow = c(3, 2), mar = c(7, 4, 3, 1))
for (i in 1:6) {
  acf(as.mcmc(gamma_gibbis_thinned[, i]), ylab = "Autocorrelation", ci = F, lag.max = 4000, main = substitute(paste(gamma[i]), list(i = i)))
}
for (i in 7:11) {
  acf(as.mcmc(gamma_gibbis_thinned[, i]), ylab = "Autocorrelation", ci = F, lag.max = 4000, main = substitute(paste(gamma[i]), list(i = i)))
}
# Posterior means of gamma
cat("Posterior means of gamma:\n")
for (i in 1:11) {
  cat(paste("gamma[", i, "]: ", mean(gamma_gibbis_thinned[, i]), "\n", sep = ""))
}
# Posterior standard deviations of gamma
cat("\nPosterior standard deviations of gamma:\n")
for (i in 1:11) {
  cat(paste("gamma[", i, "]: ", sd(gamma_gibbis_thinned[, i]), "\n", sep = ""))
}
#### Between group variation and spatial correlation 
# Load the necessary packages
library(reshape2)
library(ggplot2)
library(viridis)
library(bayesplot)
library(RColorBrewer)
# Set column names for the data
new_colnames <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
# variation plot
colnames(g) <- new_colnames
g <- data.frame(upsilon_POST_5)
data <- var(g[sapply(g, is.numeric)])
data1 <- melt(data)
colnames(data1) <- c("Var1", "Var2", "Variations")
x_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
y_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
upsilon_plot <- ggplot(data1, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() +
  scale_fill_viridis(discrete = FALSE) +
  labs(
    title = "Between Group Variation",
    x = "Regions",
    y = "Regions"
  ) +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = y_labels)
# Set column names for the data
new_colnames <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
colnames(g) <- new_colnames
#spatial correlation plot
gg <- data.frame(gamma_POST_5)
data <- cor(gg[sapply(gg, is.numeric)])
data1 <- melt(data)
colnames(data1) <- c("Var1", "Var2", "Variations")
gamma_plot <- ggplot(data1, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() +
  scale_fill_viridis(discrete = FALSE) +
  labs(
    title = "Spatial Correlation",
    x = "Regions",
    y = "Regions",
    fill = "Spatial correlation"
  ) +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = y_labels)
# Posterior predictions check
# The last observation
last_obs_5 <- array(NA, c(ncol(X), m, 1))
for (j in 1:m) {
  last_obs_5[, j, 1] <- X_init[(n_obs + 1), , j]
}

z_pred_outsample_5 <- array(NA, c(n_iter, 1, m))
for (i in 1:1) {
  for (j in 1:m) {
    z_pred_outsample_5[, i, j] <- rnorm(n_iter, BETA_POST_5 %*% last_obs_5[, j, i] + gamma_POST_5[, j] + upsilon_POST_5[, j], 1)
  }
}

# VEDI QUI
Y_init <- latent_cat
posterior_meansz1 <- matrix(NA, 1, m)
posterior_meansz1 <- t(colMeans(z_pred_outsample_5, dims = 1))
par(mfrow = c(1, 1))
plot(latent[100, ], type = "b", main = "Real data", ylab = "Dependent variable")
plot(posterior_meansz1[, 1], type = "b", main = "Predicted data", ylab = "Latent variable")
library("ggplot2")
d1 <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), y = latent[100, ]) # Real data
d2 <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), y = posterior_meansz1[, 1]) # Predicted data
par(mar = c(5, 4, 4, 4) + 0.3) # Additional space for second y-axis
plot(d1$x, d1$y, pch = 16, col = "blue", xlab = "", ylab = "", cex = 2, main = "Predictions with GLMM using Rank Likelihood") # Create first plot
par(new = TRUE) # Add new plot
plot(d2$x, d2$y, pch = 17, col = "purple", cex = 2, axes = FALSE, xlab = "Groups", ylab = "Observed variable Y") # Create second plot without axes
axis(1, at = d2$x, labels = d2$x)
###############Estimates for the Classical model##############
# Load data
mydata <- read.csv('H:\\Mythesisresearch\\latentdata.csv')
# View(mydata)  # Uncomment if you want to view the data
attach(mydata)
# Load necessary libraries
library(brms)
library(rstan)
library(StanHeaders)
# Fit classical model
model <- brm(formula = response ~ V1 + V2 + V3 + V4, data = mydata, family = cumulative)
summary(model)
# Extract posterior samples
posterior_samples <- posterior_samples(model)
head(posterior_samples)
# Compute point estimates
point_estimates <- apply(posterior_samples[, 8:11], 2, mean)
head(posterior_samples[, 8:11])
# Compute Mean Squared Error (MSE) and Root Mean Squared Error (RMSE)
x <- matrix(point_estimates, 4, 1)
y <- posterior_samples[, 8:11]
mse <- matrix(NA, 4000, 4)
for (i in 1:4000) {
  for (j in 1:4) {
    mse[i, j] <- (y[i, j] - x[j, 1])^2
  }
}
mse1 <- mean(mse)
RMSE <- sqrt(mse1)
# Compute Mean Absolute Error (MAE)
mae <- matrix(NA, 4000, 4)
for (i in 1:4000) {
  for (j in 1:4) {
    mae[i, j] <- abs(y[i, j] - x[j, 1])
  }
}
mae1 <- mean(mae)
# Compute Coverage Probability
param_intervals <- posterior_interval(model, prob = 0.95)
param_intervals_2 <- param_intervals[8:11, ]
point_estimates <- matrix(point_estimates, 4, 1)
coverage_prob <- matrix(NA, 4, 1) 
for (j in 1:4) {
  if (point_estimates[j, 1] >= param_intervals_2[j, 1] && 
      point_estimates[j, 1] <= param_intervals_2[j, 2]) {
    coverage_prob[j, 1] <- 1  
  } else {
    coverage_prob[j, 1] <- 0
  }
}
mean(coverage_prob)
