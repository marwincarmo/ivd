library(mlmRev)
library(brms)
library(nimble)
library(coda)

school_dat = mlmRev::Hsb82

## Ensure that school id is a continuous vector
school_dat$schoolid <- NA
k <- 0
for( i in unique(school_dat$school) ) {
  k <- k+1
  school_dat[school_dat$school == i, "schoolid"] <- k
}

# brms model --------------------------------------------------------------
# 
mod0 <- brm( bf( mAch ~ 1 + ( 1|c|schoolid),
                sigma ~ 1 + ( 1 |c| schoolid ) ),
            cores = 4, iter = 10000, warmup = 8000,
            data = school_dat)
summary(mod0)


# nimble model --------------------------------------------------------

school_dat = mlmRev::Hsb82

## Ensure that school id is a continuous vector
school_dat$schoolid <- NA
k <- 0
for( i in unique(school_dat$school) ) {
  k <- k+1
  school_dat[school_dat$school == i, "schoolid"] <- k
}

model_code <- nimbleCode({ 
  
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], sd = tau[i])
    mu[i] <- loc_int + u[groupid[i], 1]
    tau[i] <-  exp( scl_int + u[groupid[i], 2] )
  }
  
  ## Obtain correlated random effects
  for(j in 1:J) {
    ## normal scaling for random effects
    for( k in 1:P ){
      z[k,j] ~ dnorm(0, sd = 1)
    }
    ## Transpose U to get lower cholesky
    ## then compute the hadamard (element-wise) product with the ss vector
    u[j,1:P] <- t( sigma_rand[1:P, 1:P] %*% U[1:P, 1:P]  %*% z[1:P,j])
  }
  
  # Ranef cor prior
  zscore ~ dnorm(0, sd = 1)
  rho <- tanh(zscore)
  
  # Fixed intercept location
  loc_int ~ dnorm(0, sd = 10)
  # Fixed intercept scale
  scl_int ~ dnorm(0, sd = 10)
  ## Random effects SD
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 3), 0, )
  }
  
  ## Upper cholesky of random effects correlation 
  #U[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  #L[1:P, 1:P] <- t(U[1:P, 1:P])
  
  ## Construct L 'manually'
  L[1,1] <- 1
  L[2,1] <- rho
  L[2,2] <- sqrt(1 - rho^2)
  L[1,2] <- 0
  

  ## Construct U 'manually'
  U[1,1] <- 1
  U[2,1] <- 0
  U[2,2] <- sqrt(1 - rho^2)
  U[1,2] <- rho
  ##
  #R[1:P, 1:P] <- t(U[1:P, 1:P] ) %*% U[1:P, 1:P] # using upper cholesky
  R[1:P, 1:P] <- L[1:P, 1:P] %*% t(L[1:P, 1:P] ) # using lower cholesky
  
})


constants <- list(N = nrow(school_dat),
                   K = 1,
                   P = 2,
                   J = length(unique(school_dat$school)),
                   groupid = school_dat$schoolid)

nimble_data <- list(y = school_dat$mAch,
                 X = rep(1, nrow(school_dat)))

inits <- list(loc_int = rnorm(1, 5, 10), ## TODO: Check inits
                  scl_int =  rnorm(1, 1, 3),
                  sigma_rand = diag(rlnorm(constants$P, 0, 1))
                , L = diag(1,constants$P) 
              )

school_model <- nimbleModel(code = model_code, name = "school_model", constants = constants,
                    data = nimble_data, inits = inits)


mcmc.out <- nimbleMCMC(code = model_code, constants = constants,
                       data = nimble_data, inits = inits,
                       nchains = 4, niter = 10000, nburnin = 8000,
                       monitors = c("loc_int", "scl_int", "R", "sigma_rand", "L"#, "U"
                                    ),
                       summary = TRUE, WAIC = TRUE)

## Compute Rhats and n_eff:
# Convert nimbleMCMC output to an mcmc.list object
mcmc_chains <- mcmc.list(lapply(mcmc.out$samples, mcmc))

# Compute R-hat values
rhat_values <- gelman.diag(mcmc_chains, multivariate = FALSE)$psrf[, 1]

# Compute effective sample size (ESS)
ess_values <- effectiveSize(mcmc_chains)

summary_stats <- summary(mcmc_chains)

summary_table <- data.frame(
  Estimate = summary_stats$statistics[, "Mean"],
  Est.Error = summary_stats$statistics[, "SD"],
  Low95CI = summary_stats$quantiles[, "2.5%"],
  Upp95CI = summary_stats$quantiles[, "97.5%"],
  Rhat = rhat_values,
  ESS = ess_values
)

round(summary_table,2)

# Reconstructing R from U -------------------------------------
# > round(summary_table,2)
# Estimate Est.Error Low95CI Upp95CI Rhat    ESS
# R[1, 1]              1.00      0.00    1.00    1.00  NaN   0.00
# R[2, 1]             -0.75      0.15   -0.98   -0.41 2.13  51.19
# R[1, 2]             -0.75      0.15   -0.98   -0.41 2.13  51.19
# R[2, 2]              1.00      0.00    1.00    1.00 1.10   0.00
# U[1, 1]              1.00      0.00    1.00    1.00  NaN   0.00
# U[2, 1]              0.00      0.00    0.00    0.00  NaN   0.00
# U[1, 2]             -0.75      0.15   -0.98   -0.41 2.13  51.19
# U[2, 2]              0.61      0.19    0.21    0.91 2.16  46.22
# loc_int             12.61      0.24   12.13   13.08 1.07 221.19
# scl_int              1.83      0.01    1.81    1.85 1.01 550.65
# sigma_rand[1, 1]     2.33      0.21    1.97    2.77 1.53 124.88
# sigma_rand[2, 1]     0.00      0.00    0.00    0.00  NaN   0.00
# sigma_rand[1, 2]     0.00      0.00    0.00    0.00  NaN   0.00
# sigma_rand[2, 2]     0.16      0.08    0.09    0.39 1.60  84.82

# Reconstructing R from L -------------------------------------
# > round(summary_table,2)
# Estimate Est.Error Low95CI Upp95CI Rhat    ESS
# R[1, 1]              1.00      0.00    1.00    1.00  NaN   0.00
# R[2, 1]             -0.55      0.13   -0.80   -0.29 1.01 506.46
# R[1, 2]             -0.55      0.13   -0.80   -0.29 1.01 506.46
# R[2, 2]              1.00      0.00    1.00    1.00 1.00   0.00
# L[1, 1]              1.00      0.00    1.00    1.00  NaN   0.00
# L[2, 1]             -0.55      0.13   -0.80   -0.29 1.01 506.46
# L[1, 2]              0.00      0.00    0.00    0.00  NaN   0.00
# L[2, 2]              0.82      0.09    0.60    0.96 1.01 452.23
# loc_int             12.66      0.24   12.22   13.15 1.06 159.59
# scl_int              1.83      0.01    1.81    1.85 1.01 876.16
# sigma_rand[1, 1]     2.94      0.18    2.63    3.33 1.04 178.83
# sigma_rand[2, 1]     0.00      0.00    0.00    0.00  NaN   0.00
# sigma_rand[1, 2]     0.00      0.00    0.00    0.00  NaN   0.00
# sigma_rand[2, 2]     0.09      0.01    0.06    0.12 1.02 459.42

# brms estimates ----------------------------------------------
# ~schoolid (Number of levels: 160) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                      2.97      0.19     2.62     3.35 1.01     1501     2837
# sd(sigma_Intercept)                0.09      0.01     0.07     0.12 1.00     3036     4491
# cor(Intercept,sigma_Intercept)    -0.54      0.13    -0.77    -0.27 1.00     3394     4199
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          12.64      0.24    12.15    13.12 1.00      848     1541
# sigma_Intercept     1.83      0.01     1.81     1.85 1.00     3046     5427