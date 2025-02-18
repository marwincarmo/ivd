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
    ## Transpose L to get lower cholesky
    ## then compute the hadamard (element-wise) product with the ss vector
    u[j,1:P] <- t( sigma_rand[1:P, 1:P] %*% L[1:P, 1:P]  %*% z[1:P,j])
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
  
  ## Lower cholesky of random effects correlation 
  L[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  
  ## Construct L 'manually'
  # L[1,1] <- 1
  # L[2,1] <- rho
  # L[2,2] <- sqrt(1 - rho^2)
  # L[1,2] <- 0
  ##
  R[1:P, 1:P] <- t(L[1:P, 1:P] ) %*% L[1:P, 1:P]

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
                  sigma_rand = diag(rlnorm(constants$P, 0, 1)),
                  L = diag(1,constants$P) )

school_model <- nimbleModel(code = model_code, name = "school_model", constants = constants,
                    data = nimble_data, inits = inits)


mcmc.out <- nimbleMCMC(code = model_code, constants = constants,
                       data = nimble_data, inits = inits,
                       nchains = 4, niter = 10000, nburnin = 8000,
                       monitors = c("loc_int", "scl_int", "R", "sigma_rand"),
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

