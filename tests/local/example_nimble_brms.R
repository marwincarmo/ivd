library(mlmRev)
library(brms)
library(nimble)

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
mod0 <- brm( bf( mAch ~ 1 + ( 1|c|school),
                sigma ~ 1 + ( 1 |c| school ) ),
            cores = 4,
            data = school_dat)
summary(mod0)


# nimble LKJ model --------------------------------------------------------

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

pumpCode <- nimbleCode({ 
  
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], sd = tau[i])
    mu[i] <- loc_int + u[groupid[i], 1]
    tau[i] <-  exp( scl_int + u[groupid[i], 2] )
  }
  
  # for (j in 1:J) {
  #   for (p in 1:2) {
  #     u[j, p] ~ dnorm(0, 1)
  #   }
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
  for(i in 1:(P-1)) {
    for(j in (i+1):P) {
      zscore ~ dnorm(0, sd = 1)
      rho[i,j] <- tanh(zscore)  # correlations between effects i and j
      
    }
  }
  
  # Fixed intercept location
  loc_int ~ dnorm(0, sd = 10)
  # Fixed intercept scale
  scl_int ~ dnorm(0, sd = 10)
  ## Random effects SD
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 3), 0, )
  }
  ## Lower cholesky of random effects correlation 
  #L[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  L[1,1] <- 1
  L[2,1] <- rho[1,2]
  L[2,2] <- sqrt(1 - rho[1,2]^2)
  L[1,2] <- 0
  ##
  R[1:P, 1:P] <- t(L[1:P, 1:P] ) %*% L[1:P, 1:P]

})


pumpConsts <- list(N = nrow(school_dat),
                   K = 1,
                   P = 2,
                   J = length(unique(school_dat$school)),
                   groupid = school_dat$schoolid)

pumpData <- list(y = school_dat$mAch,
                 X = rep(1, nrow(school_dat)))

pumpInits <- list(loc_int = rnorm(1, 5, 10), ## TODO: Check inits
                  scl_int =  rnorm(1, 1, 3),
                  sigma_rand = diag(rlnorm(pumpConsts$P, 0, 1)),
                  L = diag(1,pumpConsts$P) )

pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

# pump$check()
# 
# pump$getNodeNames()
# pump$getVarNames()

mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits,
                       nchains = 4, niter = 8000, nburnin = 6000,
                       monitors = c("loc_int", "scl_int", "R", "sigma_rand"),
                      # monitors2 = c("sigma_rand"),
                       summary = TRUE, WAIC = TRUE)
#chain_summary <- mcmc.out$summary$all.chains

## Compute Rhats and n_eff:
# Convert nimbleMCMC output to an mcmc.list object
mcmc_samples <- mcmc.list(
  mcmc(mcmc.out$samples$chain1),
  mcmc(mcmc.out$samples$chain2)
)
summary(mcmc_samples)

# Compute R-hat values
rhat_values <- gelman.diag(mcmc_samples, multivariate = FALSE)$psrf[, 1]

# Compute effective sample size (ESS)
ess_values <- effectiveSize(mcmc_samples)

summary_stats <- summary(mcmc_samples)

# Create a summary table
summary_table <- data.frame(
  #Parameter = rownames(summary_stats$statistics),
  Estimate = summary_stats$statistics[, "Mean"],
  Est.Error = summary_stats$statistics[, "SD"],
  Low95CI = summary_stats$quantiles[, "2.5%"],
  Upp95CI = summary_stats$quantiles[, "97.5%"],
  Rhat = rhat_values,
  ESS = ess_values
)

summary_table
summary_stats$statistics
cn <- rownames(summary_stats$statistics )
tau_index <- grep('tau',  cn )
summary_stats$statistics[ -tau_index, ]
