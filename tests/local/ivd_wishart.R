devtools::load_all( )

library(mlmRev)
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

school_dat$mAch_s <- scale(school_dat$mAch,  center = TRUE,  scale = TRUE )

location_formula = mAch_s ~  meanses + ( 1 | schoolid)
scale_formula =  ~ meanses  + (1 | schoolid)
data = school_dat
niter = 2000
nburnin = 2000
WAIC = TRUE
workers = 4

run_MCMC_allcode <- function(seed, data, constants, code, niter, nburnin, useWAIC = WAIC, inits, thin, ...) {
  ## See Nimble cheat sheet: https://r-nimble.org/cheatsheets/NimbleCheatSheet.pdf
  
  
  ## Create model object
  myModel <- nimble::nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = inits)
  ## Compile baseline model
  cmpModel <- nimble::compileNimble(myModel)
  
  
  ## Configure MCMC
  config <- nimble::configureMCMC(myModel)
  ## Enable WAIC if useWAIC is TRUE
  if (useWAIC) {
    config$enableWAIC <- useWAIC
  }
  config$monitors <- c("beta", "zeta", "ss", "sigma_rand", "u", "R", "chol_Sigma")
  config$addMonitors(c("mu", "tau"))
  
  ## build mcmc object
  myMCMC <- nimble::buildMCMC(config)
  
  ## Recompile myMCMC linking it to cmpModel
  compMCMC <- nimble::compileNimble(myMCMC, project = cmpModel)
  
  ## Run model
  results <- nimble::runMCMC(compMCMC, niter = niter, setSeed = seed, nburnin = nburnin, WAIC = useWAIC, thin = thin, ...)
  return( results )
}

uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })


  niter <- niter + nburnin
  dat <- prepare_data_for_nimble(data = data, location_formula = location_formula, scale_formula = scale_formula)
  data <- dat[[1]]
  groups <- dat$groups
  group_id <- dat$group_id
  
  ## Nimble part:
  ## Nimble constants
  constants <- list(N = length(data$Y),
                    J = groups,
                    K = ncol(data$X),  ## number of fixed location effects
                    Kr = ncol(data$Z), ## number of random location effects
                    S = ncol(data$X_scale),  ## number of fixed scale effects
                    Sr = ncol(data$Z_scale),  ## number of random scale effects                    
                    P = ncol(data$Z) + ncol(data$Z_scale),  ## number of random effects
                    groupid = group_id,
                    bval = matrix(c(rep(1,  ncol(data$Z)), rep(0.5, ncol(data$Z_scale)) ), ncol = 1)) ## Prior probability for dbern 
  ## Nimble inits
  inits <- list(beta = rnorm(constants$K, 5, 10), ## TODO: Check inits
                zeta =  rnorm(constants$S, 1, 3),
                sigma_rand = rlnorm(constants$P, 0, 1),
                L = diag(1,constants$P),
                identityMatrix = (diag(1,constants$P)),
                df = constants$P )
  
  #identityMatrix <- (diag(1,constants$P))
  #df <- constants$P + 1
  modelCode <- nimbleCode({
    ## Likelihood components:
    for(i in 1:N) {
      Y[i] ~ dnorm(mu[i], sd = tau[i]) ## explicitly ask for SD not precision
      
      ## Location
      if(K>1) {
        if(Kr>1) {
          mu[i] <- sum(beta[1:K] * X[i, 1:K]) + sum( u[groupid[i], 1:Kr] * Z[i, 1:Kr] )
        } else {
          mu[i] <- sum(beta[1:K] * X[i, 1:K]) + u[groupid[i], 1]
        }
      } else {
        mu[i] <- beta[1] + u[groupid[i], 1] * Z[i, 1]
      }
      
      ## Scale
      if(S>1) { 
        if(Sr>1) {
          tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + sum(u[groupid[i], (Kr+1):(Kr+Sr)] * Z_scale[i, 1:Sr]) )
        } else {
          tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + u[groupid[i], (Kr+1)] ) 
        }
      } else {
        tau[i] <- exp( zeta[1] + u[groupid[i], (Kr+1)] )
      }
    }
    
    ## Obtain correlated random effects
    for(j in 1:J) {
      ## Bernoulli for Spike and Slab
      for(p in 1:P){
        ss[p,j] ~ dbern(bval[p,1])
      }
      
      ## normal scaling for random effects
      for(k in 1:P){
        z[k,j] ~ dnorm(0, sd = 1)
      }
      
      ## Cholesky decomposition of covariance matrix
      u[j,1:P] <- chol_Sigma[1:P, 1:P] %*% z[1:P,j] * ss[1:P,j]
    }
    
    ## Priors:
    ## Fixed effects: Location
    for (k in 1:K) {
      beta[k] ~ dnorm(0, sd = 1000)
    }
    ## Fixed effects: Scale
    for (s in 1:S) {
      zeta[s] ~ dnorm(0, sd = 1000)
    }
    
    ## Random effects covariance matrix
    R[1:P, 1:P] ~ dwish(identityMatrix[1:P, 1:P], df) ## Wishart prior with df degrees of freedom and R_prior scale matrix
    
    ## Cholesky decomposition of R
    chol_Sigma[1:P, 1:P] <- chol(R[1:P, 1:P])
    
    ## Random effects standard deviations (optional, can extract from Sigma if needed)
    for(p in 1:P){
      sigma_rand[p] <- sqrt(R[p, p])
    }
  })
  
future::plan(multisession, workers = workers)

results <- future_lapply(1:workers, function(x) run_MCMC_allcode(seed = x, 
                                                                   data = data, 
                                                                   constants = constants,
                                                                   code = modelCode, 
                                                                   niter = niter, 
                                                                   nburnin = nburnin,
                                                                   useWAIC = WAIC, 
                                                                   inits = inits),
                           future.seed = TRUE, 
                           future.packages = c("nimble"), 
                           future.globals = list(
                             uppertri_mult_diag = uppertri_mult_diag,  # Custom nimble function
                             run_MCMC_allcode = run_MCMC_allcode,      # Custom MCMC function
                             data = data,                              # Data object
                             constants = constants,                    # Constants list
                             modelCode = modelCode,                    # Nimble model code
                             niter = niter,                            # Number of iterations
                             nburnin = nburnin,                        # Burn-in
                             WAIC = WAIC,                              # WAIC option
                             inits = inits                             # Initial values
                           ))

out <- list()
mcmc_chains <- lapply(results, as.mcmc)
combined_chains <- mcmc.list(mcmc_chains)

## Compute logLik:
## Check that Y,  mu and tau are of same length, in case grep picks up other variables
if(length(grep("mu", colnames(combined_chains[[1]]$samples))) != length(grep("tau", colnames(combined_chains[[1]]$samples))) &
   length(grep("mu", colnames(combined_chains[[1]]$samples))) != length(data$Y)) {
  stop("mu and tau are not of same lenght -- check ivd.R")
}

## Collect mu and tau
## Get mu's across chains
mu_combined <- lapply(combined_chains, function(chain) {
  mu_indices <- grep("mu", colnames(chain$samples))
  mu_samples <- chain$samples[, mu_indices, drop = FALSE]
  return(mu_samples)
})

## Get tau's across chains
tau_combined <- lapply(combined_chains, function(chain) {
  tau_indices <- grep("tau", colnames(chain$samples))
  tau_samples <- chain$samples[, tau_indices, drop = FALSE]
  return(tau_samples)
})

N <- length( data$Y )
chains <- length(mu_combined)  # Number of chains
iterations <- nrow(mu_combined[[1]])  # Number of iterations (assuming all chains have same iterations)

## Compute R hats:
x <- mcmc.list( lapply(combined_chains, FUN = function(x) mcmc(x$samples)) )
## Extract dimensions
#iterations <- nrow(cleaned_chains[[1]])
parameters <- ncol(x[[1]])
#chains <- length(cleaned_chains)
## Initialize a 3D array
samples_array <- array(NA, dim = c(iterations, chains, parameters))

## Fill the 3D array with the data from the list
for (i in seq_along(x)) {
  samples_array[, i, ] <- x[[i]]
}
## Use the monitor function from rstan to obtain Rhat (coda's gelman.rhat does not work reliably)
print("Compiling results...")
monitor_results <- rstan::monitor(samples_array, print = FALSE)
## Extract and print R-hat values
out$rhat_values <- monitor_results[, "Rhat"]
if( any(out$rhat_values > 1.1) ) warning("Some R-hat values are greater than 1.10 -- increase warmup and/or sampling iterations." )

## Effective sample size
out$n_eff <- monitor_results[, "n_eff"]

## Save the rest to the out object
out$samples <- combined_chains
out$nimble_constants <- constants
out$X_location_names <- colnames(data$X) # save fixed effects names for summary table renaming
out$X_scale <- data$X_scale
out$Z_location_names <- colnames(data$Z) # save random effects names for summary table renaming
out$Z_scale <- data$Z_scale
out$Y <- data.frame("group_id" = group_id, "Y" = data$Y)
out$workers <- workers

class(out) <- c("ivd", "list")



# Summary table -----------------------------------------------------------

summary(out)

