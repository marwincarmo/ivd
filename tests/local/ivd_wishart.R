devtools::load_all( )

library(mlmRev)
library(nimble)
library(coda)
library(data.table)

# A. Prepare data ---------------------------------------------------------

## Grand mean center student SES
#saeb$student_ses <- c(scale(saeb$student_ses, scale = FALSE))

## Calculate school-level SES
school_ses <- saeb[, .(school_ses = mean(student_ses, na.rm = TRUE)), by = school_id]

## Join the school_ses back to the original dataset
saeb <- saeb[school_ses, on = "school_id"]

## Define student level SES as deviation from the school SES
saeb$student_ses <- saeb$student_ses - saeb$school_ses

## Grand mean center school ses
saeb$school_ses <- c(scale(saeb$school_ses, scale = FALSE))



# B. Model 1: Random intercept only ---------------------------------------


# m1 <- ivd(location_formula = math_proficiency ~ 1 + (1|school_id),
#           scale_formula =  ~ 1 + (1|school_id),
#           data = saeb,
#           niter = 3000, nburnin = 12000, WAIC = TRUE, workers = 6)

location_formula = math_proficiency ~ 1 + (1|school_id)
scale_formula =  ~ 1 + (1|school_id)
data = saeb
niter = 100
nburnin = 100
WAIC = TRUE
workers = 4
n_eff = 'local'

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
  config$monitors <- c("beta", "zeta", "R", "ss", "sigma_rand", "u")
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
              identityMatrix = diag(1,constants$P),
              scale_matrix <- diag(constants$P) * 0.1,  # Example: identity matrix scaled by 0.1
              df_value <- constants$P + 1  # Degrees of freedom, must exceed P - 1
)

modelCode <- nimbleCode({
  ## Likelihood components:
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], sd = tau[i]) ## explicitly ask for SD not precision
    ## Check if K (number of fixed location effects) an S (number of fixed scale effecs)
    ## are greater than 1, if not, use simplified computation to avoid indexing issues in nimble
    ## Location
    ## Check if we have more than just an intercept:
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
    ## Check if we have more than just an fixed intercept:
    if(S>1) { 
      if(Sr>1) {
        tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + sum(u[groupid[i], (Kr+1):(Kr+Sr)] * Z_scale[i, 1:Sr]) )  
      } else {
        tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + u[groupid[i], (Kr+1)] ) 
      }
    } else {
      ## This assumes that if there is only one fixed intercept in scale, there is also exactly one random intercept in scale,
      ## and no other effects
      tau[i] <- exp( zeta[1] + u[groupid[i], (Kr+1)] )
    }
  }
  ## Obtain correlated random effects
  for(j in 1:J) {
    ## Bernoulli for Spike and Slab
    for(p in 1:P){
      ss[p,j] ~ dbern(bval[p,1]) ## bval is a constant
    }    
    ## normal scaling for random effects
    for( k in 1:P ){
      z[k,j] ~ dnorm(0, sd = 1)
    }
    ## Transpose L to get lower cholesky
    ## then compute the hadamard (element-wise) product with the ss vector
    u[j,1:P] <- t( sigma_rand[1:P, 1:P] %*% L[1:P, 1:P]  %*% z[1:P,j] * ss[1:P,j] )
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
  ## Random effects SD
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 3), 0, )
  }
  ## Lower cholesky of random effects correlation 
  #L[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  L[1:P, 1:P] ~ dinvwish_chol(cholesky = chol(scale_matrix), df = df_value, scale_param = TRUE)
  
  ##
  R[1:P, 1:P] <- t(L[1:P, 1:P] ) %*% L[1:P, 1:P]
  #R[1:P, 1:P] <- diag(1 / sqrt(diag(L %*% t(L)))) %*% (L %*% t(L)) %*% diag(1 / sqrt(diag(L %*% t(L))))
  
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
head(combined_chains[[1]]$samples[,2])

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

## Split Rhat and split n_eff:
## Vehtari et al doi:10.1214/20-BA1221 available at
## http://www.stat.columbia.edu/~gelman/research/published/rhat.pdf

## Initialize a new array with double the chains, half the iterations
split_samples <- array(NA, dim = c(iterations / 2, chains * 2, parameters))

## Split each chain into two halves
for (c in 1:chains) {
  ## First half of the iterations for the first split chain
  split_samples[, (c * 2) - 1, ] <- samples_array[1:(iterations / 2), c, ]
  ## Second half of the iterations for the second split chain
  split_samples[, c * 2, ] <- samples_array[(iterations / 2 + 1):iterations, c, ]
}

## m - number of chains after splitting
m <- dim(split_samples )[2]
## n be the length of th chain
n <- dim(split_samples )[1]
s <- m*n

## B ingredients
tdm <- apply(split_samples, 3, function(slice) {
  apply(slice, 2, mean)
})
tdd <- colMeans(tdm)  

## Eq. 3.1
result <- tdm - matrix(tdd, nrow = nrow(tdm), ncol = ncol(tdm), byrow = TRUE)
B <- apply(result, 2, function(x) sum(x^2)*n/(m-1))

## Eq. 3.2
W <- apply(split_samples, 3, function(param_samples) {
  chain_variances <- apply(param_samples, 2, function(chain) {
    chain_mean <- mean(chain)
    sum((chain - chain_mean)^2) / (n - 1)  # s2m: Variance for each chain
  })
  mean(chain_variances)  # Average variance across all split chains
})

## Eq. 3.3
vtp <- (n-1)*W/n + B/n

## Compute R-hat
Rhat <- sqrt(vtp / W)

if (n_eff == 'local') {
  ## Compute split-chain n_eff
  ## ACF is computed using FFT as per Vehtari et al. 
  ## Compute for multiple chains, following eq 10:
  
  mn_s2m_ptm <- apply(split_samples, 3, function(param_samples) {
    
    chain_variances <- apply(param_samples, 2, function(chain) {
      chain_mean <- mean(chain)
      sum((chain - chain_mean)^2) / (n - 1)  # s2m: Variance for each chain
    })
    
    chain_rho <- apply(param_samples, 2, function(samp_per_chain) {
      acf_values <- .autocorrelation_fft( samp_per_chain )
      ## Truncate according to Geyer (1992)
      position <-  min( seq(2:length(acf_values))[acf_values[-length(acf_values)] + acf_values[-1] < 0] )
      ## position contains NA for constants, needs to be addressed here:
      
      if( !is.na(position) ) {
        ## Pad with NA's so that all vectors are of same length. Saves me storing the position object
        ## pad with NA so that mean() can be calculated over differing rho's per chains
        rho <- append(acf_values[1:position+1], rep(NA, length(acf_values)-position), after = position)
      } else {
        rho <- rep(NA, n)
      }
    })
    
    s2m_rtm <- lapply(seq_along(chain_variances), function(i) {
      chain_variances[i] * chain_rho[,i]
    })
    
    ## average across chains    
    ## Convert list to a matrix
    matrix_form <- do.call(cbind, s2m_rtm)
    avg_s2m_rtm <- rowMeans(matrix_form, na.rm = TRUE)
  })
  
  ## Eq 10: W - mn_s2m_ptm
  numerator <- matrix(W, nrow = nrow(mn_s2m_ptm), ncol = length(W), byrow = TRUE) - mn_s2m_ptm
  rho_t <- 1 - numerator / matrix(vtp, nrow = nrow(numerator), ncol = length(vtp), byrow = TRUE)
  
  n_eff <- round( n*m / ( 1+2*colSums(rho_t, na.rm = TRUE) ) )
  
}  else if(n_eff == 'stan') {
  ## Based on rstan, takes forever...
  monitor_results <- rstan::monitor(samples_array, print = FALSE)
  n_eff <- monitor_results$n_eff
}
## Extract and print R-hat values
out$rhat_values <- Rhat
if( any(out$rhat_values > 1.1) ) warning("Some R-hat values are greater than 1.10 -- increase warmup and/or sampling iterations." )

## Effective sample size
out$n_eff <- n_eff

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

