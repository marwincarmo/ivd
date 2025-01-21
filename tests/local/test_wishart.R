library(nimble)
# Now your model code should work
modelCode <- nimbleCode({
  ## Likelihood components:
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], sd = tau[i])
    if(K>1) {
      if(Kr>1) {
        mu[i] <- sum(beta[1:K] * X[i, 1:K]) + sum( u[groupid[i], 1:Kr] * Z[i, 1:Kr] )
      } else {
        mu[i] <- sum(beta[1:K] * X[i, 1:K]) + u[groupid[i], 1]
      }
    } else {
      mu[i] <- beta[1] + u[groupid[i], 1] * Z[i, 1]        
    }
    
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
  
  for (k in 1:K) {
    beta[k] ~ dnorm(0, sd = 1000)
  }
  for (s in 1:S) {
    zeta[s] ~ dnorm(0, sd = 1000)
  }  
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 3), 0, )
  }
  
  cholesky <- chol(scale_matrix)
  # Make sure to specify scale_param
  L[1:P, 1:P] ~ dinvwish(cholesky = cholesky,
                         df = 5,
                         scale_param = 1)
  
  R[1:P, 1:P] <- t(L[1:P, 1:P]) %*% L[1:P, 1:P]
  ## Lower cholesky of random effects correlation 
  
})

library(mlmRev)

school_dat = mlmRev::Hsb82

## Ensure that school id is a continuous vector
school_dat$schoolid <- NA
k <- 0
for( i in unique(school_dat$school) ) {
  k <- k+1
  school_dat[school_dat$school == i, "schoolid"] <- k
}
school_dat$mAch_s <- scale(school_dat$mAch,  center = TRUE,  scale = TRUE )
school_dat$ses_s <- scale(school_dat$ses)

dat <- prepare_data_for_nimble(data = school_dat, location_formula = mAch_s ~  meanses + ( 1 | schoolid),
                               scale_formula =  ~ meanses  + (1 | schoolid))
data = dat[[1]]
groups <- dat$groups
group_id <- dat$group_id
constants <- list(N = length(data$Y),
                  J = groups,
                  K = ncol(data$X),  ## number of fixed location effects
                  Kr = ncol(data$Z), ## number of random location effects
                  S = ncol(data$X_scale),  ## number of fixed scale effects
                  Sr = ncol(data$Z_scale),  ## number of random scale effects                    
                  P = ncol(data$Z) + ncol(data$Z_scale),  ## number of random effects
                  groupid = group_id,
                  bval = matrix(c(rep(1,  ncol(data$Z)), rep(0.5, ncol(data$Z_scale)) ), ncol = 1)) ## Prior probability for dbern 
# Initialize values
inits <- list(beta = rnorm(constants$K, 5, 10), ## TODO: Check inits
              zeta =  rnorm(constants$S, 1, 3),
              sigma_rand = diag(rlnorm(constants$P, 0, 1)),
              L = diag(1,constants$P) )

# Try to build the model
try({
  test_model <- nimbleModel(
    code = modelCode,
    data = data,
    constants = constants,
    inits = inits
  )
  
  # If successful, compile it
  c_test_model <- compileNimble(test_model)
  
  # Calculate initial log probability
  cat("Initial log probability:", c_test_model$calculate(), "\n")
})


## Compile baseline model
cmpModel <- nimble::compileNimble(test_model)

