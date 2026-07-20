##' Define data from formula
##' @param data Data object in long format
##' @param location_formula Formula for location
##' @param scale_formula Formula for scale
##' @keywords internal
prepare_data_for_nimble <- function(data, location_formula, scale_formula) {
  
  ## Collapse a possibly multi-line deparse() into one string. deparse() wraps
  ## long formulas across lines; the parsing below assumed a single line, which
  ## silently dropped predictors / the grouping term on models with many terms.
  .flatten_formula <- function(f) paste(deparse(f), collapse = " ")

  ## Helper function to prepare model parts
  prepare_model_part <- function(data, formula, is_scale_model = FALSE) {
    ## Parse the formula to get response and predictors
    response_var <- if(is_scale_model) NA else all.vars(formula)[1]
    
    fixed_effects <- strsplit(.flatten_formula(formula), split = "\\+ \\(", perl = TRUE)[[1]][1]
    ## slplit out random effects, first split contains grouping variable
    random_effects_F <- strsplit(.flatten_formula(formula), split = "\\+ \\(", perl = TRUE)[[1]][2]
    ## split at | 
    random_effects <- strsplit(random_effects_F, split = "\\|", perl = TRUE)[[1]][1]

    predictors <- all.vars(formula)[-length(all.vars(formula) )]
    if (!is_scale_model) {
      predictors <- predictors[-1]  # Exclude the response variable for location model
    }    
    
    ## Creating X matrix
    X_formula <- update.formula(formula,   fixed_effects )
      #update.formula(formula, paste("~", paste(predictors, collapse = "*")))
    X_matrix <-  model.matrix(X_formula, data)
    
    ## For Z, random effects predictors
    Z_matrix <- if( !is.na(random_effects) ) {
                  model.matrix( formula( paste("~", random_effects) ), data)
                } else {
                  stop("Random effects missing")
                }
    list(X = X_matrix, Z = Z_matrix) # Adjusting for intercept
  }
  
  ## Extracting the grouping variable from the location formula
  location_formula_string <- .flatten_formula(location_formula)
  grouping_variable_match <- regmatches(location_formula_string, regexec("\\|\\s*(\\w+)", location_formula_string))
  if (length(grouping_variable_match[[1]]) < 2) {
    stop("Grouping variable not found in the location formula.")
  }
  grouping_variable <- grouping_variable_match[[1]][2]
  ## Extracting the grouping variable from the scale formula
  ## Only support models where grouping variable is the same for location and scale
  scale_formula_string <- .flatten_formula(scale_formula)
  scl_grouping_variable_match <- regmatches(scale_formula_string, regexec("\\|\\s*(\\w+)", scale_formula_string))
  if (length(scl_grouping_variable_match[[1]]) < 2) {
    stop("Grouping variable not found in the scale formula.")
  }
  ## Check that both location and scale have same grouping variable
  if(grouping_variable != scl_grouping_variable_match[[1]][2]) {
    stop("Location and scale grouping variable needs to be the same.")
  }
    
  ## Drop rows with missing values in any model variable so that the response,
  ## design matrices and grouping index stay aligned. Otherwise model.matrix()
  ## silently drops NA rows from X/Z while Y and group_id keep their full length.
  model_vars <- intersect(unique(c(all.vars(location_formula),
                                   all.vars(scale_formula))), names(data))
  keep <- stats::complete.cases(data[, model_vars, drop = FALSE])
  if (!all(keep)) {
    message("ivd: dropping ", sum(!keep),
            " row(s) with missing values in model variables.")
    data <- data[keep, , drop = FALSE]
  }

  ## Recode the grouping variable to the gap-free 1..J integer index NIMBLE
  ## needs, keeping the original labels so user-facing output (summary rows,
  ## plot labels) can report the user's own cluster IDs. factor() orders
  ## numeric IDs numerically and everything else alphabetically. Row order of
  ## `data` is never changed -- the model indexes u[group_id[i], ] per row, so
  ## rows need not be sorted by group.
  group_factor <- factor(data[[grouping_variable]])
  group_labels <- levels(group_factor)
  data[[grouping_variable]] <- as.integer(group_factor)
  
  ## Processing location and scale models
  location_data <- prepare_model_part(data, formula = location_formula)
  scale_formula_cleaned <- gsub("sigma = ", "", .flatten_formula(scale_formula))  # Remove "sigma = " if present
  scale_data <- if(!is.null(scale_formula_cleaned) && nzchar(scale_formula_cleaned)) {
                  prepare_model_part(data, formula = as.formula(paste(scale_formula_cleaned)), TRUE)
  } else {
    list(X = NULL, Z = NULL)
  }

  ## Check if repsonse_var has attributes, due to scaling with scale()
  ## Remove attributes, if present
  if( is.null(attributes( data[[all.vars(location_formula)[1]]] ) ) ) {
    Y <- data[[all.vars(location_formula)[1]]] # Assuming the first variable is the response
  } else if ( !is.null(attributes( data[[all.vars(location_formula)[1]]] ) )) {
    Y <- c(data[[all.vars(location_formula)[1]]]) # Assuming the first variable is the response
  }
  
  
  # Assemble the data structure for NIMBLE
  list(data = list(
         Y = Y,  
         X = location_data$X, 
         Z = location_data$Z, 
         X_scale = scale_data$X, 
         Z_scale = scale_data$Z
       ), 
       groups = length(unique(data[[grouping_variable]])),
       group_id = data[[grouping_variable]],
       group_labels = group_labels,
       response_var = all.vars(location_formula)[1]
  )
}
##' Render a single-line, carriage-return progress string for parallel chains
##'
##' Builds the live status line shown by `ivd(progress = TRUE)`: a spinner, the
##' number of chains being fit, and elapsed time. Deliberately has no completion
##' bar -- chains run in parallel and finish together, so a fraction bar would
##' sit at 0 then jump to full. Leads with "\\r" so repeated prints overwrite
##' the same terminal line.
##' @param total Integer number of chains (workers) being fit.
##' @param t0 Start time (`Sys.time()`), used to compute elapsed time.
##' @param spinner Optional single-character spinner frame.
##' @return A length-1 character string.
##' @keywords internal
.progress_line <- function(total, t0, spinner = "", workers = total) {
  el <- as.integer(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  elapsed <- sprintf("%02d:%02d", el %/% 60L, el %% 60L)
  on_workers <- if (workers < total) sprintf(" (%d workers)", workers) else ""
  sprintf("\r%s ivd: fitting %d %s%s | %s elapsed ",
          spinner, total, if (total == 1) "chain" else "chains", on_workers,
          elapsed)
}

##' Extract samples to mcmc object
##' @param obj
##' @return mcmc object
##' @author Philippe Rast
##' @keywords internal
.extract_to_mcmc <- function(obj) {
  e_to_mcmc <- lapply(obj$samples, FUN = function(x) mcmc(x$samples))
  return(e_to_mcmc)
}

##' Reconstruct the posterior mean of the location predictor `mu`
##'
##' `mu` is no longer monitored (it stores O(N x iterations) values per chain).
##' Because `mu[i] = X[i, ] %*% beta + Z[i, ] %*% u[group_i, 1:Kr]` is *linear*
##' in the monitored `beta` and `u`, the posterior mean of `mu` equals the
##' linear predictor evaluated at the posterior means of `beta` and `u` -- no
##' per-iteration `mu` storage required. Used by `plot.ivd()` for the cluster
##' outcome plot.
##' @param obj An `ivd` object (must carry the location design matrices `X`/`Z`).
##' @return Numeric vector of length N: the posterior mean of `mu` per observation.
##' @keywords internal
.reconstruct_mu_means <- function(obj) {
  if (is.null(obj$X) || is.null(obj$Z)) {
    stop("Cannot reconstruct cluster means: location design matrices (X, Z) ",
         "are missing from the ivd object. Refit with the current version of ivd().",
         call. = FALSE)
  }
  Kr <- obj$nimble_constants$Kr
  J  <- obj$nimble_constants$J

  ## Pool draws across chains; only beta and u columns are needed.
  all_draws <- do.call(rbind, .extract_to_mcmc(obj))
  cn <- colnames(all_draws)

  ## Fixed location effects beta[1..K], ordered by their numeric index so they
  ## line up with the columns of X.
  beta_means <- colMeans(all_draws[, grep("^beta\\[", cn), drop = FALSE])
  beta_means <- beta_means[order(as.integer(gsub("\\D", "", names(beta_means))))]

  ## Random location effects u[j, p], p <= Kr, as a J x Kr matrix of means.
  u_means <- colMeans(all_draws[, grep("^u\\[", cn), drop = FALSE])
  idx <- regmatches(names(u_means), gregexpr("[0-9]+", names(u_means)))
  jj <- as.integer(vapply(idx, `[`, character(1), 1)) # group index
  pp <- as.integer(vapply(idx, `[`, character(1), 2)) # random-effect index
  u_loc <- matrix(0, nrow = J, ncol = Kr)
  loc <- pp <= Kr
  u_loc[cbind(jj[loc], pp[loc])] <- u_means[loc]

  group_id <- obj$Y$group_id
  as.numeric(obj$X %*% beta_means) +
    rowSums(obj$Z * u_loc[group_id, , drop = FALSE])
}


##' Fast Fourier transform algorithm to compute the ACF across the whole chain lengt.
##' As noted by Vehtari et al. (2021), Section 3.2
##' @title Use fast Fourier transform to compute ACF
##' @param chain 
##' @return acf
##' @author Philippe Rast
##' @keywords internal
##' @importFrom stats fft nextn
.autocorrelation_fft <- function(chain) {
  ## Ensure the input is a numeric vector
  ts <- as.numeric(chain)

  ## Center the time series (subtract the mean)
  ts_centered <- ts - mean(ts)

  ## Length of the chain
  n <- length(ts_centered)

  ## Zero-padding the series to avoid circular convolution. fft() has no
  ## length argument (its 2nd argument is `inverse`), so pad explicitly;
  ## nextn() rounds up to a highly composite length for FFT speed.
  padded_length <- nextn(2 * n)
  ts_padded <- c(ts_centered, rep(0, padded_length - n))

  ## Compute the FFT of the centered, zero-padded series
  fft_ts <- fft(ts_padded)

  ## Compute the inverse FFT of the product of FFT and its conjugate
  acf_raw <- Re(fft(fft_ts * Conj(fft_ts), inverse = TRUE))

  ## Extract the relevant part and normalize
  acf_raw <- acf_raw[1:n]

  ## Normalize the result to match the acf() function output
  acf <- acf_raw / acf_raw[1]

  return(acf)
}

##' Truncate an autocorrelation sequence following Geyer (1992)
##'
##' Keeps the autocorrelations up to (and including) the first lag pair whose
##' sum is negative and pads the remainder with `NA` so that rho vectors from
##' chains with different truncation points can be averaged with
##' `rowMeans(..., na.rm = TRUE)`.
##' @title Geyer (1992) truncation of an ACF sequence
##' @param acf_values Autocorrelation sequence starting at lag 0, as returned
##'   by [.autocorrelation_fft()].
##' @return Numeric vector of `length(acf_values)`: the autocorrelations at
##'   lags 1, 2, ... up to the truncation point, padded with `NA`. All-`NA`
##'   when the ACF itself is undefined (constant chain).
##' @author Philippe Rast
##' @keywords internal
.geyer_truncate <- function(acf_values) {
  n <- length(acf_values)
  if (n < 2) return(rep(NA_real_, n))
  pair_sums <- acf_values[-n] + acf_values[-1]
  ## A constant chain has an undefined (NaN) ACF; treat as no usable lags.
  if (anyNA(pair_sums)) return(rep(NA_real_, n))
  crossings <- which(pair_sums < 0)
  ## When no pair sum ever goes negative (short or strongly autocorrelated
  ## chains) keep every available lag instead of erroring: min(integer(0))
  ## would return Inf and 1:Inf downstream crashed n_eff = "local".
  position <- if (length(crossings)) crossings[1] else n - 1
  c(acf_values[2:(position + 1)], rep(NA_real_, n - position))
}

##' Per-chain samples with human-readable parameter labels
##'
##' Restricts each chain's MCMC matrix to the parameters shown by
##' \code{summary()} (dropping per-cluster \code{u}, redundant \code{R} /
##' \code{sigma_rand} entries, the always-1 location \code{ss}, and any
##' monitored \code{mu}/\code{tau}) and renames the columns to the summary
##' labels (\code{beta} -> location names, \code{zeta} -> \code{scl_*},
##' \code{ss} -> \code{pip}, \code{(Intercept)} -> \code{Intc}, ...).
##' Shared by \code{codaplot()}, \code{as.mcmc.list.ivd()} and
##' \code{as_draws.ivd()}.
##' @title Renamed per-chain MCMC samples
##' @param obj An object of class `ivd`.
##' @return A list with one element per chain, each an iterations x
##'   parameters matrix of class `mcmc` with renamed columns.
##' @author Philippe Rast
##' @keywords internal
.renamed_mcmc_list <- function(obj) {
  ## Extract to mcmc object
  extract_samples <- .extract_to_mcmc(obj)
  Kr <- obj$nimble_constants$Kr
  ## Extract relevant names with summary_table function
  mat_transposed <- .summary_table(t(extract_samples[[1]]), Kr)

  ## Exclude mu and tau indexes (absent unless return_logLik = TRUE / legacy).
  ## Guard the empty case: `x[-integer(0), ]` selects ZERO rows, not all.
  drop_idx <- c(grep('^mu\\[', rownames(mat_transposed)),
                grep('^tau\\[', rownames(mat_transposed)))
  mat_kept <- if (length(drop_idx)) mat_transposed[-drop_idx, , drop = FALSE] else mat_transposed

  raw_internal_names <- rownames(mat_kept)
  internal_names <- raw_internal_names

  ## Location fixed effects
  beta_index <- grep('^beta\\[', internal_names)
  if(length(beta_index) > 0 && length(beta_index) == length(obj$X_location_names)) {
    internal_names[beta_index] <- obj$X_location_names
  }
  ## Scale fixed effects
  zeta_index <- grep('^zeta\\[', internal_names)
  if(length(zeta_index) > 0 && length(zeta_index) == length(colnames(obj$X_scale))) {
    internal_names[zeta_index] <- paste0("scl_", colnames(obj$X_scale))
  }
  ## Random effects SD (diagonal of sigma_rand)
  sigma_rand_index <- grep('^sigma_rand\\[(\\d+)]', internal_names) # Only diagonal
  sd_names <- c(obj$Z_location_names, paste0("scl_", colnames(obj$Z_scale)))
  if(length(sigma_rand_index) > 0 && length(sigma_rand_index) == length(sd_names)) {
    internal_names[sigma_rand_index] <- paste0("sd_", sd_names)
  }

  ## Rewrite correlation variable
  ## number of random effects:
  cols <- length(sigma_rand_index)
  ## Place holder variable: R is indexed as vech
  M <- matrix(1:cols^2, ncol = cols )
  ## record positions
  vech <- M[lower.tri(M)]
  ## match variable names to position
  corrvar <- expand.grid(sd_names, sd_names)[vech,]
  R_index <-  grep('R\\[', internal_names)
  if( nrow(corrvar) != length(R_index) )stop("Check R_index in summary.R" )
  internal_names[R_index] <- paste0("R[",paste(corrvar[, 1], corrvar[, 2], sep = ", "), "]")

  ## Link PIP to actual clustering units
  ## find the positions of the scale random effects in the model.
  ## Scale random effects occupy rows (Kr+1):(Kr+Sr) of u/ss -- offset by the
  ## number of *location* random effects Kr (using Sr here only worked when
  ## Kr == Sr).
  scale_ranef <- colnames(obj$Z_scale)
  scale_indexes <- seq_len(length(scale_ranef)) + obj$nimble_constants$Kr
  ## build patterns and replacements
  patterns <- paste0("\\[", scale_indexes, ",")
  replacements <- paste0("[", scale_ranef, ",")
  ## create a vector with the new rownames
  new_rownames <- Reduce(function(x, pattern_replacement) {
    gsub(pattern_replacement[1], pattern_replacement[2], x)
  },
  mapply(c, patterns, replacements, SIMPLIFY = FALSE),
  init = internal_names)
  ## assign back
  internal_names <- new_rownames

  pip_pos <- grep("ss", internal_names)
  internal_names[pip_pos] <- sub("^ss", "pip", internal_names[pip_pos])

  ## (Intercept) is annoying long. Change to Int.
  Int_index <- grep("\\(Intercept\\)", internal_names)
  internal_names[Int_index] <- gsub("\\(Intercept\\)",  "Intc", internal_names[Int_index])

  ## Filter each chain in the list for the relevant parameters
  lapply(extract_samples, function(chain) {
    chain_names <- colnames(chain)
    matching_cols <- chain_names %in% raw_internal_names
    filtered_chain <- chain[, matching_cols, drop = FALSE]
    colnames(filtered_chain) <- internal_names
    return(filtered_chain)
  })
}

##' Resolve and validate the `priors` argument of [ivd()]
##'
##' Fills a partial user specification with the package defaults and
##' validates it. Only the hyperparameters are tunable -- the families are
##' fixed (normal for fixed effects, half-t for random-effect SDs, LKJ for
##' the correlation): this keeps the NIMBLE model code static so that a
##' user's `priors` can never inject model code.
##' @title Resolve the prior specification for ivd()
##' @param priors Named list with any of `beta_intercept`, `beta`, `zeta`
##'   (each `c(mean = , sd = )`), `sigma_rand` (`c(df = , scale = )`),
##'   `lkj_eta` (single positive number), and `nu` (`c(shape = , rate = )`,
##'   student-t df prior). Partial specifications are filled with the
##'   defaults.
##' @param mean_pred,sd_pred Empirical mean and SD of the response, used for
##'   the default location-intercept prior.
##' @return Fully resolved named list of the same shape, with an
##'   `empirical_intercept` attribute recording whether the intercept prior
##'   kept its data-dependent default.
##' @author Philippe Rast
##' @keywords internal
.ivd_priors <- function(priors, mean_pred, sd_pred) {
  defaults <- list(
    beta_intercept = c(mean = mean_pred, sd = 3 * sd_pred),
    beta = c(mean = 0, sd = 1000),
    zeta = c(mean = 0, sd = 3),
    sigma_rand = c(df = 3, scale = 1),
    lkj_eta = 1,
    ## gamma prior of the student-t df (nu = 2 + gamma); used only when
    ## ivd(family = "student")
    nu = c(shape = 2, rate = 0.1)
  )
  if (is.null(priors)) priors <- list()
  if (!is.list(priors)) {
    stop("`priors` must be a named list; see ?ivd.", call. = FALSE)
  }
  if (length(priors) && (is.null(names(priors)) || any(!nzchar(names(priors))))) {
    stop("All elements of `priors` must be named.", call. = FALSE)
  }
  unknown <- setdiff(names(priors), names(defaults))
  if (length(unknown)) {
    stop("Unknown prior component(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(defaults), collapse = ", "), ".",
         call. = FALSE)
  }

  resolved <- defaults
  for (nm in names(priors)) {
    spec <- priors[[nm]]
    if (nm == "lkj_eta") {
      if (!is.numeric(spec) || length(spec) != 1 || is.na(spec)) {
        stop("`priors$lkj_eta` must be a single number.", call. = FALSE)
      }
      resolved$lkj_eta <- as.numeric(spec)
      next
    }
    if (!is.numeric(spec) || is.null(names(spec)) || any(!nzchar(names(spec)))) {
      stop("`priors$", nm, "` must be a named numeric vector, e.g. c(",
           paste(names(defaults[[nm]]), collapse = " = , "), " = ).",
           call. = FALSE)
    }
    bad <- setdiff(names(spec), names(defaults[[nm]]))
    if (length(bad)) {
      stop("Unknown element(s) in `priors$", nm, "`: ",
           paste(bad, collapse = ", "), ". Use ",
           paste(names(defaults[[nm]]), collapse = ", "), ".", call. = FALSE)
    }
    resolved[[nm]][names(spec)] <- spec
  }

  if (anyNA(unlist(resolved))) {
    stop("`priors` must not contain missing values.", call. = FALSE)
  }
  for (nm in c("beta_intercept", "beta", "zeta")) {
    if (resolved[[nm]][["sd"]] <= 0) {
      stop("`priors$", nm, "` needs sd > 0.", call. = FALSE)
    }
  }
  if (resolved$sigma_rand[["df"]] <= 0 || resolved$sigma_rand[["scale"]] <= 0) {
    stop("`priors$sigma_rand` needs df > 0 and scale > 0.", call. = FALSE)
  }
  if (resolved$lkj_eta <= 0) {
    stop("`priors$lkj_eta` must be positive.", call. = FALSE)
  }
  if (resolved$nu[["shape"]] <= 0 || resolved$nu[["rate"]] <= 0) {
    stop("`priors$nu` needs shape > 0 and rate > 0.", call. = FALSE)
  }

  attr(resolved, "empirical_intercept") <- !("beta_intercept" %in% names(priors))
  resolved
}
