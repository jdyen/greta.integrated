# internal function: define age-abundance log-likelihood
age_abundance_loglik <- function(data, params) {

  # unpack params
  n_obs <- params$n_obs
  classes <- params$classes
  density <- params$density
  mat <- params$matrix
  inits <- params$inits
  process_class <- params$process_class
  age_to_stage <- params$age_to_stage_conversion

  # set up iterated states
  if (length(dim(matrix)) == 3) {
    iterated_states <- iterate_matrix_dynamic(matrix = mat, initial_state = inits, density = density)
  } else {
    n_iter <- ncol(data$data)
    iterated_states <- iterate_matrix(matrix = mat, initial_state = inits, niter = n_iter, density = density)
  }

  # add a conversion step if we have stage-structured data
  if (process_class != "leslie") {
    modelled_states <- age_to_stage %*% iterated_states
  } else {
    modelled_states <- iterated_states
  }
    
  # create vectors of modelled and observed data
  mu <- c(modelled_states)
  observed_tmp <- as_data(c(data$data))
  
  # add bias transform if required
  observed <- data$bias$bias(observed_tmp, data$bias$params)

  # size obs model
  distribution(observed) <- do.call(data$likelihood, list(mu))
  
  # could return mu at some point but want to keep clear for now
  NULL
  
}

# internal function: define stage-abundance log-likelihood
stage_abundance_loglik <- function(data, params) {
  
  # unpack params
  n_obs <- params$n_obs
  classes <- params$classes
  density <- params$density
  mat <- params$matrix
  inits <- params$inits
  process_class <- params$process_class
  stage_to_age <- params$stage_to_age_conversion

  # set up iterated states
  if (length(dim(matrix)) == 3) {
    iterated_states <- iterate_matrix_dynamic(matrix = mat, initial_state = inits, density = density)
  } else {
    n_iter <- ncol(data$data)
    iterated_states <- iterate_matrix(matrix = mat, initial_state = inits, niter = n_iter, density = density)
  }
  
  # add a conversion step if we have stage-structured data
  if (process_class == "leslie") {
    modelled_states <- stage_to_age %*% iterated_states
  } else {
    modelled_states <- iterated_states
  }
  
  # create vectors of modelled and observed data
  mu <- c(modelled_states$all_states)
  observed_tmp <- as_data(c(data$data))
  
  # add bias transform if required
  observed <- data$bias$bias(observed_tmp, data$bias$params)
  
  # size obs model
  distribution(observed) <- do.call(data$likelihood, list(mu))
  
  NULL
  
}

# internal function: define age-recapture log-likelihood
binned_age_recapture_loglik <- function(data, params) {
  
  ## TWO POSSIBILITIES: observed ages, modelled ages; observed ages, modelled stages
  process_class <- params$process_class
  if (process_class == "stage") {
    # we're modelling stages, need a conversion
  } else {
    # we're modelling ages, no conversion needed
  }
  
  # unpack params
  survival <- params$survival
  stage_to_age <- params$stage_to_age_conversion
  classes <- params$classes
  bias <- params$bias
  
  # calculate first and final observations
  first_obs <- apply(data$data, 1, function(x) min(which(x > 0)))
  final_obs <- apply(data$data, 1, function(x) max(which(x > 0)))
  
  # pull out first and final ages
  first_class <- apply(data$data, 1, function(x) x[min(which(x > 0))])
  final_class <- apply(data$data, 1, function(x) x[max(which(x > 0))])

  # number of years alive
  n_alive <- final_obs - first_obs
  
  # are any individuals never recaptured?
  single_obs <- first_obs == final_obs
  
  # focus on those with >1 observation
  n_alive <- n_alive[!single_obs]
  first_class <- first_class[!single_obs]
  first_seen <- first_obs[!single_obs]
  last_seen <- final_obs[!single_obs]
  
  # probabiity of stage j given age k (should be n_obs x n_stage)
  p_stage_first_capture <- stage_to_age[first_class, ]
  
  # create Phi = P(surv_n_alive | stage_at_time_1) (n_alive x n_stage)
  max_alive <- max(n_alive)
  survival <- cbind(survival,
                    do.call(cbind, lapply(seq_len(max_alive), function(x) survival[, classes])))
  
  # pre-calculate all possible survival trajectories
  first_seen_lifespan <- matrix(as.numeric(unlist(strsplit(unique(paste(first_seen, n_alive, sep = '_')), '_'))), ncol = 2, byrow = TRUE)
  id_match <- match(paste(first_seen, n_alive, sep = '_'), unique(paste(first_seen, n_alive, sep = '_')))
  idx <- apply(first_seen_lifespan, 1, function(x) x[1]:(x[1] + x[2]))
  surv_mat <- ones(nrow(first_seen_lifespan), classes)
  for (i in seq_len(nrow(first_seen_lifespan))) {
    mat_tmp <- survival[idx[[i]], ]
    out <- NULL
    for (j in seq_len(classes)) {
      out <- c(out, seq(1, length(idx[[i]]) ^ 2, by = (length(idx[[i]]) + 1)) + (j - 1) * length(idx[[i]]))
    }
    idy <- rep(seq_len(classes), each = length(idx[[i]]))
    surv_mat[i, ] <- tapply(mat_tmp[out], idy, "prod")
  }
  
  # what are survival probabilities of each stage at first sight?
  p_survival_hist_stage <- surv_mat[id_match, ]
  
  # calculate P(survival_history | age = k)
  p_survival_hist_stage <- rowSums(p_survival_hist_stage * p_stage_first_capture)
  
  # probs of binary capture history, assume detection is constant
  binary_capture_history <- apply(data$data, 1, function(x) x[min(which(x > 0)):max(which(x > 0))])
  binary_capture_history <- ifelse(do.call(c, binary_capture_history) > 0, 1, 0)
  
  # probs of final detection (for each age??)
  # p(never_observed_again | final_size, final_obs) = \Sum_ages p(final_age | final_size) p(not_observed | final_obs, final_age) 
  # p_final_obs <- SOMETHING
  # not_detected is independent of age/size
  # survival depends on year and age
  # set up recursively (but still needs one entry for each individual)
  
  distribution(binary_capture_history) <- binomial(size = 1, p = bias)
  survival_hist_ones <- ones(length(p_survival_hist_stage))
  distribution(survival_hist_ones) <- binomial(size = 1, p = p_survival_hist_stage)
  # final_obs_ones <- ones(length(p_final_obs))
  # distribution(final_obs_ones) <- binomial(size = 1, p = p_final_obs)
  
}

# internal function: define stage-recapture log-likelihood
binned_stage_recapture_loglik <- function(data, params) {
  
  # unpack params
  survival <- params$survival
  age_to_stage <- params$age_to_stage_conversion
  classes <- params$classes
  bias <- params$bias
  
  # calculate first and final observations
  first_obs <- apply(data$data, 1, function(x) min(which(x > 0)))
  final_obs <- apply(data$data, 1, function(x) max(which(x > 0)))
  
  # pull out first and final ages
  first_class <- apply(data$data, 1, function(x) x[min(which(x > 0))])
  final_class <- apply(data$data, 1, function(x) x[max(which(x > 0))])
  
  # number of years alive
  n_alive <- final_obs - first_obs
  
  # are any individuals never recaptured?
  single_obs <- first_obs == final_obs
  
  # focus on those with >1 observation
  n_alive <- n_alive[!single_obs]
  first_class <- first_class[!single_obs]
  first_seen <- first_obs[!single_obs]
  last_seen <- final_obs[!single_obs]
  
  # probabiity of stage j given age k (should be n_obs x n_stage)
  p_stage_first_capture <- age_to_stage[first_class, ]
  
  # create Phi = P(surv_n_alive | stage_at_time_1) (n_alive x n_stage)
  max_alive <- max(n_alive)
  survival <- cbind(survival,
                    do.call(cbind, lapply(seq_len(max_alive), function(x) survival[, n_age])))
  
  # pre-calculate all possible survival trajectories
  first_seen_lifespan <- matrix(as.numeric(unlist(strsplit(unique(paste(first_seen, n_alive, sep = '_')), '_'))), ncol = 2, byrow = TRUE)
  id_match <- match(paste(first_seen, n_alive, sep = '_'), unique(paste(first_seen, n_alive, sep = '_')))
  idx <- apply(first_seen_lifespan, 1, function(x) x[1]:(x[1] + x[2]))
  surv_mat <- ones(nrow(first_seen_lifespan), classes)
  for (i in seq_len(nrow(first_seen_lifespan))) {
    mat_tmp <- survival[idx[[i]], ]
    out <- NULL
    for (j in seq_len(classes)) {
      out <- c(out, seq(1, length(idx[[i]]) ^ 2, by = (length(idx[[i]]) + 1)) + (j - 1) * length(idx[[i]]))
    }
    idy <- rep(seq_len(classes), each = length(idx[[i]]))
    surv_mat[i, ] <- tapply(mat_tmp[out], idy, "prod")
  }
  
  # what are survival probabilities of each stage at first sight?
  p_survival_hist_stage <- surv_mat[id_match, ]
  
  # calculate P(survival_history | age = k)
  p_survival_hist_stage <- rowSums(p_survival_hist_stage * p_stage_first_capture)
  
  # probs of binary capture history, assume detection is constant
  binary_capture_history <- apply(data$data, 1, function(x) x[min(which(x > 0)):max(which(x > 0))])
  binary_capture_history <- ifelse(do.call(c, binary_capture_history) > 0, 1, 0)
  
  # probs of final detection (for each age??)
  # p(never_observed_again | final_size, final_obs) = \Sum_ages p(final_age | final_size) p(not_observed | final_obs, final_age) 
  # p_final_obs <- SOMETHING
  # not_detected is independent of age/size
  # survival depends on year and age
  # set up recursively (but still needs one entry for each individual)
  
  distribution(binary_capture_history) <- binomial(size = 1, p = bias)
  survival_hist_ones <- ones(length(p_survival_hist_stage))
  distribution(survival_hist_ones) <- binomial(size = 1, p = p_survival_hist_stage)
  # final_obs_ones <- ones(length(p_final_obs))
  # distribution(final_obs_ones) <- binomial(size = 1, p = p_final_obs)
  
}

# internal function: define stage-to-age log-likelihood
stage_to_age_loglik <- function(data, params) {

  # unpack parameters
  stage_to_age <- params$stage_to_age_conversion
  
  # estimate size-age distribution from data
  size_sums <- apply(data$data, 1, sum)

  # size-age model
  distribution(data$data) <- do.call(data$likelihood, list(size = size_sums, p = stage_to_age))

}

# internal function: define age-to-stage log-likelihood
age_to_stage_loglik <- function(data, params) {
  
  # unpack parameters
  age_to_stage <- params$age_to_stage_conversion
  
  # estimate size-age distribution from data
  size_sums <- apply(data$data, 1, sum)
  
  # size-age model
  distribution(data$data) <- do.call(data$likelihood, list(size = size_sums, p = age_to_stage))
  
}
