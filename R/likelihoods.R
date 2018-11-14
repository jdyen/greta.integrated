# internal function: define abundance log-likelihood
abundance_loglik <- function(data, params) {

  age_dist <- vector('list', length = n_site)
  for (i in seq_len(n_site)) {
    age_dist[[i]] <- iterate_matrix_dynamic(matrix = mat$matrix[system_list == all_systems[i], , ],
                                            initial_state = c(mat$init[, i]),
                                            dens_param = mat$dens_param[i],
                                            dens_form = dens_type)
  }
  
  ## NEED SOME WAY TO HANDLE THIS -- REQUIRES > 1 LIKELIHOOD
  # BUT PARAMS ARE PASSED, SO JUST NEED TO MAKE SURE age_size_dist EXISTS,
  ## AND MAKE SURE THIS FUNCTION KNOWS IT EXISTS
  ### BUT ABUNDANCE_DATA CLASS IS MATCHING PROCESS TO DATA-- WILL ERROR BEFORE THIS
  #  NEED A NEW CLASS FOR SIZE_AGE_DATA? 
  # MAYBE: "stage_abundance()" and "age_abundance()" options?
  # IF stage_abundance() but leslie matrix, warn but allow it
  # convert sizes to ages
  modelled_sizes <- vector('list', length = length(age_dist))
  for (i in seq_along(age_dist))
    modelled_sizes[[i]] <- age_size_dist %*% age_dist[[i]]
  
  # create vectors of fitted and observed data
  mu <- do.call(c, modelled_sizes)
  size_vec <- do.call(c, size_binned)
  
  # size obs model
  # ALLOW USERS TO PASS THEIR OWN LIKELIHOOD (WITHIN REASON?)
  distribution(size_vec) <- poisson(mu)
  
  
  
}

# internal function: define recapture log-likelihood
recapture_loglik <- function(data, params) {
  
  # calculate first and final observations
  first_obs <- apply(catch_size_class, 1, function(x) min(which(x > 0)))
  final_obs <- apply(catch_size_class, 1, function(x) max(which(x > 0)))
  
  # number of years alive
  n_alive <- final_obs - first_obs
  
  # are any individuals never recaptured?
  single_obs <- first_obs == final_obs
  
  # focus on those with >1 observation
  n_alive <- n_alive[!single_obs]
  first_sizes <- first_size_class[!single_obs]
  first_seen <- first_obs[!single_obs]
  last_seen <- final_obs[!single_obs]
  
  # probabiity of age j given size class k (should be n_obs x n_age)
  p_age_first_capture <- age_size_dist[first_sizes, ]
  
  # create Phi = P(surv_n_alive | age_at_time_1) (n_alive x n_age)
  max_alive <- max(n_alive)
  survival_lmr <- mat$surv_params[system_list == 'LOWERMURRAY', ]
  survival_lmr <- cbind(survival_lmr,
                        do.call(cbind,
                                ## INDEXING ERRORS MIGHT BE DUE TO THIS
                                lapply(seq_len(max_alive),
                                       function(x) survival_lmr[, n_age])))
  
  # pre-calculate all possible survival trajectories
  first_seen_lifespan <- matrix(as.numeric(unlist(strsplit(unique(paste(first_seen, n_alive, sep = '_')), '_'))), ncol = 2, byrow = TRUE)
  id_match <- match(paste(first_seen, n_alive, sep = '_'), unique(paste(first_seen, n_alive, sep = '_')))
  idx <- apply(first_seen_lifespan, 1, function(x) x[1]:(x[1] + x[2]))
  surv_mat <- ones(nrow(first_seen_lifespan), n_age)
  for (i in seq_len(nrow(first_seen_lifespan))) {
    mat_tmp <- survival_lmr[idx[[i]], ]
    out <- NULL
    for (j in seq_len(n_age)) {
      out <- c(out, seq(1, length(idx[[i]]) ^ 2, by = (length(idx[[i]]) + 1)) + (j - 1) * length(idx[[i]]))
    }
    idy <- rep(seq_len(n_age), each = length(idx[[i]]))
    surv_mat[i, ] <- tapply(mat_tmp[out], idy, 'prod')
  }
  
  # what are survival probabilities for each age at first sight?
  p_survival_hist_age <- surv_mat[id_match, ]
  
  # calculate P(survival_history | size = k)
  p_survival_hist_size <- rowSums(p_survival_hist_age * p_age_first_capture)
  
  # probs of binary capture history, assume detection is constant
  detection <- beta(1, 1)
  binary_capture_history <- apply(catch_size_class, 1, function(x) x[min(which(x > 0)):max(which(x > 0))])
  binary_capture_history <- ifelse(do.call(c, binary_capture_history) > 0, 1, 0)
  
  # probs of final detection (for each age??)
  # p(never_observed_again | final_size, final_obs) = \Sum_ages p(final_age | final_size) p(not_observed | final_obs, final_age) 
  # p_final_obs <- SOMETHING
  # not_detected is independent of age/size
  # survival depends on year and age
  # set up recursively (but still needs one entry for each individual)
  
  distribution(binary_capture_history) <- binomial(size = 1, p = detection)
  survival_hist_ones <- ones(length(p_survival_hist_size))
  distribution(survival_hist_ones) <- binomial(size = 1, p = p_survival_hist_size)
  # final_obs_ones <- ones(length(p_final_obs))
  # distribution(final_obs_ones) <- binomial(size = 1, p = p_final_obs)
  
}

# internal function: define size-at-age log-likelihood
size_age_loglik <- function(data, params) {

  alpha_set <- ifelse(is.na(size_age_obs), 0, 1)
  age_size_dist <- dirichlet(alpha = alpha_set)
  
  # estimate size-age distribution from data
  size_sums <- apply(size_age_obs, 1, sum)

  
  # size-age model
  ## ALLOW ALTERNATIVE LIKLEIHOODS?
  distribution(size_age_obs) <- multinomial(size = size_sums, p = age_size_dist)
  
}
