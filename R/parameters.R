# internal function: create a set of parameters based on an integrated process object
define_parameters <- function(process, classes_alt, n_fixed = NULL) {
  
  # do we need to deal with predictors?
  if (!is.null(n_fixed)) {

    # what dims do we need to get to?
    survival_dims <- c(n_fixed, process$classes)
    transition_dims <- c(n_fixed, sum(process$masks$transition))
    fecundity_dims <- c(n_fixed, sum(process$masks$fecundity))
    
    # setup priors with correct dims
    survival <- change_dims_extract_op(process$priors$survival, survival_dims)
    transition <- change_dims_extract_op(process$priors$transition, transition_dims)
    fecundity <- change_dims_extract_op(process$priors$fecundity, fecundity_dims)

    # setup initial conditions
    inits <- change_dims(process$priors$initials, c(process$classes, 1L))

  } else {
    
    # what dims do we need to get to?
    survival_dims <- process$classes
    transition_dims <- sum(process$masks$transition)
    fecundity_dims <- sum(process$masks$fecundity)
    
    # setup priors with correct dims
    survival <- change_dims(process$priors$survival, survival_dims)
    transition <- change_dims(process$priors$transition, transition_dims)
    fecundity <- change_dims(process$priors$fecundity, fecundity_dims)
    
    # setup initial conditions
    inits <- change_dims(process$priors$initials, c(process$classes, 1L))

  }
  
  age_stage <- stage_age <- NULL
  if (!is.null(classes_alt)) {
    if (process$type == "age")
      stage_age <- dirichlet(alpha = ones(classes_alt, process$classes))
    if (process$type == "stage")
      age_stage <- dirichlet(alpha = ones(classes_alt, process$classes))
  }
  
  # return outputs
  list(inits = inits,
       survival = survival, transition = transition, fecundity = fecundity,
       age_to_stage_conversion = age_stage, stage_to_age_conversion = stage_age)
  
}

construct_matrix <- function(data, parameters) {

  # pull out the masks; we need these to work out how many non-zero entries there are
  masks <- data$process$masks
  
  # extract parameters to save some characters
  surv <- parameters$survival
  tran <- parameters$transition
  fec <- parameters$fecundity
  
  # can only have lists if the parameters have predictor variables
  if (!is.null(data$predictors)) {

    # how many rows does this matrix need? (one for each set of predictors)    
    n_obs <- nrow(data$predictors$fixed)
    
    # reformat random variables
    random <- data$predictors$random
    n_cluster <- ncol(random)
    n_group <- apply(random, 2, max)
    random_vec <- c(sweep(random, 2, c(0, n_cluster[-length(n_cluster)]), "+"))
    
    # extract noise prior to save some typing
    noise <- data$process$priors$noise
    
    # define additional noise in transitions
    surv_sigma <- change_dims(noise, c(n_cluster, parameters$classes))
    tran_sigma <- change_dims(noise, c(n_cluster, sum(masks$transition)))
    fec_sigma <- change_dims(noise, c(n_cluster, sum(masks$fecundity)))

    # we can't make a linear predictor if the raw parameters come from multiple distributions
    if (length(surv$children) > 1)
      stop("transformed prior for survival has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(tran$children) > 1)
      stop("transformed prior for transition has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(fec$children) > 1)
      stop("transformed prior for fecundity has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(surv_noise$children) > 1)
      stop("transformed prior for noise has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    
    # add some random effects
    surv_noise_baseline <- normal(0.0, rep(surv_sigma, times = n_group))
    surv_noise <- rowSums(greta_array(surv_noise_baseline[random_vec], c(n_obs, n_cluster)))
    tran_noise_baseline <- normal(0.0, rep(tran_sigma, times = n_group))
    tran_noise <- rowSums(greta_array(tran_noise_baseline[random_vec], c(n_obs, n_cluster)))
    fec_noise_baseline <- normal(0.0, rep(fec_sigma, times = n_group))
    fec_noise <- rowSums(greta_array(fec_noise_baseline[random_vec], c(n_obs, n_cluster)))
    
    # setup linear predictors
    surv_link <- data$predictors %*% surv$children[[1]] + surv_noise
    tran_link <- data$predictors %*% tran$children[[1]] + tran_noise
    fec_link <- data$predictors %*% fec$children[[1]] + fec_noise
    
    # recombine link nodes with correct transformations
    if (length(surv$op_args)) {
      survival_all <- do.call(surv$op, list(surv_link), surv$op_args)
    } else {
      survival_all <- do.call(surv$op, list(surv_link))
    }
    if (length(tran$op_args)) {
      transition_all <- do.call(tran$op, list(tra_link), tran$op_args)
    } else {
      transition_all <- do.call(tran$op, list(tran_link))
    }
    if (length(fec$op_args)) {
      fecundity_all <- do.call(fec$op, list(fec_link), fec$op_args)
    } else {
      fecundity_all <- do.call(fec$op, list(fec_link))
    }
    
    # check bounds of new nodes and warn if not reasonable
    survival_bounds <- extract_bounds(survival_all)
    transition_bounds <- extract_bounds(transition_all)
    fecundity_bounds <- extract_bounds(fecundity_all)
    
    # are these reasonable?
    if (survival_bounds[1] < 0 | survival_bounds[2] > 1)
      warning("the transformed prior (including predictor effects) for survival has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
    if (transition_bounds[1] < 0 | transition_bounds[2] > 1)
      warning("the transformed prior (including predictor effects) for transition has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
    if (fecundity_bounds[1] < 0)
      warning("the transformed prior (including predictor effects) for fecundity has a lower bound less than 0; is this reasonable?", call. = FALSE)
    
    # construct population matrix
    mat1 <- mat2 <- zeros(n_obs, parameters$classes, parameters$classes)
    
    # vectorise matrix fill: transition
    idx <- which(masks$transition == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat1[idy] <- transition_all
    
    # need to standardise survival and transition matrices
    warning("internal warning: einsum sweep approximation is not working correctly")
    mat1_sums <- apply(mat1, c(1, 3), "sum")
    mat1 <- rescale_array(mat1, mat1_sums, survival_all)
    
    # vectorise matrix fill: fecundity
    idx <- which(masks$fecundity == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat2[idy] <- fecundity_all

  } else { # must not have predictors
    
    # construct population matrix
    mat1 <- mat2 <- zeros(parameters$classes, parameters$classes)
    
    # vectorise matrix fill: transition
    mat1[masks$transition == 1] <- tran
    
    # need to standardise transition matrices
    mat1 <- sweep(mat1, 2, colSums(mat1), "/")
    mat1 <- sweep(mat1, 2, surv, "*")
    
    # vectorise matrix fill: fecundity
    mat2[masks$fecundity == 1] <- fec
    
  }

  # combine transition and fecundity    
  mat <- mat1 + mat2
  
  mat
   
}