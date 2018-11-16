# internal function: create a set of parameters based on an integrated process object
define_parameters <- function(process, classes_alt, n_predictors = NULL) {
  
  # do we need to deal with predictors?
  if (!is.null(n_predictors)) {

    # what dims do we need to get to?
    survival_dims <- c(n_predictors, process$classes)
    transition_dims <- c(n_predictors, sum(process$masks$transition))
    fecundity_dims <- c(n_predictors, sum(process$masks$fecundity))
    
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
  list(classes = process$classes, classes_alt = classes_alt,
       inits = inits,
       n_predictors = n_predictors,
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
  if (!is.null(parameters$n_predictors)) {

    # how many rows does this matrix need? (one for each set of predictors)    
    n_obs <- nrow(data$predictors)

    # extract noise prior to save some typing
    noise <- data$process$priors$noise
    
    # define additional noise in transitions
    surv_noise <- change_dims_extract_op(noise, c(n_obs, parameters$classes))
    tran_noise <- change_dims_extract_op(noise, c(n_obs, sum(masks$transition)))
    fec_noise <- change_dims_extract_op(noise, c(n_obs, sum(masks$fecundity)))
    
    # we need to combine these two parameters with the same transformation
    if (surv$op != surv_noise$op)
      warning("transformation applied to survival does not match transformation of noise; using survival transformation", call. = FALSE)
    if (tran$op != tran_noise$op)
      warning("transformation applied to transition does not match transformation of noise; using survival transformation", call. = FALSE)
    if (fec$op != fec_noise$op)
      warning("transformation applied to fecundity does not match transformation of noise; using survival transformation", call. = FALSE)
    
    # we can't make a linear predictor if the raw parameters come from multiple distributions
    if (length(surv$children) > 1)
      stop("transformed prior for survival has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(tran$children) > 1)
      stop("transformed prior for transition has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(fec$children) > 1)
      stop("transformed prior for fecundity has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    if (length(surv_noise$children) > 1)
      stop("transformed prior for noise has multiple inputs, which will not work in a model that includes predictors", call. = FALSE)
    
    # setup linear predictors
    surv_link <- predictors$data %*% surv$children[[1]] + surv_noise$children[[1]]
    tran_link <- predictors$data %*% tran$children[[1]] + tran_noise$children[[1]]
    fec_link <- predictors$data %*% fec$children[[1]] + fec_noise$children[[1]]
    
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
    
    # construct population matrix
    mat1 <- mat2 <- zeros(n_obs, parameters$classes, parameters$classes)
    
    # vectorise matrix fill: transition
    idx <- which(masks$transition == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat1[idy] <- transition_all
    
    # need to standardise survival and transition matrices
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