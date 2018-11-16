# internal function: create a set of parameters based on an integrated process object
define_parameters <- function(classes, priors, masks, process_type, predictors = NULL, classes_alt = NULL) {
  
  # start with the predictors, see if we need to deal with them
  if (!is.null(predictors)) {
    
    # what do the predictors look like?
    n_obs <- nrow(predictors)
    n_pred <- ncol(predictors)

    # what dims do we need to get to?
    survival_dims <- c(n_pred, sum(masks$survival))
    transition_dims <- c(n_pred, sum(masks$transition))
    fecundity_dims <- c(n_pred, sum(masks$fecundity))

    # setup priors with correct dims
    survival <- change_dims_extract_op(priors$survival, survival_dims)
    transition <- change_dims_extract_op(priors$transition, transition_dims)
    fecundity <- change_dims_extract_op(priors$fecundity, fecundity_dims)
    
    # define additional noise in transitions
    survival_noise <- change_dims_extract_op(priors$noise, c(n_obs, survival_dims[2]))
    transition_noise <- change_dims_extract_op(priors$noise, c(n_obs, transition_dims[2]))
    fecundity_noise <- change_dims_extract_op(priors$noise, c(n_obs, fecundity_dims[2]))
    
    # check that ops match for noise and main components
    if (survival$op != survival_noise$op)
      warning("transformation applied to survival does not match transformation of noise; using survival transformation", call. = FALSE)
    if (transition$op != transition_noise$op)
      warning("transformation applied to transition does not match transformation of noise; using survival transformation", call. = FALSE)
    if (fecundity$op != fecundity_noise$op)
      warning("transformation applied to fecundity does not match transformation of noise; using survival transformation", call. = FALSE)
    
    # setup initial conditions
    inits <- change_dims(priors$initials, c(classes, 1L))
    
    # setup linear predictors
    survival_link <- predictors %*% survival + survival_noise
    transition_link <- predictors %*% transition + transition_noise
    fecundity_link <- predictors %*% fecundity + fecundity_noise

    # recombine link nodes with correct transformations
    if (length(survival$op_args)) {
      survival_all <- do.call(survival$op, list(survival_link), survival$op_args)
    } else {
      survival_all <- do.call(survival$op, list(survival_link))
    }
    if (length(transition$op_args)) {
      transition_all <- do.call(transition$op, list(transition_link), transition$op_args)
    } else {
      transition_all <- do.call(transition$op, list(transition_link))
    }
    if (length(fecundity$op_args)) {
      fecundity_all <- do.call(fecundity$op, list(fecundity_link), fecundity$op_args)
    } else {
      fecundity_all <- do.call(fecundity$op, list(fecundity_link))
    }

    # construct population matrix
    mat1 <- mat2 <- mat3 <- zeros(classes, classes)
    
    # vectorise matrix fill: survival
    idx <- which(masks$survival == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat1[idy] <- survival_all
    
    # vectorise matrix fill: transition
    idx <- which(masks$transition == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat2[idy] <- transition_all

    # need to standardise survival and transition matrices
    mat4 <- mat1 + mat2
    sweep_sums <- apply(mat4, c(1, 3), "sum")
    mat4 <- sweep(mat4, c(1, 3), sweep_sums, "/")
    sweep_survival <- beta(1, 1, c(n_obs, classes))
    mat4 <- sweep(mat4, c(1, 3), sweep_survival, "*")
    
    # vectorise matrix fill: fecundity
    idx <- which(masks$fecundity == 1) - 1
    idy <- n_obs * rep(idx, each = n_obs) + seq_len(n_obs)
    mat3[idy] <- fecundity_all

    # combine survival, transition, and fecundity    
    mat <- mat3 + mat4

  } else {
    
    # we haven't got predictors, so this isn't hard
    n_obs <- 1

    # what dims do we need to get to?
    survival_dims <- sum(masks$survival)
    transition_dims <- sum(masks$transition)
    fecundity_dims <- sum(masks$fecundity)
    
    # setup priors with correct dims
    survival <- change_dims(priors$survival, survival_dims)
    transition <- change_dims(priors$transition, transition_dims)
    fecundity <- change_dims(priors$fecundity, fecundity_dims)
    
    # setup initial conditions
    inits <- change_dims(priors$initials, c(classes, 1L))
   
    # construct population matrix
    mat1 <- mat2 <- mat3 <- zeros(classes, classes)
    
    # vectorise matrix fill: survival
    idx <- which(masks$survival == 1)
    mat1[idx] <- survival
    
    # vectorise matrix fill: transition
    idx <- which(masks$transition == 1)
    mat2[idx] <- transition
    
    # need to standardise survival and transition matrices
    mat4 <- mat1 + mat2
    sweep_sums <- apply(mat4, 2, "sum")
    mat4 <- sweep(mat4, 2, sweep_sums, "/")
    sweep_survival <- beta(1, 1, classes)
    mat4 <- sweep(mat4, 2, sweep_survival, "*")

    # vectorise matrix fill: fecundity
    idx <- which(masks$fecundity == 1)
    mat3[idx] <- fecundity

    # combine survival, transition, and fecundity    
    mat <- mat3 + mat4
    
  }
  
  age_stage <- stage_age <- NULL
  if (!is.null(classes_alt)) {
    if (process_type == "age")
      stage_age <- dirichlet(alpha = ones(classes_alt, classes))
    if (process_type == "stage")
      age_stage <- dirichlet(alpha = ones(classes_alt, classes))
  }
  
  # return outputs
  list(n_obs = n_obs, classes = classes,
       matrix = mat, inits = inits,
       survival = survival, transition = transition, fecundity = fecundity,
       sweep_survival = sweep_survival,
       age_to_stage_conversion = age_stage, stage_to_age_conversion = stage_age)
  
}
