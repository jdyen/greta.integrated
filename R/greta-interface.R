#' @name greta_interface
#' @title work with greta objects created with \link[greta.integrated]{integrated_model}
#'
#' @description A set of helper functions.
#' 
#' @param object a \link[greta]{greta_model} object
#'
#' @details something
#'
#' @return A \link[greta]{initials} object
#' 
#' @export
#' 
#' @import greta
#' 
#' @examples
#' \dontrun{
#' 
#'  library(integrated)
#' 
#' # prepare an example model
#' model <- integrated_model()
#'                         
#' # summarise fitted model
#' model
#' summary(model)
#' plot(model)
#' }

#' @export initialise
#' @rdname greta_interface
#' 
initialise <- function(object) {
  UseMethod("initialise")
}

#' @export initialize
#' @rdname greta_interface
#' 
initialize <- function(object) {
  UseMethod("initialise")
}

#' @method initialise greta_model
#' @export
# internal function: prepare greta initials from integrated_model
initialise.greta_model <- function(object) {

  # extract the underlying distribution nodes for each parameter  
  dependencies <- greta.integrated:::get_node(object$target_greta_arrays[[1]])$child_names(recursive = TRUE)
  all_nodes <- object$dag$node_list[dependencies]
  is_distrib <- sapply(all_nodes, function(x) !is.null(x$distribution))
  distrib_nodes <- all_nodes[is_distrib]

  # apply initialise function separately to each node
  inits_list <- lapply(distrib_nodes, initialise_internal)

  # turn params into a flat vector but store these in a named list
  inits_list <- list(unlist(inits_list))
  names(inits_list) <- names(object$target_greta_arrays)
  
  # return output as a greta initials object
  do.call(initials, inits_list)

}

#' @method initialise default
#' @export
# internal function: prepare greta initials from integrated_model
initialise.default <- function(object) {
  
  # can only work with greta models
  stop(paste0("cannot set initial values from an object of class ", class(object), "; have you compiled the integrated_model to a greta model?"), call. = FALSE)
  
}

# internal function: set initials based on bounds of a single greta array
initialise_internal <- function(x) {
  
  # extract bounds
  bounds <- extract_bounds(x)
  
  # work out if the bounds are finite
  finite_bounds <- bounds[is.finite(bounds)]
  
  # is either bound positive or negative inf?
  is_inf <- is.infinite(bounds)
  is_pos_inf <- is_inf & bounds > 0
  is_neg_inf <- is_inf & bounds < 0
  
  # start at zero, keep it this way if both bounds are inf
  out <- array(0, dim = x$dim)
  
  # if both bounds are finite, use the mean to set a starting condition
  if (length(finite_bounds) == 2)
    out <- array(mean(finite_bounds), dim = x$dim)
  
  # if one of the bounds is real and the other is positive infinite, set all inits to 1
  if (length(finite_bounds) & any(is_pos_inf))
    out <- array(1.0, dim = x$dim)
  
  # if one of the bounds is real and the other is negative infinite, set all inits to -1
  if (length(finite_bounds) & any(is_neg_inf))
    out <- array(-1.0, dim = x$dim)
  
  # return outputs
  out
  
}

#' @export parameters
#' @rdname greta_interface
#' 
parameters <- function(object, ...) {
  UseMethod("parameters")
}

#' @method parameters integrated_model
#' @export
# internal function: extract greta parameters from integrated_model
parameters.integrated_model <- function(object, ...) {

  # are there extra parameters?
  extras <- list(...)
  extras_names <- as.list(substitute(list(...)))[-1L]
  
  # pull out params from each process
  main_param <- lapply(object, extract_direct)
  dens_param <- lapply(object, extract_dens)
  bias_param <- lapply(object, extract_bias)
  sigma_param <- lapply(object, extract_sigmas)

  # flatten these to a single vector  
  param_vec1 <- do.call(c, lapply(main_param, function(x) x$parameters))
  param_vec2 <- do.call(c, lapply(dens_param, function(x) x$parameters))
  param_vec3 <- do.call(c, lapply(bias_param, function(x) x$parameters))
  param_vec4 <- do.call(c, lapply(sigma_param, function(x) x$parameters))
  param_vec <- c(param_vec1, param_vec2, param_vec3, param_vec4)
  
  # make sure we keep names for each parameter to help reconstruct
  names_vec1 <- do.call(c, lapply(main_param, function(x) x$names))
  names_vec2 <- do.call(c, lapply(dens_param, function(x) x$names))
  names_vec3 <- do.call(c, lapply(bias_param, function(x) x$names))
  names_vec4 <- do.call(c, lapply(sigma_param, function(x) x$names))
  names_vec <- c(names_vec1, names_vec2, names_vec3, names_vec4)

  # pull out dims of each parameter, also to help reconstruct  
  dim_mat1 <- do.call(rbind, lapply(main_param, function(x) do.call(rbind, x$dims)))
  dim_mat2 <- do.call(rbind, lapply(dens_param, function(x) do.call(rbind, x$dims)))
  dim_mat3 <- do.call(rbind, lapply(bias_param, function(x) do.call(rbind, x$dims)))
  dim_mat4 <- do.call(rbind, lapply(sigma_param, function(x) do.call(rbind, x$dims)))
  dim_mat <- rbind(dim_mat1, dim_mat2, dim_mat3, dim_mat4)
  
  # are there extra parameters to pass to greta::model?
  if (length(extras)) {
    
    # flatten extras to a list and add to parameters
    param_vec <- c(param_vec, do.call(c, extras))
    
    # what dimensions do these have?
    extras_dims <- t(sapply(extras, dim))
    rownames(extras_dims) <- extras_names
    dim_mat <- rbind(dim_mat, extras_dims)
    
    # add names for extra parameters
    names_vec <- c(names_vec, rep(as.character(extras_names), times = apply(extras_dims, 1, "prod")))

  }

  # return all
  list(parameters = param_vec, names = names_vec, dims = dim_mat)
  
}

#' @method parameters default
#' @export
# internal function: extract greta parameters from integrated_model
parameters.default <- function(object, ...) {
  
  # can only work with greta models
  stop(paste0("cannot extract parameters from an object of class ", class(object), collapse = "; "), call. = FALSE)
  
  
}

# internal function: extract visible parameters from integrated model
extract_direct <- function(x) {

  # which vars to we want to track?  
  include <- c("inits", "survival", "transition", "fecundity",
               "age_to_stage_conversion", "stage_to_age_conversion")
  
  # we want these in a flat vector
  parameters <- do.call(c, x[include])
  
  # but need to keep track of names and dims so we can reconstruct
  xnames <- rep(names(x[include]), times = sapply(x[include], length))
  xdims <- sapply(x[include], dim)
  
  # return
  list(parameters = parameters, names = xnames, dims = xdims)
  
}

# internal function: extract bias parameters from integrated model
extract_bias <- function(x) {
  
  x_tmp <- x$bias
  include <- !sapply(x_tmp, is.null)[seq(2, length(x_tmp), by = 2)]
  
  # we want these in a flat vector
  parameters <- do.call(c, x_tmp[2, include])
  
  # this is still a list, so re-apply the c() function
  parameters <- do.call(c, parameters)
  
  # but need to keep track of names and dims so we can reconstruct
  xnames <- rep(colnames(x_tmp)[include], times = sapply(x_tmp[2, include], length))
  xdims <- lapply(x_tmp[2, include], function(x) c(length(x), 1))
  xnames <- paste0("bias_", xnames)
  names(xdims) <- paste0("bias_", xnames)
  
  # return
  list(parameters = parameters, names = xnames, dims = xdims)
  
}

# internal function: extract bias parameters from integrated model
extract_dens <- function(x) {
  
  x_tmp <- x$dens

  # we want these in a flat vector
  parameters <- c(x_tmp$parameters)
  
  # but need to keep track of names and dims so we can reconstruct
  xnames <- c(x_tmp$name)
  xdims <- lapply(parameters, function(x) c(length(x), 1))
  xnames <- paste0("dens_", xnames)
  names(xdims) <- paste0("dens_", xnames)
  
  # return
  list(parameters = parameters, names = xnames, dims = xdims)
  
}

# internal function: extract bias parameters from integrated model
extract_sigmas <- function(x) {
  
  if (is.null(x$sigmas)) {
    
    parameters <- names <- dims <- NULL
    
  } else {
    
    x_tmp <- x$sigmas
    
    # we want these in a flat vector
    parameters <- do.call(c, x_tmp)
    
    # but need to keep track of names and dims so we can reconstruct
    xnames <- rep(names(x_tmp), times = sapply(x_tmp, length))
    xdims <- lapply(x_tmp, dim)
    xnames <- paste0("sigma_", xnames)
    names(xdims) <- paste0("sigma_", names(xdims))
    
  }
  
  # return
  list(parameters = parameters, names = xnames, dims = xdims)
  
}
