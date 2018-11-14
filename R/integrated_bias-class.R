#' @name integrated_bias
#' @title create integrated bias objects
#'
#' @description An \code{integrated_bias} object contains a bias model
#'   that can be shared among multiple components in an integrated population
#'   model
#' 
#' @param params named list of parameters (see details for information on setting prior distributions)
#' @param ... additional arguments to \link[base]{print}, \link[base]{summary}, and \link[graphics]{plot} methods (currently ignored)
#' @param x an \code{integrated_bias} object
#' @param object an \code{integrated_bias} object
#'
#' @details Prior distributions can be specified as single-dimensional
#'   greta distribution, e.g., \code{normal(0, 1)}. Link functions and transformations
#'   can be specified directly in-line, e.g., \code{ilogit(normal(0, 1))} specifies
#'   normal priors with a mean of zero and a standard deviation of one, transformed
#'   with an inverse-logit link.
#'
#' @return An object of class \code{integrated_bias}, which can be used to create
#'    \link[greta.integrated]{integrated_data} and \link[greta.integrated]{integrated_model} objects
#' 
#' @import greta
#' 
#' @examples
#' \dontrun {
#' 
#' library(integrated)
#' 
#' # a really basic age-structured model with five age classes
#' process <- leslie(5, density = "none")
#' 
#' # setting custom priors
#' process <- detection(params = list(abundance = ilogit(normal(0, 1))))
#' }

#' @export
#' @rdname integrated_bias
#' 
detection <- function(params = list()) {
  
  # set type
  type <- "detection"
  
  # initialise model parameters
  param_list <- list(abundance = beta(1, 1),
                     recapture = beta(1, 1),
                     linear = normal(0, 1))
  param_list[names(params)] <- params
  
  # do the parameters have reasonable bounds?
  abundance_bounds <- extract_bounds(param_list$abundance)
  recapture_bounds <- extract_bounds(param_list$recapture)
  
  # warn if not
  if (abundance_bounds[1] < 0 | abundance_bounds[2] > 1)
    warning("the prior for abundance has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (recapture_bounds[1] < 0 | recapture_bounds[2] > 1)
    warning("the prior for recapture has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  
  # collate and return outputs  
  bias <- list(type = type,
               params = param_list)
  
  # return outputs
  integrated_bias(bias)
  
}

#' @export
#' @rdname integrated_bias
#' 
is.integrated_bias <- function(object) {
  inherits(object, "integrated_bias")
}

#' @export
#' @rdname integrated_bias
#' 
print.integrated_bias <- function(x, ...) {
  cat(paste0("This is an integrated_bias object\n"))
}

#' @export
#' @rdname integrated_bias
#' 
summary.integrated_bias <- function(object, ...) {
  
  NULL
  
}

# internal function: create integrated_bias object
as.integrated_bias <- function(object) {
  as_class(object, name = "integrated_bias", type = "list")
}
