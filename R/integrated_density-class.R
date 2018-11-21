#' @name integrated_density
#' @title create integrated density objects
#'
#' @description An \code{integrated_density} object specifies the form of density
#'   dependence and associated prior distributions for an integrated population analysis
#' 
#' @param lambda parameter of Ricker and Beverton-Holt models of density dependence
#' @param ... additional arguments to \link[base]{print}, \link[base]{summary}, and \link[graphics]{plot} methods (currently ignored)
#' @param x an \code{integrated_density} object
#' @param object an \code{integrated_density} object
#'
#' @details Prior distributions can be specified as single-dimensional
#'   greta distribution, e.g., \code{normal(0, 1)}. Link functions and transformations
#'   can be specified directly in-line, e.g., \code{ilogit(normal(0, 1))} specifies
#'   normal priors with a mean of zero and a standard deviation of one, transformed
#'   with an inverse-logit link.
#'
#' @return An object of class \code{integrated_density}, which can be used to create
#'    and \link[greta.integrated]{integrated_process} object
#' 
#' @import greta
#' @import tensorflow
#' 
#' @examples
#' \dontrun{
#' 
#' library(integrated)
#' 
#' # a really basic age-structured model with five age classes
#' process <- leslie(5, density = ricker(lambda = uniform(0, 1)))
#' 
#' # setting custom priors
#' process <- leslie(5, density = beverton(lambda = uniform(0, 1)),
#'                   priors = list(survival = ilogit(normal(0, 1)),
#'                                 fecundity = exp(normal(0, 1))))
#' }

#' @export
#' @rdname integrated_density
#' 
ricker <- function(parameters = uniform(0, 1)) {

  # specify functional form
  form <- function(x, params) {
    exp(- x * params)
  }
  
  # collate and return outputs  
  density <- list(form = form,
                  parameters = parameters,
                  name = "ricker")
  
  # return outputs
  as.integrated_density(density)
  
}

#' @export
#' @rdname integrated_density
#' 
beverton <- function(parameters = uniform(0, 1)) {
  
  # specify functional form
  form <- function(x, params) {
    tf$ones(1, dtype = x$dtype) / (tf$ones(1, dtype = x$dtype) + x * params)
  }
  
  # collate and return outputs  
  density <- list(form = form,
                  parameters = parameters,
                  name = "beverton")
  
  # return outputs
  as.integrated_density(density)
  
}

#' @export
#' @rdname integrated_density
#' 
element_wise <- function(parameters = uniform(0, 1), masks = list()) {
  
  # setup masking matrices
  mask_list <- list()
  mask_list[names(masks)] <- masks
  
  # specify functional form
  form <- function(x, params) {
    tf$ones(1, dtype = x$dtype) / (tf$ones(1, dtype = x$dtype) + x * params)
  }
  
  # collate and return outputs  
  density <- list(form = form,
                  parameters = parameters,
                  name = "custom")
  
  # return outputs
  as.integrated_density(density)
  
}

#' @export
#' @rdname integrated_density
#' 
no_density <- function(parameters = uniform(0, 1)) {
  
  # specify functional form
  form <- function(x, params) {
    1.0
  }
  
  # collate and return outputs  
  density <- list(form = form,
                  parameters = parameters,
                  name = "none")
  
  # return outputs
  as.integrated_density(density)
  
}

#' @export
#' @rdname integrated_density
#' 
is.integrated_density <- function(object) {
  inherits(object, "integrated_density")
}

#' @export
#' @rdname integrated_density
#' 
print.integrated_density <- function(x, ...) {
  cat(paste0("This is an integrated_density object\n"))
}

#' @export
#' @rdname integrated_density
#' 
summary.integrated_density <- function(object, ...) {
  
  NULL
  
}

#' @export
#' @rdname integrated_density
#' 
plot.integrated_density <- function(x, ...) {
  
  # make a nice plot of scaling coefficient against range of values
  NULL
  
}

# internal function: create integrated_density object
as.integrated_density <- function(object) {
  as_class(object, name = "integrated_density", type = "list")
}
