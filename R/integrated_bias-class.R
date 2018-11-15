#' @name integrated_bias
#' @title create integrated bias objects
#'
#' @description An \code{integrated_bias} object contains a bias model
#'   that can be shared among multiple components in an integrated population
#'   model
#' 
#' @param x 
#' @param p_detect parameter of detection model (see details for information on setting prior distributions)
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
#' \dontrun{
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
no_bias <- function() {
  
  # set up no_bias function
  bias_fn <- function(x, params) x
  
  bias <- list(bias = bias_fn,
               params = NULL)
  
  # return outputs
  as.integrated_bias(bias)
    
}

#' @export
#' @rdname integrated_bias
#' 
detection <- function(p_detect = beta(1, 1)) {
  
  # is the parameter of a reasonable class?
  if ("greta_array" %in% class(p_detect)) {
    
    # check if it's a distribution
    node <- get_node(p_detect)
    
    # is it a distribution?
    if (!is.null(node$distribution)) {
      
      # if so, do the parameters have reasonable bounds?
      p_bounds <- extract_bounds(p_detect)
      
    } else {
     
      # if not, does it have reasonable dims?
      if (any(dim(p_detect) != 1))
        stop("p_detect must be a scalar numeric, scalar data greta_array, or greta distribution", call. = FALSE)
      
    }
    
  } else {
    
    if (is.numeric(p_detect)) {
      
      # check length
      if (length(p_detect) > 1)
        stop("p_detect must be a scalar numeric, scalar data greta_array, or greta distribution", call. = FALSE)
      
      # do the parameters have reasonable bounds?
      p_bounds <- range(p_detect)
      
    }
    
  }
    
  # warn if bounds are not logical
  if (p_bounds[1] < 0 | p_bounds[2] > 1)
    warning("the prior for p_detect has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)

  # set up correction function
  bias_fn <- function(x, params) params$detection * x

  bias <- list(bias = bias_fn,
               params = list(detection = p_detect))
  
  # return outputs
  as.integrated_bias(bias)
  
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
