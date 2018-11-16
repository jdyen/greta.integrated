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

#' @method initialise greta_model
#' @export
# internal function: prepare greta initials from integrated_model
initialise.greta_model <- function(object) {
  
  # extract parameters
  params <- object$target_greta_arrays
  
  inits_list <- lapply(params, initialise_internal)
  inits_list <- inits_list[!sapply(inits_list, function(x) is.null(x))]
  
  # return output as a greta initials object
  do.call(initials, inits_list)

}

#' @method initialise default
#' @export
# internal function: prepare greta initials from integrated_model
initialise.default <- function(object) {
  
  # can only work with greta models
  stop(paste0("cannot set initial values from an object of class ", class(object)), call. = FALSE)
  
}

# internal function: set initials based on bounds of a single greta array
initialise_internal <- function(x) {
  
  # extract bounds
  bounds <- extract_bounds(x)
  
  # set a default
  out <- NULL
  
  # if we actually have bounds, we can update the default
  if (!is.null(bounds)) {
    
    # work out if the bounds are finite
    finite_bounds <- bounds[is.finite(bounds)]

    # is either bound positive or negative inf?
    is_inf <- is.infinite(bounds)
    is_pos_inf <- is_inf & bounds > 0
    is_neg_inf <- is_inf & bounds < 0
    
    # start at zero, keep it this way if both bounds are inf
    out <- array(0, dim = dim(x))
    
    # if both bounds are finite, use the mean to set a starting condition
    if (length(finite_bounds) == 2)
      out <- array(mean(finite_bounds), dim = dim(x))
    
    # if one of the bounds is real and the other is positive infinite, set all inits to 1
    if (length(finite_bounds) & any(is_pos_inf))
      out <- array(1.0, dim = dim(x))
    
    # if one of the bounds is real and the other is negative infinite, set all inits to -1
    if (length(finite_bounds) & any(is_neg_inf))
      out <- array(-1.0, dim = dim(x))

  }
  
  # return outputs
  out
  
}

integrated_summary <- NULL
