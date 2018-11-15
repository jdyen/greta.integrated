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

#' @export
#' @rdname greta_interface
#' 
initialise <- function(object) {
  UseMethod("initialise")
}

# internal function: prepare greta initials from integrated_model
initialise.greta_model <- function(object) {
  
  # extract parameters
  params <- object$target_greta_arrays
  
  inits_list <- lapply(params, initialise_internal)
  inits_list <- inits_list[!sapply(inits_list, function(x) is.null(x))]
  
  # return output as a greta initials object
  do.call(initials, inits_list)

}

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

    # if not, just set them to zero    
    out <- zeros(dim(x))
    
    # if they are, use the mean of two values or the single finite value
    if (length(finite_bounds))
      out <- greta_array(mean(finite_bounds), dim(x))
    
  }
  
  # return outputs
  out
  
}

integrated_summary <- NULL
