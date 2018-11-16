#' @name integrated_summary
#' @title summarise greta models containing \link[greta.integrated]{integrated_model} components
#'
#' @description A set of helper functions.
#' 
#' @param object a \link[coda]{mcmc.list} object created with \link[greta]{mcmc}
#'
#' @details something
#'
#' @return A \link[greta.integrated]{integrated_summary} object
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
#' @rdname integrated_summary
#' 
integrated_summary <- function(object) {
  NULL
}
