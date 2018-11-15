#' @name integrated_model
#' @title create integrated model objects
#'
#' @description An \code{integrated_model} object combines multiple \link[greta.integrated]{integrated_data}
#'   objects with one or more \link[greta.integrated]{integrated_process} objects and returns all
#'    \link[greta]{greta_array} objects required to fit an integrated population model. Model fitting
#'    is handled separately with the \link[greta]{model} and \link[greta]{mcmc} functions.
#' 
#' @param process an \link[greta.integrated]{integrated_process} object
#' @param ... one or more \link[greta.integrated]{integrated_data} objects
#' @param x an \code{integrated_model} object
#' @param object an \code{integrated_model} object
#'
#' @details something
#'
#' @return A named list of \code{greta_array} objects, which can be passed to
#'    \link[greta]{model}. This list has associated `print`, `plot`, and `summary` methods.
#' 
#' @export
#' 
#' @import greta
#' 
#' @examples
#' \dontrun{
#' 
#' ### ADD EXAMPLES OF ALL FUNCTIONS HERE
#' 
#' library(integrated)
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
#' @rdname integrated_model
#' 
integrated_model <- function(...) {

  data_modules <- list(...)

  # any make-or-break errors? Kill these first
  # anything not an integrated_data object?

  # any data_modules not currently implemented ("community", what else?)
  
  
  # are any data modules predictors? Deal with this second (expand process)
  
  
  # check process models  
  
  # check process models are all the same
  ## NEED A STRING FOR THIS?
  ## USE identical()
  process_list <- sapply(data_modules, extract_process)
  if (length(unique(process_list) > 1))
    stop(paste0("data are connected to ", length(unique(process_list)), " different processes"), call. = FALSE)
  
  # check the bias layers -- should share params where needed
  
  # DEFINE PROCESS
  # expand to use multiple processes      
  process <- process_list[1]
  # have to add greta array setup here
  ###
  
  # initialise mu values
  mu_param <- NULL
  
  # create a named parameters list to pass to all likelihood calcs
  parameters <- list(survival = survival,
                     fecundity = fecundity)

  # subset data_modules to those that aren't predictors
  data_modules_response <- data_modules[not_predictors]
  
  for (i in seq_along(data_modules_response)) {
    
    data_tmp <- data_modules_response[[i]]

    # switch based on type
    loglik_fun <- switch(data_tmp$type,
                         age_abundance = age_abundance_loglik,
                         stage_abundance = stage_abundance_loglik,
                         age_recapture = age_recapture_loglik,
                         stage_recapture = stage_recapture_loglik,
                         size_to_age = size_to_age_loglik,
                         age_to_size = age_to_size_loglik)
    
    # calculate
    loglik_fun(data_tmp$data, parameters)
    
  } 
  
  # return outputs with class definition 
  as.integrated_model(parameters)
  
}

#' @export
#' @rdname integrated_model
#' 
is.integrated_model <- function(object) {
  inherits(object, 'integrated_model')
}

#' @export
#' @rdname integrated_model
#' 
print.integrated_model <- function(x, ...) {
  cat(paste0('This is an integrated_model object\n'))
}

#' @export
#' @rdname integrated_model
#' 
plot.integrated_model <- function(x, ...) {
  
  plot(x$greta_model, ...)
  
}

#' @export
#' @rdname integrated_model
#' 
summary.integrated_model <- function(object, ...) {
  
  NULL
  
}

# internal function: create integrated_model object
as.integrated_model <- function(object) {
  as_class(object, name = "integrated_model", type = "list")
}
