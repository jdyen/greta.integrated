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

  # collect inputs
  data_modules <- list(...)

  # is everything an integrated_data object?
  if(!all(sapply(data_modules, is.integrated_data)))
    stop("integrated_model only takes integrated_data objects as arguments", call. = FALSE)

  # which data types are we dealing with?
  data_types <- sapply(data_modules, function(x) x$data_type)
  
  # are all of these implemented?
  implemented <- c("age_abundance", "stage_abundance", "age_recapture",
                   "stage_recapture", "stage_to_age", "age_to_stage",
                   "predictors")
  if (!all(data_types %in% implemented)) {
    problem_types <- data_types[!(data_types %in% implemented))]
    stop(paste0("one or more data types are not currently implemented (", paste(problem_types, collapse = ", "), ")"), call. = FALSE)
  }
  
  # are any data modules predictors? Deal with this second (expand process)
  is_predictor <- data_types == "predictors"
  any_predictors <- any(is_predictor)
  
  # how many different process models are we dealing with?
  process_hash <- sapply(data_modules, function(x) x$process$hash)
  unique_process_hash <- unique(process_hash)
  n_process <- length(unique_process_hash)
  process_id <- match(process_hash, unique_process_hash)

  # check that we haven't doubled up on predictors
  if (any_predictors) {
    process_predictors <- tapply(is_predictor, process_id, sum)
    if (any(process_predictors) > 1)
      stop("cannot assign multiple sets of predictors to the same process", call. = FALSE)
  }
      
  # let's work through the processes one-by-one
  for (i in seq_len(n_process)) {
    
    # pull out any process that matches hash i
    process_list[[i]] <- data_modules[[which(process_id == i)[1]]]$process
    
    # does this have predictors?
    includes_predictors <- i %in% process_id[is_predictor]
    
    # if it includes predictors, we need to massage the structure carefully
    if (includes_predictors) {
      
      # pull out the predictor data set (can't be more than one per process)
      predictors <- data_module[[which(process_id == i & is_predictor)]]
      
      # pull together a set of parameters based on the details of this process
      parameters[[i]] <- define_parameters(process_list[[i]]$classes, process_list[[i]]$priors, process_list[[i]]$masks, predictors)
      
    } else {
      
      # pull together a set of parameters based on the details of this process
      parameters[[i]] <- define_parameters(process_list[[i]]$classes, process_list[[i]]$priors, process_list[[i]]$masks)
      
    }
    
    
  }

  ### PARAMETERS are with the processes??
  # create a named parameters list to pass to all likelihood calcs
  parameters <- list(n_obs,
                     n_iter,
                     classes,
                     matrix,
                     inits,
                     density,
                     process_class,
                     age_to_stage_conversion,
                     stage_to_age_conversion,
                     bias)

  # subset data_modules to those that aren't predictors
  if (any_predictors) {
    data_modules <- data_modules[!is_predictor]
    process_id <- process_id[!is_predictor]
  }
  
  for (i in seq_along(data_modules)) {
  
    process_tmp <- process_list[[process_id[i]]]
      
    data_tmp <- data_modules[[i]]

    # switch based on type
    loglik_fun <- switch(data_tmp$type,
                         age_abundance = age_abundance_loglik,
                         stage_abundance = stage_abundance_loglik,
                         binned_age_recapture = binned_age_recapture_loglik,
                         binned_stage_recapture = binned_stage_recapture_loglik,
                         binary_age_recapture = binary_age_recapture_loglik,
                         binary_stage_recapture = binary_stage_recapture_loglik,
                         stage_to_age = stage_to_age_loglik,
                         age_to_stage = age_to_stage_loglik)
    
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
