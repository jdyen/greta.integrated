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
#' @import greta.dynamics
#' @import tensorflow
#' @import methods
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
                   "stage_recapture", "stage_to_age", "age_to_stage")
  if (!all(data_types %in% implemented)) {
    problem_types <- data_types[!(data_types %in% implemented)]
    stop(paste0("one or more data types are not currently implemented (", paste(problem_types, collapse = ", "), ")"), call. = FALSE)
  }
  
  # how many different process models are we dealing with?
  process_hash <- sapply(data_modules, function(x) x$process$hash)
  unique_process_hash <- unique(process_hash)
  n_process <- length(unique_process_hash)
  process_id <- match(process_hash, unique_process_hash)

  # let's work through the processes one-by-one
  process_list <- vector("list", length = n_process)
  parameters <- vector("list", length = n_process)
  for (i in seq_len(n_process)) {
    
    # pull out all data_modules corresponding to process i
    data_sub <- data_modules[process_id == i]
    
    # are any data modules predictors? Deal with this second (expand process)
    includes_predictors <- sapply(data_sub, function(x) !is.null(x$predictors))

    # need to know number of predictors if they're included
    n_predictors <- NULL
    if (any(includes_predictors)) {
      
      # can't work with predictors if they're not shared by all data sets
      if (!all(includes_predictors))
        stop("predictors must be included for all data modules that share a single process", call. = FALSE)
      
      # how many predictors?
      n_predictors <- sapply(data_sub[includes_predictors], function(x) ncol(x$predictors))
      
      # do these match in number? (assuming they match in predictor sets too; perhaps warn about this?)
      if (length(unique(n_predictors)) > 1)
        stop("all data modules that share a single process must have the same number of predictors", call. = FALSE)
      
    }

    # pull out any process that matches hash i
    process_list[[i]] <- data_modules[[which(process_id == i)[1]]]$process
    
    # do we need to deal with age-stage conversions?
    classes_alt <- sapply(data_modules[which(process_id == i)], function(x) x$classes_alt)
    
    # would like this not be a list
    classes_alt <- unlist(classes_alt)
    
    # if so, we need to make sure there's only one set of conversions
    age_or_stage <- ifelse(process_list[[i]]$type == "leslie", "stage-to-age", "age-to-stage")
    process_type <- ifelse(process_list[[i]]$type == "leslie", "age", "stage")
    if (length(unique(classes_alt)) > 1)
      stop(paste0("there are multiple ", age_or_stage, " conversions with different dimensions; perhaps set up separate process models for each"), call. = FALSE)
    classes_alt <- unique(classes_alt)

    # if it includes predictors, we need to massage the structure carefully
    if (any(includes_predictors)) {
      
      # pull together a set of parameters based on the details of this process
      parameters[[i]] <- define_parameters(process = process_list[[i]],
                                           classes_alt = classes_alt,
                                           n_predictors = n_predictors)

    } else {
      
      # pull together a set of parameters based on the details of this process
      parameters[[i]] <- define_parameters(process = process_list[[i]],
                                           classes_alt = classes_alt)
      
    }
    
  }

  for (i in seq_along(data_modules)) {
  
    # pull out the parameters we need
    parameters_tmp <- parameters[[process_id[i]]]
    
    # we need to add a couple of things from the process module
    parameters_tmp$density <- process_list[[process_id[i]]]$density
    parameters_tmp$process_class <- process_list[[process_id[i]]]$type

    # prepare matrix
    parameters_tmp$matrix <- construct_matrix(data_modules[[i]], parameters_tmp)
    
    # choose appropriate likelihood based on type of data
    loglik_fun <- switch(data_modules[[i]]$data_type,
                         age_abundance = age_abundance_loglik,
                         stage_abundance = stage_abundance_loglik,
                         binned_age_recapture = binned_age_recapture_loglik,
                         binned_stage_recapture = binned_stage_recapture_loglik,
                         binary_age_recapture = binary_age_recapture_loglik,
                         binary_stage_recapture = binary_stage_recapture_loglik,
                         stage_to_age = stage_to_age_loglik,
                         age_to_stage = age_to_stage_loglik)
    
    # define likelihood (doesn't return anything)
    loglik_fun(data_modules[[i]], parameters_tmp)
    
  } 
  
  # return a clean list of parameters
  names(parameters) <- paste0("process_", number_to_word(seq_len(n_process)))
  
  # return outputs with class definition 
  as.integrated_model(parameters)
  
}

#' @export
#' @rdname integrated_model
#' 
is.integrated_model <- function(object) {
  inherits(object, "integrated_model")
}

#' @export
#' @rdname integrated_model
#' 
print.integrated_model <- function(x, ...) {
  cat(paste0("This is an integrated_model object\n"))
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
