#' @name integrated_data
#' @title create integrated data objects
#' 
#' @description \code{integrated_data} defines a data object with appropriate likelihood based on
#'  a process model defined with \link[greta.integrated]{integrated_process}
#' 
#' @param data a single data input (see details for descriptions of possible input types)
#' @param process an \link[greta.integrated]{integrated_process} object
#' @param bias a bias function that connects observed and modelled data (see details)
#' @param settings a named list of settings passed to data formatting functions (see details)
#' @param ... additional arguments to \link[base]{print} and \link[base]{summary} methods (currently ignored)
#' @param x an \code{integrated_data} object
#' @param object an \code{integrated_data} object
#'
#' @details Do something. The settings list can be used to specify how the data are binned, either with
#'   specific breaks for binning or with the number of breaks to use. If these are not provided, the 
#'   functions use the \code{classes} element of \code{process} to determine the number
#'   of bins (\code{nbreaks = classes + 1}).
#' 
#' @return An object of class \code{integrated_data}, which contains information on the data module and
#'   can be passed to \link[greta.integrated]{integrated_model}
#' 
#' @import greta
#' 
#' @examples
#' \dontrun{
#' 
#' library(integrated)
#' 
#' # prepare an example model
#' data <- add_data()
#'                         
#' # summarise data module
#' model
#' summary(model)
#' plot(model)
#' }

#' @export
#' @rdname integrated_data
#' 
abundance <- function(data,
                      process,
                      bias,
                      settings = list()) {

  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # need to treat matrix-like data and list data differently
  if (is.data.frame(data) | is.matrix(data)) {
    
    # are the data formatted as if there's an individual age or size in each column?
    if (ncol(data) != process$classes) {
      
      # edge case: data are total abundances through time
      if (nrow(data) == 1 | ncol(data) == 1)
        stop("data have one row or column; are you providing total abundance data? If so, try flattening data to a vector", call. = FALSE)
      
      # if data aren't provided as counts, we need breaks
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if data are not entered as counts", call. = FALSE)
        
      # what's the mismatch between classes and data?
      greater_less <- ifelse(ncol(data) > process$classes, "more", "fewer")
      
      # let the user know what's going on
      cat(paste0("data has ", greater_less, " columns than classes; each row will be treated as individual ", class_type, "s and not counts\n"))
      
      # use hist_fn to bin the data according to the provided breaks
      data_clean <- t(apply(data, 1, hist_fn, breaks = all_settings$breaks))
      
    } else {
      
      # most likely case: data provided as binned counts per class
      cat(paste0("data has one column for each class; each row will be treated as counts of individuals per class\n"))
      
      # data shouldn't need any work; return as is
      data_clean <- data
      
    }
    
  } else {
    
    # are data formatted in a list with multiple entries?
    if (is.list(data)) {
      
      # how long is each element?
      list_len <- sapply(data, length)
      
      # if they're all the same, might be binned data provided as a list rather than matrix-like
      if (all(list_len == classes)) {
        
        # let the user know what's going on
        cat(paste0("each element of data has an entry for each class; elements will be treated as counts of individuals per class\n"))

        # directly convert list to matrix (each element has same length)        
        data_clean <- do.call(rbind, data)
        
      } else {  # elements differ in length, probably have individual sizes or ages
        
        # if so, we need breaks to bin the data        
        if (is.null(all_settings$breaks))
          stop("breaks must be provided if data are not entered as counts", call. = FALSE)
        
        # let the user know what's going on
        cat(paste0("converting from list to matrix; are data ", class_type, "s and not counts?\n"))
        
        # bin the data using hist_fn, applied to each element of the input list
        data_clean <- t(sapply(data, hist_fn, breaks = all_settings$breaks))
        
      }
      
    } else {
      
      # edge case: data are total abundances through time
      if (is.numeric(data) | is.integer(data))
        stop("are you providing total abundance data as a single vector? These aren't supported yet but will be soon", call. = FALSE)

      stop("abundance data must be a matrix, data.frame, or list", call. = FALSE) 
      
    }
    
  }

  # want to tie everything together in a single output
  data_module <- list(data = clean_data,
                      process = process, 
                      bias = bias)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
cmr <- function(data,
                process,
                bias,
                settings = list()) {
  
  # turn data into a list and check that each element has the correct columns
  if (is.matrix(data) | is.data.frame(data)) {
    
    data <- list(data)
    
  } else {
    
    if (!is.list(data)) {
      stop("stage_recapture data must be a matrix, data.frame, or list of matrices or data.frames",
           call. = FALSE) 
    }
    
  }
  
  # check data format
  for (i in seq_along(data)) {
    if (!("size" %in% colnames(data[[i]]))) {
      stop("stage_recapture models require size measurements at each recapture",
           call. = FALSE)
    }
    
    if (!all(c("size", "id", "time") %in% colnames(data[[i]]))) { 
      stop("stage_recapture data should be in long format with size, id, and time columns",
           call. = FALSE)
    }
  }
  
  # prepare cmr histories
  # mark-recapture model to estimate detection and survival probabilities
  # calculate size-based catch history for each individual
  catch_size <- with(cmr_data, tapply(weight, list(idfish, year), mean))
  catch_size <- ifelse(is.na(catch_size), 0, catch_size)
  
  # observed at least once?
  observed <- apply(catch_size, 1, sum) > 0
  catch_size <- catch_size[observed, ]
  
  # convert to size classes
  catch_size_class <- matrix(cut(catch_size, size_breaks, labels = FALSE),
                             ncol = ncol(catch_size))
  catch_size_class <- ifelse(is.na(catch_size_class), 0, catch_size_class)
  
  # first and final size classes
  first_size_class <- apply(catch_size_class, 1, function(x) x[min(which(x > 0))])
  final_size_class <- apply(catch_size_class, 1, function(x) x[max(which(x > 0))])
  
  # has it shrunk *many* classes?
  size_errors <- final_size_class < (first_size_class - 1)
  
  # if so, remove these observations
  catch_size_class <- catch_size_class[!size_errors, ]
  first_size_class <- first_size_class[!size_errors]
  final_size_class <- final_size_class[!size_errors]
  
  # calculate first and final observations
  first_obs <- apply(catch_size_class, 1, function(x) min(which(x > 0)))
  final_obs <- apply(catch_size_class, 1, function(x) max(which(x > 0)))
  
  # number of years alive
  n_alive <- final_obs - first_obs
  
  # are any individuals never recaptured?
  single_obs <- first_obs == final_obs
  
  # focus on those with >1 observation
  n_alive <- n_alive[!single_obs]
  first_sizes <- first_size_class[!single_obs]
  first_seen <- first_obs[!single_obs]
  last_seen <- final_obs[!single_obs]
  

  
  data_module <- list(data = clean_data,
                      process = process, 
                      bias = bias)
  
  # return outputs
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
predictors <- function(data,
                       process,
                       bias,
                       settings = list()) {
  
  # prepare predictor data
  
  data_module <- list(data = clean_data,
                      process = process,
                      bias = bias)
  
  # return outputs
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
growth <- function(data,
                   process,
                   bias,
                   settings = list()) {
  
  # prepare growth trajectories
  # convert matrix or data.frame data to a list
  if (is.matrix(data) | is.data.frame(data)) {
    
    # check if data are formatted correctly
    if (ncol(data) != nrow(data)) {
      data <- make_growth_data_matrix(data = data,
                                      classes = integrated_process$classes,
                                      settings = settings)
    } 
    
    data <- list(data)
    
  }
  
  # data should be a list or matrix
  if (!is.list(data)) {
    stop("growth data must be a matrix, data.frame, or list of matrices or data.frames",
         call. = FALSE)
  }
  
  # error if the number of data classes doesn't match the process
  if (!all(integrated_process$classes == sapply(data, nrow))) {
    stop("number of classes in growth data must match number of classes in integrated_process")
  }
  
  if (integrated_process$replicates > 1) {      
    if (length(data) > 1) {
      if (length(data) != length(integrated_process$replicate_id)) {
        stop(paste0("growth data have ", length(data), " elements but ",
                    " integrated process has ", integrated_process$replicates,
                    " replicates"),
             call. = FALSE)
      }
    } else {
      stop(paste0("one growth data set should be supplied for each of the ",
                  integrated_process$replicates, " process replicates"),
           call. = FALSE)
    }
  }
  
  # create data module from growth data matrix
  data_module <- define_individual_growth_module(data = data,
                                                 integrated_process = integrated_process, 
                                                 observation_model = observation_model)
  
  data_module <- list(data = clean_data,
                      process = process, 
                      bias = bias)

  
  # return outputs
  as.integrated_data(data_module)
  
}   

#' @export
#' @rdname integrated_data
#' 
community <- function(data,
                      process,
                      bias,
                      settings = list()) {
  
  # prepare cmr histories
  
  data_module <- list(data = clean_data,
                      process = process, 
                      bias = bias)
  
  # return outputs
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
is.integrated_data <- function(object) {
  inherits(object, integrated_data)
}

#' @export
#' @rdname integrated_data
#' 
print.integrated_data <- function(x, ...) {
  cat(paste0(This is an integrated_data object\n))
}

#' @export
#' @rdname integrated_data
#' 
summary.integrated_data <- function(object, ...) {
  
  NULL
  
}

# internal function: create integrated_data object
as.integrated_data <- function(object) {
  as_class(object, name = integrated_data, type = list)
}
