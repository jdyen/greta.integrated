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
age_abundance <- function(data, process, bias = no_bias(), settings = list()) {

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
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "abundance")
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
stage_abundance <- function(data, process, bias = no_bias(), settings = list()) {
  
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
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "abundance")
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
age_recapture <- function(data, process, bias = no_bias(), settings = list()) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # check that our data look like they might be capture histories
  if (is.matrix(data) | is.data.frame(data)) {
    
    # do the entries look like binned or binary data?
    is_binary <- all(unique(data) %in% c(0, 1))
    is_binned <- all(unique(data) %in% seq_len(process$classes))
    
    # if data aren't binned we need to bin them
    if (!is_binary & !is_binned) {
      
      # this won't work if we don't have breaks
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if recapture data are not binned or binary", call. = FALSE)
      
      # the cut function can bin the data
      data_binned <- matrix(cut(data, all_settings$breaks, labels = FALSE), ncol = ncol(data))
      
      # just need to clean up some NAs afterwards
      data_binned <- ifelse(is.na(data_binned), 0, data_binned)
      
      # that's it, just return this
      data_clean <- data_binned
      
      # set the data type so we know which likelihood to use
      data_type <- "binned_recapture"
      
    }
    
    # can only do so much with binary data
    if (is_binary) {
      
      # let the user know that we don't like binary data
      cat(paste0("binary capture histories can only inform detection and total survival probabilities\n"))
      
      # return what we have
      data_clean <- data
      
      # set data type so the likelihood can be worked out quickly
      data_type <- "binary_recapture"
      
    }
    
    # return as is if already bined
    if (is_binned) {

      # return what we have
      data_clean <- data
      
      # set the data type so we know which likelihood to use
      data_type <- "binned_recapture"
      
    }
    
    # are any individuals not observed at least once?
    not_observed <- apply(data_clean, 1, sum) == 0
    
    # if not, remove unobserved individuals (with a warning)
    if (any(not_observed)) {
      warning(paste0("removing ", sum(not_observed), " individuals that were never observed"), call. = FALSE)
      data_clean <- data_clean[!not_observed, ]
    }
    
    # check that bin IDs are logical
    if (!is_binary) {
      
      # what are the first and final classes?
      first_class <- apply(data_clean, 1, function(x) x[min(which(x > 0))])
      final_class <- apply(data_clean, 1, function(x) x[max(which(x > 0))])
       
      # has any individual regressed multiple classes?
      class_errors <- final_class < (first_class - 1)

      # this is fine in some models, so just give a warning
      warning(paste0(sum(class_errors), " individuals regressed by two or more classes, is this reasonable?"), call. = FALSE)
      
    }
    
  } else {
    
    # can't handle other data types
    stop("recapture data must be a matrix of capture histories", call. = FALSE)
    
  }

  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = data_type)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
stage_recapture <- function(data, process, bias = no_bias(), settings = list()) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # check that our data look like they might be capture histories
  if (is.matrix(data) | is.data.frame(data)) {
    
    # do the entries look like binned or binary data?
    is_binary <- all(unique(data) %in% c(0, 1))
    is_binned <- all(unique(data) %in% seq_len(process$classes))
    
    # if data aren't binned we need to bin them
    if (!is_binary & !is_binned) {
      
      # this won't work if we don't have breaks
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if recapture data are not binned or binary", call. = FALSE)
      
      # the cut function can bin the data
      data_binned <- matrix(cut(data, all_settings$breaks, labels = FALSE), ncol = ncol(data))
      
      # just need to clean up some NAs afterwards
      data_binned <- ifelse(is.na(data_binned), 0, data_binned)
      
      # that's it, just return this
      data_clean <- data_binned
      
      # set the data type so we know which likelihood to use
      data_type <- "binned_recapture"
      
    }
    
    # can only do so much with binary data
    if (is_binary) {
      
      # let the user know that we don't like binary data
      cat(paste0("binary capture histories can only inform detection and total survival probabilities\n"))
      
      # return what we have
      data_clean <- data
      
      # set data type so the likelihood can be worked out quickly
      data_type <- "binary_recapture"
      
    }
    
    # return as is if already bined
    if (is_binned) {
      
      # return what we have
      data_clean <- data
      
      # set the data type so we know which likelihood to use
      data_type <- "binned_recapture"
      
    }
    
    # are any individuals not observed at least once?
    not_observed <- apply(data_clean, 1, sum) == 0
    
    # if not, remove unobserved individuals (with a warning)
    if (any(not_observed)) {
      warning(paste0("removing ", sum(not_observed), " individuals that were never observed"), call. = FALSE)
      data_clean <- data_clean[!not_observed, ]
    }
    
    # check that bin IDs are logical
    if (!is_binary) {
      
      # what are the first and final classes?
      first_class <- apply(data_clean, 1, function(x) x[min(which(x > 0))])
      final_class <- apply(data_clean, 1, function(x) x[max(which(x > 0))])
      
      # has any individual regressed multiple classes?
      class_errors <- final_class < (first_class - 1)
      
      # this is fine in some models, so just give a warning
      warning(paste0(sum(class_errors), " individuals regressed by two or more classes, is this reasonable?"), call. = FALSE)
      
    }
    
  } else {
    
    # can't handle other data types
    stop("recapture data must be a matrix of capture histories", call. = FALSE)
    
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = data_type)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
size_at_age <- function(x, ...) {
  UseMethod("size_at_age")
}

size_at_age.formula <- function(x, data, process, bias = no_bias(), settings = list()) {
  
  # this won't work if we haven't got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to match size and age data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # parse formula to give outputs
  var_names <- all.vars(x)
  response <- get(var_names[1], data)
  predictor <- get(var_names[2], data)
  
  # basic checks to see we haven't missed something
  if (length(response) != length(predictor))
    stop(paste0(var_names[1], " and ", var_names[2], " should be the same length"), call. = FALSE)
  
  # are the data likely to be binned?
  is_binned <- all(response %% 1 == 0)
  
  # if data are not binned, we need to change this
  if (!is_binned) {
    
    # let the user know what's going on
    cat("attempting to bin data because size data are not rounded\n")
    
    # if data are binned, we need breaks to bin them
    if (is.null(all_settings$breaks))
      stop("breaks must be provided if size-at-age data are not binned", call. = FALSE)
    
    # bin away
    response <- cut(response, all_settings$breaks, labels = FALSE)
    
  }
  
  # lots of ones can help us count up transition categories
  ones_vec <- rep(1, length(response))
  
  # need to truncate ages if they exceed the number of available classes
  predictor <- ifelse(predictor >= process$classes, process_classes, predictor)
  
  # need to turn vectors into a transition matrix
  n_bin <- max(response)
  data_clean <- matrix(0, nrow = n_bin, ncol = process$classes)
  for (i in seq_len(n_bin)) {
    predictor_sub <- predictor[response == i]
    classification <- tapply(ones_vec, predictor_sub, sum)
    data_clean[i, match(names(classification), seq_len(process$classes))] <- classification
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "size_at_age")
  
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

size_at_age.default <- function(x, process, bias = no_bias(), settings = list()) {

  # this won't work if we haven't got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to match size and age data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # is the input data a matrix-like data type?
  if (!is.matrix(x) & !is.data.frame(x))
    stop("size-at-age data must be a matrix or data.frame", call. = FALSE)

  # what dims does x have?
  dims <- dim(x)
  
  # if at least one of columns or rows don't line up with the number of classes, this won't work
  if (!(process$classes %in% dims))
    stop(("one dimension of size-at-age data must match the number of classes in process (", process$classes, ")"), call. = FALSE)

  # format so that ages are always in columns
  if (dims[2] != process$classes) {
    
    data_clean <- t(data)
    
  } else {
    
    # if ncol(x) == nrow(x), warn that we assume ages in columns
    if (dims[1] == dims[2])
      cat("size-at-age data must have ages in columns; is this data set formatted correctly?", call. = FALSE)
    
    data_clean <- data
    
  }
          
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "size_at_age")
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}   

#' @export
#' @rdname integrated_data
#' 
community <- function(data, process, bias = no_bias(), settings = list()) {
  
  warning("community data are not currently implemented; this integrated_data object will be ignored in subsequent models", call. = FALSE)
  
  # want to tie everything together in a single output
  data_module <- list(data = NULL,
                      process = process, 
                      bias = bias,
                      data_type = "community")
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
predictors <- function(data, process, bias = no_bias(), settings = list()) {
  
  # make sure predictor data are in a matrix or data.frame
  if (!is.matrix(data) & !is.data.frame(data))
    stop("predictor data must be a matrix or data.frame", call. = FALSE)
  
  # are there any NAs in the data?
  na_col_check <- apply(data, 2, function(x) any(is.na(x)))
  if (any(na_col_check))
    warning("there are missing values in the predictor data; these will be ignored but fitted models will be more reproducible if NAs are handled prior to model fitting", call. = FALSE)
  
  # are there are any completely missing rows?
  na_row_check <- apply(data, 1, function(x) all(is.na(x)))
  
  # if so, warn and remove
  if (any(na_row_check)) {
    warning(paste0("there are ", sum(na_row_check), " rows with completely missing data; these will be removed from all analyses"), call. = FALSE)
    data_clean <- data[!na_row_check, ]
  } else {
    data_clean <- data
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "predictor")
  
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
