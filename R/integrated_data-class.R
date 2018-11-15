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
  all_settings <- list(breaks = NULL,
                       likelihood = poisson)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  process_class <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  is_list <- is.list(data) & !is_matrix_like
  not_ok <- !is_matrix_like & !is.list
  
  if (not_ok)
    stop("age_abundance data must be a matrix, data.frame, or list", call. = FALSE) 

  # need to treat matrix-like data and list data differently
  if (is_matrix_like) {

    # default: return data as is if there's a column for each class
    data_clean <- data
    
    # is this a model with age data but a stage structure?
    if (process_class == "stage")
      warning("it looks like you're fitting a stage-structured model to age data; is this correct?", call. = FALSE)

    # potentially reformat data if there isn't one column per class
    if (ncol(data) != process$classes & process_class != "stage") {
      
      # edge case: data are total abundances through time; need to use a different function
      if (nrow(data) == 1)
        stop("data only have one row; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      if (ncol(data) == 1)
        stop("data only have one column; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      
      stop(paste0("data should have one column per age class but have ", ncol(data), " columns and there are ", process$classes, " age classes"), call. = FALSE)

    }
    
  }
  
  if (is_list) {
    
    # how long is each element?
    list_len <- sapply(data, length)
    
    # if they're all the same, might be binned data provided as a list rather than matrix-like
    if (all(list_len == classes)) {
      
      data_clean <- do.call(rbind, data)
      
    } else {  # elements differ in length, probably have individual ages
      
      # can't handle ages older than number of classes
      if (max(data) > process$classes) {
        warning(paste0("some individuals are older than the total number of classes; ages will be truncated to a maximum of ", process$classes), call. = FALSE)
        data <- ifelse(data > process$classes, process$classes, data)
      }

      # if so, we need breaks to bin the data        
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if data are not binned", call. = FALSE)
      
      # let the user know what's going on
      cat(paste0("converting data from list to matrix; this assumes data are individual ", class_type, "s and not counts\n"))
      
      # bin the data using hist_fn, applied to each element of the input list
      data_clean <- t(sapply(data, hist_fn, breaks = all_settings$breaks))
      
    }
    
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "age_abundance",
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
stage_abundance <- function(data, process, bias = no_bias(), settings = list()) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = poisson)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  process_class <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  is_list <- is.list(data) & !is_matrix_like
  not_ok <- !is_matrix_like & !is.list
  
  if (not_ok)
    stop("stage_abundance data must be a matrix, data.frame, or list", call. = FALSE) 
  
  # need to treat matrix-like data and list data differently
  if (is_matrix_like) {
    
    # default: return data as is if there's a column for each class
    data_clean <- data
    
    # is this a model with age data but a stage structure?
    if (process_class == "age")
      warning("it looks like you're fitting an age-structured model to stage data; is this correct?", call. = FALSE)
    
    # potentially reformat data if there isn't one column per class
    if (ncol(data) != process$classes & process_class != "age") {
      
      # edge case: data are total abundances through time; need to use a different function
      if (nrow(data) == 1)
        stop("data only have one row; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      if (ncol(data) == 1)
        stop("data only have one column; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      
      stop(paste0("data should have one column per stage but have ", ncol(data), " columns and there are ", process$classes, " stages"), call. = FALSE)
      
    }
    
  }
  
  if (is_list) {
    
    # how long is each element?
    list_len <- sapply(data, length)
    
    # if they're all the same, might be binned data provided as a list rather than matrix-like
    if (all(list_len == classes)) {
      
      data_clean <- do.call(rbind, data)
      
    } else {  # elements differ in length, probably have individual ages
      
      # can't handle ages older than number of classes
      if (max(data) > process$classes)
        stop(paste0("there are ", max(data), " stages and ", process$classes, " classes; the number of stages should not exceed the number of classes"), call. = FALSE)

      # if ok, we need breaks to bin the data        
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if data are not binned", call. = FALSE)
      
      # let the user know what's going on
      cat(paste0("converting data from list to matrix; this assumes data are individual ", class_type, "s and not counts\n"))
      
      # bin the data using hist_fn, applied to each element of the input list
      data_clean <- t(sapply(data, hist_fn, breaks = all_settings$breaks))
      
    }
    
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "stage_abundance",
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
age_recapture <- function(data, process, bias = no_bias(), settings = list()) {
  
  # need to unpack the settings
  all_settings <- list(likelihood = multinomial)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)

  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("age_recapture data must be a matrix or data.frame containing capture histories", call. = FALSE) 
  
  # do the entries look like binned or binary data?
  is_binary <- all(unique(data) %in% c(0, 1))
  is_binned <- all(unique(data) %in% seq_len(process$classes))
  
  # if data aren't binned we need to bin them
  if (!is_binary & !is_binned)
    stop(paste0("range of data does not match the number of classes (", process$classes, ")"), call. = FALSE)
  
  # can only do so much with binary data
  if (is_binary) {
    
    # let the user know that we don't like binary data
    cat(paste0("binary capture histories can only inform detection and total survival probabilities\n"))
    
    # return what we have
    data_clean <- data
    
    # set data type so the likelihood can be worked out quickly
    data_type <- "binary_age_recapture"
    
  }
  
  # return as is if already binned
  if (is_binned) {
    
    # return what we have
    data_clean <- data
    
    # set the data type so we know which likelihood to use
    data_type <- "binned_age_recapture"
    
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
    
    # has any individual gotten younger through time?
    class_errors <- final_class < first_class
    
    # this is fine in some models, so just give a warning
    warning(paste0(sum(class_errors), " individuals seemed to get younger through time, is this reasonable?"), call. = FALSE)
    
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = data_type,
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
stage_recapture <- function(data, process, bias = no_bias(), settings = list()) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = multinomial)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  
  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("stage_recapture data must be a matrix or data.frame containing capture histories", call. = FALSE) 
  
  # do the entries look like binned or binary data?
  is_binary <- all(unique(data) %in% c(0, 1))
  is_binned <- all(unique(data) %in% seq_len(process$classes))
  
  # if data aren't binned we need to bin them
  if (!is_binary & !is_binned) {
    
    # this won't work if we don't have breaks
    if (is.null(all_settings$breaks))
      stop("data do not appear to be binned or binary; attempted to bin data but breaks were not provided", call. = FALSE)

    # warn that this is not a great way to set up data
    warning("data do not appear to be binned or binary; data will be binned according to settings$breaks", call. = FALSE)
    
    # the cut function can bin the data
    data_binned <- matrix(cut(data, all_settings$breaks, labels = FALSE), ncol = ncol(data))
    
    # just need to clean up some NAs afterwards
    data_binned <- ifelse(is.na(data_binned), 0, data_binned)
    
    # that's it, just return this
    data_clean <- data_binned
    
    # set the data type so we know which likelihood to use
    data_type <- "binned_stage_recapture"
    
  }
  
  # can only do so much with binary data
  if (is_binary) {
    
    # let the user know that we don't like binary data
    cat(paste0("binary capture histories can only inform detection and total survival probabilities\n"))
    
    # return what we have
    data_clean <- data
    
    # set data type so the likelihood can be worked out quickly
    data_type <- "binary_stage_recapture"
    
  }
  
  # return as is if already binned
  if (is_binned) {
    
    # return what we have
    data_clean <- data
    
    # set the data type so we know which likelihood to use
    data_type <- "binned_stage_recapture"
    
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
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = data_type,
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
} 

#' @export
#' @rdname integrated_data
#' 
stage_to_age <- function(x, ...) {
  UseMethod("stage_to_age")
}

stage_to_age.formula <- function(x, data, process, bias = no_bias(), settings = list()) {
  
  # this won't work if we haven't got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to model ages from stage data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = multinomial)
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
    cat("attempting to bin data because stage data are not rounded\n")
    
    # if data are binned, we need breaks to bin them
    if (is.null(all_settings$breaks))
      stop("breaks must be provided if stage-to-age data are not binned", call. = FALSE)
    
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
                      data_type = "stage_to_age",
                      likelihood = all_settings$likelihood)
  
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

stage_to_age.default <- function(x, process, bias = no_bias(), settings = list()) {

  # this won't work if we haven't got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to model ages from stage data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # unpack settings
  all_settings <- list(likelihood = multinomial)
  all_settings[names(settings)] <- settings
  
  # is the input data a matrix-like data type?
  if (!is.matrix(x) & !is.data.frame(x))
    stop("stage-to-age data must be a matrix or data.frame", call. = FALSE)

  # what dims does x have?
  dims <- dim(x)
  
  # if at least one of columns or rows don't line up with the number of classes, this won't work
  if (!(process$classes %in% dims))
    stop(("one dimension of stage-to-age data must match the number of classes in process (", process$classes, ")"), call. = FALSE)

  # format so that ages are always in columns
  if (dims[2] != process$classes) {
    
    data_clean <- t(data)
    
  } else {
    
    # if ncol(x) == nrow(x), warn that we assume ages in columns
    if (dims[1] == dims[2])
      cat("stage-to-age data must have ages in columns; is this data set formatted correctly?", call. = FALSE)
    
    data_clean <- data
    
  }
          
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "stage_to_age",
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}   

#' @export
#' @rdname integrated_data
#' 
age_to_stage <- function(x, ...) {
  UseMethod("age_to_stage")
}

age_to_stage.formula <- function(x, data, process, bias = no_bias(), settings = list()) {
  
  # this won't work if we haven't got a stage-structured process model
  if (process$type == "leslie")
    stop("trying to model stages from age data without a stage-based model; this seems like a bad idea", call. = FALSE)
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = multinomial)
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
    cat("attempting to bin data because age data are not rounded\n")
    
    # if data are binned, we need breaks to bin them
    if (is.null(all_settings$breaks))
      stop("breaks must be provided if age-to-stage data are not binned", call. = FALSE)
    
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
                      data_type = "age_to_stage",
                      likelihood = all_settings$likelihood)
  
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

age_to_stage.default <- function(x, process, bias = no_bias(), settings = list()) {
  
  # this won't work if we haven't got a stage-structured process model
  if (process$type == "leslie")
    stop("trying to model stages from age data without a stage-based model; this seems like a bad idea", call. = FALSE)
  
  # unpack settings
  all_settings <- list(likelihood = multinomial)
  all_settings[names(settings)] <- settings
  
  # is the input data a matrix-like data type?
  if (!is.matrix(x) & !is.data.frame(x))
    stop("age-to-stage data must be a matrix or data.frame", call. = FALSE)
  
  # what dims does x have?
  dims <- dim(x)
  
  # if at least one of columns or rows don't line up with the number of classes, this won't work
  if (!(process$classes %in% dims))
    stop(("one dimension of age-to-stage data must match the number of classes in process (", process$classes, ")"), call. = FALSE)
  
  # format so that ages are always in columns
  if (dims[2] != process$classes) {
    
    data_clean <- t(data)
    
  } else {
    
    # if ncol(x) == nrow(x), warn that we assume ages in columns
    if (dims[1] == dims[2])
      cat("age-to-stage data must have stages in columns; is this data set formatted correctly?", call. = FALSE)
    
    data_clean <- data
    
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data_clean,
                      process = process, 
                      bias = bias,
                      data_type = "age_to_stage",
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}   

#' @export
#' @rdname integrated_data
#' 
community <- function(data, process, bias = no_bias(), settings = list()) {
  
  warning("community data are not currently implemented; this integrated_data object will be ignored in subsequent models", call. = FALSE)
  
  # unpack settings
  all_settings <- list(likelihood = binomial)
  all_settings[names(settings)] <- settings
  
  # want to tie everything together in a single output
  data_module <- list(data = NULL,
                      process = process, 
                      bias = bias,
                      data_type = "community",
                      likelihood = all_settings$likelihood)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
predictors <- function(data, process, bias = no_bias(), settings = list()) {
  
  # unpack settings
  all_settings <- list(likelihood = normal)
  all_settings[names(settings)] <- settings
  
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
                      data_type = "predictor",
                      likelihood = all_settings$likelihood)
  
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
