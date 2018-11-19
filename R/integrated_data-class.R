#' @name integrated_data
#' @title create integrated data objects
#' 
#' @description \code{integrated_data} defines a data object with custom likelihood based on
#'  a process model defined with \link[greta.integrated]{integrated_process}
#' 
#' @param data a single data input (see details for descriptions of possible input types)
#' @param process an \link[greta.integrated]{integrated_process} object
#' @param likelihood something
#' @param bias a bias function that connects observed and modelled data (see details)
#' @param settings a named list of settings passed to data formatting functions (see details)
#' @param predictors
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
#' @export
#' 
#' @import greta
#' @import greta.dynamics
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
add_data <- function(data, process, likelihood, bias = no_bias(), settings = list(), predictors = NULL) {
  
  # compile data for the chosen likelihood
  data <- do.call(paste0(likelihood$type, "_internal"), list(data, process, settings))

  # are predictors provided?
  if (!is.null(predictors)) {
    
    if (!("integrated_predictor" %in% class(predictors)))
      stop("predictors must be created with the `add_predictors()` function", call. = FALSE)
    
    # check that dims match and remove observations from data if any predictors are removed
    data <- match_dims(data, predictors)
  
  }
  
  # want to tie everything together in a single output
  data_module <- list(data = data$data,
                      process = process, 
                      likelihood = likelihood,
                      bias = bias,
                      classes = data$classes,
                      predictors = predictors)
  
  # return outputs with class definition
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
age_abundance <- function(distribution = poisson) {
  
  list(type = "age_abundance",
       distribution = distribution)

}

#' @export
#' @rdname integrated_data
#' 
stage_abundance <- function(distribution = poisson) {
  
  list(type = "stage_abundance",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
age_cjs <- function(distribution = binomial) {

  if (!identical(binomial, distribution))
    warning("a binomial distribution will be used for cjs model", call. = FALSE)
    
  distribution <- binomial
  
  list(type = "age_cjs",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
stage_cjs <- function(distribution = binomial) {
  
  if (!identical(binomial, distribution))
    warning("a binomial distribution will be used for cjs model", call. = FALSE)
  
  distribution <- binomial
  
  list(type = "stage_cjs",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
stage_to_age <- function(distribution = multinomial) {
  
  if (!identical(multinomial, distribution))
    warning("a multinomial distribution will be used for stage-to-age model", call. = FALSE)
  
  list(type = "stage_to_age",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
age_to_stage <- function(distribution = multinomial) {
  
  if (!identical(multinomial, distribution))
    warning("a multinomial distribution will be used for age-to-stage model", call. = FALSE)
  
  list(type = "age_to_stage",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
occupancy <- function(distribution = bernoulli) {
  
  stop("occupancy data are not supported (yet)", call. = FALSE)
  
  list(type = "population_occupancy",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
abundance <- function(distribution = poisson) {
  
  stop("population abundance data are not supported (yet)", call. = FALSE)
  
  list(type = "population_abundance",
       distribution = distribution)
  
}

#' @export
#' @rdname integrated_data
#' 
community <- function(distribution = bernoulli) {
  
  stop("community data are not supported (yet)", call. = FALSE)

  list(type = "community",
       distribution = distribution)

}

#' @export
#' @rdname integrated_data
#' 
is.integrated_data <- function(object) {
  inherits(object, "integrated_data")
}

#' @export
#' @rdname integrated_data
#' 
print.integrated_data <- function(x, ...) {
  cat(paste0("This is an integrated_data object\n"))
}

#' @export
#' @rdname integrated_data
#' 
summary.integrated_data <- function(object, ...) {
  NULL
}

# internal function: compile age-abundance data
age_abundance_internal <- function(data, process, settings) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = poisson)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  process_class <- ifelse(process$type == "leslie", "age", "stage")
  
  # how many classes are there in the process?
  classes <- c(process$classes, NA)
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  is_list <- is.list(data) & !is_matrix_like
  not_ok <- !is_matrix_like & !is_list
  
  if (not_ok)
    stop("age_abundance data must be a matrix, data.frame, or list", call. = FALSE) 
  
  # need to treat matrix-like data and list data differently
  if (is_matrix_like) {
    
    # default: return data as is if there's a column for each class
    data_clean <- data
    
    # is this a model with age data but a stage structure?
    if (process_class == "stage") {
      warning("it looks like you're fitting a stage-structured model to age data; is this correct?", call. = FALSE)
      classes[2] <- nrow(data_clean)
    }
    
    # what if there is not one column per class?
    if (nrow(data) != classes[1] & process_class != "stage") {
      
      # edge case: data are total abundances through time; need to use a different function
      if (nrow(data) == 1)
        stop("data only have one row; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      if (ncol(data) == 1)
        stop("data only have one column; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      
      # there must be something wrong
      stop(paste0("data should have one row per age class but have ", nrow(data), " rows and there are ", classes[1], " age classes"), call. = FALSE)
      
    }
    
  }
  
  if (is_list) {
    
    # how long is each element?
    list_len <- sapply(data, length)
    
    # if they're all the same, might be binned data provided as a list rather than matrix-like
    if (all(list_len == classes[1])) {
      
      data_clean <- do.call(cbind, data)
      
    } else {  # elements differ in length, probably have individual ages
      
      # cannot handle ages older than number of classes
      if (max(data) > classes[1]) {
        warning(paste0("some individuals are older than the total number of classes; ages will be truncated to a maximum of ", classes[1]), call. = FALSE)
        data <- ifelse(data > process$classes, classes[1], data)
      }
      
      # if so, we need breaks to bin the data        
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if data are not binned", call. = FALSE)
      
      # let the user know what's going on
      cat(paste0("converting data from list to matrix; this assumes data are individual ", class_type, "s and not counts\n"))
      
      # bin the data using hist_fn, applied to each element of the input list
      data_clean <- sapply(data, hist_fn, breaks = all_settings$breaks)
      
    }
    
  }
  
  # return outputs
  list(data = data_clean, classes = classes)
  
}

# internal function: compile stage-abundance data
stage_abundance_internal <- function(data, process, settings) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL,
                       likelihood = poisson)
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  process_class <- ifelse(process$type == "leslie", "age", "stage")
  
  # how many classes are there in the process?
  classes <- c(process$classes, NA)
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  is_list <- is.list(data) & !is_matrix_like
  not_ok <- !is_matrix_like & !is_list
  
  if (not_ok)
    stop("stage_abundance data must be a matrix, data.frame, or list", call. = FALSE) 
  
  # need to treat matrix-like data and list data differently
  if (is_matrix_like) {
    
    # default: return data as is if there's a column for each class
    data_clean <- data
    
    # is this a model with age data but a stage structure?
    if (process_class == "age") {
      warning("it looks like you're fitting an age-structured model to stage data; is this correct?", call. = FALSE)
      classes[2] <- nrow(data_clean)
    }
    
    # potentially reformat data if there is not one column per class
    if (nrow(data) != classes[1] & process_class != "age") {
      
      # edge case: data are total abundances through time; need to use a different function
      if (nrow(data) == 1)
        stop("data only have one row; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      if (ncol(data) == 1)
        stop("data only have one column; are you providing total abundance data? If so, try the `abundance()` function", call. = FALSE)
      
      stop(paste0("data should have one row per stage but have ", nrow(data), " rows and there are ", classes[1], " stages"), call. = FALSE)
      
    }
    
  }
  
  if (is_list) {
    
    # how long is each element?
    list_len <- sapply(data, length)
    
    # if they're all the same, might be binned data provided as a list rather than matrix-like
    if (all(list_len == classes[1])) {
      
      data_clean <- do.call(cbind, data)
      
    } else {  # elements differ in length, probably have individual ages
      
      # cannot handle ages older than number of classes
      if (max(data) > classes[1])
        stop(paste0("there are ", max(data), " stages and ", classes[1], " classes; the number of stages should not exceed the number of classes"), call. = FALSE)
      
      # if ok, we need breaks to bin the data        
      if (is.null(all_settings$breaks))
        stop("breaks must be provided if data are not binned", call. = FALSE)
      
      # let the user know what's going on
      cat(paste0("converting data from list to matrix; this assumes data are individual ", class_type, "s and not counts\n"))
      
      # bin the data using hist_fn, applied to each element of the input list
      data_clean <- sapply(data, hist_fn, breaks = all_settings$breaks)
      
    }
    
  }
  
  # return outputs
  list(data = data_clean, classes = classes)
  
}

# internal function: compile age-cjs data
age_cjs_internal <- function(data, process, settings) {
  
  # need to unpack the settings
  all_settings <- list()
  all_settings[names(settings)] <- settings
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  
  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("age_recapture data must be a matrix or data.frame containing capture histories", call. = FALSE) 
  
  # calculate number of classes
  classes <- c(process$classes, NA)
  
  # do the entries look like binned or binary data?
  is_binary <- all(unique(data) %in% c(0, 1))
  is_binned <- all(unique(data) %in% seq_len(classes[1]))
  
  # if data are not binned we need to bin them
  if (!is_binary & !is_binned)
    stop(paste0("range of data does not match the number of classes (", classes[1], ")"), call. = FALSE)
  
  # can only do so much with binary data
  if (is_binary) {
    
    # let the user know that we do not like binary data
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
  
  # return outputs
  list(data = data_clean, classes = classes)

} 

# internal function: compile stage-cjs data
stage_cjs_internal <- function(data, process, settings) {
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # warn that predictors are currently ignored for mark-recapture data
  if (!is.null(predictors))
    warning("predictors are not supported (yet) for stage_recapture data; provided predictors will be ignored", call. = FALSE)
  
  # is the process model a stage or age based model?
  class_type <- ifelse(process$type == "leslie", "age", "stage")
  
  # what format do the data take?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  
  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("stage_recapture data must be a matrix or data.frame containing capture histories", call. = FALSE) 
  
  # calculate classes
  classes <- c(process$classes, NA)
  
  # do the entries look like binned or binary data?
  is_binary <- all(unique(data) %in% c(0, 1))
  is_binned <- all(unique(data) %in% seq_len(classes[1]))
  
  # if data are not binned we need to bin them
  if (!is_binary & !is_binned) {
    
    # this will not work if we do not have breaks
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
    
    # let the user know that we do not like binary data
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

  # return outputs
  list(data = data_clean, classes = classes)
  
} 

# internal function: compile stage-to-age data
stage_to_age_internal <- function(data, process, settings) {
  
  # are data matrix_like?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  
  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("stage_to_age data must be a matrix or data.frame containing individual stages and ages or counts of stages and ages", call. = FALSE) 
  
  # send out to a function suited to the data type
  if (ncol(data) != nrow(data)) {
    out <- stage_to_age_internal_vectors
  } else {
    out <- stage_to_age_internal_array
  }
  
  # return outputs
  list(data = out$data_clean, classes = out$classes)
  
}

# internal function: compile stage-to-age data from vectors
stage_to_age_internal_vectors <- function(data, process, settings) {
  
  # this will not work if we have not got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to model ages from stage data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # parse formula to give outputs
  data_classes <- apply(data, 2, max)
  response <- data[, 1]
  predictor <- data[, 2]

  # check number of classes
  classes <- c(process$classes, data_classes[1])
  if (classes[1] < data_classes[2])
    stop(paste0("there are more age classes in the data than there are in the process model; consider specifying a new process model with ", data_classes[2], " age classes"), call. = FALSE)
  
  # basic checks to see if we have missed something
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
  predictor <- ifelse(predictor >= classes[1], classes[1], predictor)
  
  # need to turn vectors into a transition matrix
  data_clean <- matrix(0, nrow = classes[2], ncol = classes[1])
  for (i in seq_len(classes[2])) {
    predictor_sub <- predictor[response == i]
    classification <- tapply(ones_vec, predictor_sub, sum)
    data_clean[i, match(names(classification), seq_len(classes[1]))] <- classification
  }
  
  # return outputs
  list(data = data_clean, classes = classes)

}

# internal function: compile stage-to-age data from an array
stage_to_age_internal_array <- function(data, process, settings) {
  
  # this will not work if we have not got a Leslie matrix process
  if (process$type != "leslie")
    stop("trying to model ages from stage data without an age-based model; this seems like a bad idea", call. = FALSE)
  
  # unpack settings
  all_settings <- list()
  all_settings[names(settings)] <- settings
  
  # what dims does x have?
  dims <- dim(x)
  
  # if at least one of columns or rows do not line up with the number of classes, this will not work
  if (!(process$classes %in% dims))
    stop(paste0("one dimension of stage-to-age data must match the number of classes in process (", process$classes, ")"), call. = FALSE)
  
  # format so that ages are always in columns
  if (dims[2] != process$classes) {
    
    data_clean <- t(data)
    
  } else {
    
    # if ncol(x) == nrow(x), warn that we assume ages in columns
    if (dims[1] == dims[2])
      cat("stage-to-age data must have ages in columns; is this data set formatted correctly?", call. = FALSE)
    
    data_clean <- data
    
  }
  
  # set classes
  classes <- c(ncol(data_clean), nrow(data_clean))
  
  # return outputs
  list(data = data_clean, classes = classes)
  
}   

# internal function: compile age-to-stage data
age_to_stage_internal <- function(data, process, settings) {
  
  # are data matrix_like?
  is_matrix_like <- is.data.frame(data) | is.matrix(data)
  
  # can only handle matrix-like data types
  if (!is_matrix_like)
    stop("stage_to_age data must be a matrix or data.frame containing individual stages and ages or counts of stages and ages", call. = FALSE) 
  
  # send out to a function suited to the data type
  if (ncol(data) != nrow(data)) {
    out <- age_to_stage_internal_vectors
  } else {
    out <- age_to_stage_internal_array
  }
  
  # return outputs
  list(data = out$data_clean, classes = out$classes)
  
}

# internal function: compile age-to-stage data from vectors
age_to_stage_internal_vectors <- function(data, process, settings) {
  
  # this will not work if we do not have a stage-structured process model
  if (process$type == "leslie")
    stop("trying to model stages from age data without a stage-based model; this seems like a bad idea", call. = FALSE)
  
  # need to unpack the settings
  all_settings <- list(breaks = NULL)
  all_settings[names(settings)] <- settings
  
  # parse formula to give outputs
  data_classes <- apply(data, 2, max)
  response <- data[, 1]
  predictor <- data[, 2]
  
  # check number of classes
  classes <- c(process$classes, data_classes[1])
  if (classes[1] < data_classes[2])
    stop(paste0("there are more stages in the data than there are in the process model; consider specifying a new process model with ", data_classes[2], " stages"), call. = FALSE)
  
  # basic checks to see we missed something
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
  predictor <- ifelse(predictor >= classes[1], classes[1], predictor)
  
  # need to turn vectors into a transition matrix
  data_clean <- matrix(0, nrow = classes[2], ncol = classes[1[]])
  for (i in seq_len(classes[2])) {
    predictor_sub <- predictor[response == i]
    classification <- tapply(ones_vec, predictor_sub, sum)
    data_clean[i, match(names(classification), seq_len(classes[1]))] <- classification
  }
  
  # return outputs
  list(data = data_clean, classes = classes)
  
}

# internal function: compile age-to-stage data from an array
age_to_stage_internal_array <- function(data, process, settings) {
  
  # this will not work if we do not have a stage-structured process model
  if (process$type == "leslie")
    stop("trying to model stages from age data without a stage-based model; this seems like a bad idea", call. = FALSE)
  
  # unpack settings
  all_settings <- list(likelihood = multinomial)
  all_settings[names(settings)] <- settings
  
  # what dims does x have?
  dims <- dim(x)
  
  # if at least one of columns or rows do not line up with the number of classes, this will not work
  if (!(process$classes %in% dims))
    stop(paste0("one dimension of age-to-stage data must match the number of classes in process (", process$classes, ")"), call. = FALSE)
  
  # format so that ages are always in columns
  if (dims[2] != process$classes) {
    
    data_clean <- t(data)
    
  } else {
    
    # if ncol(x) == nrow(x), warn that we assume ages in columns
    if (dims[1] == dims[2])
      cat("age-to-stage data must have stages in columns; is this data set formatted correctly?", call. = FALSE)
    
    data_clean <- data
    
  }
  
  # set classes
  classes <- c(ncol(data_clean), nrow(data_clean))
  
  # return outputs
  list(data = data_clean, classes = classes)
  
}   

# internal function: create integrated_data object
as.integrated_data <- function(object) {
  as_class(object, name = "integrated_data", type = "list")
}
