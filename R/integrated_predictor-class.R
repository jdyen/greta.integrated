#' @name integrated_predictor
#' @title create integrated predictor objects
#' 
#' @description \code{integrated_predictor} defines a data object containing predictor data
#'   that can be passed to an \link[greta.integrated]{integrated_data} constructor.
#' 
#' @param formula formula specifying the predictor structure (see \link[lme4]{lmer} for details)
#' @param data a single data input (see details for descriptions of possible input types)
#' @param ... additional arguments to \link[base]{print} and \link[base]{summary} methods (currently ignored)
#' @param x an \code{integrated_data} object
#' @param object an \code{integrated_data} object
#'
#' @details Add something.
#' 
#' @return An object of class \code{integrated_predictor}, which contains information on predictors and
#'   can be passed to an \link[greta.integrated]{integrated_data} constructor.
#' 
#' @export
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
#' @rdname integrated_predictor
#' 
add_predictors <- function(formula, data) {
  
  # try and find the data if not supplied
  if (missing(data)) {
    data_env <- environment(formula)
  } else {
    if (!is.data.frame(data))
      data <- as.data.frame(data)
    data_env <- as.environment(data)
  }
  
  # parse formula
  response <- all.vars(formula)[1] 
  terms <- terms(formula, data = data)
  random <- grep("\\|", attributes(terms)$term.labels)
  var_names <- all.vars(delete.response(terms))
  full_var_list <- colnames(attributes(terms)$factors)
  if (length(random)) {
    full_var_list_fixed <- full_var_list[-grep("\\|", full_var_list)]
  } else {
    full_var_list_fixed <- full_var_list
  }
  
  # use correct var_names when random is missing
  if (length(random)) {
    
    # check there are no interactions in the random terms
    if (length(grep('\\*', full_var_list[random])))
      stop('cannot include interactions in random effects; use separate (1 | random) terms for each random variable', call. = FALSE)
    
    # separate names of fixed and random vars
    fixed_vars <- var_names[-random]
    random_vars <- var_names[random]
    
  } else {
    
    # all vars are fixed
    fixed_vars <- var_names
    random_vars <- NULL
    
  }
  
  # create x and z objects to return
  if (length(fixed_vars)) {
    x_tmp <- mget(fixed_vars, envir = data_env, inherits = TRUE)
  }
  if (length(random_vars)) {
    z_tmp <- mget(random_vars, envir = data_env, inherits = TRUE)
    z_tmp <- lapply(z_tmp, function(x) as.integer(as.factor(x)))
  }
  
  # create model matrix of fixed variables
  if (length(fixed_vars)) {
    x <- model.matrix(as.formula(paste0("~", paste(full_var_list_fixed, collapse = " + "))), data = x_tmp)
  } else {
    x <- matrix(1, nrow = nrow(z_tmp))
  }
  
  # create model matrix of random variables
  if (length(random_vars)) {
    z <- model.matrix(as.formula(paste0(" ~ ", paste(random_vars, collapse = " + "))), data = z_tmp)
    z <- z[, -1]
  } else {
    z <- NULL
  }
  
  # are there are any missing data in the fixed or random effects?
  if (!is.null(x)) {
    na_rows <- apply(x, 1, function(x) any(is.na(x)))
  } else {
    na_rows <- rep(FALSE, nrow(z))
  }
  if (!is.null(z))
    na_rows <- na_rows | apply(z, 1, function(x) any(is.na(x)))
  
  # remove any rows with missing fixed or random effects
  if (any(na_rows)) {
    warning(paste0("there are ", sum(na_rows), " rows with missing data; these will be removed from all analyses"), call. = FALSE)
    if (!is.null(x))
      x <- x[!na_rows, ]
    if (!is.null(z))
      z <- z[!na_rows]
  }
  
  # collate outputs
  data_module <- list(fixed = x,
                      random = z,
                      removed = na_rows)
  
  # set predictor class and return outputs
  as.integrated_predictor(data_module)
  
}

#' @export
#' @rdname integrated_predictor
#' 
is.integrated_predictor <- function(object) {
  inherits(object, "integrated_predictor")
}

#' @export
#' @rdname integrated_predictor
#' 
print.integrated_predictor <- function(x, ...) {
  cat(paste0("This is an integrated_predictor object\n"))
}

#' @export
#' @rdname integrated_predictor
#' 
summary.integrated_predictor <- function(object, ...) {
  NULL
}

# internal function: create integrated_predictor object
as.integrated_predictor <- function(object) {
  as_class(object, name = "integrated_predictor", type = "list")
}
