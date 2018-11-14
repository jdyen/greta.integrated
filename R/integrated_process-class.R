#' @name integrated_process
#' @title create integrated process objects
#'
#' @description An \code{integrated_process} object contains the underlying
#'   process model for an integrated population analysis
#' 
#' @param classes something
#' @param density function of class \link[greta.integrated]{integrated_density}
#' @param priors named list of prior distributions (see details for information on setting prior distributions)
#' @param masks masking of matrices
#' @param ... additional arguments to \link[base]{print}, \link[base]{summary}, and \link[graphics]{plot} methods (currently ignored)
#' @param x an \code{integrated_process} object
#' @param object an \code{integrated_process} object
#'
#' @details something. Prior distributions can be specified as single-dimensional
#'   greta distribution, e.g., \code{normal(0, 1)}. Link functions and transformations
#'   can be specified directly in-line, e.g., \code{ilogit(normal(0, 1))} specifies
#'   normal priors with a mean of zero and a standard deviation of one, transformed
#'   with an inverse-logit link.
#'
#' @return An object of class \code{integrated_process}, which can be used to create
#'    \link[greta.integrated]{integrated_data} and \link[greta.integrated]{integrated_model} objects
#' 
#' @import greta
#' 
#' @examples
#' \dontrun {
#' 
#' library(integrated)
#' 
#' # a really basic age-structured model with five age classes
#' process <- leslie(5, density = ricker(lambda = uniform(0, 1)))
#' 
#' # setting custom priors
#' process <- leslie(5, density = bh(lambda = uniform(0, 1)),
#'                   priors = list(survival = ilogit(normal(0, 1)),
#'                                 fecundity = exp(normal(0, 1))))
#' }

#' @export
#' @rdname integrated_process
#' 
leslie <- function(classes, density = no_density(), priors = list()) {
  
  # set default masking for a leslie matrix
  masks <- list(survival = ifelse(row(diag(classes)) == classes & col(diag(classes)) == classes, 1, 0),
                    transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1, 1, 0),
                    fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  
  # create an unstructured matrix with these masks
  process <- unstructured(classes = classes, density = density, priors = priors, masks = masks)
  
  # set type appropriately
  process$type <- "leslie"
  
  # return outputs
  process

}

#' @export
#' @rdname integrated_process
#' 
lefkovitch <- function(classes, density = no_density(), priors = list()) {
  
  # set default masking
  masks <- list(survival = diag(classes),
                transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1, 1, 0),
                fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  
  # create an unstructured matrix with these masks
  process <- unstructured(classes = classes, density = density, priors = priors, masks = masks)
  
  # set type appropriately
  process$type <- "lefkovitch"

  # return outputs
  process
  
}

#' @export
#' @rdname integrated_process
#' 
age <- function(classes, density = no_density(), priors = list()) {
  
  leslie(classes, density, priors)
  
}

#' @export
#' @rdname integrated_process
#' 
stage <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  # set default masking
  mask_list <- list(survival = diag(classes),
                    transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1, 1, 0),
                    fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  
  # overwrite defaults with user-specified masks
  mask_list[names(masks)] <- masks
  
  # create an unstructured matrix with these masks
  process <- unstructured(classes = classes, density = density, priors = priors, masks = mask_list)
  
  # set type appropriately
  if (length(masks)) {
    process$type <- "stage"
  } else {
    process$type <- "lefkovitch"
  }
      
  # return outputs
  process
  
}

#' @export
#' @rdname integrated_process
#' 
unstructured <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  # set type
  type <- "unstructured"
  
  # initialise model priors
  prior_list <- list(survival = beta(1, 1),
                     transition = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # set default masking
  mask_list <- list(survival = diag(classes),
                    transition = matrix(1, nrow = classes, ncol = classes) - diag(classes),
                    fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  
  # overwrite defaults with user-specified masks
  mask_list[names(masks)] <- masks
  
  # what are the bounds on the priors?
  survival_bounds <- extract_bounds(prior_list$survival)
  transition_bounds <- extract_bounds(prior_list$transition)
  fecundity_bounds <- extract_bounds(prior_list$survival)
  
  # are these reasonable?
  if (survival_bounds[1] < 0 | survival_bounds[2] > 1)
    warning("the prior for survival has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (transition_bounds[1] < 0 | transition_bounds[2] > 1)
    warning("the prior for transition has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (fecundity_bounds[1] < 0)
    warning("the prior for fecundity has a lower bound less than 0; is this reasonable?", call. = FALSE)
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  priors = prior_list,
                  masks = mask_list)
  
  # return outputs
  as.integrated_process(process)
  
}

#' @export
#' @rdname integrated_process
#' 
ipm <- function(classes, density = no_density(), priors = list()) {
  
  # set type
  type <- "ipm"
  
  # will default parameters work for ipm setup?
  stop("ipm models are not yet implemented", call. = FALSE)
  
  # warn that classes is now computational not process-based
  if (classes < 50)
    warning(paste0("classes = ", classes, " seems quite low; an ipm model is approximating continuous classes, and typically classes > 50"), call. = FALSE)
  
  # initialise model priors
  prior_list <- list(survival = beta(1, 1),
                     transition = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  priors = prior_list)
  
  # return outputs
  as.integrated_process(process)
  
}

#' @export
#' @rdname integrated_process
#' 
occupancy <- function(classes, density = no_density(), priors = list()) {
  
  # set type
  type <- "occupancy"
  
  # initialise model priors
  prior_list <- list(survival = beta(1, 1),
                     transition = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  priors = prior_list)
  
  # return outputs
  as.integrated_process(process)
  
}

#' @export
#' @rdname integrated_process
#' 
is.integrated_process <- function(object) {
  inherits(object, "integrated_process")
}

#' @export
#' @rdname integrated_process
#' 
print.integrated_process <- function(x, ...) {
  cat(paste0("This is an integrated_process object\n"))
}

#' @export
#' @rdname integrated_process
#' 
summary.integrated_process <- function(object, ...) {
  
  NULL
  
}

#' @export
#' @rdname integrated_process
#' 
plot.integrated_process <- function(x, ...) {
  
  # make a nice plot of the matrix with colours for zero/nonzero cells
  NULL
  
}

# internal function: create integrated_process object
as.integrated_process <- function(object) {
  as_class(object, name = "integrated_process", type = "list")
}
