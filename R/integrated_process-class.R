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
#' \dontrun{
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
leslie <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  # set default masking for a leslie matrix
  mask_list <- list(transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1 | row(diag(classes)) == classes & col(diag(classes)) == classes, 1, 0),
                    fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  mask_list[names(masks)] <- masks
  
  # create an unstructured matrix with these masks
  process <- unstructured(classes = classes, density = density, priors = priors, masks = mask_list)
  
  # set type appropriately
  process$type <- "leslie"
  
  # return outputs
  process

}

#' @export
#' @rdname integrated_process
#' 
lefkovitch <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  # set default masking
  mask_list <- list(transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1 | diag(classes), 1, 0),
                    fecundity = ifelse(row(diag(classes)) == 1 & col(diag(classes)) != 1, 1, 0))
  mask_list[names(masks)] <- masks
  
  # create an unstructured matrix with these masks
  process <- unstructured(classes = classes, density = density, priors = priors, masks = mask_list)
  
  # set type appropriately
  process$type <- "lefkovitch"

  # return outputs
  process
  
}

#' @export
#' @rdname integrated_process
#' 
age <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  leslie(classes, density, priors, masks)
  
}

#' @export
#' @rdname integrated_process
#' 
stage <- function(classes, density = no_density(), priors = list(), masks = list()) {
  
  # set default masking
  mask_list <- list(transition = ifelse(row(diag(classes)) == col(diag(classes)) + 1 | diag(classes), 1, 0),
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
                     noise = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # set default masking
  mask_list <- list(transition = matrix(1, nrow = classes, ncol = classes),
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
  
  # add a process hash
  process$hash <- paste(sample(c(LETTERS, letters, seq_len(10) - 1), size = 20, replace = TRUE), collapse = "")
  
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
                     noise = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  priors = prior_list)
  
  # add a process hash
  process$hash <- paste(sample(c(LETTERS, letters, seq_len(10) - 1), size = 20, replace = TRUE), collapse = "")
  
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
                     noise = normal(0, 10, truncation = c(0, Inf)))
  
  # overwrite defaults with user-specified priors
  prior_list[names(priors)] <- priors
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  priors = prior_list)
  
  # add a process hash
  process$hash <- paste(sample(c(LETTERS, letters, seq_len(10) - 1), size = 20, replace = TRUE), collapse = "")
  
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
  
  ### REDO with DiagrammeR to create a stage-transition diagram
  
  
  # what are we working with?
  tran <- x$masks$transition
  fec <- x$masks$fecundity

  # do any overlap?
  all <- tran + fec
  overlap <- any(all > 1)
  
  # define a matrix of colours
  col_mat <- ifelse(tran, 1, 0)
  col_mat <- ifelse(fec, 2, col_mat)
  
  # set a colour palette
  col_pal <- c("gray95", "#57A0D3", "#1D2951", "#0E4D92")
  
  # easy if they don't overlap
  if (!overlap) {

    # pull out old settings to restore
    old_mar <- par()$mar

    # remove ridiculously wide borders
    par(mar = c(5, 3, 0.7, 0.7))
    
    # plot with colours according to col_mat
    image(t(col_mat),
          col = col_pal[1:3],
          xaxt = "n", yaxt = "n",
          bty = "n")
    
    # add a legend
    legend(0.5, -0.25, xpd = TRUE,
           fill = col_pal,
           horiz = TRUE,
           border = c("gray70", col_pal[2:3]),
           bty = "n",
           cex = 1.25,
           xjust = 0.5,
           legend = c("none", "transition", "fecundity"))
    
    # add some labels
    type <- ifelse(x$type == "leslie", "Age", "Stage")
    mtext(paste0(type, " at time t"), side = 1, adj = 0.5, line = 1, cex = 1.5)
    mtext(paste0(type, " at time t + 1"), side = 2, adj = 0.5, line = 1.5, cex = 1.5)
    
    # mark low and high ends
    mtext("Low", side = 1, adj = 0.01, line = 0.2, cex = 1)
    mtext("High", side = 1, adj = 0.99, line = 0.2, cex = 1)
    mtext("Low", side = 2, adj = 0.01, line = 0.2, cex = 1)
    mtext("High", side = 2, adj = 0.99, line = 0.2, cex = 1)
    
    # reset plot settings
    par(mar = old_mar)
    
  } else {
    
    # add some new colours if they overlap
    col_mat[tran + fec == 2] <- 3

    # pull out old settings to restore
    old_mar <- par()$mar
    
    # remove ridiculously wide borders
    par(mar = c(5, 3, 0.7, 0.7))
    
    # plot with colours according to col_mat
    image(t(col_mat),
          col = col_pal,
          xaxt = "n", yaxt = "n",
          bty = "n")
    
    # add a legend
    legend(0.5, -0.25, xpd = TRUE,
           fill = col_pal,
           horiz = TRUE,
           border = c("gray70", col_pal[2:4]),
           bty = "n",
           cex = 1,
           xjust = 0.5,
           legend = c("none", "transition", "fecundity", "both"))
    
    # add some labels
    type <- ifelse(x$type == "leslie", "Age", "Stage")
    mtext(paste0(type, " at time t"), side = 1, adj = 0.5, line = 1, cex = 1.5)
    mtext(paste0(type, " at time t + 1"), side = 2, adj = 0.5, line = 1.5, cex = 1.5)
    
    # mark low and high ends
    mtext("Low", side = 1, adj = 0.01, line = 0.2, cex = 1)
    mtext("High", side = 1, adj = 0.99, line = 0.2, cex = 1)
    mtext("Low", side = 2, adj = 0.01, line = 0.2, cex = 1)
    mtext("High", side = 2, adj = 0.99, line = 0.2, cex = 1)
    
    # reset plot settings
    par(mar = old_mar)
    
  }

}

# internal function: create integrated_process object
as.integrated_process <- function(object) {
  as_class(object, name = "integrated_process", type = "list")
}
