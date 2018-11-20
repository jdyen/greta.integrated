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
plot.integrated_process <- function(x, y, ...) {
  
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("the DiagrammeR package must be installed to plot greta.integrated process modules",
         call. = FALSE)
  }
  
  # set up graph
  transition_mat <- t(x$masks$transition)
  fecundity_mat <- t(x$masks$fecundity)
  full_mat <- transition_mat + fecundity_mat

  # pull out process type
  type <- ifelse(x$type == "leslie", "Age", "Stage")
  
  gr <- DiagrammeR::from_adj_matrix(full_mat,
                                    mode = "directed",
                                    use_diag = TRUE)
  
  # how many nodes?
  n_nodes <- nrow(gr$nodes_df)
  
  # edge types
  to <- gr$edges_df$to
  from <- gr$edges_df$from
  
  # change type back to stage if there are multiple "age x+" terms
  if (sum(to == from) > 1)
    type <- "Stage"
  
  # identify different node types
  node_type <- rep("pre_reprod", n_nodes)
  node_type[from[to == 1 & from > 1]] <- "reprod"
  max_reprod <- max(which(node_type == "reprod"))
  node_type[seq_len(n_nodes) > max_reprod] <- "post_reprod"

  # change node types and colours based on node type
  node_shapes <- rep("square", n_nodes)
  node_shapes[node_type == "reprod"] <- "circle"
  node_shapes[node_type == "post_reprod"] <- "diamond"

  col_pal <- c("gray50", "#57A0D3", "#1D2951", "#0E4D92")
  col_pal_light <- ggplot2::alpha(col_pal, 0.5)
  node_edge_colours <- rep(col_pal[1], n_nodes)
  node_edge_colours[node_type == "reprod"] <- col_pal[4]
  node_edge_colours[node_type == "post_reprod"] <- col_pal[3]

  node_colours <- rep(col_pal_light[2], n_nodes)
  node_colours[node_type == "reprod"] <- col_pal_light[4]
  node_colours[node_type == "post_reprod"] <- col_pal_light[3]

  node_size <- rep(0.9, n_nodes)
  node_size[node_type == "reprod"] <- 0.9
  node_size[node_type == "post_reprod"] <- 1.0

  # add some labels for the nodes (age or stage depending on type of model)
  node_labels <- paste(type, seq_len(n_nodes), sep = " ")
  
  # if it's a Leslie matrix and to == from, we have an "age+" situation
  if (type == "Age")
    node_labels[from[to == from]] <- paste0(node_labels[from[to == from]], "+")
  
  edge_style <- rep("solid", length(to))
  
  # node options
  gr$nodes_df$type <- node_type
  gr$nodes_df$fontcolor <- col_pal[4]
  gr$nodes_df$fontsize <- 12
  gr$nodes_df$penwidth <- 2
  
  gr$nodes_df$shape <- node_shapes
  gr$nodes_df$color <- node_edge_colours
  gr$nodes_df$fillcolor <- node_colours
  gr$nodes_df$width <- node_size
  gr$nodes_df$height <- node_size * 0.8
  gr$nodes_df$label <- node_labels
  
  # edge options
  gr$edges_df$color <- "Gainsboro"
  gr$edges_df$fontname <- "Avenir"
  gr$edges_df$fontcolor <- "gray70"
  gr$edges_df$fontsize <- 14
  gr$edges_df$penwidth <- 4

  edge_types <- rep("transition", length(from))
  edge_types <- ifelse(from == to, "survival", edge_types)
  edge_types <- ifelse(from > to, "fecundity", edge_types)

  gr$edges_df$label <- edge_types
  gr$edges_df$style <- edge_style

  # set the layout type
  gr$global_attrs$value[gr$global_attrs$attr == "layout"] <- "dot"
  
  # make it horizontal
  gr$global_attrs <- rbind(gr$global_attrs,
                           data.frame(attr = "rankdir",
                                      value = "LR",
                                      attr_type = "graph"))
  
  grViz <- DiagrammeR::render_graph(gr)
  attr(grViz, "dgr_graph") <- gr
  grViz

}

# internal function: create integrated_process object
as.integrated_process <- function(object) {
  as_class(object, name = "integrated_process", type = "list")
}
