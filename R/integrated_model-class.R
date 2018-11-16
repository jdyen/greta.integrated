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
  
    # prepare matrix
    parameters[[process_id[i]]]$matrix <- construct_matrix(data_modules[[i]], parameters[[process_id[i]]])
    
    # we need to add a couple of things from the process module
    parameters[[process_id[i]]]$density <- process_list[[process_id[i]]]$density
    parameters[[process_id[i]]]$process_class <- process_list[[process_id[i]]]$type

    # make parameters cleaner if predictors are included (ops were kept separate for convenience)
    if (!is.null(parameters[[process_id[i]]]$n_predictors)) {
      parameters[[process_id[i]]]$survival <- parameters[[process_id[i]]]$survival$children[[1]]
      parameters[[process_id[i]]]$transition <- parameters[[process_id[i]]]$transition$children[[1]]
      parameters[[process_id[i]]]$fecundity <- parameters[[process_id[i]]]$fecundity$children[[1]]
    }
    
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
    loglik_fun(data_modules[[i]], parameters[[process_id[i]]])
    
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

#' @rdname integrated_model
#' @param y unused default argument
#' @param colour base colour used for plotting. Defaults to \code{greta} colours
#'   in violet.
#'
#' @details The plot method produces a visual representation of the defined
#'   model. It uses the \code{DiagrammeR} package, which must be installed
#'   first. Here's a key to the plots:
#'   \if{html}{\figure{plotlegend.png}{options: width="100\%"}}
#'   \if{latex}{\figure{plotlegend.pdf}{options: width=7cm}}
#'
#' @return \code{plot} - a \code{\link[DiagrammeR:grViz]{DiagrammeR::grViz}}
#'   object, with the
#'   \code{\link[DiagrammeR:create_graph]{DiagrammeR::dgr_graph}} object used to
#'   create it as an attribute \code{"dgr_graph"}.
#'
#' @export
plot.integrated_model <- function(x,
                                  y,
                                  colour = "#996bc7",
                                  ...) {
  
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("the DiagrammeR package must be installed to plot greta models",
         call. = FALSE)
  }
  
  # set up graph
  dag_mat <- x$dag$adjacency_matrix
  
  gr <- DiagrammeR::from_adj_matrix(dag_mat,
                                    mode = "directed",
                                    use_diag = FALSE)
  
  n_nodes <- nrow(gr$nodes_df)
  
  names <- names(x$dag$node_list)
  types <- x$dag$node_types
  to <- gr$edges_df$to
  from <- gr$edges_df$from
  
  node_shapes <- rep("square", n_nodes)
  node_shapes[types == "variable"] <- "circle"
  node_shapes[types == "distribution"] <- "diamond"
  node_shapes[types == "operation"] <- "circle"
  
  node_edge_colours <- rep(greta_col("lighter", colour), n_nodes)
  node_edge_colours[types == "distribution"] <- greta_col("light", colour)
  node_edge_colours[types == "operation"] <- "lightgray"
  
  node_colours <- rep(greta_col("super_light", colour), n_nodes)
  node_colours[types == "distribution"] <- greta_col("lighter", colour)
  node_colours[types == "operation"] <- "lightgray"
  node_colours[types == "data"] <- "white"
  
  node_size <- rep(1, length(types))
  node_size[types == "variable"] <- 0.6
  node_size[types == "data"] <- 0.5
  node_size[types == "operation"] <- 0.2
  
  # get node labels
  node_labels <- vapply(x$dag$node_list,
                        member,
                        "plotting_label()",
                        FUN.VALUE = "")
  
  # add greta array names where available
  visible_nodes <- lapply(x$visible_greta_arrays, get_node)
  known_nodes <- vapply(visible_nodes,
                        member,
                        "unique_name",
                        FUN.VALUE = "")
  known_nodes <- known_nodes[known_nodes %in% names]
  known_idx <- match(known_nodes, names)
  node_labels[known_idx] <- paste(names(known_nodes),
                                  node_labels[known_idx],
                                  sep = "\n")
  
  # for the operation nodes, add the operation to the edges
  op_idx <- which(types == "operation")
  op_names <- vapply(x$dag$node_list[op_idx],
                     member,
                     "operation_name",
                     FUN.VALUE = "")
  op_names <- gsub("`", "", op_names)
  
  ops <- rep("", length(types))
  ops[op_idx] <- op_names
  
  # get ops as tf operations
  edge_labels <- ops[to]
  
  # for distributions, put the parameter names on the edges
  distrib_to <- which(types == "distribution")
  
  parameter_list <- lapply(x$dag$node_list[distrib_to],
                           member,
                           "parameters")
  
  node_names <- lapply(parameter_list,
                       function(parameters) {
                         vapply(parameters,
                                member,
                                "unique_name",
                                FUN.VALUE = "")
                       })
  
  # for each distribution
  for (i in seq_along(node_names)) {
    
    from_idx <- match(node_names[[i]], names)
    to_idx <- match(names(node_names)[i], names)
    param_names <- names(node_names[[i]])
    
    # assign them
    for (j in seq_along(from_idx)) {
      idx <- from == from_idx[j] & to == to_idx
      edge_labels[idx] <- param_names[j]
    }
    
  }
  
  edge_style <- rep("solid", length(to))
  
  # put dashed line between target and distribution
  # for distributions, put the parameter names on the edges
  names <- names(x$dag$node_list)
  types <- x$dag$node_types
  distrib_idx <- which(types == "distribution")
  
  # find those with targets
  targets <- lapply(x$dag$node_list[distrib_idx],
                    member,
                    "target")
  
  keep <- !vapply(targets, is.null, TRUE)
  distrib_idx <- distrib_idx[keep]
  
  
  target_names <- vapply(x$dag$node_list[distrib_idx],
                         member,
                         "target$unique_name",
                         FUN.VALUE = "")
  distribution_names <- names(target_names)
  distribution_idx <- match(distribution_names, names)
  target_idx <- match(target_names, names)
  
  # for each distribution
  for (i in seq_along(distribution_idx)) {
    
    idx <- which(to == target_idx[i] & from == distribution_idx[i])
    edge_style[idx] <- "dashed"
    
  }
  
  # node options
  gr$nodes_df$type <- "lower"
  gr$nodes_df$fontcolor <- greta_col("dark", colour)
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
  gr$edges_df$fontname <- "Helvetica"
  gr$edges_df$fontcolor <- "gray"
  gr$edges_df$fontsize <- 11
  gr$edges_df$penwidth <- 3
  
  gr$edges_df$label <- edge_labels
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
