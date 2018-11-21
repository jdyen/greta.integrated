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
  data_names <- as.list(substitute(list(...)))[-1L]
  names(data_modules) <- data_names

  # is everything an integrated_data object?
  if(!all(sapply(data_modules, is.integrated_data)))
    stop("integrated_model only takes integrated_data objects as arguments", call. = FALSE)

  # which data types are we dealing with?
  data_types <- sapply(data_modules, function(x) x$likelihood$type)
  
  # are all of these implemented?
  implemented <- c("age_abundance", "stage_abundance", "age_cjs",
                   "stage_cjs", "stage_to_age", "age_to_stage")
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
      if (!all(includes_predictors)) {
        if (!(sapply(data_sub[!includes_predictors], function(x) x$likelihood$type) %in% c("stage_to_age", "age_to_stage")))
          stop("predictors must be included for all data modules that share a single process", call. = FALSE)
      }
      
      # how many predictors?
      n_fixed <- sapply(data_sub[includes_predictors], function(x) ncol(x$predictors$fixed))
      n_random <- sapply(data_sub[includes_predictors], function(x) ncol(x$predictors$random))
      
      # do these match in number? (assuming they match in predictor sets too; perhaps warn about this?)
      if (length(unique(n_fixed)) > 1 | length(unique(n_random)) > 1)
        stop("all data modules that share a single process must have the same number of predictors", call. = FALSE)
      
      n_fixed <- unique(n_fixed)
      n_random <- unique(n_random)
      
    }

    # pull out any process that matches hash i
    process_list[[i]] <- data_modules[[which(process_id == i)[1]]]$process
    
    # do we need to deal with age-stage conversions?
    idx <- sapply(data_sub, function(x) !(x$likelihood$type %in% c("stage_to_age", "age_to_stage")))
    classes_alt <- sapply(data_sub[idx], function(x) x$classes[2])
    
    # would like this not to be a list
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
                                           n_fixed = n_fixed,
                                           n_random = n_random)

    } else {
      
      # pull together a set of parameters based on the details of this process
      parameters[[i]] <- define_parameters(process = process_list[[i]],
                                           classes_alt = classes_alt)
      
    }
    
    # add the input data name
    parameters[[i]]$data_names <- names(data_sub)
    
  }

  for (i in seq_along(data_modules)) {
  
    # pull out a temporary subset process to work with
    parameters_tmp <- parameters[[process_id[i]]]
    
    # prepare matrix if not a conversion likelihood
    if (!(data_modules[[i]]$likelihood$type) %in% c("stage_to_age", "age_to_stage")) {
      
      parameters_tmp$matrix <- construct_matrix(data_modules[[i]], parameters[[process_id[i]]])
      
      # save the matrix with a unique name
      parameters[[process_id[i]]] <- c(parameters[[process_id[i]]], list(parameters_tmp$matrix))
      names(parameters[[process_id[i]]])[length(parameters[[process_id[i]]])] <- paste0("matrix", i)
      
    }
    
    # we need to add a couple of things from the process module
    parameters_tmp$density <- process_list[[process_id[i]]]$density
    parameters_tmp$process_class <- process_list[[process_id[i]]]$type

    # choose appropriate likelihood based on type of data
    loglik_fun <- switch(data_modules[[i]]$likelihood$type,
                         age_abundance = age_abundance_loglik,
                         stage_abundance = stage_abundance_loglik,
                         age_cjs = age_cjs_loglik,
                         stage_cjs = stage_cjs_loglik,
                         stage_to_age = stage_to_age_loglik,
                         age_to_stage = age_to_stage_loglik)
    
    # define likelihood (doesn't return anything)
    loglik_fun(data_modules[[i]], parameters_tmp)
    
  } 
  
  # make parameters cleaner if predictors are included (ops were kept separate for convenience)
  for (i in seq_len(n_process)) {
    if (!is.null(parameters[[i]]$n_fixed) | !is.null(parameters[[i]]$n_random)) {
      parameters[[i]]$survival <- parameters[[i]]$survival$children[[1]]
      parameters[[i]]$transition <- parameters[[i]]$transition$children[[1]]
      parameters[[i]]$fecundity <- parameters[[i]]$fecundity$children[[1]]
    }
    data_subset <- data_modules[process_hash == unique(process_hash)[i]]
    parameters[[i]]$data_types <- sapply(data_subset, function(x) x$likelihood$type)
    parameters[[i]]$data_details <- sapply(data_subset, function(x) dim(x$data))
    parameters[[i]]$process <- process_list[[i]]$type
    parameters[[i]]$likelihood <- sapply(data_subset, function(x) x$likelihood$distribution)
    parameters[[i]]$density <- process_list[[i]]$density
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
#' @export
#' 
plot.integrated_model <- function(x, y, ...) {
  
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("the DiagrammeR package must be installed to plot greta.integrated process modules",
         call. = FALSE)
  }

  # how many processes are we dealing with?
  n_process <- length(x)
  
  # set up a separate graph for each process
  for (i in seq_len(n_process)) {
    
    # subset to a single process
    x_sub <- x[[i]]
    
    # pass to internal plotting function
    gr_tmp <- plot_integrated_internal(x_sub)
    
    if (i == 1) {
      gr_ndfs <- gr_tmp$nodes_df
      gr_edfs <- gr_tmp$edges_df
    } else {
      gr_edges <- gr_tmp$edges_df
      gr_edges$from <- gr_edges$from + max(gr_ndfs$id)
      gr_edges$to <- gr_edges$to + max(gr_ndfs$id)
      gr_edfs <- DiagrammeR::combine_ndfs(gr_edfs, gr_edges)
      gr_ndfs <- DiagrammeR::combine_ndfs(gr_ndfs, gr_tmp$nodes_df)
    }
    
  }
  
  # create graph from nodes and edges
  gr_all <- DiagrammeR::create_graph(gr_ndfs, gr_edfs)
  
  # set the layout type
  gr_all$global_attrs$value[gr_all$global_attrs$attr == "layout"] <- "dot"
  
  # make it horizontal
  gr_all$global_attrs <- rbind(gr_all$global_attrs,
                               data.frame(attr = "rankdir",
                                          value = "LR",
                                          attr_type = "graph")) 
  
  # render the graph to the RStudio viewer
  grViz <- DiagrammeR::render_graph(gr_all)
  attr(grViz, "dgr_graph") <- gr_all
  
  # return the DiagrammeR graph defn
  grViz
  
} 

# internal function to set up diagrammer plots
plot_integrated_internal <- function(x) {

  # what are its key features?
  process_type <- x$process
  type <- x$data_types
  
  # are we dealing with any age-to-stage conversions?
  conversion <- type %in% c("age_to_stage", "stage_to_age")
  
  # are there predictors?
  n_rand <- unlist(x$n_random)
  n_fix <- unlist(x$n_fixed)
  includes_predictors <- !is.null(n_fix) | !is.null(n_rand)
  n_rand <- ifelse(is.null(n_rand), 0, n_rand)
  n_fix <- ifelse(is.null(n_fix), 0, n_fix)
  
  # is there density dependence?
  includes_density <- x$density$name != "none"
  
  # separate conversions and non-conversions
  true_data <- !conversion
  n_data <- length(type)
  
  # add rows for conversions
  if (process_type == "leslie") {
    converted <- grep("stage", type)
  } else {
    converted <- grep("^age", type)
  }
  
  # create a suitable matrix
  ndata_plus1 <- n_data + 1
  
  full_mat <- matrix(0, ndata_plus1, ndata_plus1)
  
  # setup a matrix with links to processes only where not converted
  if (any(conversion)) {
    
    is_conversion <- which(conversion) + 1
    
    full_mat[1, which(!(seq_len(n_data) %in% converted) | conversion) + 1] <- 1
    full_mat[is_conversion, which(true_data) + 1] <- 1
    full_mat <- rbind(full_mat, rep(0, ncol(full_mat)))
    full_mat <- cbind(full_mat, rep(0, nrow(full_mat)))
    full_mat[is_conversion, ncol(full_mat)] <- 1
    full_mat[nrow(full_mat), which(true_data) + 1] <- 1
    
    node_type <- c("process", rep("true_data", ncol(full_mat) - 1))
    node_type[is_conversion] <- "converter"
    node_type[length(node_type)] <- "conversion"

  } else {
    
    full_mat[1, 2:ndata_plus1] <- 1
    
    node_type <- c("process", rep("true_data", ncol(full_mat) - 1))
    
  }

  
  if (includes_predictors) {
    
    full_mat <- cbind(full_mat, rep(0, nrow(full_mat)))
    full_mat <- rbind(full_mat, c(1, rep(0, nrow(full_mat))))
    node_type <- c(node_type, "predictor")
    
  }
  
  if (includes_density) {
    
    full_mat <- cbind(full_mat, rep(0, nrow(full_mat)))
    full_mat <- rbind(full_mat, c(1, rep(0, nrow(full_mat))))
    node_type <- c(node_type, "density")

  }
  
  # build graph from adjacency matrix
  gr <- DiagrammeR::from_adj_matrix(full_mat,
                                    mode = "directed",
                                    use_diag = TRUE)
  
  # how many nodes?
  n_nodes <- nrow(gr$nodes_df)
  
  # edge types
  to <- gr$edges_df$to
  from <- gr$edges_df$from

  # change node types and colours based on node type
  node_shapes <- rep("rectangle", n_nodes)
  node_shapes[node_type == "process"] <- "circle"
  node_shapes[node_type == "converter"] <- "diamond"
  node_shapes[node_type == "density"] <- "diamond"
  
  col_pal <- c("#57A0D3", "#4F7942", "#0E4D92", "#964000")
  col_pal_light <- paste0(col_pal, "50")
  node_edge_colours <- rep(col_pal[3], n_nodes)
  node_edge_colours[node_type == "process"] <- col_pal[2]
  node_edge_colours[node_type == "converter"] <- col_pal[4]
  node_edge_colours[node_type == "density"] <- col_pal[4]
  
  node_colours <- rep(col_pal_light[3], n_nodes)
  node_colours[node_type == "process"] <- col_pal_light[2]
  node_colours[node_type == "converter"] <- col_pal_light[4]
  node_colours[node_type == "density"] <- col_pal_light[4]
  
  node_size <- rep(0.9, n_nodes)
  node_size[node_type == "converter"] <- 1.0
  node_size[node_type == "density"] <- 1.0
  
  node_labels <- rep("Node", ndata_plus1)
  node_labels[node_type == "process"] <- paste0(paste0(toupper(substr(process_type, 1, 1)), substr(process_type, 2, nchar(process_type))), "\nmatrix")
  node_labels[node_type == "conversion"] <- ifelse(process_type == "leslie", "Stage-to-age", "Age-to-stage")
  node_labels[node_type == "conversion"] <- paste(node_labels[node_type == "conversion"], paste0("(", names(type)[conversion], ")"), sep = "\n")
  node_labels[node_type == "converter"] <- paste(ifelse(process_type == "leslie", "Age", "Stage"), " conversion")
  node_labels[node_type == "predictor"] <- paste("Predictor variables", paste0("(", n_fix, " fixed, ", n_rand, " random", ")"), sep = "\n")
  node_labels[node_type == "density"] <- paste0("Density\n(", x$density$name, ")")
  
  data_type <- rep("Data", n_data)
  data_type[grep("abundance", type)] <- "Abundance"
  data_type[grep("cjs", type)] <- "Recapture"
  data_type <- paste(type, paste0("(", names(type), ")"), sep = "\n")
  
  node_labels[node_type == "true_data"] <- data_type[!conversion]
  
  edge_style <- rep("solid", length(to))
  edge_style[to == which(node_type == "process")] <- "dashed"

  font_colour <- rep(col_pal[3], n_nodes)
  font_colour[node_type == "process"] <- col_pal[2]
  font_colour[node_type == "converter"] <- col_pal[4]
  font_colour[node_type == "density"] <- col_pal[4]
  
  # node options
  gr$nodes_df$type <- node_type
  gr$nodes_df$fontname <- "Avenir"
  gr$nodes_df$fontcolor <- font_colour
  gr$nodes_df$fontsize <- 12
  gr$nodes_df$penwidth <- 2
  
  gr$nodes_df$shape <- node_shapes
  gr$nodes_df$color <- node_edge_colours
  gr$nodes_df$fillcolor <- node_colours
  gr$nodes_df$width <- node_size
  gr$nodes_df$width[node_type %in% c("true_data", "conversion", "predictor")] <- 1.5
  gr$nodes_df$width[node_type %in% c("converter")] <- 1.6
  gr$nodes_df$height[node_type %in% c("converter")] <- 1
  gr$nodes_df$width[node_type %in% c("density")] <- 1.6
  gr$nodes_df$height[node_type %in% c("density")] <- 1
  gr$nodes_df$height <- node_size * 0.8
  gr$nodes_df$label <- node_labels
  
  # edge options
  gr$edges_df$color <- "Gainsboro"
  gr$edges_df$fontname <- "Avenir"
  gr$edges_df$fontcolor <- "LightGray"
  gr$edges_df$fontsize <- 14
  
  edge_labels <- rep("", length(from))
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, poisson)), "Poisson", "Distribution")
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, negative_binomial)), "Negative binomial", ll_types)
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, binomial)), "Binomial", ll_types)
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, bernoulli)), "Bernoulli", ll_types)
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, multinomial)), "Multinomial", ll_types)
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, normal)), "Normal", ll_types)
  ll_types <- ifelse(sapply(x$likelihood, function(x) identical(x, lognormal)), "Lognormal", ll_types)

  if (any(conversion)) {
    edge_labels[from == which(node_type == "converter" & node_type != "predictor")] <- c(ll_types[!conversion], ll_types[conversion])
  } else {
    edge_labels[from == which(node_type == "process")] <- ll_types
  }
  
  gr$edges_df$label <- edge_labels
  gr$edges_df$style <- edge_style
  
  edge_weights <- rep(4, length(to))
  edge_weights[from == which(node_type == "conversion")] <- 0
  gr$edges_df$penwidth <- edge_weights
  
  gr
  
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
