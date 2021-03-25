
setOldClass("igraph")

#' @title Graph matching results class
#'
#' @description An 'S4' class for the results of a graph matching function
#'
#' @slot corr data.frame indicating the correspondence between two graphs
#' @slot  dimensions of the original two graphs
#' @slot call The call to the graph matching function
#' @slot detail List with other more detailed information
#' 
#'
#' @seealso graphMatch_methods
#'
#' @rdname gm_constructor
setClass("graphMatch",
  slots = c(
    corr = "data.frame",
    nnodes = "integer",
    call = "call"
  ),
  contains = 'list')

#' @rdname gm_constructor
#'
#' @param corr data.frame indicating the correspondence between two graphs
#' @param nnodes dimensions of the original two graphs
#' @param call The call to the graph matching function
#' @param detail List with other more detailed information
#' 
#' @return graphMatch object
#'
#' @export
setGeneric(
  name = "graphMatch",
  def = function(
    corr,
    nnodes,
    call = NULL,
    detail = list()
  ) {

    # check args ....
    gm <- new(
      "graphMatch",
      corr = corr,
      nnodes = nnodes
    )
    gm@call <- call
    for (n in names(detail)) {
      gm[[n]] <- detail[[n]]
    }
    # gm@.Data <- detail
    gm
  }
)


get_perm <- function(match){
  make_perm(match@nnodes[1], match@nnodes[2], match@corr)
}


as.character.graphMatch <- function(from) {
  paste0(
    "Call:", as.character(from@call), "\n",
    "Match (", 
    paste(as.character(from@nnodes), collapse = " x "),
    "):\n",
    as.character(from@corr)
  )
}



#' @name as
#' @title Methods for the graphMatch class
#' 
#' @description These methods provide functionality to plot,
#'  view, inspect, and convert graphMatch objects.
#' 
#' 
#' @rdname graphMatch_methods
setAs("graphMatch", "character", as.character.graphMatch)


#' @name as
#' @rdname graphMatch_methods
setAs("graphMatch", "Matrix", function(from) {
  from[]
})


#' @name as
#' @rdname graphMatch_methods
setAs("graphMatch", "data.frame", function(from) {
  from@corr
})


#' @rdname graphMatch_methods
#' 
#' @param x graphMatch object
setMethod("as.data.frame", signature("graphMatch"),
  function(x) {
    x@corr
  }
)

show.graphMatch <- function(object){
    print(object@call)
    cat(paste0("\nMatch (", 
      paste(as.character(object@nnodes), collapse = " x "),
      "):\n"
    ))
    show(object@corr)
}



######################################################################

#' @rdname graphMatch_methods
#' @param object graphMatch object
setMethod("show", signature("graphMatch"), show.graphMatch)


#' @rdname graphMatch_methods
setMethod("print", signature("graphMatch"), 
  function(x) show.graphMatch(x))


.leftmult = function(x, y){
  #x is graphMatch, y is a matrix
  
}


#' @rdname graphMatch_methods
#' 
#' @param i row index for the correspondence data.frame
#' @param j col index for the correspondence data.frame
#' @param drop ignored
setMethod("[",
  signature(x = "graphMatch",
    i = 'missing', j = 'missing', drop = 'missing') ,
  function(x, i = NULL, j = NULL, drop = NULL) {
          get_perm(x)
  })

#' @rdname graphMatch_operators
#' 
#' @title Operator methods for graphMatch objects
#' 
#' @description Methods to use graphMatch objects as operators on 
#'  igraph and matrix-like objects.
#' 
#' @param x Either graphMatch object or a matrix-like object
#' @param y Either graphMatch object or a matrix-like object
setMethod("%*%", signature(x = "graphMatch", y = "ANY"), function(x, y) {
  x[] %*% y
})


#' @rdname graphMatch_operators
setMethod("%*%", signature(x = "ANY", y = "graphMatch"),
  function(x, y){
    x %*% y[]
  })

#' @rdname graphMatch_methods
setMethod("dim", signature(x = "graphMatch"), function(x) { x@nnodes })

#' @rdname graphMatch_methods
setMethod("length", signature(x = "graphMatch"), function(x) { nrow(x@corr)})


#' @rdname graphMatch_methods
setMethod("[",
  signature(x = "graphMatch",
    i = 'ANY', j = 'ANY', drop = 'ANY') ,
  function(x, i = NULL, j = NULL, drop = NULL) {
          x@corr[i, j]
  })


#' @rdname graphMatch_methods
setMethod('str', signature(object = "graphMatch"), function(object){
  as.character(object)
})



# @todo function to get seeds, nonseeds, match/w seeds match w/o seeds


identity_match <- function(x, y) {
  if(is.igraph(x)) {
    nmin <- min(igraph::vcount(x), igraph::vcount(y))
  } else {
    nmin <- min(nrow(x), nrow(y))
  }
  
  graphMatch(
    corr = data.frame(corr_A = 1:nmin, corr_B = 1:nmin),
    nnodes = c(nmin, nmin),
    call = match.call()
  )
}


#' @title Plotting methods for visualizing matches
#'
#' @description Two functions are provided, \code{match_plot_igraph}
#' which makes a ball and stick plot from 'igraph' objects
#' and \code{match_plot_matrix} which shows an adjacency
#' matrix plot.
#'
#'
#' @param x First graph, either an igraph object or a Matrix
#' @param y second graph, either an igraph object or a Matrix
#' @param match result from a match call. Requires element
#'  \code{corr} as a data.frame with names corr_A, corr_B.
#' @param color Whether to color edges according to which
#'  graph(s) they are in.
#' @param linetype Whether to set edge linetypes according
#'  to which graph(s) they are in.
#' @param ... additional parameters passed to either the
#'  'igraph' plot function or the Matrix image function.
#'
#' @returns Both functions return values invisibly.
#' \code{match_plot_igraph} returns the union of the
#'  matched graphs as an 'igraph' object with additional
#'  edge attributes \code{edge_match, color, lty}.
#'  \code{match_plot_matrix} returns the difference between
#'  the matched graphs.
#'
#' @details
#' Grey edges/pixels indicate common edges, red
#' indicates edges only in graph A and green
#' represents edges only graph B. The corresponding
#' linetypes are solid, short dash, and long dash.
#'
#' The plots can be recreated from the output with the code \cr
#' \code{plot(g)} \cr
#' for \code{g <- match_plot_igraph(...)} and  \cr
#' \code{col <- colorRampPalette(c("#AA4444", "#888888", "#44AA44"))} \cr
#' \code{image(m, col.regions = col(256))} \cr
#' for \code{m <- match_plot_match(...)}.
#'
#' This only plots and returns the matched vertices.
#'
#' @rdname plot_graphMatch
#' 
#' 
#'
#' @examples
#' set.seed(123)
#' graphs <- sample_correlated_gnp_pair(20, .9, .3)
#' A <- graphs$graph1
#' B <- graphs$graph2
#' res <- graph_match_percolation(A, B, 1:4)
#'
#' plot(A, B, res)
#' plot(A[], B[], res)
setMethod("plot", signature(x = "igraph", y = "igraph"),
  function(x,y, match = NULL, 
    color = TRUE, linetype = TRUE,...) {
    if(is.null(match)) {
      match <- identity_match(x,y)
    }
    match_plot_igraph(x, y, match, 
      color, linetype, ...)

  })

#' @rdname plot_graphMatch
#' @param col.regions NULL for default colors, otherwise see \link[Matrix]{image-methods}
#' @param at NULL for default at values for at (ensures zero is grey), otherwise see \link[Matrix]{image-methods}
#' @param colorkey NULL for default colorkey, otherwise see \link[Matrix]{image-methods}
#'
#' @export
setMethod("plot", signature(x = "Matrix", y = "Matrix"),
  function(x,y, match = NULL,
    col.regions = NULL, at = NULL, colorkey = NULL, ...) {

    if(is.null(match)) {  
      match <- identity_match(x,y)
    }
    match_plot_matrix(x, y, match, col.regions, at, colorkey, ...)

  })

permuted_subgraph <- function(g, corr_g) {
  if(is.null(igraph::V(A)$name)){
    igraph::set.vertex.attribute(A, "name", corr_g, corr_g)
  }

  A <- igraph::permute.vertices(
    igraph::induced_subgraph(A, corr_g),
    rank(corr_g)
  )
}


#' @title Summary methods for graphMatch objects
#' 
#' @param object graphMatch object
#' @param A igraph or matrix-like object
#' @param B igraph or matrix-like object
#' @param true_label the true correspondence (if available)
#' @param directed whether to treat the graphs 
#'  as directed (TRUE) or not directed (FALSE) default is NULL
#'  which will treat the graphs as directed if either adjacency
#'  matrix is not symetric.
#' @rdname graphMatch_summary
setMethod("summary", signature("graphMatch"),
  function(object, A = NULL, B = NULL, true_label = NULL, directed = NULL) {

    # Matched nodes
    corr <- object@corr
    object$n_match <- nrow(corr) - nrow(object$seeds)
    if(!is.null(true_label)){
      object$n_true_match <-
        sum(true_label[corr$corr_A] == corr$corr_B) - nrow(object$seeds)

    }
    if(!is.null(A) & !is.null(B)){
      graph_pair <- check_graph(A, B)
      A <- graph_pair[[1]]
      B <- graph_pair[[2]]


      # Matched edges
      object$edge_match_info <-
        edge_match_info(corr, A, B, directed)
    }


    # objective value: ||A-PBP^T||_F
    # object@.Data$Permutation <- get_perm(nrow(A[[1]]), nrow(B[[1]]), corr)
    class(object) <- c("summary.graphMatch")
    object
  }
)

show.summary.graphMatch <- function(match) {
    cat("Call: \n")
    print(match@call)
    cat("\n# Matches:", match$n_match)
    if(!is.null(match$n_true_match)){
      cat("\n# True Matches: ", match$n_true_match)
    }
    if(!is.null(match$seeds)) { 
      cat(", # Seeds: ", nrow(match$seeds0))
    }
    cat("\n")
    if(!is.null(match$edge_match_info)) {
      ep <- as.data.frame(t(match$edge_match_info))
      colnames(ep) <- NULL
      print(ep)
    }
}

setClass("summary.graphMatch",
  # slots = c(
  #   corr = "data.frame",
  #    = "integer",
  #   call = "call",
  #   detail = "list"
  # ),
  contains = 'graphMatch'
)


setMethod("show", signature("summary.graphMatch"),
  function(object){
    show.summary.graphMatch(object)
  }
)

