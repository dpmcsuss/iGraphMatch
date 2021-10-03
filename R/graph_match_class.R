
setOldClass("igraph")

#'
#' @title Graph matching results class
#' @description An 'S4' class for the results of a graph matching function
#'
#' @slot corr data.frame indicating the correspondence between two graphs
#' @slot nnodes of the original two graphs
#' @slot call The call to the graph matching function
#'
#' @details graphMatch objects are returned by any of the graph matching methods
#'   implemented in the iGraphMatch package. These objects are primarily to
#'   represent the found correspondence between the two vertex sets. This is
#'   represented by a data.frame with two columns indicating the aligned
#'   vertex-pairs across the two graphs.
#'
#' @seealso \link[iGraphMatch:as.data.frame,graphMatch-method]{graphMatch_methods},
#' \link[iGraphMatch:summary,graphMatch-method]{graphMatch_summary},
#' \link[iGraphMatch:\%*\%,graphMatch,ANY-method]{graphMatch_operators},
#' \link[iGraphMatch:plot,igraph,igraph-method]{graphMatch_plot}
#'
#' @rdname graphMatch_constructor
setClass("graphMatch",
  slots = c(
    corr = "data.frame",
    nnodes = "integer",
    call = "call"
  ),
  contains = 'list'
)

#' @rdname graphMatch_constructor
#'
#' @param corr data.frame indicating the correspondence between two graphs
#' @param nnodes dimensions of the original two graphs
#' @param call The call to the graph matching function
#' @param detail List with other more detailed information
#'
#' @return graphMatch object
#'
#' @examples
#' # sample a pair of correlated random graphs from G(n,p)
#' set.seed(123)
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#'
#' # match g1 & g2 using percolation algorithm with some known node pairs as seeds
#' match <- gm(A = g1, B = g2, seeds = 1:3, method = 'indefinite')
#'
#' # graphMatch object
#' match
#' match@call # graph matching function
#' match@nnodes # sizes of two graphs
#' match@corr # matching correspondence
#'
#' match$corr_A # matching correspondence in the first graph
#' match$corr_B # matching correspondence in the second graph
#' match$seeds # vector of logicals indicating seeded nodes
#'
#' # matching details unique to the FW methodology with indefinite relaxation
#' match$iter # number of iterations
#' match$soft # doubly stochastic matrix from the last iteration, can be used to extract soft matching
#' match$lap_method # method for solving lap
#'
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

    if (!is.data.frame(corr)) {
      stop("Correspondence corr must be stored as a data.frame")
    }
    gm <- new(
      "graphMatch",
      corr = corr,
      nnodes = nnodes
    )
    if (is.null(call)) {
      #browser()
      call <- call("no call")
    }
    gm@call <- call
    for (n in names(detail)) {
      gm[[n]] <- detail[[n]]
    }
    if ("seeds" %in% names(detail)) {
      gm$seeds <- check_seeds(
        gm$seeds,
        max(nnodes),
        logical = TRUE
      )[corr$corr_A]
    }
    # gm@.Data <- detail
    gm
  }
)

#' @method as.character graphMatch
#' @export
as.character.graphMatch <- function(x, ...) {
  paste0(
    "Call:", as.character(x@call), "\n",
    "Match (",
    paste(as.character(x@nnodes), collapse = " x "),
    "):\n",
    as.character(x@corr)
  )
}



#'
# #' @rdname graphMatch_methods
setAs("graphMatch", "character",
  function(from)  as.character.graphMatch(from)
  )


# #' @rdname graphMatch_methods
setAs("graphMatch", "Matrix", function(from) {
  from[]
})


# #' @rdname graphMatch_methods
setAs("graphMatch", "data.frame", function(from) {
  from@corr
})



#' @title Methods for the graphMatch class
#'
#' @description These methods provide functionality to plot, view, inspect, and
#'   convert graphMatch objects.
#'
#' @details Methods for the graphmatch
#'
#' @return \code{dim} returns a vector of length two indicating the number of
#'   vertices in each original graph. \code{length} returns the number of found
#'   vertex-pair matches. \code{m[i,j]} will index the 2xlength data.frame of
#'   vertex-pair matches. This is true for any i,j unless both are missing. In
#'   that case, \code{m[]} returns a sparse matrix of dimension dim(m) where
#'   \code{m[][i,j]} is 0 unless m matches node i with node j. (Note this is not
#'   guaranteed to be a permutation matrix unless \code{dim(m)[1] = dim(m)[2] =
#'   length(m)}.
#'
#'
#' @param x graphMatch object
#'
#' @examples
#' # sample a pair of correlated random graphs from G(n,p)
#' set.seed(123)
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#'
#' # match g1 & g2 using FW methodology with indefinite relaxation
#' match <- gm(A = g1, B = g2, seeds = 1:3, method = 'indefinite')
#'
#' # print graphMatch object
#' match
#' print(match)
#' show(match)
#'
#' # print matching correspondence
#' match@corr
#' match$corr_A # matching correspondence in the first graph
#' match$corr_B # matching correspondence in the second graph
#'
#' # get nonseed matching correspondence
#' match[!match$seeds]
#'
#' @seealso \link[iGraphMatch:plot,igraph,igraph-method]{graphMatch_plot}, \link[iGraphMatch:\%*\%,graphMatch,ANY-method]{graphMatch_operators}
#'
#' @rdname graphMatch_methods
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


#' @rdname graphMatch_methods
#'
#' @param object graphMatch object
#'
setMethod("show", signature("graphMatch"), show.graphMatch)


#' @rdname graphMatch_methods
setMethod(
  "print",
  signature("graphMatch"),
  function(x) show.graphMatch(x)
)


#' @rdname graphMatch_methods
#'
#' @param i row index for the correspondence data.frame
#' @param j col index for the correspondence data.frame
#' @param drop ignored
#'
#' @examples
#' # get corresponding permutation matrix for the match result
#' match[] # preferred approach
#' # or equivalently
#' get_perm_mat(match)
#'
setMethod("[",
  signature(x = "graphMatch",
    i = 'missing', j = 'missing', drop = 'missing') ,
  function(x, i = NULL, j = NULL, drop = NULL) {
          get_perm_mat(x)
  }
)

#' @rdname graphMatch_operators
#'
#' @title Operator methods for graphMatch objects
#'
#' @description Methods to use graphMatch objects as operators on
#'  igraph and matrix-like objects.
#'
#'
#' @param x Either graphMatch object or a matrix-like object
#' @param y Either graphMatch object or a matrix-like object
#'
#' @return These methods return an object of the same type
#'  as the non-graphMatch object. If m is the match of g1
#'  to g2 (both igraph objects), then m %*% g2 returns g2
#'  permuted so as to match with g1. Conversely, g1 %*% m
#'  returns g1 permuted so as to match with g2.
#'
#' @examples
#'
#' set.seed(123)
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#'
#' # match g1 & g2 using FW methodology with indefinite relaxation
#' match <- gm(A = g1, B = g2, seeds = 1:3, method = 'indefinite')
#'
#' # permute the second graph according to the match result: P %*% g2 %*% P^T
#' match %*% g2 # return an igraph object
#' # equivalent to:
#' match %*% g2 %*% t(match)
#'
#' match %*% g2[] # return a matrix
#' # equivalent to:
#' P <- match[]
#' P %*% g2[] %*% Matrix::t(P)
#'
setMethod("%*%", signature(x = "graphMatch", y = "ANY"),
  function(x, y) {
    x[] %*% y %*% t(x[])
  }
)

#' @rdname graphMatch_operators
setMethod("%*%", signature(x = "ANY", y = "graphMatch"),
  function(x, y){
    t(y[]) %*% x %*% y[]
  }
)

#' @rdname graphMatch_operators
setMethod("%*%", signature(x = "graphMatch", y = "igraph"),
  function(x, y) {
    m <- min(dim(x))
    cA <- x[,1]
    cB <- x[,2]
    permuted_subgraph(y, cB[cA[cA <= m]])
  }
)

#' @rdname graphMatch_operators
setMethod("%*%", signature(x = "igraph", y = "graphMatch"),
  function(x, y) {

    m <- min(dim(y))
    cA <- y[,1]
    cB <- y[,2]
    permuted_subgraph(x, cA[order(cB[cB <= m])])
  }
)

#' @rdname graphMatch_methods
#'
#' @examples
#'
#' # sizes of two graphs
#' dim(match)
setMethod("dim", signature(x = "graphMatch"), function(x) { x@nnodes })

#' @rdname graphMatch_methods
#'
#' @examples
#'
#' # number of matched node pairs
#' length(match)
#'
setMethod("length", signature(x = "graphMatch"), function(x) { nrow(x@corr)})


reverse_match <- function(x) {
  x@corr[,c(1,2)] <- x@corr[,c(2,1)]
  names(x@corr) <- rev(names(x@corr))
  x@nnodes <- rev(x@nnodes)
  for (n in names(x)) {
    if (is.matrix(x[[n]]) || is(x[[n]], "Matrix")) {
      x[[n]] <- t(x[[n]])
    }
    # if (n == "seeds") {
    #   x$seeds <- x$seeds[, c(2,1)]
    #   names(x$seeds) <- rev(names(x$seeds))
    # }
  }
  x
}

#' @rdname graphMatch_methods
#'
#' @examples
#'
#' # reverse the matching correspondence between two graphs
#' t(match)
#' rev(match)
setMethod("t", signature(x = "graphMatch"), reverse_match)


#' @rdname graphMatch_methods
setMethod("rev", signature(x = "graphMatch"), reverse_match)



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

#'
#' @param name name of element in the list
#'
#' @rdname graphMatch_methods
setMethod("$", signature(x = "graphMatch"),
  function(x, name) {
    if (name == "corr")
      return(x@corr)

    if (name == "corr_A")
      return(x@corr$corr_A)

    if (name == "corr_B")
      return(x@corr$corr_B)

    i <- which(names(x) == name)
    if (length(i) == 0)
      return(NULL)
    x@.Data[[i]]
  }
)

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
#' @param linetype Whether to set edge line types according
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
#' Grey edges/pixels indicate common edges, blue
#' indicates edges only in graph A and red
#' represents edges only graph B. The corresponding
#' linetypes are solid, long dash, and short dash.
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
#' res <- gm(A, B, 1:4, method = "percolation")
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
  if(is.null(igraph::V(g)$name)){
    g <- igraph::set.vertex.attribute(g, "name", corr_g, corr_g)
  }

  g <- igraph::permute.vertices(
    igraph::induced_subgraph(g, corr_g),
    rank(corr_g)
  )
}


summary_graphMatch <- function(object, A = NULL, B = NULL, true_label = NULL, directed = NULL) {

  # Matched nodes
  corr <- object@corr
  object$n_match <- nrow(corr) - sum(object$seeds)
  if(!is.null(true_label)){
    object$n_true_match <-
      sum(true_label[corr$corr_A] == corr$corr_B) - sum(object$seeds)

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

#' @method summary graphMatch
#' @export
summary.graphMatch <- function(object, ...) {
  summary_graphMatch(object, ...)
}


setGeneric("summary")

#' @title Summary methods for graphMatch objects
#'
#' @param object graphMatch object
#' @param A igraph or matrix-like object
#' @param B igraph or matrix-like object
#' @param true_label the true correspondence (if available)
#' @param directed whether to treat the graphs as directed (TRUE) or not
#'   directed (FALSE) default is NULL which will treat the graphs as directed if
#'   either adjacency matrix is not symmetric.
#'
#' @return \code{summary} returns the graph matching formula, and a summary of
#'   graph matching results including the number of matches, the number of
#'   correct matches (if the true correspondence is available), and common
#'   edges, missing edges, extra edges, common non-edges and the objective
#'   function value.
#'
#' @examples
#' set.seed(123)
#' graphs <- sample_correlated_gnp_pair(20, .9, .3)
#' A <- graphs$graph1
#' B <- graphs$graph2
#' match <- gm(A, B, 1:4, method = "percolation")
#'
#' summary(match, A, B)
#' summary(match, A, B, true_label = 1:20) # also output the number of correct matches
#'
#' @rdname graphMatch_summary
setMethod("summary", signature("graphMatch"), summary_graphMatch)

show.summary.graphMatch <- function(match) {
    cat("Call: ")
    print(match@call)
    cat("\n# Matches:", match$n_match)
    if(!is.null(match$n_true_match)){
      cat("\n# True Matches: ", match$n_true_match)
    }
    if(!is.null(match$seeds)) {
      cat(", # Seeds: ", sum(match$seeds))
    }
    cat(", # Vertices: ", paste(dim(match), collapse = ", "))
    cat("\n")
    match$edge_match_info
    if(!is.null(match$edge_match_info)) {
      ep <- as.data.frame(t(match$edge_match_info))
      ep <- cbind(
        data.frame(name = rownames(ep)),
        ep)
      rownames(ep) <- NULL
      if(ep$name[1] == "layer") {
        colnames(ep) <- sapply(ep[1,], as.character)
        ep <- ep[2:nrow(ep),]
      } else {
        colnames(ep) <- NULL
      }
      print(ep, row.names = FALSE)
    }
}

setClass("summary.graphMatch",
  contains = 'graphMatch'
)


setMethod("show", signature("summary.graphMatch"),
  function(object){
    show.summary.graphMatch(object)
  }
)


setMethod("print", signature("summary.graphMatch"),
  function(x){
    show.summary.graphMatch(x)
  }
)

