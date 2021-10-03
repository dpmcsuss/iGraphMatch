#' @title Graph Matching Methods
#'
#' @description \code{gm} is used to match a pair of given graphs, with
#'   specifications of the adjacency matrices of for a pair of graphs, possible
#'   prior knowledge, and a graph matching method.
#'
#' @param A A matrix, 'igraph' object, or list of either.
#' @param B A matrix, 'igraph' object, or list of either.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be a matrix or a data frame, with the first
#'   column being the indices of \eqn{G_1} and the second column being the
#'   corresponding indices of \eqn{G_2}.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex
#'   similarities. Mandatory for the "IsoRank" method.
#' @param method Choice for graph matching methods. One of "indefinite",
#'   "convex", "PATH", "percolation", "IsoRank", "Umeyama", or a user-defined
#'   graph matching function. Please check Details and Examples sections for
#'   instructions on how to define your own function.
#' @param ... Arguments passed to graph matching methods. Please refer to
#'   Details section for more information.
#'
#' @rdname gm
#'
#' @details If \code{method} is a function, it should take two matrices or
#'   igraph objects, seeds and similarity scores as arguments for minimum.
#'   Additionally, it can also take other arguments if needed. The self-defined
#'   function should return a graphMatch class object with matching
#'   correspondence, sizes of two input graphs, matching formula, and other
#'   algorithm hyperparameter details.
#'
#'   The \code{method} argument can also take one of the implemented algorithms,
#'   including \link[=graph_match_indefinite]{"indefinite"},
#'   \link[=graph_match_convex]{"convex"}, \link[=graph_match_PATH]{"PATH"},
#'   \link[=graph_match_percolation]{"percolation"}, \link[=graph_match_IsoRank]{"IsoRank"},
#'   and \link[=graph_match_Umeyama]{"Umeyama"}.
#'   In this case, one can pass additional arguments to the \code{gm} function
#'   according to the specified method.
#'   For a detailed list of additional arguments for each one of the implemented method,
#'   please click on the corresponding method name for its help page.
#'
#'
#' @return \code{gm} returns an object of class "\code{\link{graphMatch}}".
#' The functions \code{summary} and \code{plot} can be used to get a summary of
#' graph matching results and visualization of matches.
#' Implemented operations for the returned \code{graphMatch} objects enables obtaining
#' flexible representations of matching results,
#' such as extracting a subset of matching correspondence, get the corresponding
#' permutation matrix, and get the permuted graphs using matching results.
#'
#' An object of class "\code{\link{graphMatch}}" contains the following slots,
#' and each can be called by using @
#'
#'   \describe{
#'     \item{corr}{data.frame indicating the correspondence between two graphs}
#'     \item{nnodes}{of the original two graphs}
#'     \item{call}{The call to the graph matching function}
#'   }
#'
#' Additionally, \code{gm} also returns a list of matching details of the specified method.
#' Please refer to the help page for each implemented method, i.e.
#' \link[=graph_match_indefinite]{"indefinite"},
#' \link[=graph_match_convex]{"convex"}, \link[=graph_match_PATH]{"PATH"},
#' \link[=graph_match_percolation]{"percolation"}, \link[=graph_match_IsoRank]{"IsoRank"},
#' and \link[=graph_match_Umeyama]{"Umeyama"} for details on the corresponding returned list.
#'
#'
#' @examples
#' # match G_1 & G_2 with some known node pairs as seeds
#' set.seed(123)
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.5, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 4
#'
#' m_rds <- gm(g1, g2, seeds, method = "indefinite", start = "rds", max_iter = 20)
#' summary(m_rds, g1, g2, true_label = 1:10)
#' m_rds@call
#' m_rds@corr
#' m_rds[!m_rds$seeds] # matching correspondence for nonseed vertices
#'
#' m_rds$soft # doubly stochastic matrix from the last iteration
#' m_rds[] # P: permutation matrix
#' m_rds$iter
#'
#' m_rds %*% g2 # permute the second graph: P %*% g2 %*% P^T
#'
#' plot(g1, g2, m_rds)
#' plot(g1[], g2[], m_rds)
#'
#' m_path <- gm(g1, g2, method = "PATH", lap_method = "lapmod")
#'
#'
#' # match two multi-layer graphs
#' set.seed(123)
#' gp_list <- replicate(3, sample_correlated_gnp_pair(20, .3, .5), simplify = FALSE)
#' A <- lapply(gp_list, function(gp)gp[[1]])
#' B <- lapply(gp_list, function(gp)gp[[2]])
#'
#' m_perco <- gm(A, B, seeds, method = "percolation", ExpandWhenStuck = FALSE)
#' summary(m_perco, A, B)
#'
#' sim <- as.matrix(init_start(start = "bari", nns = 20, soft_seeds = 1:5))
#' m_Iso <- gm(A, B, similarity = sim, method = "IsoRank", lap_method = "greedy")
#' summary(m_Iso, A, B)
#'
#' m_Umeyama <- gm(A, B, method = "Umeyama")
#'
#'
#' # customized graph matching algorithm
#' graph_match_rand <- function(A, B, seeds = NULL, similarity = NULL, rand_seed){
#'   totv1 <- nrow(A[[1]])
#'   totv2 <- nrow(B[[1]])
#'   nv <- max(totv1, totv2)
#'
#'   corr_A <- 1:nv
#'   set.seed(rand_seed)
#'   corr_B <- c(1:nv)[sample(nv)]
#'   corr <- data.frame(corr_A, corr_B)
#'
#'   graphMatch(
#'     corr = corr,
#'     nnodes = c(totv1, totv2),
#'     detail = list(
#'       rand_seed = rand_seed
#'     )
#'   )
#' }
#'
#' m_self <- gm(
#'   g1, g2,
#'   method = graph_match_rand,
#'   rand_seed = 123 # pass additional argument 'rand_seed' to input
#' )
#' summary(m_self, g1, g2)
#' m_self$rand_seed # graph_match_rand method hyperparameter
#' m_self@call
#' m_self@nnodes
#' m_self@corr
#'
#'
#'
#' @export
#'
#'
gm <- function(A, B, seeds = NULL, similarity = NULL, method = "indefinite", ...){

  methods <- c("indefinite", "convex", "PATH", "percolation", "IsoRank", "Umeyama")
  if (!is.function(method) && !(method %in% methods)) {
    stop("Method must be one of: ", paste0(paste(methods, collapse = ", "),
                                           " or a function that takes a pair of graphs and other prior knowledge if applicable."))
  }

  # A, B argument checks
  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2

  # seeds argument check
  seed_check <- check_seeds(seeds, nv = max(totv1, totv2))
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds

  # similarity score matrix argument check
  similarity_raw <- similarity
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)

  if(is.function(method)){
    m <- method(A, B, seeds, similarity, ...)
  } else if(method == "indefinite"){
    m <- graph_match_indefinite(A, B, seeds, similarity, ...)
  } else if(method == "convex"){
    m <- graph_match_convex(A, B, seeds, similarity, ...)
  } else if(method == "PATH"){
    m <- graph_match_PATH(A, B, seeds, similarity, ...)
  } else if(method == "percolation"){
    if(nrow(seeds) == 0 & is.null(similarity_raw)){
      stop("At least one of seeds and similarity score should be known for this method.")
    }
    m <- graph_match_percolation(A, B, seeds, similarity, ...)
  } else if(method == "IsoRank"){
    if(is.null(similarity_raw)){
      stop("Similarity scores are mandatory for this method. Please input a value for the 'similarity' argument.")
    }
    similarity <- check_sim(similarity_raw, seeds, nonseeds, totv1, totv2, for_nonseeds = FALSE)
    m <- graph_match_IsoRank(A, B, seeds, similarity, ...)
  } else if(method == "Umeyama"){
    m <- graph_match_Umeyama(A, B, seeds, similarity)
  }
  tryCatch(
    {
      m@nnodes <- c(totv1, totv2)
      m@call <- match.call()
      m
    },
    error = function(e) {
      stop(
        'Customized graph matching method function must return a graphMatch class object. See ?gm for examples.',
        e
      )
    }
  )
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
"Thanks for using iGraphMatch!
We'd love to get feedback on what you like, what you don't like,
and how you are using the package.
See ?iGraphMatch for contact information.")
}
