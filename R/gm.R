#' @title Graph Matching Methods
#'
#' @description \code{gm} is used to match a pair of given graphs, with specifications
#'   of the adjacency matrices of for a pair of graphs, possible prior knowledge, and
#'   a graph matching method.
#'
#' @param A A matrix, 'igraph' object, or list of either.
#' @param B A matrix, 'igraph' object, or list of either.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix or a data frame, with the first
#'   column being the indices of \eqn{G_1} and the second column being the
#'   corresponding indices of \eqn{G_2}.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param method Choice for graph matching methods.
#' @param ... Arguments passed to graph matching methods. Please refer to Details section
#'  for more information.
#'
#' @rdname gm
#'
#' @details If \code{method} is a function, it should take two matrices or igraph objects as
#' arguments for minimum. Additionally, it can also take prior information in the form of
#' \code{seeds} or \code{similarity} or both, and other arguments if needed. The self-defined
#' function should return a graphMatch class object with matching correspondence, sizes of two
#' input graphs, matching formula, and other algorithm hyperparameter details.
#'
#'
#' @return \code{graph_match_indefinite}, \code{graph_match_convex} and \code{graph_match_PATH}
#'   return an object of class "gm" which is a list containing the following
#'   components:
#'
#'   \describe{
#'     \item{corr_A}{matching correspondence in \eqn{G_1}}
#'     \item{corr_B}{matching correspondence in \eqn{G_2}}
#'     \item{soft}{the doubly stochastic matrix from the last iteration with which one can
#'           extract more than one matching candidates}
#'     \item{iter}{number of iterations until convergence or reaches the \code{max_iter}}
#'     \item{max_iter}{Maximum number of replacing matches}
#'     \item{lap_method}{Choice for solving the LAP}
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'   }
#'
#'
#'
#' @examples
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 3
#'
#' # customized graph matching algorithm
#' graph_match_rand <- function(A, B, rand_seed){
#'   totv1 <- nrow(A[[1]])
#'   totv2 <- nrow(B[[1]])
#'   nv <- max(totv1, totv2)
#'
#'   corr_A <- 1:nv
#'   set.seed(rand_seed)
#'   corr_B <- c(1:nv)[sample(nv)]
#'   corr <- data.frame(corr_A, corr_B)
#'
#'   cl <- match.call()
#'   graphMatch(
#'     corr = corr,
#'     nnodes = c(totv1, totv2),
#'     call = cl,
#'     detail = list(
#'       rand_seed = rand_seed
#'     )
#'   )
#' }
#'
#' m_self <- gm(g1, g2, method = graph_match_rand, rand_seed = 123) # pass additional argument 'rand_seed' to input
#' summary(m_self, g1, g2)
#' m_self$rand_seed # graph_match_rand method hyperparameter
#' m_self@call
#' m_self@nnodes
#' m_self@corr
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
    m <- method(A, B, ...)
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
