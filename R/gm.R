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
#'   vector. If not, seeds must be  a matrix or a data frame, with the first
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
#'   The \code{method} argument can also take one of the implemented algorithms.
#'   In this case, one can pass additional arguments to the \code{gm} function
#'   according to the specified method.
#'
#'   For "indefinite", "convex" and "PATH" methods, additional arguments with their
#'   default values contain the followings:
#'
#'   \describe{
#'     \item{\code{start = "bari"}}{A matrix or a character. Any \code{nns-by-nns}
#'       matrix or character value like "bari", "rds" or "convex" to
#'       initialize the starting matrix.}
#'     \item{\code{max_iter = 20}}{A number. Maximum number of replacing matches.
#'       When the algorithm reaches \code{max_iter} before convergence, it would
#'       stop and give a warning message.
#'       Note that the default value is 100 for "convex".}
#'     \item{\code{tol = 1e-05}}{A number. Tolerance of edge disagreements for "convex" and
#'       "PATH" methods. The default value is 1e-05.}
#'     \item{\code{lap_method = NULL}}{Choice for lap method. One of "lapjv", "lapmod",
#'       or "clue".
#'       When \code{lap_method} takes default value, the function automatically selects
#'       a preferred lap method based on graph sizes.}
#'  }
#'
#'  The "IsoRank" method also has \code{max_iter} and \code{lap_method} as additional
#'  arguments. The default for \code{max_iter} is 50. \code{lap_method} takes either "greedy"
#'  (default) or "LAP", where the former extracts mapping using the greedy algorithm and
#'  not necessarily grows a mapping to the entire graphs, whereas the latter solves a LAP.
#'
#'  The "percolation" method has two additional arguments:
#'
#'  \describe{
#'    \item{\code{r = 2}}{A number. Threshold of neighboring pair scores.}
#'    \item{\code{ExpandWhenStuck = FALSE}}{A logical. TRUE if expand the seed
#'      set when Percolation algorithm stops before matching all the vertices.}
#'  }
#'
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
#' Additionally, \code{gm} also returns a list of matching details corresponding
#' to the specified method. The returned list for the "percolation" method contains
#'
#'   \describe{
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'     \item{match_order}{the order of vertices getting matched}
#'   }
#'
#' Notably, \code{seeds} information is returned for all the methods. In addition to
#' \code{seeds}, "convex", "indefinite", "PATH", "IsoRank" and "Umeyama" also return
#'
#'   \describe{
#'     \item{soft}{the doubly stochastic matrix from the last iteration or
#'       functional similarity score matrix with which one can extract more than
#'       one matching candidates}
#'     \item{lap_method}{Choice for extracting matches or solving the LAP}
#'   }
#'
#' "IsoRank" method also returns match_order when \code{lap_method="greedy"}. "PATH",
#' "convex" and "indefinite" methods return two more components, max_iter and iter
#' to indicate the maximum number of replacing matches and the actual number of iterations.
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
#'     call = match.call(),
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
#'
#' @references V. Lyzinski and D. E. Fishkind and M. Fiori and J. T. Vogelstein and C. E. Priebe
#' and G. Sapiro (2016), \emph{Graph Matching: Relax at Your Own Risk}. IEEE TPAMI, pages 60-73.
#' @references V. Lyzinski and D. E. Fishkind and C. E. Priebe (2014), \emph{Seeded Graph Matching
#' for Correlated Erdos-Renyi Graphs}.J. Mach. Learn. Res., pages 3513-3540.
#'
#' @references Y. Aflalo and A. Bronstein and R. Kimmel (2014), \emph{On convex
#' relaxation of graph isomorphism}. Proceedings of the National Academy of Sciences,
#' pages 2942-2947.
#'
#' @references M. Zaslavskiy, F. Bach and J. Vert (2009), \emph{A Path following
#' algorithm for the graph matching problem}. IEEE Trans Pattern Anal Mach Intell,
#' pages 2227-2242.
#'
#' @references L. Yartseva and M. Grossglauser (2013), \emph{On the performance
#'   of percolation graph matching}. COSN, Boston, MA, USA, pages 119–130.
#' @references E. Kazemi, S. H. Hassani, and M. Grossglauser (2015),
#' \emph{Growing a graph matching from a handful of seeds}. Proc. of the VLDB
#' Endowment, 8(10):1010–1021.
#'
#' @references R. Singh, J. Xu, B. Berger (2008), \emph{Global alignment of
#' multiple protein interaction networks with application to functional
#' orthology detection}. Proc Natl Acad Sci. USA, pages 12763-12768.
#'
#' @references S. Umeyama (1988), \emph{An eigendecomposition approach to weighted
#'   graph matching problems}. IEEE TPAMI. USA, pages 695-703.
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
