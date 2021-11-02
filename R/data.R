#' Chemical synapses and electrical synapses networks of roundworm
#'
#' C.Elegans networks consist of the chemical synapses network and the
#' electrical synapses network of the roundworm, where each of 279 nodes
#' represents a neuron and each edge represents the intensity of synapses
#' connections between two neurons.
#'
#' Two networks are weighted and directed graphs with self-loops.
#' There are 2194 and 1031 edges in two graphs respectively and the empirical
#' Pearson's correlation between two graphs is 0.17.
#' Two networks are stored in a list in the form of igraph objects, where
#' the first network in the list is the chemical synapses network and the
#' other one is the electrical synapses network.
#'
#' @usage data(C.Elegans)
#'
#' @references Chen, L., Vogelstein, J. T., Lyzinski, V., & Priebe, C. E.
#' (2016). \emph{A joint graph inference case study: the C. elegans chemical
#' and electrical connectomes.} Worm, 5(2), e1142041.
#'
#' Sulston, J. E., Schierenberg, E., White, J. G., & Thomson, J.N. (1983).
#' \emph{The embryonic cell lineage of the nematode caenorhabditis
#'  elegans.} Developmental biology, 100(1):64–119.
#'
#' @examples
#' data(C.Elegans)
#' g1 <- C.Elegans[[1]]
#' g2 <- C.Elegans[[2]]
#' plot(g1, g2)
#'
"C.Elegans"


#' Email communication networks of Enron Corporation
#'
#' The Enron network data consists of email messages between 184 employees
#' of the Enron Corporation where each graph represents one week of emails
#' and each edge indicates whether there is email sent from one employee to
#' the other.
#'
#' Two networks are unweighted and directed with self-loops.
#' There are 488 and 482 edges in two networks respectively and the empirical
#' Pearson's correlation between two graphs is 0.85.
#' Two email communication networks for two different weeks are stored in a
#' list in the form of igraph objects.
#'
#' @usage data(Enron)
#'
#' @references Originally released by William Cohen at CMU. \href{http://www.cs.cmu.edu/~enron/}{More details}
#' on the origins and research uses of the dataset.
#'
#' @examples
#' data(Enron)
#' g1 <- Enron[[1]]
#' g2 <- Enron[[2]]
#' plot(g1, g2)
#'
"Enron"



#' Britain Transportation Network
#'
#' The Britain Transportation Network reflects the transportation connections in
#' the UK, with five layers representing ferry, rail, metro, coach, and bus.
#'
#' The data consists of a smaller template graph with 53 nodes and 56 connections
#' across five layers, a larger world graph with candidates of the template graph
#' with 2075 nodes and 8368 connections, and a list of candidate matches for each
#' template node, where the true correspondence is guaranteed to be among the candidates.
#'
#' The template graph was constructed based on a random walk starting from a randomly
#' chosen hub node, a node that has connections in all the layers.
#' All edges in the template are common edges shared by two graphs, where 40\%, 24.1\%,
#' 37.5\%, 31.7\% and 25.6\% of edges in the world graph are in template for each layer.
#' All graphs are unweighted, directed, and do not have self-loops.
#'
#' @format A list of length 3, corresponding to the template graph, world graph, and
#' candidate data frame with first column indicating template node ID's and second column
#' indicating world node ID's.
#' The template graph and world graph are stored as lists of five adjacency matrices,
#' representing ferry, rail, metro, coach, and bus transportation connections respectively.
#'
#' @usage data(Transportation)
#'
#' @seealso The original Britain Transportation Network
#' \href{https://math.bu.edu/people/sussman/data/Transportation.rda}{data}.
#' The template graph and world graph in the `Transportation` data are
#' induced subgraphs of the original graphs , keeping only the candidate nodes.
#'
#' @references Gallotti, R., Barthelemy, M. (2015). \emph{The multilayer temporal
#' network of public transport in Great Britain.} Sci Data 2, 140056 .
#' https://doi.org/10.1038/sdata.2014.56.
#'
#' J. D. Moorman, Q. Chen, T. K. Tu, Z. M. Boyd and A. L. Bertozzi, (2018).
#' \emph{Filtering Methods for Subgraph Matching on Multiplex Networks.} 2018 IEEE
#' International Conference on Big Data (Big Data), pp. 3980-3985,
#' doi: 10.1109/BigData.2018.8622566.
#'
#'
#' @examples
#' tm <- Transportation[[1]]
#' cm <- Transportation[[2]]
#' candidate <- Transportation[[3]]
#' tn <- nrow(tm[[1]])
#' wn <- nrow(cm[[1]])
#' similarity <- with(candidate, Matrix::sparseMatrix(i = tem, j = wor, x = 1,
#'                             dims = c(tn,wn)))
"Transportation"
