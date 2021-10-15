#' Chemical synapses and electrical synapses networks of roundworm
#'
#' C.Elegans networks consist of the chemical synapses network and the
#' electrical synapses network of the roundworm, where each of 279 nodes
#' represents a neuron and each edge represents the intensity of synapses
#' connections between two neurons.
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



