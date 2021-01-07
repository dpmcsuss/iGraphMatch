#' Chemical synapses and electrical synapses networks of roundworm
#'
#' C.Elegans networks consist of the chemical synapses network and the
#' electrical synapses network of the roundworm, where each of 279 nodes
#' represents a neuron and each edge represents the intensity of synapses
#' connections between two neurons.
#' Two networks are weighted and directed graphs with self-loops.
#' There are 2194 and 1031 edges in two graphs respectively and the empirical
#' Pearson's correlation bewteen two graphs is 0.17.
#' Two networks are stored in a list in the form of igraph objects, where
#' the first network in the list is the chemical synapses network and the
#' other one is the electrical synapses network.
#'
#' @usage data(C.Elegans)
#' 
#' @examples 
#' data(C.Elegans)
#' g1 <- C.Elegans[[1]]
#' g2 <- C.Elegans[[2]]
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
#' Pearson's correlation bewteen two graphs is 0.85.
#' Two email communication networks for two different weeks are stored in a
#' list in the form of igraph objects.
#'
#' @usage data(Enron)
#' 
#' #' @examples 
#' data(Enron)
#' g1 <- Enron[[1]]
#' g2 <- Enron[[2]]
#' 
"Enron"



#' Britain Transportation network
#' 
#' The network reflects the transportation connections in the UK, with five 
#' layers representing ferry, rail, metro, coach, and bus.
#' The original graph has 262385 nodes and 461555 connections in five layers in 
#' total.
#' A smaller graph was constructed based on a random walk starting from a 
#' randomly chosen hub node, a node that has connections in all the layers.
#' The template graph has 53 nodes and 56 connections in total and is an induced 
#' subgraph of the original graph.
#' In addition, the dataset also contains a list of candidate matches for each 
#' template node, where the true correspondence is guaranteed to among the candidates.
#' The world graph, the template graph and the candidate matches information are stored 
#' in a list where two graphs are stored in two lists of igraph objects representing
#' different layers, and the candidate information is stored in a data frame where the 
#' ids of the template nodes and the world nodes are in columns named "tem" and "wor".
#'
#' @references Gallotti Riccardo and Barthelemy Marc. (2015), \emph{The Multilayer 
#' Temporal Network of Public Transport in Great Britain}. Scientific Data. 
#' 
#' @usage 
#' data(Transportation)
#' 
#' @examples 
#' data(Transportation)
#' template <- Transportation[[1]]
#' world <- Transportation[[2]]
#' candidate <- Transportation[[3]]
#' 
"Transportation"
