#' @title Split an igraph object into aligned graphs by attribute
#'
#' @description Given an igraph object and an edge attribute, this function
#' finds all unique values of the edge attribute in the graph
#' and returns a list of igraph objects on the same vertex set
#' where each element of the list has a graph containing only
#' those edges with specified attributed.
#'
#' @param g An igraph object
#' @param e_attr the name of an edge attribute in g
#' @param strip_vertex_attr Whether to remove all vertex
#'  attribute from the new graphs
#'
#' @returns A named list of igraph objects with names corresponding to the values of
#'  the edge attributes.
#'
#' @examples
#' g <- igraph::sample_gnm(20, 60)
#' igraph::E(g)$color <-
#'   sample(c("red", "green"), 60, replace = TRUE)
#' split_igraph(g, "color")
#'
#' @export
split_igraph <- function(g, e_attr, strip_vertex_attr = FALSE) {
  if (!igraph::is.igraph(g)) {
    stop("g must be an igraph object")
  }

  all_attr <- igraph::get.edge.attribute(g, e_attr)
  u_attr <- unique(all_attr)
  try({
    u_attr <- sort(u_attr)
  }, silent = TRUE)

  if (strip_vertex_attr) {
    for (v_attr in igraph::vertex_attr_names(g)){
      g <- igraph::delete_vertex_attr(g, v_attr)
    }
  }

  sapply(u_attr, function(u) {
    igraph::subgraph.edges(
      g, igraph::E(g)[all_attr == u], delete.vertices = FALSE
    )
  }, USE.NAMES = TRUE, simplify = FALSE)
}
