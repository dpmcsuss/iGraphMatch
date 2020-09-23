g <- igraph::sample_gnm(20, 60)
igraph::E(g)$color <- sample(c("red", "green"), 60, replace = TRUE)

test_that("Correct number of layers", {
  l <- split_igraph(g, "color")
  expect_length(l, 2)
})

igraph::V(g)$weight <- runif(igraph::vcount(g))

test_that("Spripped attrs", {
  l <- split_igraph(g, "color", strip_vertex_attr = TRUE)
  expect_length(igraph::vertex_attr_names(l[[1]]), 0)
  expect_length(igraph::vertex_attr_names(l[[2]]), 0)
})


test_that("Error on non igraph", {
  expect_error(
    split_igraph(1, "a"),
    "g must be an igraph object"
  )
})



# test_that("Warning on many attributes", {
#   g <- igraph::sample_gnm(20, 60)
#   E(g)$color <- as.character(sample(seq(40), 60, replace = TRUE))
#   expect_warning(split_igraph(g, "color"))
# })
