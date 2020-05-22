test_that("Correct number of layers", {
  g <- igraph::sample_gnm(20, 60)
  igraph::E(g)$color <- sample(c("red", "green"), 60, replace = TRUE)
  l <- split_igraph(g, "color")
  expect_length(l, 2)
})

test_that("Error on non igraph", {
  expect_error(split_igraph(1, "a"))
})

# test_that("Warning on many attributes", {
#   g <- igraph::sample_gnm(20, 60)
#   E(g)$color <- as.character(sample(seq(40), 60, replace = TRUE))
#   expect_warning(split_igraph(g, "color"))
# })
