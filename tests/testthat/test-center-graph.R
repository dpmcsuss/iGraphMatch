

set.seed(123)
A <- sample_correlated_gnp_pair(n = 3, corr = .5, p = .5)$graph1
test_that("naive scheme", {
  expect_equal(
    as.matrix(center_graph(A, scheme = "naive")),
    matrix(c(0,0,1,0,0,0,1,0,0),nrow=3),
    ignore_attr = TRUE
  )
})
test_that("center scheme", {
  expect_equal(
    as.matrix(center_graph(A, scheme = "center")),
    matrix(c(-1,-1,1,-1,-1,-1,1,-1,-1),nrow=3),
    ignore_attr = TRUE
  )
})
test_that("scheme equal to 1", {
  expect_equal(
    as.matrix(center_graph(A, scheme = 1)),
    matrix(c(0,0,1,0,0,0,0,0,0),nrow=3),
    ignore_attr = TRUE
  )
})
test_that("input scheme check", {
  expect_error(
    center_graph(A, scheme = c(-1, 1, 1)),
    "scheme must be either 'center', 'naive', a positive integer, or a pair of scalars.")
  expect_error(
    center_graph(A, scheme = "other scheme"), 
    "scheme must be either 'center', 'naive', a positive integer, or a pair of scalars.")
})
test_that("use splr", {
  expect_equal(
    as.matrix(center_graph(A, scheme = 1, use_splr = FALSE)),
    matrix(c(rep(0, 6), 1, 0, 0),byrow = TRUE, nrow=3),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(center_graph(A, scheme = c(-1, 1), use_splr = FALSE)),
    matrix(c(-1,-1,1,-1,-1,-1,1,-1,-1),nrow=3),
    ignore_attr = TRUE
  )
})

A <- sample_correlated_gnp_pair(100, .5, .5)$graph1
test_that("use splr", {
  expect_equal(
    nrow(center_graph(A, scheme = 2, use_splr = FALSE)),
    100
  )
})


test_that("padding", {
  expect_equal(
    nrow(pad(A[], 2)),
    102
  )
  expect_equal(
    nrow(pad(center_graph(A, scheme = 2, use_splr = TRUE), 2)),
    102
  )
})

gp_list <- replicate(3, sample_correlated_gnp_pair(20, .3, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
A <- as(A, "matrix_list")

test_that("padding for multi-layer graphs", {
  expect_equal(
    nrow(pad(A, 2)),
    22
  )
})
