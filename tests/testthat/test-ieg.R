

#Sample graphs from edge probability matrix and correlation matrix
set.seed(123)
n <- 10
p_mat <- matrix(runif(n^2),n)
c_mat <- matrix(runif(n^2),n)
g <- sample_correlated_ieg_pair(n, p_mat, c_mat, directed = TRUE,loops = TRUE, permutation = 1:n)

test_that("number of vertices", {
  expect_equal(igraph::vcount(g$graph1),10) 
  expect_equal(igraph::vcount(g$graph2),10)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g$graph1),56)
  expect_equal(igraph::ecount(g$graph2),56)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g$graph1),c(10,10,13,13,11,9,13,13,12,8))
  expect_equal(igraph::degree(g$graph2),c(11,10,10,11,9,13,13,13,11,11)) 
})

# Sample undirected ieg pairs
set.seed(123)
n <- 10
p_mat <- matrix(runif(n^2),n)
c_mat <- matrix(runif(n^2),n)
g <- sample_correlated_ieg_pair(n, p_mat, c_mat, directed = FALSE,loops = TRUE, permutation = 1:n)

test_that("sample undirected ieg pairs", {
  expect_equal(igraph::vcount(g$graph1),10)
  expect_equal(igraph::vcount(g$graph2),10)
})

# Sample random dot product graphs

n <- 10
xdim <- 3
scale <- 8
X <- matrix(rgamma(n*(xdim+1),scale,1),n,xdim+1)
X <- X/rowSums(X)
X <- X[,1:xdim]
g <- sample_correlated_rdpg(X,rho=0.5)

test_that("sample undirected rdpg pairs", {
  expect_equal(igraph::vcount(g$graph1),10)
  expect_equal(igraph::vcount(g$graph2),10)
})


g <- sample_correlated_rdpg(X,rho=c_mat)

test_that("sample undirected rdpg pairs w. rho mat", {
  expect_equal(igraph::vcount(g$graph1),10)
  expect_equal(igraph::vcount(g$graph2),10)
})
