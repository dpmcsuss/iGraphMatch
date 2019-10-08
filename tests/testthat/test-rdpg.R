context("rdpg")

#Sample graphs from edge probability matrix and correlation matrix
set.seed(123)
n <- 50
xdim <- 3
scale <- 8
X <- matrix(rgamma(n*(xdim+1),scale,1),n,xdim+1)
X <- X/rowSums(X)
X <- X[,1:xdim]
g<-sample_correlated_rdpg(X,rho=0.5,nc=nrow(X),directed=TRUE,loops=TRUE,permutation=c(sample(1:n)))

test_that("number of vertices", {
  expect_equal(igraph::vcount(g$graph1),50) 
  expect_equal(igraph::vcount(g$graph2),50)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g$graph1),465)
  expect_equal(igraph::ecount(g$graph2),478)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g$graph1),c(19,22,13,17,14,16,14,17,15,16,15,18,14,21,17,16,14,18,20,22,17,3,24,18,32,14,21,23,20,
                                  21,20,19,16,18,24,29,15,18,15,14,22,23,15,23,23,22,19,25,19,20))
  expect_equal(igraph::degree(g$graph2),c(18,16,17,21,17,19,25,22,28,15,16,23,20,21,15,14,13,15,19,16,16,26,21,24,20,22,22,26,22,
                                  19,17,22,17,5,16,16,20,20,21,20,21,20,24,18,17,24,17,18,11,24)) 
})