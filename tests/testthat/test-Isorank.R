

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1
seeds<-1:4



test_that("isorank match with greedy LAP", {
  tt <- gm(A, B, seeds, startm, method = "IsoRank", lap_method = "greedy")
  expect_snapshot_output(print(tt))
  expect_snapshot_output(print(round(as.matrix(tt$soft), 4)))

})

test_that("isorank match with hungarian lap", {
  tt <- gm(A, B, seeds, startm, method = "IsoRank", lap_method = "LAP")
  expect_snapshot_output(print(tt))
  expect_snapshot_output(print(round(as.matrix(tt$soft), 4)))
})



# test_that("order of nodes getting matched", {
#   expect_equal(gm(A, B, seeds, startm, method = "IsoRank", lap_method = "greedy")$order,
#                c(1,2,3,4,6,10,8,7,5,9))
# })
# test_that("test LAP method", {
#   expect_equal(gm(A, B, seeds, startm, method = "IsoRank", lap_method = "LAP")$seeds,
#                data.frame(A = 1:4, B = 1:4))
# })



set.seed(12)
gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])
seeds <- 1:3

test_that("IsoRank multi-layer", {
  tt <- gm(A, B, seeds, startm, method = "IsoRank", lap_method = "LAP")
  expect_snapshot_output(print(tt))
  expect_snapshot_output(print(round(as.matrix(tt$soft), 4)))
})
