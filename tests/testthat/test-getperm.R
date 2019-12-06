context("get_perm")


corr <- data.frame(corr_A = c(1,2,3,4), corr_B = c(1,4,2,3))

test_that("permutation matrix", {
  expect_equivalent(as.matrix(get_perm(4, 4, corr)),
               matrix(c(1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0),nrow=4))
})

