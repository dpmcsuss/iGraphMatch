

test_that("vector type input", {
  expect_equal(check_seeds(seeds = c(1,3,5), nv = 5), 
               list(seeds = data.frame(A = c(1, 3, 5), B = c(1, 3, 5)),
                    nonseeds = data.frame(A = c(2, 4), B = c(2, 4))))
}) 
  
test_that("matrix input containing incorrect seeds", {
  expect_equal(check_seeds(seeds = matrix(c(1, 2, 1, 3), nrow = 2), nv = 3),
               list(seeds = data.frame(A = c(1, 2), B = c(1, 3)),
                    nonseeds = data.frame(A = 3, B = 2)))
})

test_that("no seeds", {
  expect_equal(check_seeds(seeds = NULL, nv = 3),
               list(seeds = data.frame(A = numeric(), B = numeric()), 
                    nonseeds = data.frame(A = 1:3, B = 1:3)))
})

test_that("example in documentation", {
  expect_equal(check_seeds(1:10 <= 3, nv = 10)$seeds, data.frame(A = 1:3, B = 1:3))
  expect_equal(check_seeds(c(1,4,2,7,3), nv = 10)$seeds, data.frame(A = c(1,4,2,7,3), B = c(1,4,2,7,3)))
  expect_equal(check_seeds(matrix(1:4,2), nv = 10)$seeds, data.frame(A = c(1,2), B = c(3,4)))
  expect_equal(check_seeds(matrix(1:4,2), nv = 10, logical = TRUE), 1:10 <= 2)
  expect_equal(check_seeds(as.data.frame(matrix(1:4,2)), nv = 10)$seeds, data.frame(A = c(1,2), B = c(3,4)))
})