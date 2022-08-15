
set.seed(1234)

x <- Matrix::rsparsematrix(10, 10, .5)
a <- matrix(runif(30), 10)
b <- matrix(runif(30), 10)
ns <- matrix(runif(100), 10)

s <- splr(x, a, b)


test_that("error on no b or rank",
  {
    expect_error(
      splr(x, ns),
      "please provide an already factorized low-rank matrix.*"
    )
  }
)

test_that("error on no b and a is wrong dim",
  {
    expect_error(
      splr(x, a, rank = 2),
      "b is not provided and a is not the same dimension as x"
    )
  }
)

s <- splr(x, ns, rank = 2)


test_that("various operations work as expected",{

  expect_snapshot_output(s <- splr(x, Matrix(a), Matrix(b)))
  expect_snapshot_output(sp <- splr_to_sparse(s))
  expect_snapshot_output(m <- as.matrix(s))

  expect_snapshot_output(ss <- splr_sparse_plus_constant(x, 1))
  expect_snapshot_output(dm <- as(s, "dMatrix"))
  # try({
  #   sink("/dev/null")
  # })
  expect_snapshot_output(print(s))
  expect_snapshot_output(show(s))
  expect_snapshot_output(as.character(s))

  expect_snapshot_output(str(s))
  # sink()


  expect_snapshot_output(length(s))


  expect_snapshot_output(s %*% s)
  expect_snapshot_output(s * s)

  expect_snapshot_output(ml <- matrix_list(list(ns, x)))

  expect_snapshot_output(s %*% ml)
  expect_snapshot_output(x %*% s)
  expect_snapshot_output(ns %*% s)
  expect_snapshot_output(3 * s)

  expect_snapshot_output(ml %*% s)
  expect_snapshot_output(s %*% x)
  expect_snapshot_output(s %*% ns)
  expect_snapshot_output(s * 3)

  expect_snapshot_output(ml * s)
  expect_snapshot_output(s * x)
  expect_snapshot_output(s * ns)


  expect_snapshot_output(s / x)
  expect_snapshot_output(s / ns)

  expect_snapshot_output(s - ml)
  expect_snapshot_output(x - s)
  expect_snapshot_output(ns - s)
  expect_snapshot_output(3 - s)

  expect_snapshot_output(ml - s)
  expect_snapshot_output(s - x)
  expect_snapshot_output(s - ns)
  expect_snapshot_output(s - 3)
  expect_snapshot_output(s - s)

  expect_snapshot_output(s + ml)
  expect_snapshot_output(x + s)
  expect_snapshot_output(ns + s)
  expect_snapshot_output(3 + s)

  expect_snapshot_output(ml + s)
  expect_snapshot_output(s + x)
  expect_snapshot_output(s + ns)
  expect_snapshot_output(s + 3)
  expect_snapshot_output(s + s)





  expect_snapshot_output(s * Diagonal(10))
  expect_snapshot_output(Diagonal(10) *  s)
  expect_snapshot_output(s - Diagonal(10))
  expect_snapshot_output(s + Diagonal(10))


  expect_snapshot_output(norm(s, "f"))
  expect_snapshot_output(norm(s, "F"))
  expect_snapshot_output(norm(s, "2"))


  expect_snapshot_output(rowSums(s))
  expect_snapshot_output(colSums(s))
  expect_snapshot_output(colMeans(s))
  expect_snapshot_output(rowMeans(s))
  expect_snapshot_output(sum(s))
  expect_snapshot_output(mean(s))

  expect_snapshot_output(s[1:3, 1:3])



  expect_snapshot_output(s[1:3, 1:3, drop = FALSE])

  expect_snapshot_output(s[1:3, ])
  expect_snapshot_output(s[ , 1:3])

  expect_snapshot_output(s0 <- s[numeric(),])

  expect_snapshot_output(s0[, 1:3])

  expect_snapshot_output(s0 <- s[,numeric()])

  expect_snapshot_output(s0[1:3, ])


  expect_snapshot_output(s[1:3, 1:10 < 5])
  expect_snapshot_output(s[1:10 < 5, 1:3])

  expect_snapshot_output(s[, 1:10 < 5])
  expect_snapshot_output(s[1:10 < 5, ])
  expect_snapshot_output(s[1:10 < 5, 1:10 < 3])


  expect_snapshot_output(s[matrix(c(1,2,3,4), 2)])

  expect_snapshot_output(s[])

  expect_snapshot_output(dim(s))
  expect_snapshot_output(t(s))
  expect_snapshot_output(diag(s))
  expect_snapshot_output(as(s, "dgeMatrix"))
})

test_that("error on no b",
  {
    expect_error(
      splr(x, ns),
      "please provide an already factorized low-rank matrix.*"
    )
  }
)


test_that("error on wrong dim for a",
  {
    expect_error(
      splr(x, a[1:4,], b),
      "number of rows of x not equal to number of rows of a.*"
    )
  }
)

test_that("error on wrong dim for b",
  {
    expect_error(
      splr(x, a, b[1:4,]),
      "number of columns of x not equal to number of rows of b.*"
    )
  }
)

test_that("error on wrong dim",
  {
    expect_error(
      splr(x, a[,1:2], b),
      "number of columns of a not equal to number of columns of b.*"
    )
  }
)


test_that("error on non-atomic numeric",
  {
    expect_error(
      s + c(1, 3),
      "Can only add length 1 numerics to splrmatrix"
    )
  }
)



test_that("warning on drop",
  {
    expect_warning(
      s[1:3, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another class first"
    )
  }
)


test_that("warning on drop",
  {
    expect_warning(
      s[, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another class first"
    )
  }
)


test_that("warning on drop",
  {
    expect_warning(
      s[1:3, , drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another class first"
    )
  }
)

test_that("warning on drop",
  {
    expect_warning(
      s[1:3, 1:10 < 3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another class first"
    )
  }
)

test_that("warning on drop",
  {
    expect_warning(
      s[1:10 < 3, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another class first"
    )
  }
)

test_that(
  "subtraction/negation works",
  {
    t <- s
    expect_lt(sum(abs(as.matrix(t - s))), 1e-12)
  }
)