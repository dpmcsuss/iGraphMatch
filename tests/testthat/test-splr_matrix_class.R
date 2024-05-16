
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

  expect_no_error(s <- splr(x, Matrix(a), Matrix(b)))
  expect_no_error(sp <- splr_to_sparse(s))
  expect_no_error(m <- as.matrix(s))

  expect_no_error(ss <- splr_sparse_plus_constant(x, 1))
  expect_no_error(dm <- as(s, "dMatrix"))
  # try({
  #   sink("/dev/null")
  # })
  expect_no_error(print(s))
  expect_no_error(show(s))
  # expect_no_error(as.character(s))

  expect_no_error(str(s))
  # sink()


  expect_no_error(length(s))


  expect_no_error(s %*% s)
  expect_no_error(s * s)
  expect_no_condition(s * s, message = "ambiguousMethodSelection")

  expect_no_error(ml <- matrix_list(list(ns, x)))

  expect_no_error(s %*% ml)
  expect_no_error(x %*% s)
  expect_no_error(ns %*% s)
  expect_no_error(3 * s)

  expect_no_error(ml %*% s)
  expect_no_error(s %*% x)
  expect_no_error(s %*% ns)
  expect_no_error(s * 3)

  expect_no_error(ml * s)
  expect_no_error(s * x)
  expect_no_error(s * ns)


  expect_no_error(s / x)
  expect_no_error(s / ns)

  expect_no_error(s - ml)
  expect_no_error(x - s)
  expect_no_error(ns - s)
  expect_no_error(3 - s)

  expect_no_error(ml - s)
  expect_no_error(s - x)
  expect_no_error(s - ns)
  expect_no_error(s - 3)
  expect_no_error(s - s)

  expect_no_error(s + ml)
  expect_no_error(x + s)
  expect_no_error(ns + s)
  expect_no_error(3 + s)

  expect_no_error(ml + s)
  expect_no_error(s + x)
  expect_no_error(s + ns)
  expect_no_error(s + 3)
  expect_no_error(s + s)





  expect_no_error(s * Diagonal(10))
  expect_no_error(Diagonal(10) *  s)
  expect_no_error(s - Diagonal(10))
  expect_no_error(s + Diagonal(10))


  expect_no_error(norm(s, "f"))
  expect_no_error(norm(s, "F"))
  expect_no_error(norm(s, "2"))


  expect_no_error(rowSums(s))
  expect_no_error(colSums(s))
  expect_no_error(colMeans(s))
  expect_no_error(rowMeans(s))
  expect_no_error(sum(s))
  expect_no_error(mean(s))

  expect_no_error(s[1:3, 1:3])



  expect_no_error(s[1:3, 1:3, drop = FALSE])

  expect_no_error(s[1:3, ])
  expect_no_error(s[ , 1:3])

  expect_no_error(s0 <- s[numeric(),])

  expect_no_error(s0[, 1:3])

  expect_no_error(s0 <- s[,numeric()])

  expect_no_error(s0[1:3, ])


  expect_no_error(s[1:3, 1:10 < 5])
  expect_no_error(s[1:10 < 5, 1:3])

  expect_no_error(s[, 1:10 < 5])
  expect_no_error(s[1:10 < 5, ])
  expect_no_error(s[1:10 < 5, 1:10 < 3])


  expect_no_error(s[matrix(c(1,2,3,4), 2)])

  expect_no_error(s[])

  expect_no_error(dim(s))
  expect_no_error(t(s))
  expect_no_error(diag(s))
  expect_no_error(as(s, "dgeMatrix"))
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