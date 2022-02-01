


test_that("graphMatch class functions",
  {
  set.seed(123)

  g1 <- igraph::make_lattice(10)
  g2 <- igraph::sample_gnm(200, 40)


  suppressWarnings(m <- gm(g1, g2, start = "bari", method = "indefinite"))

  as.character(m)
  str(m)

  m

  as(m, "Matrix")
  as(m, "data.frame")

  expect_snapshot_output(print(summary(m, g1, g2, true_label = 1:10)))

  m %*% g2
  m %*% g2[]


  norm(g1[] - m %*% g2[], "f")

  g1 %*% m
  g1[] %*% m


  m[1:3]
  m[,1]
  m[,2]

  m$seeds

  m$corr
  m$corr_A

  m[]

  as.data.frame(m)

  plot(g1, g2)
  plot(g1[], g2[])

  plot(g1, g2, m)
  plot(g1[], g2[], m)

  t(m)

  m$corr_A
  m$corr_B
  m$not_in_there
})

test_that(
  "as.graphMatch works",
  {
    expect_snapshot(as.graphMatch(data.frame(a = 1:5, b = 1:5)))
    expect_snapshot(as.graphMatch(data.frame(corr_A = 1:5, corr_B = 1:5)))
    expect_snapshot(as.graphMatch(1:10))
  }
)