

# sample pair of graphs w. 10 vertices
set.seed(123)

cgnp_pair <- sample_correlated_gnp_pair(n = 1000, corr =  0.3, p =  0.5)
g1 <- cgnp_pair$graph1[][1:20, 1:20]
g2 <- cgnp_pair$graph2
m <- gm(g1, g2, method = "indefinite", lap_method = "rect")

test_that(
  "rect lap works with indefinite",
  {
    expect_equal(
      m$corr,
      structure(list(
        corr_A = 1:20, 
        corr_B = c(805, 107, 393, 565, 378, 657, 968, 934, 
          661, 454, 212, 119, 424, 797, 440, 592, 811, 
          984, 556, 585)),
        class = "data.frame", row.names = c(NA, -20L)
      )
    )
  }
)