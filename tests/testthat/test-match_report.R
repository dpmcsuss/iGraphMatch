context("hard-seeding Frank Wolfe with incorrect seeds")

# sample pair of graphs w. 10 vertices
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- matrix_list(list(g$graph1[], g$graph1[]))
B <- matrix_list(list(g$graph2[], g$graph2[]))

mm <- graph_match_FW(A, B, seeds = NULL, start = "bari")
match_report(mm, A, B)

match_report(mm, A, B, 1:10)
match_plot_igraph(g$graph1, g$graph2, mm)
match_plot_matrix(A[[1]], B[[1]], mm)

ma <- matched_adjs(mm, A, B)
test_that(
  "matched_adjs permutes",
  {
    expect_equal(
      ma$B_m,
      B[mm$corr$corr_B, mm$corr$corr_B]
    )
  }
)
