

set.seed(123)

g1 <- igraph::make_lattice(10)
g2 <- igraph::sample_gnm(20, 40)

m <- gm_indefinite(g1, g2, start = "bari")


plot(g1, g2, m)