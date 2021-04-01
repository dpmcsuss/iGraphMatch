

set.seed(123)

g1 <- igraph::make_lattice(10)
g2 <- igraph::sample_gnm(20, 40)

m <- gm_indefinite(g1, g2, start = "bari")


m

summary(m, g1, g2)

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


plot(g1, g2, m)
plot(g1[], g2[], m)
