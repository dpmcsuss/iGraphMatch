library(igraph)
library(iGraphMatch)

n1 <- 10
n2 <- 15

g1 <- sample_gnm(10, 20)
g2 <- sample_gnm(15, 30)

# square matrices
sim <- matrix(runif(n1^2), n1)
m <- graph_match_FW(g1, g1, similarity = sim)
m <- graph_match_FW(g1, g1, seeds = 1:3, similarity = sim)
m <- graph_match_FW(g1, g1, seeds = 1:3, similarity = sim[4:10,4:10])


# non square
sim <- matrix(runif(n1*n2), n1)
m <- graph_match_FW(g1, g2, similarity = sim)
m <- graph_match_FW(g1, g2, seeds = 1:3, similarity = sim)
m <- graph_match_FW(g1, g2, seeds = 1:3, similarity = sim[4:n1, 4:n2])

# the other way
m <- graph_match_FW(g2, g1, similarity = t(sim))
m <- graph_match_FW(g2, g1, seeds = 1:3, similarity = t(sim))
m <- graph_match_FW(g2, g1, seeds = 1:3, similarity = t(sim[4:n1, 4:n2]))

# square sim but non square gradient
sim <- matrix(runif(n2*n2), n2)
m <- graph_match_FW(g1, g2, similarity = sim)
m <- graph_match_FW(g1, g2, seeds = 1:3, similarity = sim)
m <- graph_match_FW(g1, g2, seeds = 1:3, similarity = sim[4:n2, 4:n2])


