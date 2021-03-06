% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start.R
\name{bari_start}
\alias{bari_start}
\alias{rds_sinkhorn_start}
\alias{rds_perm_bari_start}
\alias{rds_from_sim_start}
\title{Start matrix initialization}
\usage{
bari_start(nns, ns = 0, soft_seeds = NULL)

rds_sinkhorn_start(nns, ns = 0, soft_seeds = NULL, distribution = "runif")

rds_perm_bari_start(nns, ns = 0, soft_seeds = NULL, g = 1, is_splr = TRUE)

rds_from_sim_start(nns, ns = 0, soft_seeds = NULL, sim)
}
\arguments{
\item{nns}{An integer. Number of non-seeds.}

\item{ns}{An integer. Number of hard seeds.}

\item{soft_seeds}{A vector, a matrix or a data frame. If there is no error in soft
seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
seeds, input in the form of a matrix or a data frame, with the first column being the
indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
that if there are seeds in graphs, seeds should be put before non-seeds.}

\item{distribution}{A character. Specify the distribution from which the random doubly stochastic
matrix is sampled. Should input the name of the function for generating random deviates from that
distribution.}

\item{g}{A number. Specified in the range of [0, 1] to set weights to random permutation matrix and
barycenter matrix.}

\item{is_splr}{should we return a splr matrix? (default = TRUE)}

\item{sim}{nns x nns non-negative matrix.}
}
\value{
\code{bari_start} returns a \code{nns-by-nns} matrix with 1's corresponding to the
adaptive seeds and being bari-centered at other places.

\code{rds_sinkhorn_start} returns a \code{nns-by-nns} doubly stochastic matrix
with 1's corresponding to adaptive seeds.

\code{rds_perm_bari} returns a \code{nns-by-nns} doubly stochastic matrix
with 1's corresponding to adaptive seeds.

\code{rds_from_sim_start} returns a doubly
 stochastic Matrix given by sinkhorn algorithm applied to
 a matrix of iid log-normal with mu=sim. Note,
 this ignores soft seeds.
}
\description{
initialize the start matrix for graph matching iteration.
}
\examples{
## Case without soft seeds
bari_start(3)

## Case with correct soft seeds and input is a vector
bari_start(nns=5, ns=3, soft_seeds=c(5, 7, 8))

## Case with erroneous soft seeds and the input is a matrix
bari_start(nns=5, soft_seeds=matrix(c(2, 4, 2, 3), nrow=2))

## Case without soft seeds
rds_sinkhorn_start(5)

## Case with soft seeds and the input is a data frame
rds_sinkhorn_start(nns=5,
   soft_seeds = as.data.frame(matrix(c(2, 4, 2, 3), nrow=2)),
   distribution = "rnorm")

## Case without soft seeds
rds_perm_bari_start(nns=5)

## Case with soft seeds and the input is a data frame
rds_perm_bari_start(nns=5, ns=0, soft_seeds=as.data.frame(matrix(c(2, 4, 2, 3), nrow=2)))

sim <- Matrix::rsparsematrix(10, 10, .4,
 rand.x = function(n) rep(1,n))
start_sparse <- rds_from_sim_start(10, sim = sim)
start_dense <- rds_from_sim_start(10, sim = as.matrix(sim))

}
