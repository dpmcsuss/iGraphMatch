# Function to compute the seed to non-seed portion of the gradient
get_s_to_ns <- function(A, B, seeds, nonseeds,
    perm = seq_along(nonseeds$B)) {

  nns <- sapply(nonseeds, length)
  ns <- nrow(seeds)

  # permute if needed
  pmat <- Matrix::Diagonal(nns[2], nns[2])[perm, ]

  Asn <- A[seeds$A,nonseeds$A]
  Ans <- A[nonseeds$A,seeds$A]

  Bsn <- B[seeds$B,nonseeds$B] %*% t(pmat)
  Bns <- pmat %*% B[nonseeds$B,seeds$B]

  ml_sum(tcrossprod(Ans, Bns) + crossprod(Asn, Bsn))
}