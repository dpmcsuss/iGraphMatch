# Function to compute the seed to non-seed portion of the gradient
get_s_to_ns <- function(A, B, seeds, nonseeds,
    perm = seq_along(nonseeds$B)) {

  nns <- nrow(nonseeds)
  ns <- nrow(seeds)

  # permute if needed
  pmat <- Matrix::Diagonal(nns, )[perm, ]

  Asn <- A[seeds$A,nonseeds$A]
  Ans <- A[nonseeds$A,seeds$A]

  Bsn <- B[seeds$B,nonseeds$B] %*% t(pmat)
  Bns <- pmat %*% B[nonseeds$B,seeds$B]

  tcrossprod(Ans, Bns) + crossprod(Asn, Bsn)
}