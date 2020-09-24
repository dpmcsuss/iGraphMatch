# Function to compute the seed to non-seed portion of the gradient
get_s_to_ns <- function(Alist, Blist, seeds, nonseeds,
    perm = seq(sum(seeds))) {

  # NEED TO CHANGE ???
  nns <- nrow(nonseeds)
  ns <- nrow(seeds)

  # permute if needed
  pmat <- Matrix::Diagonal(nns, )[perm, ]

  s_to_ns <- function(A,B){
    Asn <- A[seeds$A,nonseeds$A]
    Ans <- A[nonseeds$A,seeds$A]

    Bsn <- B[seeds$B,nonseeds$B] %*% t(pmat)
    Bns <- pmat %*% B[nonseeds$B,seeds$B]

    tcrossprod(Ans, Bns) + crossprod(Asn, Bsn)
  }

  if (!is(Alist, "list")){
    return(s_to_ns(Alist, Blist))
  }

  nc <- length(Alist)
  s2ns <- Matrix(0, nrow = nns, ncol = nns)
  for (ch in 1:nc){
    s2ns <- s2ns + s_to_ns(Alist[[ch]], Blist[[ch]])
    gc()
  }
  s2ns
}