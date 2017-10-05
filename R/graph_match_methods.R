#' @title Graph Match Methods
#'
#' @description Match two given graphs, returns a list of graph matching results,
#' including match correspondence vector of \eqn{G_2} with respect to \eqn{G_1},
#' doubly stochastic matrix and permutation matrix.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param seeds A logical vector. \code{TRUE} indicates the corresponding
#' vertex is a seed.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#' character value like "bari" or "convex" to initialize the starting matrix.
#' @param max_iter An integer. Maximum iteration time.
#' @param tol A number. Tolerance of edge disagreements.
#'
#' @rdname graph_match_methods
#'
#' @return \code{graph_match_FW} returns a list of graph matching results,
#' including match correspondence vector of \eqn{G_2} with respect to \eqn{G_1}
#' named \code{corr}, doubly stochastic matrix named \code{D} and permutation
#' matrix named \code{P} based on Frank-Wolfe methodology.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, rho =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 with no seeds
#' graph_match_FW(g1, g2)
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' seeds <- 1:50 <= 10
#' graph_match_FW(g1, g2, seeds, start = "bari")
#'
#' @export
#'
graph_match_FW <- function(A, B, seeds = NULL, start = "convex", max_iter = 100){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  # Add support for graphs with different orders
  totv1<-ncol(A)
  totv2<-ncol(B)
  if(totv1>totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv1-totv2
    for (j in 1:diff){B<-cbind(rbind(B,0),0)}
  }else if(totv1<totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv2-totv1
    for (j in 1:diff){A<-cbind(rbind(A,0),0)}
  }
  nv <- nrow(A)

  if(length(seeds)==1){
    # using this option is not recommended
    seeds <- 1:seeds
  }

  if(length(seeds)<nv){
    temp <- seeds
    seeds <- rep(FALSE,nv)
    seeds[temp]<- TRUE
  }else{
    seeds <- (seeds>0)
  }
  nonseeds <- !seeds

  ns <- sum(seeds)
  nn <- nv-ns

  Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  Ans <- A[nonseeds,seeds]

  Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  Bns <- B[nonseeds,seeds]

  if(grepl("atrix",class(start))){
    P <- start
  } else if(start=="bari"){
    P <- matrix(1/nn,nn,nn)
  } else if(start=="convex"){
    P <- graph_match_convex(A,B,seeds)$D[nonseeds,nonseeds]

  } else{ # Assuming start is an nn x nn doubly stochastic matrix
    P <- start
  }

  iter <- 0
  toggle <- TRUE
  tol <- 0.05

  # seed to non-seed info
  s_to_ns <- Ans %*% Matrix::t(Bns) + Matrix::t(Asn) %*% Bsn

  while(toggle && iter < max_iter){
    iter <- iter + 1

    # non-seed to non-seed info
    tAnn_P_Bnn <- Matrix::t(Ann) %*% P %*% Bnn

    Grad <- s_to_ns + Ann %*% P %*% Matrix::t(Bnn) + tAnn_P_Bnn
    Grad <- round(nn^2*(Grad-min(Grad)))

    ind <- as.vector(clue::solve_LSAP(as.matrix(Grad), maximum = TRUE))
    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)
    Pdir <- Pdir[ind, ]
    ns_Pdir_ns <- Matrix::t(Ann)[, order(ind)] %*% Bnn
    c <- sum(tAnn_P_Bnn * P)
    d <- sum(ns_Pdir_ns * P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(ns_Pdir_ns[ind2])
    u <- sum(P * (s_to_ns))
    v <- sum((s_to_ns)[ind2])
    if (c - d + e == 0 && d - 2 * e + u - v == 0) {
      alpha <- 0
    }
    else {
      alpha <- -(d - 2 * e + u - v)/(2 * (c - d + e))
    }
    f0 <- 0
    f1 <- c - e + u - v
    falpha <- (c - d + e) * alpha^2 + (d - 2 * e + u - v) *
      alpha
    if (alpha < tol && alpha > 0 && falpha > f0 && falpha >
        f1) {
      P <- alpha * P + (1 - alpha) * Pdir
    }
    else if (f0 > f1) {
      P <- Pdir
    }
    else {
      toggle <- F
    }
  }

  D_ns <- P
  corr_ns <- as.vector(clue::solve_LSAP(as.matrix(round(P*nn^2)), maximum = TRUE))
  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr,]
  D <- P
  D[nonseeds,nonseeds] <- D_ns
  list(corr = corr, P = P, D = D)
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_convex} returns a list of graph matching results,
#' including match correspondence vector of \eqn{G_2} with respect to \eqn{G_1}
#' named \code{corr}, doubly stochastic matrix named \code{D} and permutation
#' matrix named \code{P} based on convex relaxation method for undirected graphs.
#'
#' @examples
#' graph_match_convex(g1, g2, seeds)
#'
#' @export
#'
#'
graph_match_convex <- function(A, B, seeds = NULL, start = "bari", max_iter = 100, tol = 1e-5){

  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  # Add support for graphs with different orders ?
  nv <- nrow(A)
  if(length(seeds)==1){
    seeds <- 1:seeds
  }
  if(length(seeds)<nv){
    temp <- seeds
    seeds <- rep(FALSE,nv)
    seeds[temp]<- TRUE
  }else{
    seeds <- (seeds>0)
  }
  nonseeds <- !seeds

  ns <- sum(seeds)
  nn <- nv-ns

  Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  Ans <- A[nonseeds,seeds]

  Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  Bns <- B[nonseeds,seeds]

  tol0<-1
  if(start=="bari"){
    P <- matrix(1/nn,nn,nn)
  } else{ # Assuming start is an nn x nn doubly stochastic matrix
    P <- start
  }
  iter<-0
  toggle <- TRUE

  AtA <- t(Asn)%*%Asn + t(Ann)%*%Ann
  BBt <- Bns%*%t(Bns)+Bnn%*%t(Bnn)

  ABns_sn <- Ans%*%t(Bns) + t(Asn)%*%Bsn
  f <- sum((Ann %*% P - P%*% Bnn)^2)

  while(toggle && iter<max_iter){
    f_old <- f
    iter<-iter+1

    Grad<- AtA%*%P + P%*%BBt - ABns_sn - t(Ann)%*%P%*%Bnn - Ann%*%P%*%t(Bnn)
    # print("asdf")

    Grad <- round(as.matrix(nn^2*(Grad-min(Grad))))
    corr <- as.vector(clue::solve_LSAP(Grad))
    Pdir <- Matrix::Diagonal(nn)[corr,]


    # C <- rbind(Ann,Asn) %*% (P-Pdir) + t((Pdir-P) %*% cbind(Bnn,Bns))
    Cnn <- Ann %*% (P-Pdir) - (P-Pdir) %*% Bnn
    Dnn <- Ann %*% Pdir - Pdir %*% Bnn

    if(ns > 0){
      Cns <- -(P-Pdir) %*% Bns
      Csn <- Asn %*% (P-Pdir)

      Dns <- Ans - Pdir %*% Bns
      Dsn <- Asn %*% Pdir - Bsn
    }else{
      Dns <- Dsn <-Cns <- Csn <- 0
    }

    aq <- sum(Cnn^2)+sum(Cns^2)+sum(Csn^2)
    bq <- sum(Cnn*Dnn)+sum(Cns*Dns)+sum(Csn*Dsn)
    aopt <- -bq/aq

    P_new <- aopt*P+(1-aopt)*Pdir;
    f <- sum((Ann %*% P_new - P_new %*% Bnn)^2)

    f_diff <- abs(f-f_old)
    P_diff <- sum(abs(P-P_new))
    P <- P_new

    toggle <- f_diff > tol0 && f > tol && P_diff > tol0
  }

  D_ns <- P
  corr_ns <- as.vector(clue::solve_LSAP(as.matrix(round(P*nn^2)), maximum = TRUE))
  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr,]
  D <- P
  D[nonseeds,nonseeds] <- D_ns
  list(corr = corr, P = P, D = D)
}
#'
#' @return \code{graph_match_convex_directed} returns graph matching results based
#' on convex relaxation method for directed graphs.
#'
#' @examples
#' graph_match_convex_directed(g1, g2, seeds)
#'
#'
#'
graph_match_convex_directed <- function(A,B,seeds=NULL,start="bari",max_iter=100, tol2=1e-5){

  print("Warning, this doesn't work as expected. Need to think more.")
  A <- A[]
  B <- B[]

  # Add support for graphs with different orders ?
  nv <- nrow(A)
  if(length(seeds)==1){
    seeds <- 1:seeds
  }
  if(length(seeds)<nv){
    temp <- seeds
    seeds <- rep(FALSE,nv)
    seeds[temp]<- TRUE
  }else{
    seeds <- (seeds>0)
  }
  nonseeds <- !seeds

  ns <- sum(seeds)
  nn <- nv-ns

  Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  Ans <- A[nonseeds,seeds]

  Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  Bns <- B[nonseeds,seeds]

  tol<-1
  if(start=="bari"){
    P <- matrix(1/nn,nn,nn)
  } else{ # Assuming start is an nn x nn doubly stochastic matrix
    P <- start
  }
  iter<-0
  toggle <- TRUE

  AtA <- t(Asn)%*%Asn + t(Ann)%*%Ann
  BtB <- t(Bsn)%*%Bsn + t(Bnn)%*%Bnn
  AAt <- Ans%*%t(Ans) + Ann%*%t(Ann)
  BBt <- Bns%*%t(Bns) + Bnn%*%t(Bnn)

  ABns_sn <- Ans%*%t(Bns) + t(Asn)%*%Bsn
  f <- sum((Ann%*%P - P%*%Bnn)^2)+sum((t(Ann)%*%P - P%*%t(Bnn))^2)



  while(toggle && iter<max_iter){
    f_old <- f
    iter <- iter+1

    tAnn_P_Bnn <- t(Ann)%*%P%*%Bnn + Ann%*%P%*%t(Bnn)
    Grad<- AtA%*%P + P%*%BBt - ABns_sn - tAnn_P_Bnn +
      BtB%*%P + P%*%AAt - t(ABns_sn) - t(tAnn_P_Bnn);

    Grad <- round(as.matrix(nn^2*(Grad-min(Grad))))
    corr <- as.vector(solve_LSAP(Grad))
    Pdir <- Matrix::Diagonal(nn)[corr,]


    # C <- rbind(Ann,Asn) %*% (P-Pdir) + t((Pdir-P) %*% cbind(Bnn,Bns))
    Cnn <- Ann %*% (P-Pdir) - (P-Pdir) %*% Bnn
    tCnn <- t(Ann) %*% t(P-Pdir) - t(P-Pdir) %*% t(Bnn)
    Dnn <- Ann %*% Pdir - Pdir %*% Bnn
    tDnn <- t(Ann) %*% t(Pdir) - t(Pdir) %*% t(Bnn)

    if(ns > 0){
      Cns <- -(P-Pdir) %*% Bns
      tCns <- -t(P-Pdir) %*% t(Bsn)
      Csn <- Asn %*% (P-Pdir)
      tCsn <- t(Ans) %*% t(P-Pdir)

      Dns <- Ans - Pdir %*% Bns
      tDns <- t(Asn) - t(Pdir) %*% t(Bsn)
      Dsn <- Asn %*% Pdir - Bsn
      tDsn <- t(Ans) %*% t(Pdir) - t(Bns)
    }else{
      Dns <- Dsn <-Cns <- Csn <- 0
    }

    aq <- sum(Cnn^2+tCnn^2)+sum(Cns^2+tCns^2)+sum(Csn^2+tCsn^2)
    bq <- sum(Cnn*Dnn+tCnn*tDnn)+sum(Cns*Dns+tCns*tDns)+sum(Csn*Dsn+tCsn*tDsn)
    aopt <- -bq/aq

    P_new <- aopt*P+(1-aopt)*Pdir;
    f <- sum((Ann %*% P_new - P_new %*% Bnn)^2)+sum((t(Ann) %*% P_new - P_new %*% t(Bnn))^2)

    f_diff <- abs(f-f_old)
    P_diff <- sum(abs(P-P_new))
    P <- P_new

    toggle <- f_diff > tol && f > tol2 && P_diff > tol
  }

  D_ns <- P
  corr_ns <- unclass(solve_LSAP(as.matrix(round(P*nn^2)), maximum = TRUE))
  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr,]
  D <- P
  D[nonseeds,nonseeds] <- D_ns
  list(corr = corr, P = P, D = D)
}
