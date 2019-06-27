#' @title Graph Match Methods
#'
#' @description Match two given graphs, returns a list of graph matching
#'   results, including matching correspondence vector of \eqn{G_2} with respect
#'   to \eqn{G_1}, doubly stochastic matrix and permutation matrix.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   there is no error in seeds input can be a vector of seed indices in
#'   \eqn{G_1}. Or if there exists error in seeds, input in the form of a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#'   character value like "bari" or "convex" to initialize the starting matrix.
#' @param max_iter An integer. Maximum iteration time.
#' @param tol A number. Tolerance of edge disagreements.
#' @param r A number. Threshold of neighboring pair scores.
#' @param max_iter A number. Maximum number of replacing matches equals to
#'   max_iter times number of total vertices of \eqn{G_1}.
#' @param alpha A number betwen 0 and 1. Bigger alpha means putting more importance
#'   on the information in network topology over other information such as
#'   similarity scores
#' @param method A character. Choice of method to extract mapping from score matrix,
#'   including greedy method and the Hungarian algorithm.
#'
#' @rdname graph_match_methods
#'   
#' @return \code{graph_match_FW} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B} and the number of seeds. 
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 with no seeds
#' graph_match_FW(g1, g2)
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' seeds <- 1:10 <= 3
#' graph_match_FW(g1, g2, seeds, start = "bari")
#'
#' # match G_1 & G_2 with some incorrect seeds
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds)),hard_seeds)
#' graph_match_FW(g1, g2, seeds, start = "convex")
#'
#' @export
#'
graph_match_FW <- function(A, B, seeds = NULL, start = "convex", max_iter = 20){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]

  # Add support for graphs with different orders
  totv1<-ncol(A)
  totv2<-ncol(B)
  if(totv1>totv2){
    diff<-totv1-totv2
    B <- Matrix::bdiag(B[], Matrix(0,diff,diff))
  }else if(totv1<totv2){
    diff<-totv2-totv1
    A <- Matrix::bdiag(A[], Matrix(0,diff,diff))
  }
  nv <- nrow(A)

  if(is.null(seeds)){
    seeds <- rep(FALSE,nv)
    aseeds_err <- FALSE
    ns <- sum(seeds)
  } else{
    seeds_pair <- check_seeds(seeds)
    ns <- nrow(seeds_pair)

    seeds <- rep(FALSE,nv)
    seeds[seeds_pair$seed_A] <- TRUE

    # detect incorrect seeds
    seed_A <- seeds_pair$seed_A
    seed_B <- seeds_pair$seed_B
    aseeds_err <- ifelse(seed_A!=seed_B,TRUE,FALSE)
    seed_A_err <- seed_A[aseeds_err]
    seed_B_err <- seed_B[aseeds_err]

    if(sum(aseeds_err)!=0){
      B <- g2_hard_seeding(seed_A_err,seed_B_err,B)
    }
  }

  nn <- nv-ns
  nonseeds <- !seeds

  Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  Ans <- A[nonseeds,seeds]

  Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  Bns <- B[nonseeds,seeds]

  P <- init_start(start = start, nns = nn,
                  A = A, B = B, seeds = seeds)

  iter <- 0
  toggle <- TRUE

  # seed to non-seed info
  s_to_ns <- Ans %*% Matrix::t(Bns) + Matrix::t(Asn) %*% Bsn

  while(toggle && iter < max_iter){
    iter <- iter + 1

    # non-seed to non-seed info
    tAnn_P_Bnn <- Matrix::t(Ann) %*% P %*% Bnn

    Grad <- s_to_ns + Ann %*% P %*% Matrix::t(Bnn) + tAnn_P_Bnn
    Grad <- Grad-min(Grad)

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
    } else {
      alpha <- -(d - 2 * e + u - v)/(2 * (c - d + e))
    }
    f0 <- 0
    f1 <- c - e + u - v
    falpha <- (c - d + e) * alpha^2 + (d - 2 * e + u - v) *
      alpha

    if (alpha < 1 && alpha > 0 &&
        falpha > f0 && falpha > f1) {
      P <- alpha * P + (1 - alpha) * Pdir
    } else if (f0 > f1) {
      P <- Pdir
    } else {
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

  # fix match results if there are incorrect seeds
  if(sum(aseeds_err)!=0){
    corr <- fix_hard_corr(seed_A_err,seed_B_err,corr)
    P <- Matrix::Diagonal(nv)[corr,]
    D <- fix_hard_D(seed_A_err,seed_B_err,D)
  }

  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns,
            P = P, D = D)
  z
}
# correct the order of swapping graph2 according to new seeds
swap_order <- function(aseeds_matrix){
  #aseeds_matrix: first row:added seeds index in g1, second row added seeds match
  naseeds_err <- dim(aseeds_matrix)[2]
  ninter <- 0
  ninter_new <- naseeds_err
  aseeds_match_order <- matrix( ,2, )
  aseeds_matrix_T <- aseeds_matrix

  while(ninter_new!=ninter & ninter_new>1){
    aseeds_matrix <- aseeds_matrix_T
    naseeds_err <- ninter_new
    inter_match <- rep("FALSE",times = naseeds_err)
    ninter <- ninter_new
    ninter_new <- 0
    circle_index <- 0
    k <- 1

    for(i in 1:naseeds_err){
      # eliminate circle of two vertices
      if(aseeds_matrix[2,i] %in% aseeds_matrix[1,]){
        index <- which(aseeds_matrix[1,]==aseeds_matrix[2,i])
        if(aseeds_matrix[1,i]==aseeds_matrix[2,index]){
          aseeds_matrix[1,i] <- 0
          circle_index[k] <- i
          k <- k+1
        } else{
          inter_match[i] <- "TRUE"
          ninter_new <- ninter_new+1
        }
      }
    }

    if(circle_index[1]!=0){
      aseeds_matrix <- aseeds_matrix[,-circle_index]
      inter_match <- inter_match[-circle_index]
    }

    if(length(which(inter_match=="TRUE"))>=1){
      aseeds_matrix <- as.matrix(aseeds_matrix)
      aseeds_matrix_T <- aseeds_matrix[,which(inter_match=="TRUE")]
    }
    if(length(which(inter_match=="FALSE"))>=1){
      aseeds_matrix <- as.matrix(aseeds_matrix)
      aseeds_matrix_F <- aseeds_matrix[,which(inter_match=="FALSE")]
      aseeds_match_order <- cbind(aseeds_matrix_F, aseeds_match_order)
    }

  }

  #end with circle: only consider one circle (circle with more than three vertices) case
  if(length(which(inter_match=="TRUE"))>1){
    aseeds_matrix_T <- aseeds_matrix_T[,-1]
    aseeds_match_order <- cbind(aseeds_matrix_T,aseeds_match_order)
  } else if(length(which(inter_match=="TRUE"))==1){
    aseeds_match_order <- cbind(aseeds_matrix_T,aseeds_match_order)
  }

  aseeds_match_order[,-dim(aseeds_match_order)[2]]
}
# swap columns and rows of G_2 according to hard seeds
g2_hard_seeding <- function(seed_g1_err, seed_g2_err, g2){
  aseeds_matrix <- matrix(c(seed_g1_err,seed_g2_err),nrow=2,byrow = TRUE)
  if(length(seed_g1_err)>1)
  {
    swap <- swap_order(aseeds_matrix)
    swap <- as.matrix(swap)
    seed_g1_err <- swap[1,]
    seed_g2_err <- swap[2,]
  }

  # swap columns of g2
  nv <- nrow(g2)
  g2_new_real <- 1:nv
  for (i in 1:length(seed_g1_err)) {
    g2_new_real[c(seed_g1_err[i],seed_g2_err[i])] <-
      g2_new_real[c(seed_g2_err[i],seed_g1_err[i])]
  }

  g2_new <- g2[g2_new_real,g2_new_real]
  g2_new
}
# returns the true correspondence between G_1 and G_2 for hard seeding
fix_hard_corr <- function(seed_g1_err, seed_g2_err, corr_hard){
  aseeds_matrix <- matrix(c(seed_g1_err,seed_g2_err),nrow=2,byrow = TRUE)
  if(length(seed_g1_err)>1)
  {
    swap <- swap_order(aseeds_matrix)
    swap <- as.matrix(swap)
    seed_g1_err <- swap[1,]
    seed_g2_err <- swap[2,]
  }

  nv <- length(corr_hard)
  g2_new_real <- 1:nv
  for (i in 1:length(seed_g1_err)) {
    g2_new_real[c(seed_g1_err[i],seed_g2_err[i])] <-
      g2_new_real[c(seed_g2_err[i],seed_g1_err[i])]
  }

  corr_hard <- g2_new_real[corr_hard]
  corr_hard
}
# returns the true doubly stochastic matrix D and true permutation matrix for hard seeding
fix_hard_D <- function(seed_g1_err, seed_g2_err, D){
  aseeds_matrix <- matrix(c(seed_g1_err,seed_g2_err),nrow=2,byrow = TRUE)
  if(length(seed_g1_err)>1)
  {
    swap <- swap_order(aseeds_matrix)
    swap <- as.matrix(swap)
    seed_g1_err <- swap[1,]
    seed_g2_err <- swap[2,]
  }

  nv <- nrow(D)
  g2_new_real <- 1:nv
  for (i in 1:length(seed_g1_err)) {
    g2_new_real[c(seed_g1_err[i],seed_g2_err[i])] <-
      g2_new_real[c(seed_g2_err[i],seed_g1_err[i])]
  }

  D <- D[,g2_new_real]
  D
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_convex} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B} and the number of seeds. 
#'
#' @examples
#' seeds <- 1:10 <= 3
#' graph_match_convex(g1, g2, seeds)
#'
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds)),hard_seeds)
#' graph_match_convex(g1, g2, seeds)
#'
#' @export
#'
#'
graph_match_convex <- function(A, B, seeds = NULL, start = "bari", max_iter = 100, tol = 1e-5){

  A <- A[]
  B <- B[]

  # Add support for graphs with different orders ?
  nv <- nrow(A)
  if(is.null(seeds)){
    seeds <- rep(FALSE,nv)
    aseeds_err <- FALSE
    ns <- sum(seeds)
  } else{
    seeds_pair <- check_seeds(seeds)
    ns <- nrow(seeds_pair)

    seeds <- rep(FALSE,nv)
    seeds[seeds_pair$seed_A] <- TRUE

    # detect incorrect seeds
    seed_A <- seeds_pair$seed_A
    seed_B <- seeds_pair$seed_B
    aseeds_err <- ifelse(seed_A!=seed_B,TRUE,FALSE)
    seed_A_err <- seed_A[aseeds_err]
    seed_B_err <- seed_B[aseeds_err]

    if(sum(aseeds_err)!=0){
      B <- g2_hard_seeding(seed_A_err,seed_B_err,B)
    }
  }

  nn <- nv-ns
  nonseeds <- !seeds

  Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  Ans <- A[nonseeds,seeds]

  Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  Bns <- B[nonseeds,seeds]

  tol0 <- 1
  P <- init_start(start = start, nns = nn)
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

  # fix match results if there are incorrect seeds
  if(sum(aseeds_err)!=0){
    corr <- fix_hard_corr(seed_A_err,seed_B_err,corr)
    P <- Matrix::Diagonal(nv)[corr,]
    D <- fix_hard_D(seed_A_err,seed_B_err,D)
  }

  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns, 
            P = P, D = D)
  z
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
  
  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns)
  z
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_percolation} returns a list consists of
#'   matching correspondence of matched pairs with index of nodes in
#'   \eqn{G_1} named \code{corr_A} and index of nodes in \eqn{G_2} named
#'   \code{corr_B} returns and the order of matching for matched nodes in
#'   \eqn{G_1}.
#' @return \code{graph_match_percolation} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B}, the number of seeds and the order of nodes getting matched.
#'
#' @references L. Yartseva and M. Grossglauser (2013), \emph{On the performance
#'   of percolation graph matching}. COSN, Boston, MA, USA, pages 119–130.
#'
#' @examples
#' # match G_1 & G_2 using percolation graph matching method
#' seeds <- 1:5
#' graph_match_percolation(g1, g2, seeds, r = 2)
#'
#' @export
#'
#'
graph_match_percolation <- function(A, B, seeds, r = 2){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  n <- nrow(A)
  m <- nrow(B)
  seeds <- check_seeds(seeds) #unused seeds
  ns <- nrow(seeds)
  Z <- seeds #matched nodes
  M <- matrix(0,n,m) #marks matrix

  # mark neighbors
  for(i in 1:nrow(seeds)){
    A_adj <- which(A[seeds$seed_A[i],]>0)
    B_adj <- which(B[seeds$seed_B[i],]>0)
    mark <- outer(A[seeds$seed_A[i],A_adj], B[seeds$seed_B[i],B_adj], cal_mark)
    M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
  }
  M[seeds$seed_A,] <- -n
  M[,seeds$seed_B] <- -n

  # choose pairs with marks ge r
  while(max(M) >= r){
    max_ind <- which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind),1),]

    Z <- rbind(Z,max_ind)

    # update mark matrix
    A_adj <- which(A[max_ind[1],]>0)
    B_adj <- which(B[max_ind[2],]>0)
    mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
    M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
    M[max_ind[1],] <- -n
    M[,max_ind[2]] <- -n
  }

  if(nrow(Z) == n-1){
    all <- 1:n
    seed_A <- all[!(all %in% Z$seed_A)]
    seed_B <- all[!(all %in% Z$seed_B)]
    Z <- rbind(Z,cbind(seed_A,seed_B))
  }

  # matching result
  order <- order(Z$seed_A)
  corr <- Z[order(Z$seed_A),]
  names(corr) <- c("corr_A","corr_B")
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns, order = order)
  z
}
cal_mark <- function(x,y){
  1 - abs(x - y) / max(x, y)
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_ExpandWhenStuck} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B}, the number of seeds and the order of nodes getting matched.
#'
#' @references E. Kazemi, S. H. Hassani, and M. Grossglauser (2015),
#' \emph{Growing a graph matching from a handful of seeds}. Proc. of the VLDB
#' Endowment, 8(10):1010–1021.
#'
#' @examples
#' # match G_1 & G_2 using Expand When Stuck graph matching method
#' seeds <- 1:5
#' graph_match_ExpandWhenStuck(g1, g2, seeds, r = 2)
#'
#' @export
#'
#'
graph_match_ExpandWhenStuck <- function(A, B, seeds, r = 2){
  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  n <- nrow(A)
  m <- nrow(B)
  seeds <- check_seeds(seeds) #unused seeds
  ns <- nrow(seeds)
  Z <- seeds #matched nodes
  M <- matrix(0,n,m) #marks matrix
  M[seeds$seed_A,] <- -n*n
  M[,seeds$seed_B] <- -n*n

  # deferred percolation graph matching
  while(nrow(seeds) != 0){
    # mark neighbors
    for(i in 1:nrow(seeds)){
      A_adj <- which(A[seeds$seed_A[i],]>0)
      B_adj <- which(B[seeds$seed_B[i],]>0)
      mark <- outer(A[seeds$seed_A[i],A_adj], B[seeds$seed_B[i],B_adj], cal_mark)
      M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
    }

    # choose pairs with marks ge r
    while(max(M) >= r){
      max_ind <- which(M == max(M), arr.ind = TRUE)
      if(nrow(max_ind) != 1){
        degree_diff <- abs(rowSums(A)[max_ind[,1]]-rowSums(B)[max_ind[,2]])
        max_ind <- max_ind[which(degree_diff == min(degree_diff)),]
        if(is.vector(max_ind) == FALSE){
          max_ind <- max_ind[sample(nrow(max_ind),1),]
        }
      }

      Z <- rbind(Z, as.vector(max_ind))

      # update mark matrix
      A_adj <- which(A[max_ind[1],]>0)
      B_adj <- which(B[max_ind[2],]>0)
      mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
      M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      M[max_ind[1],] <- -n*n
      M[,max_ind[2]] <- -n*n
    }

    if(nrow(Z) == n-1){
      all <- 1:n
      seed_A <- all[!(all %in% Z$seed_A)]
      seed_B <- all[!(all %in% Z$seed_B)]
      Z <- rbind(Z,cbind(seed_A,seed_B))
      break
    }

    seeds_old <- seeds
    seeds <- which(M > 0, arr.ind = TRUE)
    seeds <- data.frame(seed_A=seeds[,1], seed_B=seeds[,2])

    if(nrow(seeds) == nrow(seeds_old)){
      if(sum(seeds == seeds_old)==2*nrow(seeds)){
        break
      }
    }

  }

  # matching result
  order <- order(Z$seed_A)
  corr <- Z[order(Z$seed_A),]
  names(corr) <- c("corr_A","corr_B")
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns, order = order)
  z
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_soft_percolation} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B}, the number of seeds and the order of nodes getting matched.
#'
#' @examples
#' # match G_1 & G_2 using soft percolation graph matching method
#' seeds <- 1:5
#' graph_match_soft_percolation(g1, g2, seeds, r = 2, max_iter = 2)
#'
#' @export
#'
#'
graph_match_soft_percolation <- function(A, B, seeds, r = 2, max_iter = 2){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  # initialization of score matrix M & MM
  n <- nrow(A)
  max_iter <- max_iter * n
  minusinf <- -max_iter
  M <- matrix(0,n,n) # actual score matrix
  seeds <- check_seeds(seeds)
  ns <- nrow(seeds)
  for(i in 1:nrow(seeds)){ # mark neighbors
    A_adj <- which(A[seeds$seed_A[i],]==1)
    B_adj <- which(B[seeds$seed_B[i],]==1)
    M[A_adj, B_adj] <- M[A_adj, B_adj] + 1
  }
  M[seeds$seed_A,] <- minusinf # hard seeds
  M[,seeds$seed_B] <- minusinf # hard seeds
  MM <- M # score matrix w. socres to matched pairs be -infinity

  # initialization for checking replacing cycle
  num <- 0 # # of removal seeds
  remove <- rbind(c(0,0),c(0,0)) # list of removed seeds
  remove_by <- rbind(c(0,0),c(0,0))
  colnames(remove) <- paste0(c("seed_A","seed_B"))
  cyc <- FALSE

  # percolate
  Z <- seeds # set of matched pairs
  ZZ <- c(0,0) # set of dominant conflict existing matches
  num <- 0 # count the number of replacing matches
  while(max(MM)>=r & num<=max_iter){
    # locate best match
    max_ind <- which(MM==max(MM), arr.ind = TRUE)
    conflict_log <- conflict_check(Z, max_ind, logical = TRUE)
    sum_conf <- sum(conflict_log)
    if(sum_conf>0 & sum_conf<length(conflict_log)){
      max_ind <- max_ind[!conflict_log,] # subset of non-conflict matches
    }
    if(!is.null(nrow(max_ind))){
      rnum <- sample(nrow(max_ind),1)
      max_ind <- max_ind[rnum,] # solve tie: give priority to non-conflict matches
    }

    # non-conflict new match
    if(sum_conf != length(conflict_log)){

      Z <- rbind(Z,max_ind)

      # correct MM caused by ZZ
      if(!is.null(nrow(ZZ))){
        MM[ZZ$seed_A,] <- M[ZZ$seed_A,]
        MM[,ZZ$seed_B] <- M[,ZZ$seed_B]
        ZZ <- c(0,0)
      }

      # update mark matrix M & MM
      A_adj <- which(A[unlist(max_ind[1]),]>0)
      B_adj <- which(B[unlist(max_ind[2]),]>0)
      M[A_adj, B_adj] <- M[A_adj, B_adj] + 1

      MM[A_adj, B_adj] <- MM[A_adj, B_adj] + 1
      MM[max_ind[1], max_ind[2]] <- minusinf

    } else{ # conflict new match: only when all the ties correspond to conflict match

      conf_row_ind <- conflict_check(Z, matrix(max_ind,1), logical = FALSE)
      conf_ind <- Z[conf_row_ind,]
      if(nrow(conf_ind==2)==2){ # conflict with two existing matches
        score1 <- M[conf_ind[1,1], conf_ind[1,2]]
        score2 <- M[conf_ind[2,1], conf_ind[2,2]]
        score <- max(score1,score2)
      } else{
        score <- M[conf_ind$seed_A,conf_ind$seed_B]
      }

      if(M[max_ind[1], max_ind[2]]>score){ #replace

        num <- num + 1

        # check cycle
        if(length(conf_row_ind)==2){
          remove <- rbind(remove,Z[conf_row_ind[1],])
          cyc_remove <- check_cycle(remove, Z[conf_row_ind[2],])
          remove <- rbind(remove, Z[conf_row_ind[2],])
          cyc_remove_by <- check_cycle(remove_by,max_ind)
          remove_by <- rbind(remove_by, max_ind)
        } else{
          cyc_remove <- check_cycle(remove, Z[conf_row_ind,])
          remove <- rbind(remove, Z[conf_row_ind,])
          cyc_remove_by <- check_cycle(remove_by,max_ind)
          remove_by <- rbind(remove_by, max_ind)
        }
        cyc <- cyc_remove & cyc_remove_by

        Z <- Z[-conf_row_ind,] # remove conflict match
        Z <- rbind(Z,max_ind) # add new match

        # correct MM caused by ZZ
        if(!is.null(nrow(ZZ))){
          MM[ZZ$seed_A,] <- M[ZZ$seed_A,]
          MM[,ZZ$seed_B] <- M[,ZZ$seed_B]
          ZZ <- c(0,0)
        }

        # update mark matrix: subtract removed seed's effect
        for (i in 1:length(conf_row_ind)) {
          A_adj <- which(A[Z$seed_A[conf_row_ind[i]],]>0)
          B_adj <- which(B[Z$seed_B[conf_row_ind[i]],]>0)
          M[A_adj, B_adj] <- M[A_adj, B_adj] - 1
          MM[A_adj, B_adj] <- MM[A_adj, B_adj] - 1
          MM[conf_ind[i,1], conf_ind[i,2]] <- M[conf_ind[i,1], conf_ind[i,2]]
        }

        # update mark matrix M & MM: add new match's effect
        A_adj <- which(A[unlist(max_ind[1]),]>0)
        B_adj <- which(B[unlist(max_ind[2]),]>0)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + 1
        MM[A_adj, B_adj] <- MM[A_adj, B_adj] + 1
        MM[max_ind[1], max_ind[2]] <- minusinf

      } else{ # choose another qualified match
        ZZ <- rbind(ZZ, conf_ind)
        MM[conf_ind[,1],] <- minusinf
        MM[,conf_ind[,2]] <- minusinf
      }
    }

  }# end while: percolate

  if(nrow(Z) == n-1){
    all <- 1:n
    seed_A <- all[!(all %in% Z$seed_A)]
    seed_B <- all[!(all %in% Z$seed_B)]
    Z <- rbind(Z,cbind(seed_A,seed_B))
  }

  # matching result
  order <- order(Z$seed_A)
  corr <- Z[order,]
  names(corr) <- c("corr_A","corr_B")
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns, order = order)
  z
}
conflict_check <- function(Matches, ind, logical = TRUE){

  if(logical == TRUE){
    rconflict <- ind[,1] %in% Matches$seed_A
    cconflict <- ind[,2] %in% Matches$seed_B
    conflict <- rconflict | cconflict
  } else{
    rconflict <- ind[1] == Matches$seed_A
    rind <- which(rconflict==TRUE)
    cconflict <- ind[2] == Matches$seed_B
    cind <- which(cconflict==TRUE)
    conflict <- c(rind,cind)
  }
  conflict
}
check_cycle <- function(rem, new){
  row <- which(rem[,1]==unlist(new[1]))
  col <- which(rem[,2]==unlist(new[2]))
  occ <- intersect(row,col)
  occ <- c(occ,length(rem[,1])+1)
  if(length(occ)==3){
    if(occ[3]-occ[2]==occ[2]-occ[1]){
      cycle1 <- rem[occ[1]:(occ[2]-1),]
      cycle2 <- rem[occ[2]:(occ[3]-1),]
      if(sum(cycle1[,1]==cycle2[,1])+sum(cycle1[,2]==cycle2[,2])==occ[3]-occ[1]){
        result <- TRUE
      } else{
        result <- FALSE
      }
    } else{
      result <- FALSE
    }
  } else{
    result <- FALSE
  }

  result
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_IsoRank} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B} and the number of seeds. If choose the greedy method to
#'   extract mapping, the order of nodes getting matched will also be returned.
#'
#' @references R. Singh, J. Xu, B. Berger (2008), \emph{Global alignment of
#' multiple protein interaction networks with application to functional
#' orthology detection}. Proc Natl Acad Sci. USA, pages 12763-12768.
#'
#' @examples
#' # match G_1 & G_2 using IsoRank algorithm
#' startm <- matrix(0, 10, 10)
#' diag(startm)[1:4] <- 1
#' GM_IsoRank <- graph_match_IsoRank(g1, g2, startm, alpha = .3, method = "greedy")
#'
#' @export
#'
graph_match_IsoRank <- function(A, B, start, alpha, max_iter=1000, method = "greedy"){
  A <- A[]
  B <- B[]

  # computing transition matrix A
  A <- A %*% Matrix::Diagonal(nrow(A), 1/Matrix::colSums(A))
  B <- B %*% Matrix::Diagonal(nrow(B), 1/Matrix::colSums(B))
  mat_A <- Matrix::kronecker(A, B)
  ns <- sum(diag(start)==1)
  start <- c(t(start)) # sparsify if poss
  E <- start/sum(abs(start))

  # computing R by power method
  R_new <- E
  tol <- 1e-5
  iter <- 1
  diff <- 1
  while(diff > tol & iter <= max_iter){

    R <- R_new
    if(alpha>0){
      AR <- mat_A %*% R
      AR <- alpha * AR + (1-alpha) * E
    } else{
      AR <- mat_A %*% R
    }
    R_new <- AR / sum(abs(AR))
    diff <- sum(abs(R-R_new))
    iter <- iter + 1
  }
  R <- matrix(R, byrow = TRUE, nrow = nrow(A))

  # find GNA
  if(method == "greedy"){
    corr <- NULL
    while (max(R)>0) {
      max_ind <- which(R == max(R), arr.ind = TRUE)
      max_ind <- max_ind[sample(nrow(max_ind), 1), ]
      corr <- rbind(corr, max_ind)
      R[max_ind[1],] <- -1
      R[,max_ind[2]] <- -1
    }
    corr <- data.frame(corr_A = corr[,1], corr_B = corr[,2])
    
    cl <- match.call()
    z <- list(call = cl, corr = corr, ns = ns, order = order(corr$corr_A))
    z
  } else if(method == "LAP"){
    # Hungarian alg.
    corr <- as.vector(clue::solve_LSAP(R, maximum = TRUE))
    corr <- data.frame(corr_A = 1:nrow(A), corr_B = corr)
    
    cl <- match.call()
    z <- list(call = cl, corr = corr, ns = ns)
    z
  }
}
#'
#' @rdname graph_match_methods
#' @return \code{graph_match_IsoRank} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B} and the number of seeds. 
#'
#' @references S. Umeyama (1988), \emph{An eigendecomposition approach to weighted
#'   graph matching problems}. IEEE TPAMI. USA, pages 695-703.
#'
#' @examples
#' # match G_1 & G_2 using Umeyama algorithm
#' G <- sample_correlated_gnp_pair(10, .9, .5)
#' G1 <- G$graph1
#' G2 <- G$graph2
#' GM_U <- graph_match_Umeyama(G1, G2, startm, alpha = .3)
#'
#' @export
#'
graph_match_Umeyama <- function(A, B, start, alpha = 0){
  A <- A[]
  B <- B[]
  ns <- sum(diag(start)==1)

  if(!isSymmetric(as.matrix(A)) | !isSymmetric(as.matrix(B))){
    # construct Hermitian matrices by adjacency matrices
    A <- as.matrix((A + Matrix::t(A))/2) + as.matrix((A - Matrix::t(A))/2)*1i
    B <- as.matrix((B + Matrix::t(B))/2) + as.matrix((B - Matrix::t(B))/2)*1i
  }

  U_A <- eigen(A)$vectors
  U_B <- eigen(B)$vectors
  AB <- Matrix::tcrossprod(abs(U_B), abs(U_A))
  Grad <- (1-alpha) * AB + alpha * Matrix::t(start)
  ind <- as.vector(clue::solve_LSAP(as.matrix(Grad), maximum = TRUE))

  corr <- data.frame(corr_A = 1:nrow(A), corr_B = ind)
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns)
  z
}
