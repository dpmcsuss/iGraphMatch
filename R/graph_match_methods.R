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
#' @param max_iter A number. Maximum number of replacing matches.
#'
#' @rdname graph_match_methods
#'
#' @return \code{graph_match_FW} returns a list of graph matching results,
#'   including match correspondence vector of \eqn{G_2} with respect to
#'   \eqn{G_1} named \code{corr}, doubly stochastic matrix named \code{D},
#'   permutation matrix named \code{P} based on Frank-Wolfe methodology and
#'   iteration time of the algorithm named \code{iter}.
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
graph_match_FW <- function(A, B, seeds = NULL,
  start = "convex", max_iter = 20,
  similarity = NULL, return_big = TRUE, usejv = FALSE){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]

  # Add support for graphs with different orders
  totv1<-ncol(A)
  totv2<-ncol(B)
  if(totv1>totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1<totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
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

  # permute B matrices by a random perm to avoid bias
  # call the below Q
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]
  # B <- QBQ^T
  Bnn <- rpmat %*% Bnn %*% Matrix::t(rpmat)
  Bns <- rpmat %*% Bns
  Bsn <- Bsn %*% Matrix::t(rpmat)
  # change P <- PQ^T
  P <- P %*% Matrix::t(rpmat)
  if(is.null(similarity)){
    similarity <- Matrix::Matrix(0,nn,nn)
  } else {
    similarity <- similarity %*% Matrix::t(rpmat)
  }

  rm(A,B)
  gc()

  iter <- 0
  toggle <- TRUE

  # seed to non-seed info
  if(ns > 1){
    s_to_ns <- Ans %*% Matrix::t(Bns) + Matrix::t(Asn) %*% Bsn
  } else if( ns == 1){
    s_to_ns <- outer(Ans,Bns) + outer(Asn, Bsn)
  } else {
    s_to_ns <- Matrix(0, nv, nv)
  }
  
  usejvmod <- FALSE
  if("rlapjv" %in% rownames(installed.packages()) ){
    library(rlapjv)
    # usejv <- TRUE
    if( totv1 / totv2 < 0.5 ){
      usejvmod <- TRUE
      usejv <- FALSE
    }
  } else {
    usejv <- FALSE
  }

  while(toggle && iter < max_iter){
    iter <- iter + 1

    # non-seed to non-seed info
    tAnn_P_Bnn <- Matrix::t(Ann) %*% P %*% Bnn

    Grad <- s_to_ns + Ann %*% P %*% Matrix::t(Bnn) + tAnn_P_Bnn + similarity
    if ( usejv ){
      Grad <- as.matrix(Grad)
      ind <- rlapjv::lapjv(Grad+1000,
        # round(Grad * nn ^ 2 * max(Grad)),
        maximize = TRUE)
    } else if ( usejvmod ) {
      if( class(Grad) == "splrMatrix" ){
        ind <- rlapjv::lapmod(splr.to.sparse(Grad),
          maximize = TRUE)
      } else {
        ind <- rlapjv::lapmod(Grad, maximize = TRUE)
      }
    } else {
      Grad <- as.matrix(Grad)
      Grad <- (Grad - min(Grad))
      ind <- as.vector(clue::solve_LSAP(Grad,
        maximum = TRUE))
    }

    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)
    Pdir <- Pdir[ind, ]
    ns_Pdir_ns <- Matrix::t(Ann)[, order(ind)] %*% Bnn

    cc <- innerproduct(tAnn_P_Bnn, P)
    d <- innerproduct(ns_Pdir_ns, P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(ns_Pdir_ns[ind2])
    u <- innerproduct(P, s_to_ns + similarity)
    v <- sum((s_to_ns + similarity)[ind2])


    if (cc - d + e == 0 && d - 2 * e + u - v == 0) {
      alpha <- 0
    } else {
      alpha <- -(d - 2 * e + u - v)/(2 * (cc - d + e))
    }

    f0 <- 0
    f1 <- cc - e + u - v
    falpha <- (cc - d + e) * alpha^2 + (d - 2 * e + u - v) *
      alpha

    Pold <- P
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

  if ( usejv ){
    corr_ns <- rlapjv::lapjv(P, maximize = TRUE)
  } else if ( usejvmod ) {
    if (class(P) == "splrMatrix") {
      # this is pretty hacky but it breaks it otherwise
      corr_ns <- rlapjv::lapmod(P@x, maximize = TRUE)
    } else {
      corr_ns <- rlapjv::lapmod(P, maximize = TRUE)
    }
  } else {
    corr_ns <- as.vector(clue::solve_LSAP(as.matrix(P), 
      maximum = TRUE))
  }


  # undo rand perm here
  corr_ns <- rp[corr_ns]
  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr, ]
  D <- P
  # and undo it right quick here too
  if ( class(D_ns) == "splrMatrix"){
    if ( nn < nv){
      warning("Only returning non-seed D.")
    }
    D <- D_ns
  } else {
    D[nonseeds, nonseeds] <- D_ns %*% rpmat
  }
  # and we should be home clear


  # fix match results if there are incorrect seeds
  if(sum(aseeds_err)!=0){
    corr <- fix_hard_corr(seed_A_err,seed_B_err,corr)
    P <- Matrix::Diagonal(nv)[corr,]
    D <- fix_hard_D(seed_A_err,seed_B_err,D)
  }
  if( return_big ){
    list(corr = corr, P = P, D = D, iter = iter)
  } else { 
    list(corr = corr, iter = iter)
  }
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


#' @rdname graph_match_methods
#' @title bla
#' @return \code{graph_match_convex} returns a list of graph matching results,
#'   including matching correspondence vector of \eqn{G_2} with respect to
#'   \eqn{G_1} named \code{corr}, doubly stochastic matrix named \code{D} and
#'   permutation matrix named \code{P} based on convex relaxation method for
#'   undirected graphs.
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
graph_match_convex <- function(A, B, seeds = NULL, start = "bari", max_iter = 100, tol = 1e-5){
 
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)
  
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
  tol0 <- 1e-5
  P <- init_start(start = start, nns = nn)
  iter<-0
  toggle <- TRUE
  
  AtA <- t(Asn)%*%Asn + t(Ann)%*%Ann
  BBt <- Bns%*%t(Bns)+Bnn%*%t(Bnn)
  ABns_sn <- Ans%*%t(Bns) + t(Asn)%*%Bsn
  f <- sum((Asn%*%P-Bsn)^2)+sum((Ans-P%*%Bns)^2)+sum((Ann%*%P-P%*%Bnn)^2)

  if ( "rlapjv" %in% rownames(installed.packages()) ){
    library(rlapjv)
    usejv <- TRUE
  } else {
    usejv <- FALSE
  }

  while(toggle && iter<max_iter){
    f_old <- f
    iter<-iter+1
    Grad<- 2*(AtA%*%P + P%*%BBt - ABns_sn - t(Ann)%*%P%*%Bnn - Ann%*%P%*%t(Bnn))

    Grad <- as.matrix(nn^2*(Grad-min(Grad)))
    if ( usejv ){
      corr <- rlapjv::lapjv(Grad, maximize = TRUE)
    } else {
      corr <- as.vector(clue::solve_LSAP(Grad,
        maximum = TRUE))
    }
    Pdir <- Matrix::Diagonal(nn)[corr,]
    difPPdir<-Pdir-P
    aq <- sum((Asn%*%difPPdir)^2)+sum((difPPdir%*%Bns)^2)+sum((Ann %*% difPPdir - difPPdir %*% Bnn)^2)
    bq <- sum(diag(t(Asn%*%P-Bsn)%*%Asn%*%difPPdir))+sum(diag(t(P%*%Bns-Ans)%*%difPPdir%*%Bns))+sum(diag((Ann %*% difPPdir - difPPdir %*% Bnn)%*%t(Ann %*%P - P%*% Bnn)))
    aopt <- -bq/aq
    aopt<-aopt*(aopt>0)*(aopt<1)+(aopt>1) #restrict the step size in (0,1)
    P_new <- (1-aopt)*P+aopt*Pdir;
    f <- sum((Asn%*%P_new-Bsn)^2)+sum((Ans-P_new%*%Bns)^2)+sum((Ann%*%P_new-P_new%*%Bnn)^2)
    f_diff <- abs(f-f_old)
    P_diff <- sum(abs(P-P_new))
    P <- P_new
    
    toggle <- f_diff > tol0 && f > tol && P_diff > tol0
  }
  
  D_ns <- P
  corr_ns <- as.vector(clue::solve_LSAP(as.matrix(P),maximum = TRUE))
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
  
  list(corr = corr, P = P, D = D)
}


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


#' @rdname graph_match_methods
#' @return \code{graph_match_percolation} returns matching correspondence of
#'   matched pairs with index of nodes in \eqn{G_1} named \code{corr_A} and
#'   index of nodes in \eqn{G_2} named \code{corr_B}.
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
graph_match_percolation <- function (A, B, seeds, r = 2, center = FALSE) 
{
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  n_A <- dim(A)[1]
  n_B <- dim(B)[1]
  n <- max(n_A, n_B)
  P <- matrix(0, nrow=n_A, ncol = n_B)
  seeds <- check_seeds(seeds)
  for (i in 1:nrow(seeds)){
    si <- seeds[i,]
    P[si$seed_A,si$seed_B] <- 1
  }
  Z <- seeds
  M <- A %*% P %*% B
  M[seeds$seed_A,] <- -n
  M[,seeds$seed_B] <- -n
  while (max(M) >= r) {
    max_ind <- which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind), 1), ]
    Pi <- matrix(0, nrow=n_A, ncol = n_B)
    Pi[max_ind[1], max_ind[2]] <- 1 
    delta <- A %*% Pi %*% B
    delta_norm <- (delta - min(delta)) / (max(delta) - min(delta))
    if(center){
      delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
    }
    M <- M + delta_norm
    M[max_ind[1], ] <- -n
    M[, max_ind[2]] <- -n
    Z <- rbind(Z, max_ind)
  }
  if (n_A == n_B & nrow(Z) == n - 1) {
    all <- 1:n
    seed_A <- all[!(all %in% Z$seed_A)]
    seed_B <- all[!(all %in% Z$seed_B)]
    Z <- rbind(Z, cbind(seed_A, seed_B))
  }
  corr <- Z[order(Z$seed_A), ]
  names(corr) <- c("corr_A", "corr_B")
  corr
}
cal_mark <- function(x,y){
  1 - abs(x - y) / max(x, y)
}

#' @rdname graph_match_methods
#' @return \code{graph_match_ExpandWhenStuck} returns matching correspondence of
#'   matched pairs with index of nodes in \eqn{G_1} named \code{corr_A} and
#'   index of nodes in \eqn{G_2} named \code{corr_B}.
#'
#' @references E. Kazemi, S. H. Hassani, and M. Grossglauser (2015),
#' \emph{Growing a graph matching from a handful of seeds}. Proc. of the VLDB
#' Endowment, 8(10):1010-1021.
#'
#' @examples
#' # match G_1 & G_2 using Expand When Stuck graph matching method
#' seeds <- 1:5
#' graph_match_ExpandWhenStuck(g1, g2, seeds, r = 2)
#'
#' @export
graph_match_ExpandWhenStuck <- function(A, B, seeds, r = 2){
  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  n <- nrow(A)
  m <- nrow(B)
  seeds <- check_seeds(seeds) #unused seeds
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
  corr <- Z[order(Z$seed_A),]
  names(corr) <- c("corr_A","corr_B")
  corr
}

#' @rdname graph_match_methods
#' @return \code{graph_match_soft_percolation} returns matching correspondence
#'   of matched pairs with index of nodes in \eqn{G_1} named \code{corr_A} and
#'   index of nodes in \eqn{G_2} named \code{corr_B}.
#'
#' @examples
#' # match G_1 & G_2 using soft percolation graph matching method
#' seeds <- 1:5
#' graph_match_soft_percolation(g1, g2, seeds, r = 2, max_iter = 2)
#'
#' @export
graph_match_soft_percolation <- function(A, B, seeds, r = 2, max_iter = 50, center = FALSE){
  
  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  # initialization of score matrix M & MM
  n_A <- dim(A)[1]
  n_B <- dim(B)[1]
  P <- matrix(0, nrow=n_A, ncol = n_B)
  seeds <- check_seeds(seeds)
  ns <- nrow(seeds)
  for (i in 1:nrow(seeds)){
    si <- seeds[i,]
    P[si$seed_A,si$seed_B] <- 1
  }
  M <- A %*% P %*% B
  minusinf <- -n
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
  while(max(MM)>=r & num<=max_iter & cyc==FALSE){
    # locate one best match
    max_ind <- which(MM==max(MM), arr.ind = TRUE)
    rnum <- sample(nrow(max_ind),1)
    max_ind <- max_ind[rnum,]
    
    conflict_log <- conflict_check(Z, matrix(max_ind,1), logical = TRUE)
    # non-conflict new match 
    if(conflict_log==FALSE){
      
      # correct MM caused by ZZ
      if(!is.null(nrow(ZZ))){
        MM[ZZ$seed_A,] <- M[ZZ$seed_A,]
        MM[,ZZ$seed_B] <- M[,ZZ$seed_B]
        for (i in 1:nrow(ZZ)) {
          ZZi <- ZZ[i,]
          MM[ZZi$seed_A, ZZi$seed_B] <- minusinf
        }
        ZZ <- c(0,0)
      }
      
      # update mark matrix M & MM
      Pi <- matrix(0, nrow=n_A, ncol = n_B)
      Pi[unlist(max_ind[1]), unlist(max_ind[2])] <- 1 
      delta <- A %*% Pi %*% B
      delta_norm <- (delta - min(delta)) / (max(delta) - min(delta))
      if(center){
        delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
      }
      M <- M + delta_norm
      
      MM <- MM + delta_norm
      MM[max_ind[1], max_ind[2]] <- minusinf
      
      Z <- rbind(Z,max_ind)
      
    } else{ # conflict new match
      
      conf_row_ind <- conflict_check(Z, matrix(max_ind,1), logical = FALSE)
      if(length(conf_row_ind)==2){ # conflict with two existing matches
        conf_ind <- Z[conf_row_ind,]
        score1 <- M[conf_ind[1,1], conf_ind[1,2]]
        score2 <- M[conf_ind[2,1], conf_ind[2,2]]
        score <- max(score1,score2)
      } else{
        conf_ind <- Z[conf_row_ind[1],]
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
        
        # correct MM caused by ZZ
        if(!is.null(nrow(ZZ))){
          MM[ZZ$seed_A,] <- M[ZZ$seed_A,]
          MM[,ZZ$seed_B] <- M[,ZZ$seed_B]
          for (i in 1:nrow(ZZ)) {
            ZZi <- ZZ[i,]
            MM[ZZi$seed_A, ZZi$seed_B] <- minusinf
          }
          ZZ <- c(0,0)
        }
        
        # update mark matrix: subtract removed seed's effect
        Pi <- matrix(0, nrow=n_A, ncol = n_B)
        for (i in 1:length(conf_row_ind)) {
          Pi[Z$seed_A[conf_row_ind[i]], Z$seed_B[conf_row_ind[i]]] <- 1 
          delta <- A %*% Pi %*% B
          delta_norm <- (delta - min(delta)) / (max(delta) - min(delta))
          if(center){
            delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
          }
          M <- M - delta_norm
          MM <- MM - delta_norm
          MM[conf_ind[i,1], conf_ind[i,2]] <- M[conf_ind[i,1], conf_ind[i,2]]
        }
        
        # update mark matrix M & MM: add new match's effect
        Pi <- matrix(0, nrow=n_A, ncol = n_B)
        Pi[unlist(max_ind[1]), unlist(max_ind[2])] <- 1 
        delta <- A %*% Pi %*% B
        delta_norm <- (delta - min(delta)) / (max(delta) - min(delta))
        if(center){
          delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
        }
        M <- M + delta_norm
        MM <- MM + delta_norm
        MM[max_ind[1], max_ind[2]] <- minusinf
        
        Z <- Z[-conf_row_ind,] # remove conflict match
        Z <- rbind(Z,max_ind) # add new match
        
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
  corr <- Z[order(Z$seed_A),]
  names(corr) <- c("corr_A","corr_B")
  corr
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
