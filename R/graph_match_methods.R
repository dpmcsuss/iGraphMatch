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
#' @param similaity A matrix. An \code{n-by-n} matrix containing vertex similaities.
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
graph_match_FW <- function(A, B, seeds = NULL,
  start = "convex", max_iter = 20,
  similarity = NULL, usejv = TRUE){

  # this will make the graphs be matrices if they are igraph objects
  A <- A[]
  B <- B[]

  # Add support for graphs with different orders
  totv1 <- ncol(A)
  totv2 <- ncol(B)

  # Check for square
  if(nrow(A) != totv1){
    stop("A is not square. iGraphMatch only supports ",
      "square matrices for matching.")
  }
  if(nrow(B) != totv2){
    stop("B is not square. iGraphMatch only supports ",
      "square matrices for matching.")
  }

  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
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

  nn <- nv - ns
  nonseeds <- !seeds

  # Asn <- A[seeds,nonseeds]
  Ann <- A[nonseeds,nonseeds]
  # Ans <- A[nonseeds,seeds]

  # Bsn <- B[seeds,nonseeds]
  Bnn <- B[nonseeds,nonseeds]
  # Bns <- B[nonseeds,seeds]

  P <- init_start(start = start, nns = nn,
                  A = A, B = B, seeds = seeds)

  iter <- 0
  toggle <- TRUE


  # make a random permutation to permute B
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]


  # seed to non-seed info
  s_to_ns <- get_s_to_ns(A, B, seeds, rp)
  # Ans %*% Matrix::t(Bns) + Matrix::t(Asn) %*% Bsn

  Bnn <- Bnn[rp, rp]

  P <- P[, rp]

  lap_method <- set_lap_method(usejv, totv1, totv2)

  if (is.null(similarity)){
    similarity <- Matrix::Matrix(0, nn, nn)
  } else {
    similarity <- similarity %*% Matrix::t(rpmat)
  }

  while(toggle && iter < max_iter){
    iter <- iter + 1

    # non-seed to non-seed info
    tAnn_P_Bnn <- Matrix::t(Ann) %*% P %*% Bnn

    Grad <- s_to_ns + Ann %*% P %*% Matrix::t(Bnn) + tAnn_P_Bnn + similarity

    ind <- do_lap(Grad, lap_method)

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

  corr_ns <- do_lap(P, lap_method)
  corr_ns <- rp[corr_ns]

  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr, ]
  D <- P
  D[nonseeds, nonseeds] <- D_ns %*% rpmat

  # fix match results if there are incorrect seeds
  if(sum(aseeds_err)!=0){
    corr <- fix_hard_corr(seed_A_err, seed_B_err,corr)
    P <- Matrix::Diagonal(nv)[corr,]
    D <- fix_hard_D(seed_A_err, seed_B_err,D)
  }

  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns,
            P = P, D = D)
  z
}

# correct the order of swapping graph2 according to new seeds
swap_order <- function(aseeds_matrix){

  # aseeds_matrix: first row:added seeds index in g1, second row added seeds match
  naseeds_err <- dim(aseeds_matrix)[2]
  ninter <- 0
  ninter_new <- naseeds_err
  aseeds_match_order <- matrix( , 2, )
  aseeds_matrix_T <- aseeds_matrix

  while(ninter_new != ninter & ninter_new > 1){
    aseeds_matrix <- aseeds_matrix_T
    naseeds_err <- ninter_new
    inter_match <- rep("FALSE",times = naseeds_err)
    ninter <- ninter_new
    ninter_new <- 0
    circle_index <- 0
    k <- 1

    for(i in 1:naseeds_err){
      # eliminate circle of two vertices
      if (aseeds_matrix[2, i] %in% aseeds_matrix[1, ]){
        index <- which(aseeds_matrix[1, ] == aseeds_matrix[2, i])
        if(aseeds_matrix[1, i] == aseeds_matrix[2, index]){
          aseeds_matrix[1, i] <- 0
          circle_index[k] <- i
          k <- k + 1
        } else{
          inter_match[i] <- "TRUE"
          ninter_new <- ninter_new+1
        }
      }
    }

    if (circle_index[1] != 0){
      aseeds_matrix <- aseeds_matrix[, -circle_index]
      inter_match <- inter_match[-circle_index]
    }

    if (length(which(inter_match == "TRUE")) >= 1){
      aseeds_matrix <- as.matrix(aseeds_matrix)
      aseeds_matrix_T <- aseeds_matrix[,which(inter_match=="TRUE")]
    }
    if (length(which(inter_match == "FALSE")) >= 1){
      aseeds_matrix <- as.matrix(aseeds_matrix)
      aseeds_matrix_F <- aseeds_matrix[,which(inter_match=="FALSE")]
      aseeds_match_order <- cbind(aseeds_matrix_F, aseeds_match_order)
    }

  }

  # end with circle: only consider one circle
  # (circle with more than three vertices) case
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
  aseeds_matrix <- matrix(c(seed_g1_err, seed_g2_err),
    nrow = 2, byrow = TRUE)
  
  if(length(seed_g1_err) > 1)
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

# combine these next two functions into one?

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
graph_match_convex <- function(A, B, seeds = NULL, start = "bari", 
                               max_iter = 100, similarity = NULL,
                               tol = 1e-5, usejv = TRUE){

  A <- A[]
  B <- B[]

  totv1 <- nrow(A)
  totv2 <- norw(B)

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

  similarity <- similarity[nonseeds,nonseeds]
  tol0 <- 1
  P <- init_start(start = start, nns = nn)
  iter<-0
  toggle <- TRUE

  AtA <- t(Asn)%*%Asn + t(Ann)%*%Ann
  BBt <- Bns%*%t(Bns)+Bnn%*%t(Bnn)

  ABns_sn <- Ans%*%t(Bns) + t(Asn)%*%Bsn
  f <- sum((Ann %*% P - P%*% Bnn)^2)



  lap_method <- set_lap_method(usejv, totv1, totv2)

  while(toggle && iter < max_iter){
    f_old <- f
    iter<-iter+1

    if(is.null(similarity)){
      similarity <- Matrix::Matrix(0, nn, nn)
    }
    Grad <- AtA%*%P + P%*%BBt - ABns_sn - t(Ann)%*%P%*%Bnn - Ann%*%P%*%t(Bnn) + similarity
   

    corr <- do_lap(Grad, lap_method)
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

    P_new <- aopt*P+(1-aopt)*Pdir
    f <- sum((Ann %*% P_new - P_new %*% Bnn)^2)

    f_diff <- abs(f-f_old)
    P_diff <- sum(abs(P-P_new))
    P <- P_new

    toggle <- f_diff > tol && f > tol && P_diff > tol
  }

  D_ns <- P
  corr_ns <- do_lap(P, lap_method)
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


# #'
# #' @return \code{graph_match_convex_directed} returns graph matching results based
# #' on convex relaxation method for directed graphs.
# #'
# #' @examples
# #' graph_match_convex_directed(g1, g2, seeds)
# #'
# #'
# #'
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
#' @return \code{graph_match_PATH} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B}, the number of seeds if any, the permutation matrix and the
#'   doubly stochastic matrix before projection onto the permutation set. 
#'
#' @references M. Zaslavskiy, F. Bach and J. Vert (2009), \emph{A Path following
#' algorithm for the graph matching problem}. IEEE Trans Pattern Anal Mach Intell,
#' pages 2227-2242.
#'
#' @examples
#' # match G_1 & G_2 using PATH algorithm
#' graph_match_PATH(g1, g2)
#'
#' @export
#'
#'
graph_match_PATH <- function(A, B, similarity = NULL, seeds = NULL, alpha = .5, epsilon = 1){
  totv1 <- vcount(A)
  totv2 <- vcount(B)
  
  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
  }
  
  D_A <- Matrix::Diagonal(length(degree(A)), degree(A))
  D_B <- Matrix::Diagonal(length(degree(B)), degree(B))
  A <- A[]
  B <- B[]
  L_A <- D_A - A
  L_B <- D_B - B
  n <- nrow(A)
  
  # alpha=0, convex relaxation
  convex_m <- graph_match_convex(A, B, similarity = similarity, seeds = seeds, tol = 1e-10)
  P <- convex_m$P
  lambda <- 0
  dlambda <- dlambda_min <-  1e-5
  #toggle <- TRUE
  delta_cal <- function(x, y){
    (y - x) ^ 2
  }
  delta <- outer(diag(D_A), diag(D_B), delta_cal)
  iter <- 0
  
  lap_method <- set_lap_method(FALSE, totv1, totv2)

  while (lambda < 1) {
    iter <- iter + 1
    # dlambda-adaptation
    F_cv <- (Matrix::norm(A %*% P - P %*% B, type = "F")) ^ 2
    L <- Matrix::kronecker(Matrix::t(L_B), Matrix::t(L_A))
    F_cc <- - sum(t(delta) %*% P) - 2 * t(Matrix::c.sparseVector(P)) %*% 
      L %*% Matrix::c.sparseVector(P)
    if(!is.null(similarity)){
      F_sim <- sum(similarity * P)
      F <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    } else{
      F <- (1 - lambda) * F_cv + lambda * F_cc
    }
    
    lambda <- lambda + dlambda
    if(!is.null(similarity)){
      F_sim <- sum(similarity * P)
      F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    } else{
      F_new <- (1 - lambda) * F_cv + lambda * F_cc
    }
    
    while (sum(abs(F - F_new)) < epsilon && lambda < 1) {
      dlambda <- 2 * dlambda
      lambda <- lambda + dlambda
      if(lambda > 1){
        lambda <- 1
        break
      }
      if(!is.null(similarity)){
        F_sim <- sum(similarity * P)
        F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
      } else{
        F_new <- (1 - lambda) * F_cv + lambda * F_cc
      }    
    }
    while (sum(abs(F - F_new)) > epsilon && dlambda != dlambda_min) {
      if(lambda > 1){
        lambda <- 1
        break
      } else{
        dlambda <- dlambda / 2
        lambda <- lambda - dlambda
      }
      if(!is.null(similarity)){
        F_sim <- sum(similarity * P)
        F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
      } else{
        F_new <- (1 - lambda) * F_cv + lambda * F_cc
      }
    }
    
    # Frank-Wolfe 
    AtA <- t(A) %*% A
    BBt <- B %*% t(B)
    tA_P_B <- Matrix::t(A) %*% P %*% B
    Grad_cv <- 2 * (AtA %*% P + P %*% BBt - tA_P_B - A %*% P %*% t(B))
    Grad_cc <- - t(delta) - 2 * Matrix::t(L_A) %*% P %*% L_B
    Grad <- (1 - lambda) * Grad_cv + lambda * Grad_cc
    if(!is.null(similarity)){
      Grad <- alpha * Grad + (1 - alpha) * similarity
    }
    ind <- do_lap(Grad, lap_method)
    ind2 <- cbind(1:n, ind)
    Pdir <- Matrix::Diagonal(n)[ind, ]
    
    delta_P <- P - Pdir
    C <- A %*% delta_P - delta_P %*% B
    D <- A %*% Pdir - Pdir %*% B
    aq <- sum(C^2)
    bq <- sum(C*D)
    vec_delta_P <- Matrix::c.sparseVector(delta_P)
    vec_Pdir <- Matrix::c.sparseVector(Pdir)
    c <- sum(t(delta) * delta_P)
    e <- Matrix::t(vec_delta_P) %*% L %*% vec_Pdir
    u <- Matrix::t(vec_Pdir) %*% L %*% vec_delta_P
    v <- Matrix::t(vec_delta_P) %*% L %*% vec_delta_P
    a <- 2 * (lambda - 1) * bq + lambda * (c - e + u)
    b <- 2 * (1 - lambda) * aq - 4 * lambda * v
    if(a[1,1] == 0 && b[1,1] == 0){
      alpha <- 0
    } else{
      alpha <- (2 * (lambda - 1) * bq + lambda * (c - e + u)) / 
        (2 * (1 - lambda) * aq - 4 * lambda * v)
      alpha <- alpha[1,1]
    }
    if(alpha > 1){
      alpha <- 1
    } else if(alpha < 0){
      alpha <- 0
    }
    P <- alpha * P + (1 - alpha) * Pdir
  }
  
  D <- P
  corr <- do_lap(P, lap_method)
  P <- Matrix::Diagonal(n)[corr,]
  
  if(!is.null(seeds)){
    ns <- nrow(check_seeds(seeds))
  } else{
    ns <- 0
  }
  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns, 
            P = P, D = D, iter = iter, lambda = lambda)
  z
}
#'
#' @rdname graph_match_methods
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
graph_match_percolation <- function (A, B, seeds, r = 2) 
{
  if(is.igraph(A)){
    weighted <- is.weighted(A)
  } else{
    if(min(A) < 0){
      weighted <- TRUE
    } else{
      weighted <- max(A) > 1
    }
  }
  A <- A[]
  B <- B[] 
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  n <- max(totv1, totv2)
  P <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
  seeds <- check_seeds(seeds)
  P[as.matrix(seeds)] <- 1
  Z <- seeds
  
  if(weighted){
    M <- Matrix::Matrix(0, totv1, totv2)
    for(i in 1:nrow(seeds)){
      A_adj <- which(A[seeds$seed_A[i],]>0)
      B_adj <- which(B[seeds$seed_B[i],]>0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[seeds$seed_A[i],A_adj], B[seeds$seed_B[i],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
    }
  } else{
    M <- (Matrix::t(A) %*% P %*% B + A %*% P %*% Matrix::t(B)) / 2
  }
  M[seeds$seed_A,] <- -n
  M[,seeds$seed_B] <- -n
  
  while (max(M) >= r) {
    max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind), 1), ]
    if(weighted){
      A_adj <- which(A[max_ind[1],]>0)
      B_adj <- which(B[max_ind[2],]>0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
    } else{
      Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
      Pi[max_ind[1], max_ind[2]] <- 1 
      delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
      M <- M + delta
    }
    M[max_ind[1], ] <- -n 
    M[, max_ind[2]] <- -n
    Z <- rbind(Z, max_ind)
  }
  
  order <- order(Z$seed_A)
  corr <- Z[order,]
  names(corr) <- c("corr_A","corr_B")
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = nrow(seeds), order = order)
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
  if(is.igraph(A)){
    weighted <- is.weighted(A)
  } else{
    if(min(A) < 0){
      weighted <- TRUE
    } else{
      weighted <- max(A) > 1
    }
  }
  A <- A[]
  B <- B[]
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  n <- max(totv1, totv2)
  P <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
  seeds <- check_seeds(seeds)
  seeds_ori <- seeds
  P[as.matrix(seeds)] <- 1
  M <- Matrix::Matrix(0, totv1, totv2)
  M[seeds_ori$seed_A,] <- -n
  M[,seeds_ori$seed_B] <- -n
  Z <- seeds

  # deferred percolation graph matching
  while(nrow(seeds) != 0){
    # mark neighbors
    if(weighted){
      for(i in 1:nrow(seeds)){
        A_adj <- which(A[seeds$seed_A[i],]>0)
        B_adj <- which(B[seeds$seed_B[i],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[seeds$seed_A[i],A_adj], B[seeds$seed_B[i],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
      }
    } else{
      Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
      Pi[as.matrix(seeds)] <- 1
      M <- M + (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
    }

    # choose pairs with marks ge r
    while(max(M) >= r){
      max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
      if(nrow(max_ind) != 1){
        degree_diff <- abs(rowSums(A)[max_ind[,1]]-rowSums(B)[max_ind[,2]])
        max_ind <- max_ind[which(degree_diff == min(degree_diff)),]
        if(is.vector(max_ind) == FALSE){
          max_ind <- max_ind[sample(nrow(max_ind),1),]
        }
      }

      # update mark matrix
      if(weighted){
        A_adj <- which(A[max_ind[1],]>0)
        B_adj <- which(B[max_ind[2],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
      } else{
        Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
        Pi[max_ind[1], max_ind[2]] <- 1 
        delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
        M <- M + delta
      }
      M[max_ind[1],] <- -n
      M[,max_ind[2]] <- -n
      max_ind <- data.frame(seed_A = max_ind[1], seed_B = max_ind[2])
      Z <- rbind(Z, max_ind)
    }

    seeds_old <- seeds
    seeds <- which(M > 0 & M < r, arr.ind = TRUE)
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
  z <- list(call = cl, corr = corr, ns = nrow(seeds), order = order)
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
graph_match_soft_percolation <- function(A, B, seeds, r = 2, max_iter = 100){

  # this will make the graphs be matrices if they are igraph objects
  if(is.igraph(A)){
    weighted <- is.weighted(A)
  } else{
    if(min(A) < 0){
      weighted <- TRUE
    } else{
      weighted <- max(A) > 1
    }
  }
  A <- A[]
  B <- B[]
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  n <- max(totv1, totv2)
  P <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
  seeds <- check_seeds(seeds)
  ns <- nrow(seeds)
  seeds_ori <- seeds
  P[as.matrix(seeds)] <- 1
  
  # initialization of score matrix M & MM
  if(weighted){
    M <- Matrix::Matrix(0, totv1, totv2)
    for(i in 1:nrow(seeds)){
      A_adj <- which(A[seeds$seed_A[i],]>0)
      B_adj <- which(B[seeds$seed_B[i],]>0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[seeds$seed_A[i],A_adj], B[seeds$seed_B[i],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
    }
  } else{
    M <- (Matrix::t(A) %*% P %*% B + A %*% P %*% Matrix::t(B)) / 2
  }
  M[seeds$seed_A,] <- -n
  M[,seeds$seed_B] <- -n
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
  while(max(MM)>=r & num<=max_iter){
    # locate best match
    max_ind <- which(MM==max(MM), arr.ind = TRUE)
    conflict_log <- conflict_check(Z, max_ind, logical = TRUE)
    sum_conf <- sum(conflict_log)
    if(sum_conf>0 && sum_conf<length(conflict_log)){
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
      if(weighted){
        A_adj <- which(A[max_ind[1],]>0)
        B_adj <- which(B[max_ind[2],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
          MM[A_adj, B_adj] <- MM[A_adj, B_adj] + mark
        }
      } else{
        Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
        Pi[max_ind[1], max_ind[2]] <- 1 
        delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
        M <- M + delta
        MM <- MM + delta
      }
      MM[max_ind[1], max_ind[2]] <- -n

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
        if(weighted){
          for (i in 1:length(conf_row_ind)) {
            A_adj <- which(A[Z$seed_A[conf_row_ind[i]],]>0)
            B_adj <- which(B[Z$seed_B[conf_row_ind[i]],]>0)
            if(length(A_adj) != 0 && length(B_adj) != 0){
              mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
              M[A_adj, B_adj] <- M[A_adj, B_adj] - mark
              MM[A_adj, B_adj] <- MM[A_adj, B_adj] - mark
              MM[conf_ind[i,1], conf_ind[i,2]] <- M[conf_ind[i,1], conf_ind[i,2]]
            }
          }
        } else{
          Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
          Pi[as.matrix(conf_ind)] <- 1 
          delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
          M <- M - delta
          MM <- MM - delta
          MM[as.matrix(conf_ind)] <- M[as.matrix(conf_ind)]
        }

        # update mark matrix M & MM: add new match's effect
        if(weighted){
          A_adj <- which(A[max_ind[1],]>0)
          B_adj <- which(B[max_ind[2],]>0)
          if(length(A_adj) != 0 && length(B_adj) != 0){
            mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
            M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
            MM[A_adj, B_adj] <- MM[A_adj, B_adj] + mark
          }
        } else{
          Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
          Pi[max_ind[1], max_ind[2]] <- 1 
          delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
          M <- M + delta
          MM <- MM + delta
        }
        MM[max_ind[1], max_ind[2]] <- -n

      } else{ # choose another qualified match
        ZZ <- rbind(ZZ, conf_ind)
        MM[conf_ind[,1],] <- -n
        MM[,conf_ind[,2]] <- -n
      }
    }

  }# end while: percolate

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
graph_match_IsoRank <- function(A, B, similarity, alpha = .5, max_iter = 1000, method = "greedy"){
  A <- A[]
  B <- B[]
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  
  # padding if two graphs different sizes
  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
  }
  
  # computing transition matrix A
  colS_A <- Matrix::colSums(A)
  colS_B <- Matrix::colSums(B)
  A <- A %*% Matrix::Diagonal(nrow(A), ifelse(colS_A == 0, 0, 1/colS_A))
  B <- B %*% Matrix::Diagonal(nrow(B), ifelse(colS_B == 0, 0, 1/colS_B))
  mat_A <- Matrix::kronecker(A, B)
  #start <- Matrix::c.sparseVector(similarity)
  start <- Matrix::c.sparseVector(Matrix::t(similarity)) 
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
  R <- ramify::resize(R, nrow = totv1, ncol = totv1, byrow = FALSE)
  
  # find GNA
  if(method == "greedy"){
    corr <- NULL
    while (max(R)>0) {
      max_ind <- Matrix::which(R == max(R), arr.ind = TRUE)
      max_ind <- max_ind[sample(nrow(max_ind), 1), ]
      corr <- rbind(corr, max_ind)
      R[max_ind[1],] <- -1
      R[,max_ind[2]] <- -1
    }
    corr <- data.frame(corr_A = corr[,1], corr_B = corr[,2])
    
    cl <- match.call()
    z <- list(call = cl, corr = corr, ns = 0, order = order(corr$corr_A))
    z
  } else if(method == "LAP"){
    # Hungarian alg.
    lap_method <- set_lap_method(FALSE, totv1, totv2)
    corr <- do_lap(R, lap_method)
    corr <- data.frame(corr_A = 1:nrow(A), corr_B = corr)
    
    cl <- match.call()
    z <- list(call = cl, corr = corr, ns = 0)
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
graph_match_Umeyama <- function(A, B, similarity = NULL, alpha = .5){
  A <- A[]
  B <- B[]
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  
  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
  }

  if(!isSymmetric(as.matrix(A)) | !isSymmetric(as.matrix(B))){
    # construct Hermitian matrices by adjacency matrices
    A <- as.matrix((A + Matrix::t(A))/2) + as.matrix((A - Matrix::t(A))/2)*1i
    B <- as.matrix((B + Matrix::t(B))/2) + as.matrix((B - Matrix::t(B))/2)*1i
  }

  U_A <- eigen(A)$vectors
  U_B <- eigen(B)$vectors
  AB <- Matrix::tcrossprod(abs(U_B), abs(U_A))
  if(!is.null(similarity)){
    Grad <- alpha * AB + (1-alpha) * Matrix::t(similarity)
  } else{
    Grad <- AB
  }
  Grad <- Grad - min(Grad)
  lap_method <- set_lap_method(FALSE, totv1, totv2)
  ind <- do_lap(Grad, lap_method)

  corr <- data.frame(corr_A = 1:nrow(A), corr_B = ind)
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = 0)
  z
}
