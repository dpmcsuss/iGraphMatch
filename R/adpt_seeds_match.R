#' @title Graph matching with adpative seeds or soft seeds
#'
#' @description Match two given graphs with a set of selected best matched nodes as additional seeds.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param seeds A logical vector. \code{TRUE} indicates the corresponding
#' vertex is a seed.
#' @param non_seed_core A logical vector. \code{TRUE} indicates the corresponding
#' vertex is a core vertex and is not a seed. Note in the case without junk vertices,
#' \code{non_seed_core=!seeds}.
#' @param select_seeds, A vector. A vector, a matrix or a data frame. If there is no error in soft
#' seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
#' seeds, input in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
#' that if there are seeds in graphs, seeds should be put before non-seeds.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#' character value like "bari", "convex" or "rds" to initialize the starting matrix.
#'
#' @rdname adpt_seeds_match
#' @return Returns a list of graph matching results, including core matching error, ultimate matching
#' error in selected seeds, and value for the objective function.
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, p =  0.5, rho = 0.3)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:50 <= 10
#' match <- graph_match_FW(g1, g2, seeds, start = "bari")
#' select_seeds <- select_seeds(g1, g2, "row_cor", nseeds=5, x=!seeds, match_corr=match$corr)
#' match_adpt <- graph_match_adpt_seeds(g1, g2, seeds, select_seeds=select_seeds)
#' match_adpt
#'
#' @export
graph_match_adpt_seeds <- function(A, B, seeds, non_seed_core = !seeds, select_seeds, start="convex"){
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  select_seeds <- check_soft_seeds(select_seeds)
  seed_A <- select_seeds$seed_A
  seed_B <- select_seeds$seed_B
  aseeds_err <- ifelse(seed_A!=seed_B,TRUE,FALSE)
  seed_A_err <- seed_A[aseeds_err]
  seed_B_err <- seed_B[aseeds_err]
  aseeds <- length(seed_A)

  # swap columns corresponding to error seeds in B
  if(sum(aseeds_err)!=0){
    B_hard <- g2_hard_seeding(seed_A_err,seed_B_err,B)
  }
  else{
    B_hard <- B
  }

  # update seeds
  seeds[seed_A] <- TRUE

  # graph_match with updated seeds and calculate errors
  match_hard <- graph_match_FW(A, B_hard, seeds = seeds, start = start)

  # fix match results
  if(sum(aseeds_err)!=0){
    match_hard$corr <- fix_hard_corr(seed_A_err,seed_B_err,match_hard$corr)
    corr <- match_hard$corr
    nv <- nrow(A)
    match_hard$P <- Matrix::Diagonal(nv)[corr,]
    match_hard$D <- fix_hard_D(seed_A_err,seed_B_err,match_hard$D)
  }

  core_err_hard <- mean(corr[non_seed_core]!=which(non_seed_core)) # total error, including errors in added seeds
  Bm_hard <- B[corr,corr]
  objective_hard <- sum(abs(A-Bm_hard))

  # errors of new added seeds
  new_seeds_err_hard <- sum(aseeds_err)/aseeds

  result_hard <- tibble(match_error_hard=core_err_hard, added_seeds_error_hard=new_seeds_err_hard,
                        objective_hard=objective_hard, match_hard=match_hard)
  result_hard
}
#' @rdname adpt_seeds_match
#' @examples
#' match_soft <- graph_match_soft_seeds(g1, g2, seeds, select_seeds=select_seeds)
#' match_soft
#' @export
#'
graph_match_soft_seeds <- function(A, B, seeds, non_seed_core = !seeds, select_seeds, start = "bari"){
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)

  select_seeds <- check_soft_seeds(select_seeds)
  seed_A <- select_seeds$seed_A
  seed_B <- select_seeds$seed_B
  ns <- sum(seeds)
  nv <- nrow(A)
  aseeds <- length(seed_A)
  nns <- nv-ns

  if(start == "bari"){
    start <- bari_start(nns,ns,select_seeds)
  } else if(start =="rds"){
    start <- rds_sinkhorn_start(nns,ns,select_seeds)
  } else if(start == "convex"){
    match <- graph_match_adpt_seeds(A, B, seeds, non_seed_core, select_seeds, start="convex")$match_hard
    start <- match$D[!seeds,!seeds]
  }

  # graph_match with updated seeds and calculate errors
  match_soft <- graph_match_FW(A, B, seeds = seeds, start = start)
  core_err_soft <- mean(match_soft$corr[non_seed_core]!=which(non_seed_core)) # total error, including errors in added seeds
  Bm_soft <- B[match_soft$corr,match_soft$corr]
  objective_soft <- sum(abs(A-Bm_soft))

  # errors of new added seeds
  seeds_add <- rep(FALSE,nv)
  seeds_add[seed_A] <- TRUE
  new_seeds_err_soft <- mean(match_soft$corr[seeds_add]!=which(seeds_add))

  result_soft <- tibble(match_error_soft=core_err_soft,added_seeds_error_soft=new_seeds_err_soft,
                        objective_soft=objective_soft)
  result_soft
}
#'
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
        }
        else{
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
#'
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
#'
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

  g2_sig <- ifelse(g2_new_real!=1:nv, TRUE, FALSE)
  g2_index <- which(g2_sig==TRUE)
  g2_value <- g2_new_real[g2_sig]
  index_corr <- sapply(g2_index, function(x) which(corr_hard==x))
  corr_hard[index_corr] <- g2_value

  corr_hard
}
#'
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

  corr_hard <- apply(D,1,function(v) which(v==1))
  g2_sig <- ifelse(g2_new_real!=1:nv, TRUE, FALSE)
  g2_index <- which(g2_sig==TRUE)
  g2_value <- g2_new_real[g2_sig]
  index_corr <- sapply(g2_index, function(x) which(corr_hard==x))
  corr_hard[index_corr] <- g2_value

  D <- Matrix::Diagonal(nv)[corr_hard,]
  D
}
