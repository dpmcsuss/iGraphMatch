do_lap <- function(Grad, method, ...){
  n <- nrow(Grad)
  switch(method,
    lapjv = { 
      Grad <- as.matrix(Grad)
      rlapjv::lapjv(Grad, # round(Grad * n ^ 2 * max(Grad)),
        maximize = TRUE)
    },
    lapmod = {
      if( class(Grad) == "splrMatrix" ){
        rlapjv::lapmod(splr_to_sparse(Grad),
          maximize = TRUE)
      } else {
        rlapjv::lapmod(Grad, maximize = TRUE)
      }
    },
    clue = {
      Grad <- as.matrix(Grad)
      Grad <- Grad - min(Grad)
      as.vector(clue::solve_LSAP(Grad,
        maximum = TRUE))
    },
    # sinkhorn = {
    #   lambda <- 10
    #   n_iter <- 20
    #   sinkhorn(exp(lambda * Grad), n_iter)      
    # },
    stop(paste0("The LAP method ", method,
        " is not implemented. Please use one of lapjv, lapmod, or clue."))
  )    
}

set_lap_method <- function(lap_method, totv1, totv2){
  methods <- c("lapmod", "lapjv", "clue") #, "sinkhorn")
  if (!is.null(lap_method) && !(lap_method %in% methods)){
    stop(paste("Unrecognized LAP method:", lap_method,
      "Please use one of:", methods))
  }
  if (is.null(lap_method)){
    if("rlapjv" %in% rownames(utils::installed.packages()) &&
        totv1 / totv2 < 0.5 || totv2 / totv1 < 0.5){
        lap_method <- "lapmod"
    } else {
      lap_method <- "clue"
    }
  }
  lap_method
}