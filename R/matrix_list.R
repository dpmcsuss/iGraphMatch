matrix_list <- setClass("matrix_list", contains = "list")

# ml <- matrix_list(list(matrix(1:9, 3), matrix(0:8, 3)))
# ml %*% ml

setMethod("%*%", signature(x = "matrix_list", y = "matrix_list"), 
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y[[i]]))
  })

setMethod("t", signature(x = "matrix_list"),
  function(x){
    matrix_list(lapply(x, t))
  })



setMethod("[",
  signature(x="matrix_list",i = 'ANY', j = 'ANY', drop = 'ANY'),
  function(x, i = 1:nrow(x[[1]]), j = 1:ncol(x[[1]]), drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[i, j, drop]))
  })

setMethod("[",
  signature(x="matrix_list",i = 'ANY', j = 'missing', drop = 'ANY'),
  function(x, i, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[i, , drop]))
  })

setMethod("[",
  signature(x="matrix_list", i = 'missing', j = 'ANY', drop = 'ANY'),
  function(x, j, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[, j, drop]))
  })

setMethod("[",
  signature(x="matrix_list",i = 'missing', j = 'missing', drop = 'ANY'),
  function(x, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[, , drop]))
  })

setMethod("%*%", signature(x = "matrix_list", y = "ANY"), 
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y))
  })

setMethod("%*%", signature(x = "ANY", y = "matrix_list"), 
  function(x, y){
    matrix_list(lapply(seq_along(y), function(i) x %*% y[[i]]))
  })

setMethod("sum", signature(x = "matrix_list", na.rm = "logical"),
  function(x, na.rm = FALSE){
    sum(sapply(x, sum, na.rm = na.rm), na.rm = na.rm)
  })

setMethod("^", signature(e1 = "matrix_list", e2 = "ANY"), 
  function(e1, e2){
    matrix_list(lapply(e1, function(x) x^e2))
  })

setMethod("-", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] - e2[[i]]))
  })


setMethod("+", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] + e2[[i]]))
  })


setMethod("*", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] * e2[[i]]))
  })


setMethod("/", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] / e2[[i]]))
  })

#########################
setMethod("-", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] - e2))
  })


setMethod("+", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] + e2))
  })


setMethod("*", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] * e2))
  })


setMethod("/", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] / e2))
  })

#########################
setMethod("-", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 - e2[[i]]))
  })


setMethod("+", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 + e2[[i]]))
  })

setMethod("*", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 * e2[[i]]))
  })


setMethod("/", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 / e2[[i]]))
  })
#########################

setMethod("-", signature(e1 = "matrix_list", e2 = "missing"),
  function(e1,e2){
    matrix_list(lapply(e1, function(x) - x))
  })




ml_sum <- function(x){
  if (is(x, "matrix_list")) {
    Reduce("+", x)
  } else {
    x
  }
}


# setMethod("%*%", signature(x = "Matrix", y = "splrMatrix"), .leftmult)

# setMethod("%*%", signature(x = "matrix", y = "splrMatrix"), .leftmult)
# setMethod("%*%", signature(x = "numeric", y = "splrMatrix"), .leftmult)

# setMethod("%*%", signature(x = "numLike", y = "splrMatrix"), .leftmult)


# setMethod("%*%",signature(x="ANY",y="splrMatrix"),.leftmult)
