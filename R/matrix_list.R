#' Lists of Matrices
#' 
#' Basically a list with matrix components that are all 
#' the same dimension. Mostly for internal igraphmatch use.
#' It can do various things like arithmetic and
#' indexing, all of which are done component-wise.
#'
#' @param x As in Matrix
#' @param y As in Matrix
#' @param i As in Matrix
#' @param j As in Matrix
#' @param drop As in Matrix
#' @param na.rm As in Matrix
#' @param e1 As in Matrix
#' @param e2 As in Matrix
#' @param ml A list of matrices
#' 
#' @keywords internal
#' 
#' @rdname matrix_list
#' @importClassesFrom Matrix sparseMatrix
#' 
setClass("matrix_list", contains = "list")

#' @rdname matrix_list
#' @export
matrix_list <- matrix_list <- function(ml)
  new("matrix_list", ml)

# ml <- matrix_list(list(matrix(1:9, 3), matrix(0:8, 3)))
# ml %*% ml


#' @rdname matrix_list
setMethod("%*%", signature(x = "matrix_list", y = "matrix_list"), 
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y[[i]]))
  })


#' @rdname matrix_list
setMethod("t", signature(x = "matrix_list"),
  function(x){
    matrix_list(lapply(x, t))
  })

#' @rdname matrix_list
setMethod("dim", signature(x = "matrix_list"),
  function(x) {
    dim(x[[1]])
  }
)


#' @rdname matrix_list
setMethod("[",
  signature(x="matrix_list",i = 'ANY', j = 'ANY', drop = 'ANY'),
  function(x, i = 1:nrow(x[[1]]), j = 1:ncol(x[[1]]), drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[i, j, drop = drop]))
  })


#' @rdname matrix_list
setMethod("[",
  signature(x="matrix_list",i = 'ANY', j = 'missing', drop = 'ANY'),
  function(x, i, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[i, , drop = drop]))
  })


#' @rdname matrix_list
setMethod("[",
  signature(x="matrix_list", i = 'missing', j = 'ANY', drop = 'ANY'),
  function(x, j, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[, j, drop = drop]))
  })


#' @rdname matrix_list
setMethod("[",
  signature(x="matrix_list",i = 'missing', j = 'missing', drop = 'ANY'),
  function(x, drop = FALSE) {
    matrix_list(lapply(x, function(xl) xl[, , drop = drop]))
  })


#' @rdname matrix_list
setMethod("%*%", signature(x = "matrix_list", y = "ANY"), 
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y))
  })


#' @rdname matrix_list
setMethod("%*%", signature(x = "ANY", y = "matrix_list"), 
  function(x, y){
    matrix_list(lapply(seq_along(y), function(i) x %*% y[[i]]))
  })

#' @rdname matrix_list
setMethod("%*%", signature(x = "matrix_list", y = "Matrix"), 
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y))
  })


#' @rdname matrix_list
setMethod("%*%", signature(x = "Matrix", y = "matrix_list"), 
  function(x, y){
    matrix_list(lapply(seq_along(y), function(i) x %*% y[[i]]))
  })


#' @rdname matrix_list
setMethod("sum", signature(x = "matrix_list", na.rm = "logical"),
  function(x, na.rm = FALSE){
    sum(sapply(x, sum, na.rm = na.rm), na.rm = na.rm)
  })


#' @rdname matrix_list
setMethod("^", signature(e1 = "matrix_list", e2 = "ANY"), 
  function(e1, e2){
    matrix_list(lapply(e1, function(x) x^e2))
  })


#' @rdname matrix_list
setMethod("-", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] - e2[[i]]))
  })



#' @rdname matrix_list
setMethod("+", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] + e2[[i]]))
  })



#' @rdname matrix_list
setMethod("*", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] * e2[[i]]))
  })



#' @rdname matrix_list
setMethod("/", signature(e1 = "matrix_list", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] / e2[[i]]))
  })


#########################
#' @rdname matrix_list
setMethod("-", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] - e2))
  })



#' @rdname matrix_list
setMethod("+", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] + e2))
  })



#' @rdname matrix_list
setMethod("*", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] * e2))
  })



#' @rdname matrix_list
setMethod("/", signature(e1 = "matrix_list", e2 = "ANY"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e1), function(i) e1[[i]] / e2))
  })


#########################
#' @rdname matrix_list
setMethod("-", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 - e2[[i]]))
  })



#' @rdname matrix_list
setMethod("+", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 + e2[[i]]))
  })


#' @rdname matrix_list
setMethod("*", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 * e2[[i]]))
  })



#' @rdname matrix_list
setMethod("/", signature(e1 = "ANY", e2 = "matrix_list"),
  function(e1, e2){
    matrix_list(lapply(seq_along(e2), function(i) e1 / e2[[i]]))
  })
#########################


#' @rdname matrix_list
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


#' @rdname matrix_list
setMethod("names<-", signature(x = "matrix_list", value = "ANY"),
  function(x, value) {
    x <- as.list(x)
    names(x) <- value
    matrix_list(x)
  }
)

# setMethod("%*%", signature(x = "Matrix", y = "splrMatrix"), .leftmult)

# setMethod("%*%", signature(x = "matrix", y = "splrMatrix"), .leftmult)
# setMethod("%*%", signature(x = "numeric", y = "splrMatrix"), .leftmult)

# setMethod("%*%", signature(x = "numLike", y = "splrMatrix"), .leftmult)


# setMethod("%*%",signature(x="ANY",y="splrMatrix"),.leftmult)
