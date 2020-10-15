context("splr matrix class")

x <- Matrix::rsparsematrix(10, 10, .5)
a <- matrix(runif(30), 10)
b <- matrix(runif(30), 10)
ns <- matrix(runif(100), 10)

s <- splr(x, a, b)


test_that("error on no b or rank",
  {
    expect_error(
      splr(x, ns),
      "please provide an already factorized low-rank matrix.*"
    )
  }
)

test_that("error on no b and a is wrong dim",
  {
    expect_error(
      splr(x, a, rank = 2),
      "b is not provided and a is not the same dimension as x"
    )
  }
)

s <- splr(x, ns, rank = 2)


test_that("error on no b",
  {
    expect_error(
      splr(x, ns),
      "please provide an already factorized low-rank matrix.*"
    )
  }
)


test_that("error on wrong dim for a",
  {
    expect_error(
      splr(x, a[1:4,], b),
      "number of rows of x not equal to number of rows of a.*"
    )
  }
)

test_that("error on wrong dim for b",
  {
    expect_error(
      splr(x, a, b[1:4,]),
      "number of columns of x not equal to number of rows of b.*"
    )
  }
)

test_that("error on wrong dim",
  {
    expect_error(
      splr(x, a[,1:2], b),
      "number of columns of a not equal to number of columns of b.*"
    )
  }
)



s <- splr(x, Matrix(a), Matrix(b))
sp <- splr_to_sparse(s)
m <- as.matrix(s)

ss <- splr_sparse_plus_constant(x, 1)
dm <- as(s, "dMatrix")
try({
  sink("/dev/null")
})
print(s)
show(s)
as.character(s)

str(s)
sink()


length(s)


s %*% s
s * s

ml <- matrix_list(list(ns, x))

s %*% ml
x %*% s
ns %*% s
3 * s

ml %*% s
s %*% x
s %*% ns
s * 3

ml * s
s * x
s * ns


s / x
s / ns

s - ml
x - s
ns - s
3 - s

ml - s
s - x
s - ns
s - 3
s - s

s + ml
x + s
ns + s
3 + s

ml + s
s + x
s + ns
s + 3
s + s





s * Diagonal(10)
Diagonal(10) *  s
s - Diagonal(10)
s + Diagonal(10)






test_that("error on non-atomic numeric",
  {
    expect_error(
      s + c(1, 3),
      "Can only add length 1 numerics to splrmatrix"
    )
  }
)


norm(s, "f")
norm(s, "F")
norm(s, "2")


rowSums(s)
colSums(s)
colMeans(s)
rowMeans(s)
sum(s)
mean(s)

s[1:3, 1:3]

test_that("warning on drop",
  {
    expect_warning(
      s[1:3, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another clas first"
    )
  }
)


test_that("warning on drop",
  {
    expect_warning(
      s[, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another clas first"
    )
  }
)


test_that("warning on drop",
  {
    expect_warning(
      s[1:3, , drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another clas first"
    )
  }
)

test_that("warning on drop",
  {
    expect_warning(
      s[1:3, 1:10 < 3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another clas first"
    )
  }
)
test_that("warning on drop",
  {
    expect_warning(
      s[1:10 < 3, 1:3, drop = TRUE],
      "drop = TRUE is ignored for the splrMatrix class. cast to another clas first"
    )
  }
)


s[1:3, 1:3, drop = FALSE]

s[1:3, ]
s[ , 1:3]

s0 <- s[numeric(),]

s0[, 1:3]

s0 <- s[,numeric()]

s0[1:3, ]


s[1:3, 1:10 < 5]
s[1:10 < 5, 1:3]

s[, 1:10 < 5]
s[1:10 < 5, ]
s[1:10 < 5, 1:10 < 3]


s[matrix(c(1,2,3,4), 2)]

s[]

dim(s)
t(s)
diag(s)
as(s, "dgeMatrix")



# #' @rdname splr_constructor
# #' 
# #' 
# #' 
# #' 
# #' @export
# setMethod(
#   f = "splr",
#   signature = signature(x = "Matrix", a = "Matrix", b = "Matrix"),
#   definition = function(x, a, b, ...){
#     new("splrMatrix", x = x, a = a, b = b, Dim = dim(x), Dimnames = list(NULL, NULL))
#   }
# )

# #' Convert splr Matrix to Sparse
# #' 
# #' @param data splr Matrix
# #' 
# #' @return sparse Matrix equal to x + a %*% t(b)
# #' 
# #' See \code{\link[Matrix]{Matrix}}.
# #' 
# #' @export
# splr_to_sparse <- function(data){
#     data@x + Matrix(data@a, sparse = TRUE) %*% Matrix(t(data@b), sparse = TRUE)
# }

# as.matrix.splrMatrix <- function(from,...)  {
#   as.matrix(as.matrix(from@x,...)+from@a%*%t(from@b),...)
  
  
# }

# #' Add a constant to a splr object
# #' 
# #' @param x splr object
# #' @param a scalar
# #' 
# #' @returns new splr object x + a
# #' 
# #' @export
# splr_sparse_plus_constant <- function(x, a){
#   d <- dim(x)
#   splr(x = x, a = rep(a, d[1]), b = rep(1, d[2]), Dim = dim(x), Dimnames = list(NULL, NULL))
# }

# setAs('splrMatrix','dMatrix', function(from) {
#   as(from@x + from@a %*% t(from@b),'Matrix')
# })

# as.character.splrMatrix <- function(from) {
#   paste0("Sparse\n", as.character(from@x), "\n",
#     "Left factor\n", as.character(from@a), "\n",
#     "Right factor\n", as.character(from@b))
# }

# # setMethod("as.character", "splrMatrix", as.character.splrMatrix)

# setAs("splrMatrix", "character", as.character.splrMatrix)

# setAs("splrMatrix", "matrix", function(from) as.matrix.splrMatrix(from))

# #' SPLR Methods
# #' 
# #' Methods for the splr Matrix class. Most behave like
# #' Matrix methods though things like output show the 
# #' decomposition. Use as.matrix to see the computed
# #' dense matrix.
# #' 
# #' @param x As in Matrix
# #' @param ... As in Matrix
# #' @param object As in Matrix
# #' @param e1 As in Matrix
# #' @param y As in Matrix
# #' @param e2 As in Matrix
# #' @param type As in Matrix
# #' @param na.rm As in Matrix
# #' @param dims As in Matrix
# #' @param i As in Matrix
# #' @param j As in Matrix
# #' @param drop As in Matrix
# #' @param value As in Matrix
# #' 
# #' 
# #' @keywords internal
# #' 
# #' @rdname splr
# setMethod("show", signature("splrMatrix"),
#   function(object){
#     cat("Sparse part\n")
#     show(object@x)
#     cat("plus left factor\n")
#     show(object@a)
#     cat("times right factor transpose\n")
#     show(object@b)
# })

# # setMethod("print", signature("splrMatrix"),
# #   function(x){
# #     print(x@x)
# #     print(x@a)
# #     print(x@b)
# # })



# # setMethod("coerce", signature("splrMatrix", "character"),
# #   function(from, to){
# #     if(class == "dMatrix"){
# #       splr_to_sparse(x)
# #     }
# # })



# # #' @rdname splr
# # setMethod('-',
# #   signature(e1 = 'splrMatrix', e2 = 'missing'),
# #   function(e1, e2 = NULL) {
# #     new("splrMatrix", x = -e1@x, a = -e1@a, b = e1@b, Dim = dim(e1@x))
# #   })

# .leftmult = function(x, y){
#   #y is splr, x is a matrix
  
#   a <- y@a
#   sx <- y@x
#   b <- y@b
 
  
#   if(is.null(a)) {
#     x %*% sx
  
#   } else if (is(x,"sparseMatrix")) {
#     p <- as(x%*%sx,"sparseMatrix")
#     anew <- as(x%*%a, "Matrix")
#     b <- as(b,"Matrix")
#     new("splrMatrix", x = p, a = anew, b = b, Dim = dim(x))
#   } else{
#     part1 = x %*% sx
#     part2 = x %*% a
#     part2 = part2 %*% t(b)
#     part1+part2
#   }
# }

# #' @rdname splr
# setMethod("%*%", signature(x = "splrMatrix", y = "splrMatrix"), function(x, y) {
#   new("splrMatrix",
#     x = x@x %*% y@x,
#     a = cbind2(x@a %*% (t(x@b) %*% y@a) + x@x %*% y@a, x@a) ,
#     b = cbind2(y@b, t(y@x) %*% x@b),
#     Dim = dim(x))
# })


# #' @rdname splr
# setMethod("%*%", signature(x = "splrMatrix", y = "matrix_list"), 
#   function(x, y){
#     matrix_list(lapply(seq_along(y), function(i) x %*% y[[i]]))
#   })

# #' @rdname splr
# setMethod("%*%", signature(x = "matrix_list", y = "splrMatrix"), 
#   function(x, y){
#     matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y))
#   })

# #' @rdname splr
# setMethod("%*%", signature(x = "Matrix", y = "splrMatrix"), .leftmult)

# #' @rdname splr
# setMethod("%*%", signature(x = "matrix", y = "splrMatrix"), .leftmult)
# #' @rdname splr
# setMethod("%*%", signature(x = "numeric", y = "splrMatrix"), .leftmult)

# #' @rdname splr
# setMethod("%*%", signature(x = "numLike", y = "splrMatrix"), .leftmult)


# #' @rdname splr
# setMethod("%*%", signature(x ="ANY", y ="splrMatrix"),.leftmult)

# .rightmult = function(x, y){
#   #x is splr, y is matrix
#   a = x@a
#   b = x@b
#   sx = x@x
  
#   if(is.null(a)) {
#     return(sx %*% y)
#   }
  
#   if (is(y,"sparseMatrix")) {
#     newx <- sx%*%y
#     newx <- as(newx,"sparseMatrix")
#     newB <-  t(t(b)%*%y)
#     newB <- as(newB,"Matrix")
#     new("splrMatrix", x = newx, a = a, b = newB, Dim = dim(newx))
#   } else {
    
  
    
#     part1 <- sx %*% y
#     part2 <- t(b) %*% y
#     part2 <- a %*% part2
#     part1+part2
#   }
  
  
# }
# #' @rdname splr
# setMethod("dim", signature(x = "splrMatrix"), function(x) { dim(x@x)})
# #' @rdname splr
# setMethod("length", signature(x = "splrMatrix"), function(x) { length(x@x)})
# #' @rdname splr
# setMethod("%*%", signature(x ="splrMatrix", y ="Matrix"),.rightmult)
# #' @rdname splr
# setMethod("%*%", signature(x ="splrMatrix", y ="matrix"),.rightmult)
# #' @rdname splr
# setMethod("%*%", signature(x ="splrMatrix", y ="numeric"),.rightmult)

# #' @rdname splr
# setMethod("%*%", signature(x ="splrMatrix", y ="numLike"),.rightmult)
# #' @rdname splr
# setMethod("%*%", signature(x ="splrMatrix", y ="ANY"),.rightmult)

# #doesn't return an splr
# #' @rdname splr
# setMethod('*', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'), function(e1, e2) {
#   x <- as(e1@x * e2@x,'sparseMatrix')
#   new("splrMatrix", x = x, a = e1@a %*% t(e1@b) * e2@a, b = e2@b, Dim = dim(x))
# })

# #return sparse
# .multiply <- function(e1, e2) {
#   if (length(e2) == 1) {
#     new('splrMatrix', x = e1@x*e2, a = e2*e1@a, b = e1@b, Dim = dim(e1@x))
#   } else {
#     # can we speed this up for sparse e2?
#     # right now it constructs a fully dense matrix
#     # by calling (e1@a %*% t(e1@b)) which could be bad
#     # if e2 itself is sparse
#     return(e1@x * e2 + (e1@a %*% t(e1@b)) * e2)
#     # the following should be faster
#     # need to test
#     # rank <- ncol(e1@a)
#     # return(e1@x + Reduce("+", lapply(1:rank, function(r){
#     #   Diagonal(e1@a[, r]) %*% e2 %*% Diagonal(e1@b[, r])
#     # }))
#   }
# }

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'Matrix', e2 = 'splrMatrix'), function(e1, e2) {
#   .multiply(e2, e1)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'splrMatrix', e2 = 'ddiMatrix'), function(e1, e2) {
#   .multiply(e1, e2)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'matrix', e2 = 'splrMatrix'), function(e1, e2) {
#   .multiply(e2, e1)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'numeric', e2 = 'splrMatrix'), function(e1, e2) {
#   .multiply(e2, e1)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'ANY', e2 = 'splrMatrix'), function(e1, e2) {
#   .multiply(e2, e1)
# })



# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'splrMatrix', e2 = 'matrix'), function(e1, e2) {
#   .multiply(e1, e2)
# })
# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'splrMatrix', e2 = 'Matrix'), function(e1, e2) {
#   .multiply(e1, e2)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'splrMatrix', e2 = 'numeric'), function(e1, e2) {
#   .multiply(e1, e2)
# })

# #' @rdname splr
# setMethod("*",
#   signature (e1 = 'splrMatrix', e2 = 'ANY'), function(e1, e2) {
#   .multiply(e1, e2)
# })
# #' @rdname splr
# setMethod("/",
#   signature (e1 = 'splrMatrix', e2 = 'matrix'), function(e1, e2) {
#   .multiply(e1, 1/e2)
# })
# #' @rdname splr
# setMethod("/",
#   signature (e1 = 'splrMatrix', e2 = 'Matrix'), function(e1, e2) {
#   .multiply(e1, 1/e2)
# })
# #' @rdname splr
# setMethod("/",
#   signature (e1 = 'splrMatrix', e2 = 'ANY'), function(e1, e2) {
  
#   .multiply(e1, 1/e2)
# })

# .addSplr <- function(e1, e2) {
#   new("splrMatrix", x = as(e1@x + e2@x,"sparseMatrix"), a = cbind2(e1@a, e2@a), b = cbind2(e1@b, e2@b), Dim = dim(e1))
# }

# #' @rdname splr
# setMethod('+', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'),.addSplr)
# #' @rdname splr
# setMethod('-', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'), function(e1, e2) {
#   .addSplr(e1,-e2)
# })



# .leftadd <- function(e1, e2) {
#   #e1 is splr
#   if (is(e2,"sparseMatrix")) {
#     new("splrMatrix", x = as(e1@x + e2,"sparseMatrix"), a = e1@a, b = e1@b, Dim = dim(e2))
#   } else if (is(e2, "splrMatrix")){
#     new("splrMatrix", x = as(e1@x + e2@x,"sparseMatrix"), a = cbind2(e1@a, e2@a), b = cbind(e1@b, e2@b), Dim = dim(e2))
#   } else if( is.numeric(e2) && is.atomic(e2) ){
#     new("splrMatrix", x = as(e1@x, "sparseMatrix"),
#       a = cbind2(e1@a, rep(e2, nrow(e1@a))),
#       b = cbind2(e1@b, rep(1, nrow(e1@b))), Dim = dim(e2))
#   } else {
#     e1@x + e1@a %*% t(e1@b) + e2
#   }
# }

# #' @rdname splr
# setMethod("+", signature(e1 ="splrMatrix", e2 ="Matrix"), function(e1, e2) {
#   .leftadd(e1 = e1, e2 = e2)
# })
# #' @rdname splr
# setMethod("+", signature(e1 ="splrMatrix", e2 ="numeric"), function(e1, e2) {
#   if(!is.atomic(e2)){
#     stop("Can only add atomic numerics to splrmatrix")
#   }
#   .leftadd(e1 = e1, e2 = e2)
# })
# #' @rdname splr
# setMethod("+", signature(e1 ="splrMatrix", e2 ="ANY"), function(e1, e2) {
#   .leftadd(e1 = e1, e2 = e2)
# })

# #' @rdname splr
# setMethod("-", signature(e1 = "splrMatrix", e2 = "missing"),
#   function(e1, e2 = NULL){
#     splr(-e1@x, a = -e1@a, b = e1@a)
#   })


# #' @rdname splr
# setMethod("-", signature(e1 ="splrMatrix", e2 ="Matrix"), 
#   function(e1, e2) {
#     .leftadd(e1 = e1, e2 = -e2)
#   })

# #' @rdname splr
# setMethod("-", signature(e1 ="splrMatrix", e2 ="ddiMatrix"), 
#   function(e1, e2) {
#     .leftadd(e1 = e1, e2 = -e2)
#   })

# #' @rdname splr
# setMethod("-", signature(e1 ="splrMatrix", e2 ="numeric"), function(e1, e2) {
#   if(!is.atomic(e2)){
#     stop("Can only add atomic numerics to splrmatrix")
#   }
#   .leftadd(e1 = e1, e2 = -e2)
# })
# #' @rdname splr
# setMethod("-", signature(e1 ="splrMatrix", e2 ="ANY"), function(e1, e2) {
#   .leftadd(e1 = e1, e2 = -e2)
# })

# .rightadd = function(e1, e2) {

#     if (is(e1, "sparseMatrix")) {
#       splr(as(e2@x + e1,"sparseMatrix"), a = e2@a, b = e2@b)
#     } else if( is.numeric(e1) ){
#       new("splrMatrix", x = as(e1@x,"sparseMatrix"),
#         a = cbind2(e1@a, rep(1, nrow(e1@a))),
#         b = cbind2(e1@b, rep(e2, nrow(e1@b))),
#         Dim = dim(e2))
#     }  else{ #e1 is not sparse
#       e2@x + e2@a %*% t(e2@b) + e1
#     }
# }

# #' @rdname splr
# setMethod("+", signature("Matrix","splrMatrix"), function(e1, e2) {
#   .rightadd(e1, e2)
# })
# #' @rdname splr
# setMethod("+", signature("numeric","splrMatrix"), function(e1, e2) {
#   if(!is.atomic(e1)){
#     stop("Can only add atomic numerics to splrmatrix")
#   }
#   .rightadd(e1, e2)
# })
# #' @rdname splr
# setMethod("+", signature("ANY","splrMatrix"), function(e1, e2) {
#   .rightadd(e1, e2)
# })
# #' @rdname splr
# setMethod("-", signature("Matrix","splrMatrix"), function(e1, e2) {
#   .rightadd(e1,-e2)
# })
# #' @rdname splr
# setMethod("-", signature("numeric","splrMatrix"), function(e1, e2) {
#   if(!is.atomic(e1)){
#     stop("Can only add atomic numerics to splrmatrix")
#   }
#   .rightadd(e1,-e2)
# })
# #' @rdname splr
# setMethod("-", signature("ANY","splrMatrix"), function(e1, e2) {
#   .rightadd(e1,-e2)
# })




# #frobenius norm
# Frobsmlr = function(x, a, b){
 
#     #expansion due to trevor hastie
#     xnorm <- norm(x, type = 'f')
#     xab = as.matrix(x%*%b)
#     xab = sum(xab*a)
#     aa = t(a)%*%a
#     bb = t(b)%*%b
#     ab = sum(aa*bb)
#     sqrt(pmax(0, xnorm^2+2*xab+ab))

  
  
# }


# #' @rdname splr
# setMethod("norm", signature(x ="splrMatrix", type ="character"),
#   function(x, type,...){
#     switch(type,
#       "F" = Frobsmlr(x = x@x, a = x@a, b = x@b),
#       "f" = Frobsmlr(x = x@x, a = x@a, b = x@b),
#       norm(as.matrix(x), type = type,...)
#     )
#   }, valueClass ="numeric")

# #' Matrix inner products
# #' 
# #' @param x matrix like object
# #' @param y matrix like object
# #' 
# #' @returns inner product <x, y>
# #' 
# #' @rdname innerproduct
# #' @export
# setGeneric("innerproduct", function(x, y){
#   sum(x * y)
# })

# #' @rdname innerproduct
# setMethod("innerproduct", signature(x = "splrMatrix", y = "splrMatrix"),
#   function(x, y){
#     sum(diag(t(y@a) %*% x@x %*% y@b)) + 
#       sum(diag(t(x@b) %*% t(y@x) %*% x@a)) +
#       sum(diag( (t(x@b) %*% y@b) %*% (t(y@a) %*% x@a) )) +
#       sum(x@x * y@x)
#   })


# .innerproduct_Matrix <- function(x, y){
#     sum(diag(t(x@b) %*% t(y) %*% x@a)) +
#       sum(x@x * y)
# }

# #' @rdname innerproduct
# setMethod("innerproduct", signature(x = "splrMatrix", y = "Matrix"),
#   function(x, y){ .innerproduct_Matrix(x, y)})

# #' @rdname innerproduct
# setMethod("innerproduct", signature(x = "Matrix", y = "splrMatrix"),
#   function(x, y){ .innerproduct_Matrix(y, x)})

# #' @rdname innerproduct
# setMethod("innerproduct", signature(x = "matrix_list", y = "matrix_list"),
#   function(x, y){
#   sapply(seq_along(x), function(i) innerproduct(x[[i]], y[[i]]))
# })


# #complete
# .rsum <- function(x, ...){
#   #x is splrMatrix matrix
  
#     rx = rowSums(x@x, ...)
#     cb = colSums(x@b, ...)
#     drop(rx+x@a%*%cb) 
  
 
# }
# #' @rdname splr
# setMethod("rowSums",
#   signature(x = "splrMatrix"),
#   .rsum)



# .csum = function(x, ...){
#   #x is splrMatrix matrix
   
#     cx <- colSums(x@x, ...)
#     ca <- colSums(x@a, ...)
#     drop( cx + x@b %*% ca)
  
# }

# #' @rdname splr
# setMethod("colSums",
#   signature(x = "splrMatrix"),
#   .csum)


# #toDo
# .rmean = function(x, ...){
#   #x is splrMatrix matrix
 
#     rx = rowMeans(x@x, ...)
#     cb = colMeans(x@b, ...)
#     drop(rx+x@a%*%cb)
  
# }

# #' @rdname splr
# setMethod("rowMeans",
#   signature(x = "splrMatrix"),
#   .rmean)


# #toDo
# .cmean = function(x, ...){
#   #x is splrMatrix matrix
  
#     cx = colMeans(x@x, ...)
#     ca = colMeans(x@a, ...)
#     drop(cx+x@b%*%ca)
  
# }

# #' @rdname splr
# setMethod("colMeans",
#   signature(x = "splrMatrix"),
#   .cmean)


# .sum <- function(x, ..., na.rm = FALSE){
#   sum(.csum(x), na.rm = na.rm) 
# }

# #' @rdname splr
# setMethod("sum", signature(x = "splrMatrix", na.rm = "ANY"), .sum)

# #' @rdname splr
# setMethod("mean", signature(x = "splrMatrix"), function(x, ...){
#   sum(x, ...) / x@Dim[1] / x@Dim[2]
# })

# #' @rdname splr
# setMethod("[",
#   signature(x = "splrMatrix", 
#     i = 'missing', j = 'missing', drop = 'missing') ,
#   function(x, i = NULL, j = NULL, drop = NULL) {
#           x
#   })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'numeric', drop = 'logical') 
#           , function(x, i, j, ..., drop) {
#             if (drop) {
#               return(drop(x@x[i, j,...] + (x@a[i,] %*% t(x@b)[, j]) ))
#             } else {
#               return(splr(x@x[i, j,..., drop = FALSE], x@a[i,, drop = FALSE], x@b[j,, drop = FALSE]
#                          , Dim = dim(x@x[i, j,..., drop = FALSE])) )
#             }
#     })


# col_index <- function(x, j, ..., drop) {
#   row_dim <- dim(x@x)[1]
#   i <- seq(row_dim)
#   if(row_dim == 0) {
#     i <- numeric()
#   }

#   if (drop) {
#     return(x@x[i, j,...] + (x@a[i, , drop = FALSE] %*% t(x@b)[, j, drop = FALSE]) )
#   } else {
#     return(new("splrMatrix",
#       x = x@x[i, j,..., drop = FALSE],
#       a = x@a[i, , drop = FALSE],
#       b = x@b[j, , drop = FALSE],
#       Dim = dim(x@x[i, j,..., drop = FALSE])) )
#   }
  
  
  
# }

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'numeric', drop = 'logical'),
#             col_index)

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'numeric', drop = 'missing'),
#             function(x, j, ...) col_index(x, j, drop = FALSE))

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'logical', drop = 'logical'),
#             col_index)

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'logical', drop = 'missing'),
#             function(x, j, ...) col_index(x, j, drop = FALSE))

# row_index <- function(x, i, ..., drop) {

#   row_dim <- dim(x@x)[2]
#   j <- seq(row_dim)
#   if(row_dim == 0) {
#     j <- numeric()
#   }

#   if (drop) {
#     return(drop(x@x[i, j,...] + (x@a[i,, drop = FALSE] %*% t(x@b)[, j, drop = FALSE]) ))
#   } else {
#     return( new("splrMatrix", 
#       x = x@x[i, j,...],
#       a = x@a[i, , drop = FALSE], 
#       b = x@b[j, , drop = FALSE],
#       Dim = dim(x = x@x[i, j,...])))
#   }
# }

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'missing', drop = 'logical'),
#             row_index)

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'missing', drop = 'missing'),
#             function(x, i, ...) row_index(x, i, drop = FALSE))

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'missing', drop = 'logical'),
#             row_index)

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'missing', drop = 'missing'),
#             function(x, i, ...) row_index(x, i, drop = FALSE))


# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'ANY', drop ='logical') 
#           , function(x, i, j, ..., drop) {
#             return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
#                        , Dim = dim(x@x[i, j,..., drop = FALSE])))
            
#           })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'logical', drop ='logical') 
#           , function(x, i, j, ..., drop) {
#             return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
#                        , Dim = dim(x@x[i, j,..., drop = FALSE])))
            
#           })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'ANY', drop ='missing') 
#           , function(x, i, j, ..., drop = TRUE) {
            
            
#             return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE],
#                        Dim = dim(x = x@x[i, j,..., drop = FALSE])))
#   })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'ANY', drop = 'ANY') 
#           , function(x, i, j, ..., drop) {
            
#             if (drop) {
#               return(new('splrMatrix', x = x@x[i, j, drop = FALSE], a = x@a[i,], b = x@b[j,, drop = FALSE] ))
#             } else {
#               new('splrMatrix', x = x@x[i, j, drop = FALSE]
#                   , a = as(x@a[i,, drop = FALSE],'Matrix'), b = as(x@b[j,, drop = FALSE],'Matrix')
#                   , Dim = dim(x = x@x[i, j,..., drop = FALSE]))
              
#             }
#  })


# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'ANY', drop = 'missing') 
#           , function(x, i, j, ...) {
            
            
#             new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
#                 , Dim = dim(as(x@x[i, j], "Matrix")))
#     })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'numeric', drop = 'missing') 
#           , function(x, i, j, ...) {
            
            
#             new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
#                 , Dim = dim(as(x@x[i, j], "Matrix")))
#     })

# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'numeric', drop = 'logical') 
#           , function(x, i, j, ..., drop) {
            
#             ret <- new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
#                 , Dim = dim(as(x@x[i, j], "Matrix")))
#             ret
#     })



# #' @rdname splr
# setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'missing', drop = 'missing'),
#   function(x, i = NULL, j = NULL, drop = NULL) {
#     x
#   })


# #can fix this
# #' @rdname splr
# setMethod("[",
#     signature(x = "splrMatrix", i = "matrix", j = "missing", drop = "missing"),
#     function(x, i , ...) {
#           x@x[i] +
#             (x@a[i[, 1], ] * x@b[i[, 2], ]) %*%
#               Matrix(1, dim(x@a)[2], 1)
#           #as.matrix(x@x + x@a%*%t(x@b))[i]
#  })


# #document the issues with doing this

# #' @rdname splr
# setMethod("[<-",
#   signature(x ="splrMatrix", i = 'numeric', j = 'numeric', value = 'ANY'),
#   function(x, i, j, ..., value) {
#     y <- x@x
#     y[i, j] <- value
#     a <- x@a
#     if (length(i) > 0) {
#       a[i,] <- 0
#     }
#     b <- x@b
#     if (length(j) > 0) {
#       b[j,] <- 0
#     }
#     new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
    
#    })

# #' @rdname splr
# setMethod("[<-",
#   signature(x ="splrMatrix", i = 'numeric', j = 'missing', value = 'ANY') ,
#   function(x, i, ..., value) {
#     j <- c(1:dim(x@x)[2])
#     y <- x@x
#     y[i, j] <- value
#     a <- x@a
#     a[i,] <- 0
#     b <- x@b
#     b[j,] <- 0
#     new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
#  })

# #' @rdname splr
# setMethod("[<-",
#   signature(x ="splrMatrix", i = 'missing', j = 'numeric', value = 'ANY') ,
#   function(x, j, ..., value) {
#     i <- c(1:dim(x@x)[1])
#     y <- x@x
#     y[, j] <- value
#     a <- x@a
#     a[i,] <- 0
#     b <- x@b
#     b[j,] <- 0
#     new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
#  })

 
# # Implementing this would be nice
# #' @rdname splr
# setMethod("[<-", signature(x ="Matrix", i = 'ANY', j = 'ANY', value = 'splrMatrix'),
#           function(x, i, j, ..., value) {
#             y <- x
#             y[i, j] <- value@x
#             a <- Matrix(0, dim(x)[1], dim(value@a)[2])
#             b <- Matrix(0, dim(x)[2], dim(value@b)[2])
#             a[i,] <- value@a
#             b[j,] <- value@b
#             new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
#  })



# #' @rdname splr
# setMethod("dim", signature(x = "splrMatrix"),
#           function(x) x@Dim, valueClass = "integer")

# #' @rdname splr
# setMethod('str', signature(object = "splrMatrix"), function(object){
#   cat('splrMatrix')
#   cat('\nDimension: ', dim(object@x))
#   cat('\nLower rank matrix is rank: ', min(dim(object@a)) )
# })


# #' @rdname splr
# setMethod("t", signature = signature(x ="splrMatrix") , function(x) {
#   #splr(t(x@x), x@b, x@a)
 
#   new("splrMatrix", x = t(x@x), a = x@b, b = x@a, Dim = dim(t(x@x)))
# })


# #' @rdname splr
# setMethod("diag", signature = signature(x = "splrMatrix"), function(x) {
#   drop(diag(x@x) + rowSums(x@a * x@b))
# })

# setAs(
#   "splrMatrix", "dgeMatrix",
#   function(from) from@x + Matrix::tcrossprod(from@a, from@b)
# )
