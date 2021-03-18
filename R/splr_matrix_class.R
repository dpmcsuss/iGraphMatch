#' @title Sparse Plus Low-Rank Matrices
#'
#' @description An 'S4' class for efficient computation with sparse plus
#' low-rank matrices. Stores sparse plus low-rank matrices
#' (e.g. from matrix factorization or centering graphs)
#' of the form \code{x + a \%*\% t(b)} for faster
#' computation.
#'
#' @slot x a sparse matrix
#' @slot a a low-rank factor or a matrix
#' @slot b optional. a low-rank factor for \code{a \%*\% t(b)}. if \code{b} is not provided, a will be factorized using
#' \code{\link[irlba]{irlba}} provided \code{factorize = TRUE}
#'
#' @param x as in 'Matrix'
#' @param a as in 'Matrix'
#' @param b as in 'Matrix'
#' @param ... as in 'Matrix'
#' 
#' @return splrMatrix object
#'
#' @seealso Methods are documented in \code{\link{splr}}.
#'
#' @rdname splr_constructor
#'
#' @import Matrix
#' @import methods
#' @importClassesFrom Matrix sparseMatrix
setClass("splrMatrix",
  slots = c(x = "sparseMatrix", a = "Matrix", b = "Matrix"),
  contains = 'sparseMatrix')

#' @rdname splr_constructor
#'
#' @param rank rank of the matrix to be factorized.
#' @param dimnames optional - the list of names for the matrix
#' 
#' @return splrMatrix object
#'
#' @export
setGeneric(
  name = "splr",
  def = function(x, a = NULL, b = NULL, rank = NULL,
  dimnames = list(NULL, NULL), ...) {

  standardGeneric("splr")
  #
  # x+ab'
  # x is mxn
  # and is an mxr and b is a rxm matrix
  dx <- dim(x)

  if (is.null(b)) {
    if (is.null(rank)) {
      stop("please provide an already factorized low-rank matrix (or specify the rank parameter)")
    }
    da <- dim(a)
    if (da[1] != dx[1] || da[2] != dx[2]) {
      stop('b is not provided and a is not the same dimension as x')
    }

    temp <- irlba::irlba(a, rank)

    #take the truncated svd
    newA <- as.matrix(temp$u[, 1:rank]) %*% diag(temp$d)[1:rank, 1:rank]
    newB <- as.matrix(temp$v[, 1:rank])
    new("splrMatrix", x = x, a = as(newA,"Matrix"), b = as(newB,'Matrix'), Dim = dim(x), Dimnames = dimnames)

  } else { #b is not null
    a <- as(a,'Matrix')
    b <- as(b,'Matrix')
    da = dim(a)
    db = dim(b)
    if(da[1]!= dx[1]) {
      stop("number of rows of x not equal to number of rows of a")
    }
    if(db[1]!= dx[2]) {
      stop("number of columns of x not equal to number of rows of b\nhint: b needs to be nxr when a is mxr and x is mxn")
    }
    if(da[2]!= db[2]) {
      stop("number of columns of a not equal to number of columns of b")
    }

    new("splrMatrix", x = x, a = a, b = b, Dim = dim(x), Dimnames = dimnames)
  }
})


#' @rdname splr_constructor
#'
#'
#'
#'
#' @export
setMethod(
  f = "splr",
  signature = signature(x = "Matrix", a = "Matrix", b = "Matrix"),
  definition = function(x, a, b, ...){
    new("splrMatrix", x = x, a = a, b = b, Dim = dim(x), Dimnames = list(NULL, NULL))
  }
)

#' Convert splr 'Matrix' to Sparse
#'
#' @param data splrMatrix
#'
#' @return sparse Matrix equal to x + a %*% t(b)
#'
#' See \code{\link[Matrix]{Matrix}}.
#'
#' @export
splr_to_sparse <- function(data){
    data@x + Matrix(data@a, sparse = TRUE) %*% Matrix(t(data@b), sparse = TRUE)
}

as.matrix.splrMatrix <- function(from,...)  {
  as.matrix(as.matrix(from@x,...)+from@a%*%t(from@b),...)


}

#' Add a constant to a splrMatrix object
#'
#' @param x splrMatrix object
#' @param a scalar
#'
#' @return new splrMatrix object x + a
#'
#' @export
splr_sparse_plus_constant <- function(x, a){
  d <- dim(x)
  splr(x = x, a = rep(a, d[1]), b = rep(1, d[2]), Dim = dim(x), Dimnames = list(NULL, NULL))
}

setAs('splrMatrix','dMatrix', function(from) {
  as(from@x + from@a %*% t(from@b),'Matrix')
})

as.character.splrMatrix <- function(from) {
  paste0("Sparse\n", as.character(from@x), "\n",
    "Left factor\n", as.character(from@a), "\n",
    "Right factor\n", as.character(from@b))
}

# setMethod("as.character", "splrMatrix", as.character.splrMatrix)

setAs("splrMatrix", "character", as.character.splrMatrix)

setAs("splrMatrix", "matrix", function(from) as.matrix.splrMatrix(from))

#' @title 'SPLR' Methods
#'
#' @description Methods for the splrMatrix class. Most behave like
#' Matrix methods though things like output show the
#' decomposition. Use as.matrix to see the computed
#' dense matrix.
#'
#' @param x As in 'Matrix'
#' @param ... As in 'Matrix'
#' @param object As in 'Matrix'
#' @param e1 As in 'Matrix'
#' @param y As in 'Matrix'
#' @param e2 As in 'Matrix'
#' @param type As in 'Matrix'
#' @param na.rm As in 'Matrix'
#' @param dims As in 'Matrix'
#' @param i As in 'Matrix'
#' @param j As in 'Matrix'
#' @param drop As in 'Matrix'
#' @param value As in 'Matrix'
#'
#'
#' @keywords internal
#'
#' @return Results of matrix operations for splrMatrix objects.
#'  Attempts are made such that the returned object is stored efficiently,
#'  either as a splrMatrix or sparse Matrix.
#'
#' @rdname splr
setMethod("show", signature("splrMatrix"),
  function(object){
    cat("Sparse part\n")
    show(object@x)
    cat("plus left factor\n")
    show(object@a)
    cat("times right factor transpose\n")
    show(object@b)
})

#' @rdname splr
setMethod("print", signature("splrMatrix"),
  function(x){
    cat("Sparse part\n")
    show(x@x)
    cat("plus left factor\n")
    show(x@a)
    cat("times right factor transpose\n")
    show(x@b)
})

# setMethod("print", signature("splrMatrix"),
#   function(x){
#     print(x@x)
#     print(x@a)
#     print(x@b)
# })



# setMethod("coerce", signature("splrMatrix", "character"),
#   function(from, to){
#     if(class == "dMatrix"){
#       splr_to_sparse(x)
#     }
# })



# #' @rdname splr
# setMethod('-',
#   signature(e1 = 'splrMatrix', e2 = 'missing'),
#   function(e1, e2 = NULL) {
#     new("splrMatrix", x = -e1@x, a = -e1@a, b = e1@b, Dim = dim(e1@x))
#   })

.leftmult = function(x, y){
  #y is splr, x is a matrix

  a <- y@a
  sx <- y@x
  b <- y@b
  if (is(x,"sparseMatrix")) {
    p <- as(x%*%sx,"sparseMatrix")
    anew <- as(x%*%a, "Matrix")
    b <- as(b,"Matrix")
    new("splrMatrix", x = p, a = anew, b = b, Dim = dim(x))
  } else{
    part1 = x %*% sx
    part2 = x %*% a
    part2 = part2 %*% t(b)
    part1+part2
  }
}

#' @rdname splr
setMethod("%*%", signature(x = "splrMatrix", y = "splrMatrix"), function(x, y) {
  new("splrMatrix",
    x = x@x %*% y@x,
    a = cbind2(x@a %*% (t(x@b) %*% y@a) + x@x %*% y@a, x@a) ,
    b = cbind2(y@b, t(y@x) %*% x@b),
    Dim = dim(x))
})


#' @rdname splr
setMethod("%*%", signature(x = "splrMatrix", y = "matrix_list"),
  function(x, y){
    matrix_list(lapply(seq_along(y), function(i) x %*% y[[i]]))
  })

#' @rdname splr
setMethod("%*%", signature(x = "matrix_list", y = "splrMatrix"),
  function(x, y){
    matrix_list(lapply(seq_along(x), function(i) x[[i]] %*% y))
  })

#' @rdname splr
setMethod("%*%", signature(x = "Matrix", y = "splrMatrix"), .leftmult)

#' @rdname splr
setMethod("%*%", signature(x = "matrix", y = "splrMatrix"), .leftmult)
#' @rdname splr
setMethod("%*%", signature(x = "numeric", y = "splrMatrix"), .leftmult)

#' @rdname splr
setMethod("%*%", signature(x = "numLike", y = "splrMatrix"), .leftmult)


#' @rdname splr
setMethod("%*%", signature(x ="ANY", y ="splrMatrix"),.leftmult)

.rightmult = function(x, y){
  #x is splr, y is matrix
  a <- x@a
  b <- x@b
  sx <- x@x

  if (is(y, "sparseMatrix")) {
    newx <- sx%*%y
    newx <- as(newx, "sparseMatrix")
    newB <-  t(t(b)%*%y)
    newB <- as(newB,"Matrix")
    new("splrMatrix", x = newx, a = a, b = newB, Dim = dim(newx))
  } else {



    part1 <- sx %*% y
    part2 <- t(b) %*% y
    part2 <- a %*% part2
    part1+part2
  }


}
#' @rdname splr
setMethod("dim", signature(x = "splrMatrix"), function(x) { dim(x@x)})
#' @rdname splr
setMethod("length", signature(x = "splrMatrix"), function(x) { length(x@x)})
#' @rdname splr
setMethod("%*%", signature(x ="splrMatrix", y ="Matrix"),.rightmult)
#' @rdname splr
setMethod("%*%", signature(x ="splrMatrix", y ="matrix"),.rightmult)
#' @rdname splr
setMethod("%*%", signature(x ="splrMatrix", y ="numeric"),.rightmult)

#' @rdname splr
setMethod("%*%", signature(x ="splrMatrix", y ="numLike"),.rightmult)
#' @rdname splr
setMethod("%*%", signature(x ="splrMatrix", y ="ANY"),.rightmult)

#doesn't return an splr
#' @rdname splr
setMethod('*', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e1, e2)
})

#return sparse
.multiply <- function(e1, e2) {
  if (length(e2) == 1) {
    new('splrMatrix', x = e1@x * e2, a = e2 * e1@a, b = e1@b, Dim = dim(e1@x))
  } else {
    # can we speed this up for sparse e2?
    # right now it constructs a fully dense matrix
    # by calling (e1@a %*% t(e1@b)) which could be bad
    # if e2 itself is sparse
    # return(e1@x * e2 + (e1@a %*% t(e1@b)) * e2)
    # the following should be faster
    # need to test
    rank <- ncol(e1@a)
    return(e1@x * e2 + Reduce("+", lapply(1:rank, function(r){
      Diagonal(x = e1@a[, r]) %*% e2 %*% Diagonal(x = e1@b[, r])
    })))
  }
}

#' @rdname splr
setMethod("*",
  signature (e1 = 'Matrix', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e2, e1)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'splrMatrix', e2 = 'ddiMatrix'), function(e1, e2) {
  .multiply(e1, e2)
})


#' @rdname splr
setMethod("*",
  signature (e1 = 'ddiMatrix', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e2, e1)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'matrix', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e2, e1)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'numeric', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e2, e1)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'ANY', e2 = 'splrMatrix'), function(e1, e2) {
  .multiply(e2, e1)
})



#' @rdname splr
setMethod("*",
  signature (e1 = 'splrMatrix', e2 = 'matrix'), function(e1, e2) {
  .multiply(e1, e2)
})
#' @rdname splr
setMethod("*",
  signature (e1 = 'splrMatrix', e2 = 'Matrix'), function(e1, e2) {
  .multiply(e1, e2)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'splrMatrix', e2 = 'numeric'), function(e1, e2) {
  .multiply(e1, e2)
})

#' @rdname splr
setMethod("*",
  signature (e1 = 'splrMatrix', e2 = 'ANY'), function(e1, e2) {
  .multiply(e1, e2)
})
#' @rdname splr
setMethod("/",
  signature (e1 = 'splrMatrix', e2 = 'matrix'), function(e1, e2) {
  .multiply(e1, 1/e2)
})
#' @rdname splr
setMethod("/",
  signature (e1 = 'splrMatrix', e2 = 'Matrix'), function(e1, e2) {
  .multiply(e1, 1/e2)
})
#' @rdname splr
setMethod("/",
  signature (e1 = 'splrMatrix', e2 = 'ANY'), function(e1, e2) {

  .multiply(e1, 1/e2)
})

.addSplr <- function(e1, e2) {
  new("splrMatrix", x = as(e1@x + e2@x,"sparseMatrix"), a = cbind2(e1@a, e2@a), b = cbind2(e1@b, e2@b), Dim = dim(e1))
}

#' @rdname splr
setMethod('+', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'),.addSplr)
#' @rdname splr
setMethod('-', signature = signature(e1 = 'splrMatrix', e2 = 'splrMatrix'), function(e1, e2) {
  .addSplr(e1,-e2)
})



.leftadd <- function(e1, e2) {
  #e1 is splr
  if (is(e2,"sparseMatrix")) {
    new("splrMatrix", x = as(e1@x + e2,"sparseMatrix"), a = e1@a, b = e1@b, Dim = dim(e2))
  } else if( is.numeric(e2) && length(e2) == 1 ){
    new("splrMatrix", x = as(e1@x, "sparseMatrix"),
      a = cbind2(e1@a, rep(e2, nrow(e1@a))),
      b = cbind2(e1@b, rep(1, nrow(e1@b))), Dim = dim(e1))
  } else if( is.numeric(e2) ) {
    stop("Can only add length 1 numerics to splrmatrix")
  } else {
    e1@x + e1@a %*% t(e1@b) + e2
  }
}

#' @rdname splr
setMethod("+", signature(e1 ="splrMatrix", e2 ="Matrix"), function(e1, e2) {
  .leftadd(e1 = e1, e2 = e2)
})
#' @rdname splr
setMethod("+", signature(e1 ="splrMatrix", e2 ="numeric"), function(e1, e2) {
  .leftadd(e1 = e1, e2 = e2)
})
#' @rdname splr
setMethod("+", signature(e1 ="splrMatrix", e2 ="ANY"), function(e1, e2) {
  .leftadd(e1 = e1, e2 = e2)
})

#' @rdname splr
setMethod("-", signature(e1 = "splrMatrix", e2 = "missing"),
  function(e1, e2 = NULL){
    splr(-e1@x, a = -e1@a, b = e1@a)
  })


#' @rdname splr
setMethod("-", signature(e1 ="splrMatrix", e2 ="Matrix"),
  function(e1, e2) {
    .leftadd(e1 = e1, e2 = -e2)
  })

#' @rdname splr
setMethod("-", signature(e1 ="splrMatrix", e2 ="ddiMatrix"),
  function(e1, e2) {
    .leftadd(e1 = e1, e2 = -e2)
  })

#' @rdname splr
setMethod("-", signature(e1 ="splrMatrix", e2 ="numeric"), function(e1, e2) {
  .leftadd(e1 = e1, e2 = -e2)
})
#' @rdname splr
setMethod("-", signature(e1 ="splrMatrix", e2 ="ANY"), function(e1, e2) {
  .leftadd(e1 = e1, e2 = -e2)
})


#' @rdname splr
setMethod("+", signature("Matrix","splrMatrix"), function(e1, e2) {
  .leftadd(e2, e1)
})
#' @rdname splr
setMethod("+", signature("numeric","splrMatrix"), function(e1, e2) {
  .leftadd(e2, e1)
})
#' @rdname splr
setMethod("+", signature("ANY","splrMatrix"), function(e1, e2) {
  .leftadd(e2, e1)
})
#' @rdname splr
setMethod("-", signature("Matrix","splrMatrix"), function(e1, e2) {
  .leftadd(e2, e1)
})
#' @rdname splr
setMethod("-", signature("numeric", "splrMatrix"), function(e1, e2) {
  .leftadd(-e2, e1)
})
#' @rdname splr
setMethod("-", signature("ANY","splrMatrix"), function(e1, e2) {
  .leftadd(-e2, e1)
})




#frobenius norm
Frobsmlr = function(x, a, b){

    #expansion due to trevor hastie
    xnorm <- norm(x, type = 'f')
    xab = as.matrix(x%*%b)
    xab = sum(xab*a)
    aa = t(a)%*%a
    bb = t(b)%*%b
    ab = sum(aa*bb)
    sqrt(pmax(0, xnorm^2+2*xab+ab))



}


#' @rdname splr
setMethod("norm", signature(x ="splrMatrix", type ="character"),
  function(x, type,...){
    switch(type,
      "F" = Frobsmlr(x = x@x, a = x@a, b = x@b),
      "f" = Frobsmlr(x = x@x, a = x@a, b = x@b),
      norm(as.matrix(x), type = type,...)
    )
  }, valueClass ="numeric")

#' Matrix inner products
#'
#' @param x matrix like object
#' @param y matrix like object
#'
#' @returns inner product <x, y>
#'
#' @rdname innerproduct
#' @export
setGeneric("innerproduct", function(x, y){
  standardGeneric("innerproduct")
  sum(x * y)
})

#' @rdname innerproduct
setMethod("innerproduct", signature(x = "splrMatrix", y = "splrMatrix"),
  function(x, y){
    sum(diag(t(y@a) %*% x@x %*% y@b)) +
      sum(diag(t(x@b) %*% t(y@x) %*% x@a)) +
      sum(diag( (t(x@b) %*% y@b) %*% (t(y@a) %*% x@a) )) +
      sum(x@x * y@x)
  })


.innerproduct_Matrix <- function(x, y){
    sum(diag(t(x@b) %*% t(y) %*% x@a)) +
      sum(x@x * y)
}

#' @rdname innerproduct
setMethod("innerproduct", signature(x = "splrMatrix", y = "Matrix"),
  function(x, y){ .innerproduct_Matrix(x, y)})

#' @rdname innerproduct
setMethod("innerproduct", signature(x = "Matrix", y = "splrMatrix"),
  function(x, y){ .innerproduct_Matrix(y, x)})

#' @rdname innerproduct
setMethod(
  "innerproduct",
  signature(x = "matrix_list", y = "matrix_list"),
  function(x, y){
    sum(sapply(
      seq_along(x),
      function(i) innerproduct(x[[i]], y[[i]])
    ))
  }
)


#complete
.rsum <- function(x, ...){
  #x is splrMatrix matrix

    rx = rowSums(x@x, ...)
    cb = colSums(x@b, ...)
    drop(rx+x@a%*%cb)


}
#' @rdname splr
setMethod("rowSums",
  signature(x = "splrMatrix"),
  .rsum)



.csum = function(x, ...){
  #x is splrMatrix matrix

    cx <- colSums(x@x, ...)
    ca <- colSums(x@a, ...)
    drop( cx + x@b %*% ca)

}

#' @rdname splr
setMethod("colSums",
  signature(x = "splrMatrix"),
  .csum)


#toDo
.rmean = function(x, ...){
  #x is splrMatrix matrix

    rx = rowMeans(x@x, ...)
    cb = colMeans(x@b, ...)
    drop(rx+x@a%*%cb)

}

#' @rdname splr
setMethod("rowMeans",
  signature(x = "splrMatrix"),
  .rmean)


#toDo
.cmean = function(x, ...){
  #x is splrMatrix matrix

    cx = colMeans(x@x, ...)
    ca = colMeans(x@a, ...)
    drop(cx+x@b%*%ca)

}

#' @rdname splr
setMethod("colMeans",
  signature(x = "splrMatrix"),
  .cmean)


.sum <- function(x, ..., na.rm = FALSE){
  sum(.csum(x), na.rm = na.rm)
}

#' @rdname splr
setMethod("sum", signature(x = "splrMatrix", na.rm = "ANY"), .sum)

#' @rdname splr
setMethod("mean", signature(x = "splrMatrix"), function(x, ...){
  sum(x, ...) / x@Dim[1] / x@Dim[2]
})

#' @rdname splr
setMethod("[",
  signature(x = "splrMatrix",
    i = 'missing', j = 'missing', drop = 'missing') ,
  function(x, i = NULL, j = NULL, drop = NULL) {
          x
  })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'numeric', drop = 'logical')
          , function(x, i, j, ..., drop) {
            if (drop) {
              warning("drop = TRUE is ignored for the splrMatrix class. cast to another clas first")
            }
            return(splr(x@x[i, j], x@a[i, ], x@b[j,]
                         , Dim = dim(x@x[i, j])) )
    })


col_index <- function(x, j, ..., drop) {
  row_dim <- dim(x@x)[1]
  i <- seq(row_dim)
  if(row_dim == 0) {
    i <- numeric()
  }

  if (drop) {
   warning("drop = TRUE is ignored for the splrMatrix class. cast to another clas first")
  }
  return(new(
    "splrMatrix",
    x = x@x[i, j],
    a = x@a[i, , drop = FALSE],
    b = x@b[j, , drop = FALSE],
    Dim = dim(x@x[i, j,..., drop = FALSE])
  ))


}

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'numeric', drop = 'logical'),
            col_index)

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'numeric', drop = 'missing'),
            function(x, j, ...) col_index(x, j, drop = FALSE))

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'logical', drop = 'logical'),
            col_index)

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'logical', drop = 'missing'),
            function(x, j, ...) col_index(x, j, drop = FALSE))

row_index <- function(x, i, ..., drop) {

  row_dim <- dim(x@x)[2]
  j <- seq(row_dim)
  if(row_dim == 0) {
    j <- numeric()
  }


  if (drop) {
   warning("drop = TRUE is ignored for the splrMatrix class. cast to another clas first")
  }
  return( new("splrMatrix",
    x = x@x[i, j,...],
    a = x@a[i, , drop = FALSE],
    b = x@b[j, , drop = FALSE],
    Dim = dim(x = x@x[i, j,...])))
}

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'missing', drop = 'logical'),
            row_index)

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'missing', drop = 'missing'),
            function(x, i, ...) row_index(x, i, drop = FALSE))

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'missing', drop = 'logical'),
            row_index)

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'missing', drop = 'missing'),
            function(x, i, ...) row_index(x, i, drop = FALSE))


#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'ANY', drop ='logical')
          , function(x, i, j, ..., drop) {
            return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
                       , Dim = dim(x@x[i, j,..., drop = FALSE])))

          })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'logical', drop ='logical')
          , function(x, i, j, ..., drop) {
            if (drop) {
              warning("drop = TRUE is ignored for the splrMatrix class. cast to another clas first")
            }
            return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
                       , Dim = dim(x@x[i, j,..., drop = FALSE])))

          })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'numeric', j = 'ANY', drop ='missing')
          , function(x, i, j, ..., drop = FALSE) {
            return(splr(x = x@x[i, j,..., drop = FALSE], a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE],
                       Dim = dim(x = x@x[i, j,..., drop = FALSE])))
  })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'ANY', drop = 'ANY')
          , function(x, i, j, ..., drop) {
            new('splrMatrix', x = x@x[i, j, drop = FALSE]
                , a = as(x@a[i,, drop = FALSE],'Matrix'), b = as(x@b[j,, drop = FALSE],'Matrix')
                , Dim = dim(x = x@x[i, j,..., drop = FALSE]))


 })


#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'ANY', drop = 'missing')
          , function(x, i, j, ...) {


            new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
                , Dim = dim(as(x@x[i, j], "Matrix")))
    })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'numeric', drop = 'missing')
          , function(x, i, j, ...) {


            new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
                , Dim = dim(as(x@x[i, j], "Matrix")))
    })

#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'logical', j = 'numeric', drop = 'logical')
          , function(x, i, j, ..., drop) {
            if (drop) {
              warning("drop = TRUE is ignored for the splrMatrix class. cast to another clas first")
            }
            ret <- new('splrMatrix', x = as(x@x[i, j], "Matrix"), a = x@a[i,, drop = FALSE], b = x@b[j,, drop = FALSE]
                , Dim = dim(as(x@x[i, j], "Matrix")))
            ret
    })



#' @rdname splr
setMethod("[", signature(x ="splrMatrix", i = 'missing', j = 'missing', drop = 'missing'),
  function(x, i = NULL, j = NULL, drop = NULL) {
    x
  })


#can fix this
#' @rdname splr
setMethod("[",
    signature(x = "splrMatrix", i = "matrix", j = "missing", drop = "missing"),
    function(x, i , ...) {
          x@x[i] +
            (x@a[i[, 1], ] * x@b[i[, 2], ]) %*%
              Matrix(1, dim(x@a)[2], 1)
          #as.matrix(x@x + x@a%*%t(x@b))[i]
 })


#document the issues with doing this

#' @rdname splr
setMethod("[<-",
  signature(x ="splrMatrix", i = 'numeric', j = 'numeric', value = 'ANY'),
  function(x, i, j, ..., value) {
    y <- x@x
    y[i, j] <- value
    a <- x@a
    if (length(i) > 0) {
      a[i,] <- 0
    }
    b <- x@b
    if (length(j) > 0) {
      b[j,] <- 0
    }
    new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))

   })

#' @rdname splr
setMethod("[<-",
  signature(x ="splrMatrix", i = 'numeric', j = 'missing', value = 'ANY') ,
  function(x, i, ..., value) {
    j <- c(1:dim(x@x)[2])
    y <- x@x
    y[i, j] <- value
    a <- x@a
    a[i,] <- 0
    b <- x@b
    b[j,] <- 0
    new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
 })

#' @rdname splr
setMethod("[<-",
  signature(x ="splrMatrix", i = 'missing', j = 'numeric', value = 'ANY') ,
  function(x, j, ..., value) {
    i <- c(1:dim(x@x)[1])
    y <- x@x
    y[, j] <- value
    a <- x@a
    a[i,] <- 0
    b <- x@b
    b[j,] <- 0
    new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
 })


# Implementing this would be nice
#' @rdname splr
setMethod("[<-", signature(x ="Matrix", i = 'ANY', j = 'ANY', value = 'splrMatrix'),
          function(x, i, j, ..., value) {
            y <- x
            y[i, j] <- value@x
            a <- Matrix(0, dim(x)[1], dim(value@a)[2])
            b <- Matrix(0, dim(x)[2], dim(value@b)[2])
            a[i,] <- value@a
            b[j,] <- value@b
            new("splrMatrix", x = y, a = a, b = b, Dim = dim(y))
 })



#' @rdname splr
setMethod("dim", signature(x = "splrMatrix"),
          function(x) x@Dim, valueClass = "integer")

#' @rdname splr
setMethod('str', signature(object = "splrMatrix"), function(object){
  cat('splrMatrix')
  cat('\nDimension: ', dim(object@x))
  cat('\nLower rank matrix is rank: ', min(dim(object@a)) )
})


#' @rdname splr
setMethod("t", signature = signature(x ="splrMatrix") , function(x) {
  #splr(t(x@x), x@b, x@a)

  new("splrMatrix", x = t(x@x), a = x@b, b = x@a, Dim = dim(t(x@x)))
})


#' @rdname splr
setMethod("diag", signature = signature(x = "splrMatrix"), function(x) {
  drop(diag(x@x) + rowSums(x@a * x@b))
})

setAs(
  "splrMatrix", "dgeMatrix",
  function(from) from@x + Matrix::tcrossprod(from@a, from@b)
)
