#' @title splrMatrix class definition
#'
#' @description stores sparse plus low-rank matrices (e.g. from matrix factorization or centering graphs) of the form
#' \code{x + a \%*\% t(b)} for faster computation
#'
#' @param x a sparse matrix
#' @param a a low-rank factor or a matrix
#' @param b optional. a low-rank factor for \code{a \%*\% t(b)}. if \code{b} is not provided, a will be factorized using 
#' \code{irlba} provided \code{factorize = TRUE}
#' @param rank the estimated rank of the matrix to be factorized.  If \code{rank} is not provided, a guess is made
#' @param tol the tolerance for the eigenvalues if \code{rank} is not provided
#' @param dimnames optional - the list of names for the matrix
#' @param factorize a safeguard to prevent unnecessary computation (e.g. if the matrix is very large).  If \code{a} 
#' is provided and \code{factorize = FALSE}, an error will result.
#' 
#'@import Matrix
#'@import methods
#'@importClassesFrom Matrix sparseMatrix
setClass("splrMatrix",
  slots = c(x="sparseMatrix",a="Matrix",b="Matrix"),
  contains = 'sparseMatrix')


setGeneric(
  name = "splr",
  def = function(x,a=NULL,b=NULL, ...){
  
  # rank = NULL,factorize = FALSE,tol = .00001,dimnames = list(NULL,NULL)

  # x+ab' 
  # x is mxn
  # and is an mxr and b is a rxm matrix

  
  dx=dim(x)
  
  #neither a nor b should be null
  if(is.null(a) & is.null(b)) {
    stop("\nEither specify \n1.) a \n or \n2.) a and b \n otherwise you can just save x as a sparse Matrix object")
  }
  
  if(is.null(b)) {
    da = dim(a)
    if (!factorize) {
      stop("please provide an already factorized low-rank matrix (or include factorize = TRUE)")
    }
    #dimensions of a can't be bigger than dimensions of x
    if (da[1] > dx[1] | da[2] > dx[2]) {
        stop("please provide a such that dim(a) <= dim(x)")
    }
    
    
    if (da[1] == dx[1] & da[2] == dx[2]) { #then we factorize (maybe)
     
      
      #if we know rank, we factorize on that
      if (!is.null(rank) ) {
        
        temp <- irlba(a,rank)
      } else {
        warning(paste0('rank is not provided, using tolerance = ',tol,' to calculate the rank of the factorization'))
        temp <- irlba(a)
        rank <- max(which(abs((temp$d)) > tol))
        if (rank == 5) {
          warning("Rank is estimated to be 5 -- providing a rank estimate may help in the factorization")
        }
      }
      
      #take the truncated svd
      newA <- as.matrix(temp$u[,1:rank]) %*% diag(temp$d)[1:rank,1:rank]
      newB <- as.matrix(temp$v[,1:rank])
      new("splrMatrix",x=x,a=as(newA,"Matrix"),b=as(newB,'Matrix'), Dim = dim(x),Dimnames = dimnames)
      
    } else {
      stop('b is not provided and a is not the same dimension as x')
    }
    
    
  } else { #b is not null
    a <- as(a,'Matrix')
    b <- as(b,'Matrix')
    da=dim(a)
    db=dim(b)
    if(da[1]!=dx[1]) {
      stop("number of rows of x not equal to number of rows of a")
    }
    if(db[1]!=dx[2]) {
      stop("number of columns of x not equal to number of rows of b\nhint: b needs to be nxr when a is mxr and x is mxn")
    }
    if(da[2]!=db[2]) {
      stop("number of columns of a not equal to number of columns of n")
    }
    
    new("splrMatrix",x=x,a=a,b=b,Dim = dim(x),Dimnames = dimnames)
  }
  
 
  
})


setMethod(
  f = "splr",
  signature = signature(x = "Matrix", a = "Matrix", b = "Matrix"),
  definition = function(x, a, b, ...){
    new("splrMatrix", x=x, a=a,b=b, Dim = dim(x),Dimnames = list(NULL,NULL))
  }
)

as.matrix.splr=function(x,...)  {
  
  as.matrix(as.matrix(x@x,...)+x@a%*%t(x@b),...)
  
  
}

setMethod("as.matrix","splrMatrix",as.matrix.splr)


# setAs('splrMatrix','Matrix',function(from) {
#   as(from@x + from@a%*%t(from@b),'Matrix')
# })

as.character.splr <- function(x) {
  as.character(as.matrix(x))
}
setMethod("as.character","splrMatrix",as.character.splr)

setAs("splrMatrix","character",function(from) {
  as(as.matrix(from),'character')
})





setMethod('-',signature = signature(e1 = 'splrMatrix',e2 = 'missing'),function(e1) {
  new("splrMatrix",x=-e1@x,a=-e1@a,b=e1@b,Dim = dim(e1@x))
})

.leftmult=function(x, y){
  #y is splr, x is a matrix
  
  a <- y@a
  sx <- y@x
  b <- y@b
 
  
  if(is.null(a)) {
    x %*% sx
  
  } else if (is(x,"sparseMatrix")) {
    p <- as(x%*%sx,"sparseMatrix")
    anew <- as(x%*%a, "Matrix")
    b <- as(b,"Matrix")
    new("splrMatrix",x=p,a=anew,b=b,Dim = dim(x))
  } else{
    part1=x %*% sx
    part2=x %*% a
    part2=part2 %*% t(b)
    part1+part2
  }
}

setMethod("%*%",signature(x="splrMatrix",y="splrMatrix"), function(x,y) {
  p <- as(x@x %*% y@x + x@x %*% y@a %*% t(y@b) + x@a %*% t(x@b) %*% y@x,'sparseMatrix')
  new('splrMatrix'
      ,x= p
      ,a =  x@a %*% t(x@b) %*% y@a
 
      ,b=y@b, Dim = dim(x))
})

setMethod("%*%",signature(x="Matrix",y="splrMatrix"),.leftmult)

setMethod("%*%",signature(x="matrix",y="splrMatrix"),.leftmult)
setMethod("%*%",signature(x="numeric",y="splrMatrix"),.leftmult)

setMethod("%*%",signature(x="numLike",y="splrMatrix"),.leftmult)


setMethod("%*%",signature(x="ANY",y="splrMatrix"),.leftmult)

.rightmult=function(x,y){
  #x is splr, y is matrix
  a=x@a
  b=x@b
  sx=x@x
  
  if(is.null(a)) {
    sx%*%y
  }
  
  if (is(y,"sparseMatrix")) {
    newx <- sx%*%y
    newx <- as(newx,"sparseMatrix")
    newB <-  t(t(b)%*%y)
    newB <- as(newB,"Matrix")
    new("splrMatrix",x=newx,a=a,b=newB,Dim = dim(newx))
  } else {
    
  
    
    part1 <- sx %*% y
    part2 <- t(b) %*% y
    part2 <- a%*% part2
    part1+part2
  }
  
  
}
setMethod("dim", signature(x= "splrMatrix"),function(x) { dim(x@x)})
setMethod("%*%",signature(x="splrMatrix",y="Matrix"),.rightmult)
setMethod("%*%",signature(x="splrMatrix",y="matrix"),.rightmult)
setMethod("%*%",signature(x="splrMatrix",y="numeric"),.rightmult)

setMethod("%*%",signature(x="splrMatrix",y="numLike"),.rightmult)
setMethod("%*%",signature(x="splrMatrix",y="ANY"),.rightmult)

#doesn't return an splr
setMethod('*',signature = signature(e1 = 'splrMatrix',e2 = 'splrMatrix'),function(e1,e2) {
  x <- as(e1@x * e2@x,'sparseMatrix')
  new("splrMatrix", x = x, a= e1@a %*% t(e1@b) * e2@a ,b =e2@b,Dim = dim(x))
})

#return sparse
.multiply <- function(e1,e2) {
  if (length(e2) == 1) {
    new('splrMatrix',x=e1@x*e2,a=e2*e1@a, b=e1@b,Dim = dim(e1@x*e2))
  } else {
    return(e1@x*e2 + e1@a %*% t(e1@b)*e2) 
  }
}

setMethod("*",signature (e1 = 'Matrix',e2 = 'splrMatrix'), function(e1,e2) {
  
  .multiply(e2,e1)
})

setMethod("*",signature (e1 = 'splrMatrix',e2 = 'ddiMatrix'), function(e1,e2) {
  
  .multiply(e1,e2)
})

setMethod("*",signature (e1 = 'matrix',e2 = 'splrMatrix'), function(e1,e2) {
  
  .multiply(e2,e1)
})
setMethod("*",signature (e1 = 'ANY',e2 = 'splrMatrix'), function(e1,e2) {
  
  .multiply(e2,e1)
})

setMethod("*",signature (e1 = 'splrMatrix',e2 = 'matrix'), function(e1,e2) {
  
  .multiply(e1,e2)
})
setMethod("*",signature (e1 = 'splrMatrix',e2 = 'Matrix'), function(e1,e2) {
  
  .multiply(e1,e2)
})
setMethod("*",signature (e1 = 'splrMatrix',e2 = 'ANY'), function(e1,e2) {
  
  .multiply(e1,e2)
})
setMethod("/",signature (e1 = 'splrMatrix',e2 = 'matrix'), function(e1,e2) {
  
  .multiply(e1,1/e2)
})
setMethod("/",signature (e1 = 'splrMatrix',e2 = 'Matrix'), function(e1,e2) {
  
  .multiply(e1,1/e2)
})
setMethod("/",signature (e1 = 'splrMatrix',e2 = 'ANY'), function(e1,e2) {
  
  .multiply(e1,1/e2)
})


#doesn't create another SPLR...
.addSplr <- function(e1,e2) {
  e1@x + e2@x + (e1@a %*% t(e1@b)) + e2@a %*% t(e2@b)
}

setMethod('+',signature = signature(e1 = 'splrMatrix',e2 = 'splrMatrix'),.addSplr)
setMethod('-',signature = signature(e1 = 'splrMatrix',e2 = 'splrMatrix'),function(e1,e2) {
  .addSplr(e1,-e2)
})



.leftadd=function(e1,e2) {
  #e1 is splr
  
  
  if (is(e2,"sparseMatrix")) {
    new("splrMatrix",x = as(e1@x + e2,"sparseMatrix"),a=e1@a,b=e1@b,Dim = dim(e2))
  } else {
    e1@x + e1@a %*% t(e1@b) + e2
  }
    
  
}
setMethod("+", signature(e1="splrMatrix",e2="Matrix"), function(e1,e2) {
  .leftadd(e1 = e1,e2=e2)
})
setMethod("+", signature(e1="splrMatrix",e2="ANY"), function(e1,e2) {
  .leftadd(e1 = e1,e2=e2)
})
setMethod("-", signature(e1="splrMatrix",e2="Matrix"), function(e1,e2) {
  .leftadd(e1 = e1,e2=-e2)
})
setMethod("-", signature(e1="splrMatrix",e2="ANY"), function(e1,e2) {
  .leftadd(e1 = e1,e2=-e2)
})

.rightadd=function(e1,e2) {
  

    if (is(e1,"sparseMatrix")) {
      splr(as(e2@x + e1,"sparseMatrix"),a=e2@a,b=e2@b)
    } else { #e1 is not sparse
      e2@x + e2@a%*%t(e2@b) + e1
    }
  
  
  
}

setMethod("+", signature("Matrix","splrMatrix"), function(e1,e2) {
  .rightadd(e1,e2)
})
setMethod("+", signature("ANY","splrMatrix"), function(e1,e2) {
  .rightadd(e1,e2)
})
setMethod("-", signature("Matrix","splrMatrix"), function(e1,e2) {
  .rightadd(e1,-e2)
})
setMethod("-", signature("ANY","splrMatrix"), function(e1,e2) {
  .rightadd(e1,-e2)
})




#frobenius norm
Frobsmlr=function(x,a,b){
 
  
    #expansion due to trevor hastie
    xnorm <- norm(x,type = 'f')
    xab=as.matrix(x%*%b)
    xab=sum(xab*a)
    aa=t(a)%*%a
    bb=t(b)%*%b
    ab=sum(aa*bb)
    sqrt(pmax(0,xnorm^2+2*xab+ab))

    
  
  
  
}


setMethod("norm",signature(x="splrMatrix",type="character"),
          function(x,type,...){
            switch(type,
                   "F"=Frobsmlr(x=x@x,a=x@a,b=x@b),
                   "f"=Frobsmlr(x=x@x,a=x@a,b=x@b),
                   norm(as.matrix(x),type = type,...)
            )
          },valueClass="numeric")




#complete
.rsum=function(x,...){
  #x is splrMatrix matrix
  
    rx=rowSums(x@x)
    cb=colSums(x@b)
    drop(rx+x@a%*%cb) 
  
 
}
setMethod("rowSums","splrMatrix",.rsum)



.csum=function(x,...){
  #x is splrMatrix matrix
   
    cx=colSums(x@x)
    ca=colSums(x@a)
    drop( cx+x@b%*%ca)
  
}

setMethod("colSums","splrMatrix",.csum)


#toDo
.rmean=function(x,...){
  #x is splrMatrix matrix
 
    rx=rowMeans(x@x)
    cb=colMeans(x@b)
    drop(rx+x@a%*%cb)
  
}


setMethod("rowMeans","splrMatrix",.rmean)


#toDo
.cmean=function(x,...){
  #x is splrMatrix matrix
  
    cx=colMeans(x@x)
    ca=colMeans(x@a)
    drop(cx+x@b%*%ca)
  
}

setMethod("colMeans","splrMatrix",.cmean)

setMethod("[",signature(x="splrMatrix",i = 'missing',j = 'missing',drop = 'missing') ,function(x) {
          x
})

setMethod("[",signature(x="splrMatrix",i = 'numeric',j = 'numeric',drop = 'logical') 
          ,function(x,i,j, ..., drop) {
            if (drop) {
              return(drop(x@x[i,j,...] + (x@a[i,]%*%t(x@b)[,j]) ))
            } else {
              return(new("splrMatrix",x@x[i,j,...,drop = FALSE],x@a[i,,drop =FALSE],x@b[j,,drop =FALSE]
                         , Dim = dim(x@x[i,j,...,drop = FALSE])) )
            }
    })

setMethod("[",signature(x="splrMatrix",i = 'missing',j = 'numeric',drop = 'logical') 
          ,function(x,j, ..., drop) {
            i <- c(1:dim(x@x)[1])
            if (drop) {
              return(x@x[i,j,...] + (x@a[i,,drop =FALSE]%*%t(x@b)[,j,drop =FALSE]) )
            } else {
              return(new("splrMatrix",x@x[i,j,...,drop = FALSE],x@a[i,,drop =FALSE],x@b[j,,drop =FALSE]
                         , Dim = dim(x@x[i,j,...,drop = FALSE])) )
            }
            
            
            
          })

setMethod("[",signature(x="splrMatrix",i = 'numeric',j = 'missing',drop = 'logical') 
          ,function(x,i, ..., drop) {
            j <- c(1:dim(x@x)[2])
            if (drop) {
              return(drop(x@x[i,j,...] + (x@a[i,,drop =FALSE]%*%t(x@b)[,j,drop =FALSE]) ))
            } else {
              return( new("splrMatrix",x=x@x[i,j,...], a = x@a[i,,drop = FALSE],b = x@b[j,,drop=FALSE],
                          Dim = dim(x=x@x[i,j,...])))
            }
        
 })


setMethod("[",signature(x="splrMatrix",i = 'numeric',j = 'numeric',drop='missing') 
          ,function(x,i,j, ...) {
            
            
            return(new("splrMatrix",x=x@x[i,j,...,drop = FALSE],a=x@a[i,,drop =FALSE],b=x@b[j,,drop =FALSE]
                       ,Dim = dim(x@x[i,j,...,drop=FALSE])))
            
          })

setMethod("[",signature(x="splrMatrix",i = 'numeric',j = 'numeric') 
          ,function(x,i,j, ..., drop = TRUE) {
            
            
            return(new("splrMatrix",x=x@x[i,j,...,drop = FALSE],a=x@a[i,,drop =FALSE],b=x@b[j,,drop =FALSE],
                       Dim = (x=x@x[i,j,...,drop = FALSE])))
  })

setMethod("[",signature(x="splrMatrix",i = 'logical',j = 'logical',drop = 'ANY') 
          ,function(x,i,j, ..., drop) {
            
            if (drop) {
              return(new('splrMatrix',x=x@x[i,j,drop =FALSE],a= x@a[i,],b=x@b[j,,drop =FALSE] ))
            } else {
              new('splrMatrix',x= x@x[i,j,drop = FALSE]
                  , a= as(x@a[i,,drop = FALSE],'Matrix'),b= as(x@b[j,,drop = FALSE],'Matrix')
                  , Dim = (x=x@x[i,j,...,drop = FALSE]))
              
            }
 })


setMethod("[",signature(x="splrMatrix",i = 'logical',j = 'logical',drop = 'missing') 
          ,function(x,i,j, ...) {
            
            
            new('splrMatrix',x= x@x[i,j], a= x@a[i,,drop = FALSE],b= x@b[j,,drop = FALSE]
                ,Dim = dim(x@x[i,j]))
    })



setMethod("[",signature(x="splrMatrix") 
          ,function(x) {
              x
 
})

setMethod("[",signature(x="splrMatrix",i = 'missing',j = 'numeric',drop = 'missing') 
          ,function(x,j, ..., drop = TRUE) {
            
            
            new('splrMatrix',x = x@x[,j],a= x@a,b = x@b[j,,drop = FALSE],
                Dim = dim(x@x[,j]))
           
            
          })
setMethod("[",signature(x="splrMatrix",i = 'numeric',j = 'missing',drop='missing') 
          ,function(x,i , ..., drop = TRUE) {
            j = c(1:dim(x@x)[2])
            
            new('splrMatrix',x = x@x[i,],a= as(x@a[i,,drop = FALSE],'Matrix'),b = x@b
                , Dim= dim(x@x[i,]))
            
          })

#can fix this
setMethod("[",signature(x="splrMatrix",i = 'matrix',j = 'missing',drop='missing') 
          ,function(x,i , ...) {
          x@x[i] + x@a[i[,2],]%*% t(x@b[i[,2],])
           #as.matrix(x@x + x@a%*%t(x@b))[i]
 })

#document the issues with doing this
setMethod("[<-",signature(x="splrMatrix",i = 'numeric',j = 'numeric',value= 'ANY') 
          ,function(x,i,j, ..., value) {
            y <- x@x
            y[i,j] <- value
            a <- x@a
            a[i,] <- 0
            b <- x@b
            b[j,] <- 0
            new("splrMatrix",x=y,a=a,b=b,Dim = dim(y))
            
   })

setMethod("[<-",signature(x="splrMatrix",i = 'numeric',j = 'missing',value= 'ANY') 
          ,function(x,i, ..., value) {
            j <- c(1:dim(x@x)[2])
            y <- x@x
            y[i,j] <- value
            a <- x@a
            a[i,] <- 0
            b <- x@b
            b[j,] <- 0
            new("splrMatrix",x=y,a=a,b=b,Dim = dim(y))
 })

setMethod("[<-",signature(x="splrMatrix",i = 'missing',j = 'numeric',value= 'ANY') 
          ,function(x,j, ..., value) {
            i <- c(1:dim(x@x)[1])
            y <- x@x
            y[,j] <- value
            a <- x@a
            a[i,] <- 0
            b <- x@b
            b[j,] <- 0
            new("splrMatrix",x=y,a=a,b=b,Dim = dim(y))
 })


setMethod("dim", signature(x = "splrMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod('str',signature(object = "splrMatrix"),function(object){
  cat('splrMatrix')
  cat('\nDimension: ', dim(object@x))
  cat('\nLower rank matrix is rank: ',min(dim(object@a)) )
})


setMethod("t",signature = signature(x="splrMatrix") ,function(x) {
  #splr(t(x@x),x@b,x@a)
 
  new("splrMatrix",x=t(x@x),a=x@b,b=x@a,Dim = dim(t(x@x)))
})


setMethod("diag",signature=signature(x = "splrMatrix"),function(x) {
  drop(diag(x@x) + x@a*x@b)
})






