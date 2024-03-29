% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splr_matrix_class.R
\docType{class}
\name{splrMatrix-class}
\alias{splrMatrix-class}
\alias{splr}
\alias{splr,Matrix,Matrix,Matrix-method}
\title{Sparse Plus Low-Rank Matrices}
\usage{
splr(x, a = NULL, b = NULL, rank = NULL, dimnames = list(NULL, NULL), ...)

\S4method{splr}{Matrix,Matrix,Matrix}(x, a = NULL, b = NULL, rank = NULL, dimnames = list(NULL, NULL), ...)
}
\arguments{
\item{x}{as in "Matrix"}

\item{a}{as in "Matrix"}

\item{b}{as in "Matrix"}

\item{rank}{rank of the matrix to be factorized.}

\item{dimnames}{optional - the list of names for the matrix}

\item{...}{as in "Matrix"}
}
\value{
splrMatrix object

splrMatrix object
}
\description{
An "S4" class for efficient computation with sparse plus
low-rank matrices. Stores sparse plus low-rank matrices
(e.g. from matrix factorization or centering graphs)
of the form \code{x + a \%*\% t(b)} for faster
computation.
}
\section{Slots}{

\describe{
\item{\code{x}}{a sparse matrix}

\item{\code{a}}{a low-rank factor or a matrix}

\item{\code{b}}{optional. a low-rank factor for \code{a \%*\% t(b)}. if \code{b} is not provided, a will be factorized using
\code{\link[irlba]{irlba}} provided \code{factorize = TRUE}}
}}

\seealso{
Methods are documented in \link{splrMatrix_method}.
 Other relevant methods are \link{splr_sparse_plus_constant} and
}
