
rperm <- function(n){
  Matrix::Diagonal(n)[sample(n), ]
}
