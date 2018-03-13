#' @name    matrix_k
#' @aliases matrix_k
#' @title   k Matrix
#' @description Gives k matrix
#'
#' @param n  Number of columns
#'
#' @return Matrix
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@gmail.com})
#'    }
#'
#' @references
#'  Perez-Elizalde, S., Jarquin, D., and Crossa, J. (2011)
#'  A General Bayesian Estimation Method of Linear–Bilinear Models
#'  Applied to Plant Breeding Trials With Genotype × Environment Interaction.
#'  \emph{Journal of Agricultural, Biological, and Environmental Statistics},
#'   17, 15–37.  (\href{https://link.springer.com/article/10.1007/s13253-011-0063-9}{doi:10.1007/s13253-011-0063-9})
#'
#'
#' @export


matrix_k <- function(n){
  UseMethod("matrix_k")
}


#' @export
#' @rdname matrix_k


matrix_k.default <- function(n){
  m <- matrix(0, nrow = (n-1), ncol = n)
  long <- sqrt((n - (0:(n-2)) - 1) * (n - (0:(n-2))))
  diag(m) <- (n - (0:(n-2)) - 1)
  m[upper.tri(m)] <- -1
  m / long[row(m)]
}
