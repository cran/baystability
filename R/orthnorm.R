#' @name    orthnorm
#' @aliases orthnorm
#' @title   Orthogonal Normalization
#' @description Perform Orthogonal Normalization of a matrix
#'
#' @param u  Matrix
#' @param basis Logical argument by default TRUE
#' @param norm  Logical argument by default TRUE
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
#' @import MASS
#' @import rstiefel
#'
#' @export


orthnorm <- function(u = NULL, basis = TRUE, norm = TRUE){
  UseMethod("orthnorm")
}


#' @export
#' @rdname orthnorm


orthnorm.default <- function(u = NULL, basis = TRUE, norm = TRUE){
  if (is.null(u)) return(NULL)
  if (!(is.matrix(u))) u <- as.matrix(u)
  p <- nrow(u)
  n <- ncol(u)
  if (prod(abs(La.svd(u)$d) > 1e-100) == 0)
    if (p < n){
      warning("too much vectors to orthogonalize.")
      u <- as.matrix(u[, 1:p])
      n <- p
    }
  if (basis & (p > n)){
    base <- diag(p)
    coef.proj <- crossprod(u, base)/diag(crossprod(u))
    base2 <- base - u %*% matrix(coef.proj, nrow = n, ncol = p)
    norm.base2 <- diag(crossprod(base2))
    base <- as.matrix(base[, order(norm.base2) > n])
    u <- cbind(u, base)
    n <- p
  }
  v <- u
  if (n > 1){
    for (i in 2:n){
      coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 1:(i - 1)]))
      v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = p) %*% matrix(coef.proj, nrow = i - 1)
    }
  }
  if (norm){
    coef.proj <- 1/sqrt(diag(crossprod(v)))
    v <- t(t(v) * coef.proj)
  }
  return(v)
}
