#' Genotype-by-environment variances
#'
#' @description
#' \code{bys_ls_var()} computes the variance of the replicates of every genotype-by-environment cell.
#'
#' @details
#' \deqn{s^2_{ij} = \frac{1}{n_{ij} - 1}\sum_l (y_{ijl} - \bar{y}_{ij})^2,}{s2_ij = sum_l (y_ijl - y_ij)^2 / (n_ij -
#'   1),}
#' \code{NA} for a cell with one replicate. The within-cell sum of squares of the Bayesian AMMI model is
#' \eqn{\sum_{ij}(n_{ij} - 1)s^2_{ij}}{sum_ij (n_ij - 1) s2_ij}.
#'
#' \strong{Change from the earlier versions.} \code{bayes_ammi()} used the sum of these variances in place of the
#' within-cell sum of squares.
#'
#' @inheritParams bys_ls_means
#' @return A list (class \code{c("bys_ge_var", "bys_result", "list")}) with \code{ge_variances}, the \eqn{G \times
#'   E}{G x E} matrix of variances.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bys_ls_means}}
#' @examples
#' data(cultivo2009)
#' bys_ls_var(
#'     .data = cultivo2009
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
bys_ls_var <- function(.data, .y, .gen, .env) UseMethod("bys_ls_var")

#' @rdname bys_ls_var
#' @export
bys_ls_var.default <- function(.data, .y, .gen, .env) {
  cd <- .bys_prep("bys_ls_var", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  d  <- cd$data
  G  <- length(cd$G)
  cell <- as.integer(d$Gen) + (as.integer(d$Env) - 1L) * G
  ss <- collapse::fsum((d$y - cd$Y[cell])^2, cell)
  V  <- matrix(unname(ss) / (as.vector(cd$N) - 1), G, length(cd$E), dimnames = dimnames(cd$Y))
  V[cd$N < 2L] <- NA_real_
  .bys_result(list(ge_variances = V), "ge_var", label = "Genotype-by-environment variances", call = match.call(),
              meta = list(response = cd$names[["y"]]))
}
