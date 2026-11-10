#' Genotype-by-environment interaction effects
#'
#' @description
#' \code{bys_ls_ge_effects()} computes the interaction effects (the cell means centred on genotype and environment means) and
#' their singular value decomposition, the least-squares estimate of the multiplicative part of the AMMI model.
#'
#' @details
#' \deqn{\mathbf{Z} = (z_{ij}), \quad z_{ij} = \bar{y}_{ij} - \bar{y}_{i\cdot} - \bar{y}_{\cdot j} + \bar{y}, \qquad
#'   \mathbf{Z} = \sum_k \lambda_k\boldsymbol{\alpha}_k\boldsymbol{\gamma}_k^\top,}{z_ij = y_ij - y_i. - y_.j + ybar;
#'   Z = sum_k lambda_k alpha_k gamma_k',}
#' with singular values \eqn{\lambda_1 \ge \lambda_2 \ge \cdots}{lambda_1 >= lambda_2 >= ...} and orthonormal
#' singular vectors orthogonal to \eqn{\mathbf{1}}{1}.
#'
#' @inheritParams bys_ls_means
#' @return A list (class \code{c("bys_ge_eff", "bys_result", "list")}) with \code{ge_effects} (the \eqn{G \times
#'   E}{G x E} matrix \eqn{\mathbf{Z}}{Z}) and \code{ge_svd} (\code{svd()} of it: \code{d}, \code{u}, \code{v}).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Gauch, H. G. (1992). \emph{Statistical Analysis of Regional Yield Trials: AMMI Analysis of Factorial Designs}.
#' Elsevier, Amsterdam.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bys_ls_ammi}}
#' @examples
#' data(cultivo2009)
#' fit <-
#'   bys_ls_ge_effects(
#'       .data = cultivo2009
#'     , .y    = y
#'     , .gen  = entry
#'     , .env  = site
#'     )
#' fit$ge_svd$d
#' @export
bys_ls_ge_effects <- function(.data, .y, .gen, .env) UseMethod("bys_ls_ge_effects")

#' @rdname bys_ls_ge_effects
#' @export
bys_ls_ge_effects.default <- function(.data, .y, .gen, .env) {
  cd <- .bys_prep("bys_ls_ge_effects", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  Y  <- cd$Y
  Z  <- Y - outer(rowMeans(Y), colMeans(Y), `+`) + mean(Y)
  .bys_result(list(ge_effects = Z, ge_svd = svd(Z)), "ge_eff", label = "Interaction effects", call = match.call(),
              meta = list(response = cd$names[["y"]]))
}
