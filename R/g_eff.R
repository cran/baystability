#' Genotype effects
#'
#' @description
#' \code{bys_ls_gen_effects()} computes the genotype main effects, the deviations of the genotype means from the grand mean.
#'
#' @details
#' \deqn{\hat{\tau}_i = \bar{y}_{i\cdot} - \bar{y}, \qquad \bar{y}_{i\cdot} = \frac{1}{E}\sum_j \bar{y}_{ij},}{tau_i
#'   = y_i. - ybar, y_i. = sum_j y_ij / E,}
#' the least-squares estimates of the AMMI model (equal to the deviations of the plot means when the data are
#' balanced).
#'
#' @inheritParams bys_ls_means
#' @return A \eqn{G \times 1}{G x 1} matrix with column \code{GEffs} (row names: genotypes).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bys_ls_env_effects}}, \code{\link{bys_ls_ge_effects}}
#' @examples
#' data(cultivo2009)
#' bys_ls_gen_effects(
#'     .data = cultivo2009
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
bys_ls_gen_effects <- function(.data, .y, .gen, .env) UseMethod("bys_ls_gen_effects")

#' @rdname bys_ls_gen_effects
#' @export
bys_ls_gen_effects.default <- function(.data, .y, .gen, .env) {
  cd <- .bys_prep("bys_ls_gen_effects", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  matrix(rowMeans(cd$Y) - mean(cd$Y), ncol = 1L, dimnames = list(cd$G, "GEffs"))
}
