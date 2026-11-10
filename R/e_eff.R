#' Environment effects
#'
#' @description
#' \code{bys_ls_env_effects()} computes the environment main effects, the deviations of the environment means from the grand mean.
#'
#' @details
#' \deqn{\hat{\delta}_j = \bar{y}_{\cdot j} - \bar{y}, \qquad \bar{y}_{\cdot j} = \frac{1}{G}\sum_i
#'   \bar{y}_{ij}.}{delta_j = y_.j - ybar, y_.j = sum_i y_ij / G.}
#'
#' @inheritParams bys_ls_means
#' @return An \eqn{E \times 1}{E x 1} matrix with column \code{EEffs} (row names: environments).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bys_ls_gen_effects}}, \code{\link{bys_ls_ge_effects}}
#' @examples
#' data(cultivo2009)
#' bys_ls_env_effects(
#'     .data = cultivo2009
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
bys_ls_env_effects <- function(.data, .y, .gen, .env) UseMethod("bys_ls_env_effects")

#' @rdname bys_ls_env_effects
#' @export
bys_ls_env_effects.default <- function(.data, .y, .gen, .env) {
  cd <- .bys_prep("bys_ls_env_effects", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  matrix(colMeans(cd$Y) - mean(cd$Y), ncol = 1L, dimnames = list(cd$E, "EEffs"))
}
