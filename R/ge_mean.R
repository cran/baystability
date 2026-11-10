#' Genotype-by-environment means
#'
#' @description
#' \code{bys_ls_means()} computes the table of genotype-by-environment cell means and the grand mean, the data of the
#' AMMI model.
#'
#' @details
#' With \eqn{y_{ijl}}{y_ijl} the response of genotype \eqn{i}{i} in replicate \eqn{l}{l} of environment \eqn{j}{j}
#' and \eqn{n_{ij}}{n_ij} the number of replicates of the cell,
#' \deqn{\bar{y}_{ij} = \frac{1}{n_{ij}}\sum_l y_{ijl}, \qquad \bar{y} = \frac{1}{GE}\sum_{ij}\bar{y}_{ij}.}{y_ij =
#'   sum_l y_ijl / n_ij; ybar = sum_ij y_ij / (GE).}
#' Genotypes and environments are in the order of their factor levels (sorted values when the columns are not
#' factors).
#'
#' @param .data a data frame with one row per plot.
#' @param .y,.gen,.env the columns of the response, genotypes and environments: bare names, strings or vectors.
#' @return A list (class \code{c("bys_ge_mean", "bys_result", "list")}) with \code{ge_means}, the \eqn{G \times
#'   E}{G x E} matrix of means (row and column names: genotypes and environments), and \code{grand_mean}.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bys_ls_var}}, \code{\link{bys_ls_ge_effects}}, \code{\link{bys_ls_ammi}}
#' @examples
#' data(cultivo2009)
#' bys_ls_means(
#'     .data = cultivo2009
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )
#' @export
bys_ls_means <- function(.data, .y, .gen, .env) UseMethod("bys_ls_means")

#' @rdname bys_ls_means
#' @export
bys_ls_means.default <- function(.data, .y, .gen, .env) {
  cd <- .bys_prep("bys_ls_means", .data, substitute(.y), substitute(.gen), substitute(.env), NULL, parent.frame(),
                  need_rep = FALSE)
  .bys_result(list(ge_means = cd$Y, grand_mean = mean(cd$Y)), "ge_mean", label = "Genotype-by-environment means",
              call = match.call(), meta = list(response = cd$names[["y"]]))
}
