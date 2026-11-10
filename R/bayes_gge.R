#' Bayesian GGE (sites regression) model
#'
#' @description
#' \code{bayes_gge()} fits the sites regression (SREG) model behind the GGE biplot by the method of Perez-Elizalde et
#' al. (2012): it is \code{\link{bayes_ammi}(.model = "SREG")}, under the name used by the package
#' \pkg{bayesstability}. The genotype main effect is part of the multiplicative terms; an earlier trial enters as
#' prior through \code{.prior_data}.
#'
#' @inheritParams bayes_ammi
#' @return An object of class \code{c("bys_bayes_ammi", "bys_result", "list")}; see \code{\link{bayes_ammi}()}.
#' @references
#' Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear
#' models applied to plant breeding trials with genotype x environment interaction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, \strong{17}(1), 15--37. \doi{10.1007/s13253-011-0063-9}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{bayes_gge_biplot}}
#' @examples
#' data(cultivo2008)
#' data(cultivo2009)
#' fg <-
#'   bayes_gge(
#'       .data       = cultivo2009
#'     , .y          = y
#'     , .gen        = entry
#'     , .env        = site
#'     , .rep        = rep
#'     , .nIter      = 400
#'     , .prior_data = cultivo2008
#'     , .n0         = NULL
#'     , .k          = 2
#'     , .burnin     = 100
#'     , .thin       = 1
#'     , .chains     = 1
#'     , .priors     = bys_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level      = 0.95
#'     , .seed       = 1
#'     , .align      = TRUE
#'     , .waic       = TRUE
#'     )
#' fg
#' bp <-
#'   bayes_gge_biplot(
#'       fit           = fg
#'     , burnin        = 0
#'     , thin          = 1
#'     , pb            = 0.05
#'     , plot_stable   = TRUE
#'     , plot_unstable = TRUE
#'     , ncolors       = 5
#'     , .which        = "both"
#'     , .grid         = 100
#'     )
#' bp$biplot
#' @export
bayes_gge <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .prior_data = NULL, .n0 = NULL, .k = NULL,
                      .burnin = 0L, .thin = 1L, .chains = 1L, .priors = bys_priors(), .level = 0.95, .seed = NULL,
                      .align = TRUE, .waic = TRUE) {
  UseMethod("bayes_gge")
}

#' @rdname bayes_gge
#' @export
bayes_gge.default <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .prior_data = NULL, .n0 = NULL,
                              .k = NULL, .burnin = 0L, .thin = 1L, .chains = 1L, .priors = bys_priors(),
                              .level = 0.95, .seed = NULL, .align = TRUE, .waic = TRUE) {
  cl <- match.call()
  cl[[1L]] <- quote(bayes_ammi)
  cl$.model <- "SREG"
  out <- eval(cl, parent.frame())
  attr(out, "call") <- .bys_clean_call(match.call())
  out
}
