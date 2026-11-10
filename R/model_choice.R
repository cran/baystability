#' Compare Bayesian AMMI and GGE models
#'
#' @description
#' \code{bys_compare()} compares fits of \code{\link{bayes_ammi}()} that differ in the number of multiplicative
#' terms, the model (AMMI or GGE), the error variances or the priors, by the deviance information criterion (DIC)
#' and the widely applicable information criterion (WAIC). Smaller values predict new data better.
#'
#' @details
#' With \eqn{D(\boldsymbol{\theta}) = -2\log p(\mathbf{y} \mid \boldsymbol{\theta})}{D(theta) = -2 log p(y |
#' theta)} over the observed plots,
#' \deqn{\mathrm{DIC} = D(\bar{\boldsymbol{\theta}}) + 2p_D, \qquad p_D = \overline{D(\boldsymbol{\theta})} -
#'   D(\bar{\boldsymbol{\theta}})}{DIC = D(theta_bar) + 2 pD, pD = mean D(theta) - D(theta_bar)}
#' (Spiegelhalter et al. 2002), the deviance evaluated at the posterior means of the cell means and variances; and
#' \deqn{\mathrm{WAIC} = -2\left(\sum_o \log\frac{1}{S}\sum_s p(y_o \mid \boldsymbol{\theta}^{(s)}) -
#'   \sum_o \mathrm{Var}_s\log p(y_o \mid \boldsymbol{\theta}^{(s)})\right)}{WAIC = -2 (lppd - p_waic)}
#' (Watanabe 2010; Gelman et al. 2014), accumulated over all draws during sampling. All fits must use the same data.
#'
#' The standard error of a difference in WAIC comes from the plot-by-plot differences of the contributions
#' \eqn{w_o}{w_o} to WAIC: \eqn{\mathrm{se} = \sqrt{n\,\mathrm{Var}_o(w_o^{(m)} - w_o^{(\mathrm{best})})}}{se = sqrt(n
#' var(w_m - w_best))} (Vehtari et al. 2017). Because the models are fitted to the same plots, it is smaller
#' than the standard error of either WAIC. A difference of less than about two of its standard errors does not
#' separate two models.
#'
#' @param ... fits of \code{\link{bayes_ammi}()} to the same data, named or not
#'   (names label the rows).
#' @return An object of class \code{c("bys_compare", "bys_result", "list")} with \code{table}
#'   (\code{data.table} sorted by WAIC: \code{model}, \code{terms}, \code{variances}, \code{DIC}, \code{pD},
#'   \code{WAIC}, \code{p_waic}, \code{se} (standard error of WAIC), \code{lambda_last_lower} (lower HPD limit of the
#'   last singular value), \code{lambda_last_mean}, \code{dWAIC} (difference from the best) and \code{dWAIC_se} (its
#'   standard error)).
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bys_interpret}()}.
#' @references
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures of model
#' complexity and fit. \emph{Journal of the Royal Statistical Society B}, \strong{64}(4), 583--639.
#' \doi{10.1111/1467-9868.00353}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion
#' in singular learning theory. \emph{Journal of Machine Learning Research}, \strong{11}, 3571--3594.
#'
#' Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models.
#' \emph{Statistics and Computing}, \strong{24}(6), 997--1016. \doi{10.1007/s11222-013-9416-2}
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out
#' cross-validation and WAIC. \emph{Statistics and Computing}, \strong{27}(5), 1413--1432.
#' \doi{10.1007/s11222-016-9696-4}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @family model choice
#' @examples
#' data(cultivo2009)
#' fits <-
#'   lapply(
#'       X   = 1:2
#'     , FUN = function(K) {
#'         bayes_ammi(
#'             .data       = cultivo2009
#'           , .y          = y
#'           , .gen        = entry
#'           , .env        = site
#'           , .rep        = rep
#'           , .nIter      = 600
#'           , .prior_data = NULL
#'           , .n0         = NULL
#'           , .model      = "AMMI"
#'           , .k          = K
#'           , .burnin     = 100
#'           , .thin       = 1
#'           , .chains     = 1
#'           , .priors     = bys_priors(
#'                 mu_mean     = 0
#'               , mu_var      = 1e+8
#'               , tau_mean    = 0
#'               , tau_var     = 1e+8
#'               , delta_mean  = 0
#'               , delta_var   = 1e+8
#'               , lambda_mean = 0
#'               , lambda_var  = 1e+8
#'               , nu          = 1
#'               , s2          = NULL
#'               )
#'           , .level      = 0.95
#'           , .seed       = 1
#'           , .align      = TRUE
#'           , .waic       = TRUE
#'           )
#'       }
#'     )
#' bys_compare(
#'     AMMI1 = fits[[1]]
#'   , AMMI2 = fits[[2]]
#'   )
#' @export
bys_compare <- function(...) UseMethod("bys_compare")

#' @rdname bys_compare
#' @export
bys_compare.default <- function(...) {
  fn <- "bys_compare"
  fits <- list(...)
  if (!length(fits)) .bys_stop(fn, "give at least one fit.")
  for (f in fits) .bys_check_class(f, "bys_bayes_ammi", "...", fn)
  nm <- names(fits)
  if (is.null(nm)) nm <- rep("", length(fits))
  nobs <- vapply(fits, function(f) sum(attr(f, "meta")$cells$N), numeric(1L))
  if (length(unique(nobs)) > 1L) .bys_stop(fn, "the fits use different data (", paste(nobs, collapse = ", "), " plots).")
  tab <- data.table::rbindlist(lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    s <- f$settings
    lam <- f$summary[block == "lambda"]
    data.table::data.table(
        model = if (nzchar(nm[i])) nm[i] else sprintf("%s%d%s", s$model, s$k, if (isTRUE(s$prior_trial)) "-prior" else "")
      , terms = s$k
      , prior = if (isTRUE(s$prior_trial)) sprintf("earlier trial (n0 = %g)", s$n0) else "vague"
      , DIC = f$fit$DIC, pD = f$fit$pD
      , WAIC = if (is.null(f$fit$WAIC)) NA_real_ else f$fit$WAIC
      , p_waic = if (is.null(f$fit$p_waic)) NA_real_ else f$fit$p_waic
      , se = if (is.null(f$fit$waic_se)) NA_real_ else f$fit$waic_se
      , lambda_last_lower = lam$lower[nrow(lam)], lambda_last_mean = lam$mean[nrow(lam)])
  }))
  key <- if (all(is.na(tab$WAIC))) tab$DIC else tab$WAIC
  tab[, dWAIC := key - min(key)]
  # standard error of the difference from the best model, from the plot-by-plot contributions to WAIC
  wi <- lapply(fits, function(f) f$fit$waic_i)
  best <- which.min(key)
  tab[, dWAIC_se := if (any(vapply(wi, is.null, NA))) NA_real_ else
    vapply(wi, function(w) sqrt(length(w) * stats::var(w - wi[[best]])), numeric(1L))]
  tab <- tab[order(key)]
  .bys_result(list(table = tab), "compare", label = "Comparison of Bayesian models (smaller is better)",
              call = match.call())
}

#-------------------------------------------------------------------------------------------------------------------
# draws for the packages posterior and bayesplot, and coda
#-------------------------------------------------------------------------------------------------------------------

#' The draws of a fit for posterior, bayesplot and coda
#'
#' @description
#' The draws of \code{\link{bayes_ammi}()} are a \code{coda::mcmc.list} (component \code{draws}). These methods hand
#' them to the packages \pkg{posterior} and \pkg{coda}; every plot of \pkg{bayesplot} (\code{mcmc_trace()},
#' \code{mcmc_areas()}, \code{mcmc_rank_overlay()}, ...) then works on them, as does \code{posterior::summarise_draws()}
#' (rank-normalised \eqn{\hat{R}}{R-hat}, bulk and tail effective sample sizes).
#'
#' @param x a result of \code{\link{bayes_ammi}()}.
#' @param ... unused.
#' @return \code{as_draws_df()}: a \code{posterior::draws_df}; \code{as_draws()}: a \code{posterior::draws_list};
#'   \code{as.mcmc.list()}: the \code{coda::mcmc.list}.
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in Bayesian workflow.
#' \emph{Journal of the Royal Statistical Society A}, \strong{182}(2), 389--402. \doi{10.1111/rssa.12378}
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021). Rank-normalization, folding, and
#' localization: an improved R-hat for assessing convergence of MCMC. \emph{Bayesian Analysis}, \strong{16}(2),
#' 667--718. \doi{10.1214/20-BA1221}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @examples
#' data(cultivo2009)
#' fit <-
#'   bayes_ammi(
#'       .data       = cultivo2009
#'     , .y          = y
#'     , .gen        = entry
#'     , .env        = site
#'     , .rep        = rep
#'     , .nIter      = 500
#'     , .prior_data = NULL
#'     , .n0         = NULL
#'     , .model      = "AMMI"
#'     , .k          = 2
#'     , .burnin     = 100
#'     , .thin       = 1
#'     , .chains     = 2
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
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   dr <- posterior::as_draws_df(x = fit)
#'   posterior::summarise_draws(
#'       .x = posterior::subset_draws(
#'           x        = dr
#'         , variable = c("mu", "lambda[1]", "lambda[2]")
#'         )
#'     )
#' }
#' if (requireNamespace("bayesplot", quietly = TRUE)) {
#'   bayesplot::mcmc_trace(x = fit$draws, pars = c("lambda[1]", "lambda[2]"))
#' }
#' @name bys-draws
NULL

#' @rdname bys-draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bys_bayes_ammi <- function(x, ...) posterior::as_draws_df(x$draws)

#' @rdname bys-draws
#' @exportS3Method posterior::as_draws
as_draws.bys_bayes_ammi <- function(x, ...) posterior::as_draws(x$draws)

#' @rdname bys-draws
#' @exportS3Method coda::as.mcmc.list
as.mcmc.list.bys_bayes_ammi <- function(x, ...) x$draws

.bys_print_compare <- function(x, digits = 1L, ...) {
  cat(attr(x, "label"), "\n\n")
  t <- data.table::copy(x$table)
  num <- c("DIC", "pD", "WAIC", "p_waic", "se", "dWAIC", "dWAIC_se", "lambda_last_lower", "lambda_last_mean")
  t[, (num) := lapply(.SD, round, digits = digits), .SDcols = num]
  print(t)
  invisible(x)
}

.bys_plot_compare <- function(x, ...) {
  d <- data.table::copy(x$table)
  d[, model := factor(model, levels = rev(model))]
  d[, dDIC := DIC - min(DIC)]
  d[, half := 2 * ifelse(is.na(dWAIC_se), se, dWAIC_se)]
  ggplot2::ggplot(d, ggplot2::aes(x = dWAIC, y = model)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = dWAIC - half, xmax = dWAIC + half), height = 0.25,
                            colour = .bys_palette[1L]) +
    ggplot2::geom_point(size = 2.2, colour = .bys_palette[1L]) +
    ggplot2::geom_point(ggplot2::aes(x = dDIC), shape = 4, size = 2.2, colour = .bys_palette[2L]) +
    .bys_x_cont() +
    ggplot2::labs(x = "Difference from the best model: WAIC (dot, with 2 standard errors) and DIC (cross)", y = NULL,
                  title = "Model comparison (smaller is better)") +
    bys_theme()
}

utils::globalVariables(c("dWAIC", "dWAIC_se", "dDIC", "half", "model", "WAIC", "se", "DIC", "statistic", "observed", "p_value", "lab"))
