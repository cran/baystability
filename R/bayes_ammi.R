#' Bayesian estimation of linear-bilinear models: AMMI and SREG (GGE)
#'
#' @description
#' \code{bayes_ammi()} estimates the linear-bilinear models of a multi-environment trial by the general Bayesian method
#' of Perez-Elizalde et al. (2012): the additive main effects and multiplicative interaction (AMMI) model and the
#' sites regression (SREG, GGE) model. The data of an earlier trial with the same genotypes and environments can be
#' the prior (\code{.prior_data}): its cell means, effects and singular value decomposition give conjugate priors whose
#' weight is its number of replicates, as in the paper's analysis of two years of a maize trial. Without an earlier
#' trial the priors are vague but proper. The sampler is compiled (C++); every parameter has a posterior distribution,
#' and \code{\link{bayes_ammi_biplot}()} draws the biplot with highest posterior density (HPD) regions.
#'
#' @details
#' \strong{Models.} For genotype \eqn{i}{i}, replicate \eqn{l}{l} and environment \eqn{j}{j},
#' \deqn{\textrm{AMMI: } y_{ijl} = \mu + \tau_i + \delta_j + \sum_{k=1}^{K} \lambda_k \alpha_{ik} \gamma_{jk} +
#'   e_{ijl}, \qquad \textrm{SREG: } y_{ijl} = \mu + \delta_j + \sum_{k=1}^{K} \lambda_k \alpha_{ik} \gamma_{jk} +
#'   e_{ijl},}{AMMI: y_ijl = mu + tau_i + delta_j + sum_k lambda_k alpha_ik gamma_jk + e_ijl; SREG: y_ijl = mu +
#'   delta_j + sum_k lambda_k alpha_ik gamma_jk + e_ijl,}
#' \eqn{e_{ijl} \sim N(0, \sigma^2)}{e_ijl ~ N(0, sigma^2)}, with \eqn{\sum_i\tau_i = \sum_j\delta_j = 0}{sum tau =
#' sum delta = 0}, \eqn{\lambda_1 \ge \cdots \ge \lambda_K > 0}{lambda_1 >= ... >= lambda_K > 0}, and orthonormal
#' score matrices \eqn{\mathbf{U} = (\boldsymbol{\alpha}_k)}{U = (alpha_k)} and \eqn{\mathbf{V} =
#' (\boldsymbol{\gamma}_k)}{V = (gamma_k)}, \eqn{\mathbf{U}^\top\mathbf{1} = \mathbf{0}}{U'1 = 0}, and for AMMI also
#' \eqn{\mathbf{V}^\top\mathbf{1} = \mathbf{0}}{V'1 = 0} (paper, eq. 1.2-1.3). By default all
#' \eqn{\min(G, E) - 1}{min(G, E) - 1} terms are fitted (\eqn{\min(G - 1, E)}{min(G - 1, E)} for SREG), as in the
#' paper and the earlier versions; \code{.k} fits fewer.
#'
#' \strong{An earlier trial as prior} (\code{.prior_data}). With \eqn{n_0}{n0} its number of replicates
#' (\code{.n0}), \eqn{\mathbf{Y}_0}{Y0} its table of cell means and \eqn{\mu_0}{mu0},
#' \eqn{\boldsymbol{\tau}_0}{tau0}, \eqn{\boldsymbol{\delta}_0}{delta0}, \eqn{\mathbf{U}_0\mathbf{D}_0
#' \mathbf{V}_0^\top}{U0 D0 V0'} its least-squares fit (paper, Section 2.2):
#' \deqn{\mu \sim N\left(\mu_0, \frac{\sigma^2}{n_0 GE}\right), \quad \mathbf{K}_G\boldsymbol{\tau} \sim
#'   N\left(\mathbf{K}_G\boldsymbol{\tau}_0, \frac{\sigma^2}{n_0 E}\mathbf{I}\right), \quad
#'   \mathbf{K}_E\boldsymbol{\delta} \sim N\left(\mathbf{K}_E\boldsymbol{\delta}_0, \frac{\sigma^2}{n_0
#'   G}\mathbf{I}\right),}{mu ~ N(mu0, sigma^2/(n0 G E)), K_G tau ~ N(K_G tau0, sigma^2/(n0 E) I), K_E delta ~ N(K_E
#'   delta0, sigma^2/(n0 G) I),}
#' \eqn{\lambda_k \sim N(\lambda_{0k}, \sigma^2/n_0)}{lambda_k ~ N(lambda0_k, sigma^2/n0)} truncated to the order
#' (eq. 2.9), von Mises-Fisher priors \eqn{p(\mathbf{U}) \propto \mathrm{etr}(n_0\mathbf{Y}_0\mathbf{V}_0
#' \mathbf{D}_0\mathbf{U}^\top/\sigma^2)}{p(U) ~ etr(n0 Y0 V0 D0 U' / sigma^2)} and \eqn{p(\mathbf{V}) \propto
#' \mathrm{etr}(n_0\mathbf{Y}_0^\top\mathbf{U}_0\mathbf{D}_0\mathbf{V}^\top/\sigma^2)}{p(V) ~ etr(n0 Y0' U0 D0 V' /
#' sigma^2)} (eq. 2.7-2.8), and \eqn{a S^2_0/\sigma^2 \sim \chi^2_a}{a S0^2 / sigma^2 ~ chi^2_a} with \eqn{a = n_0 -
#' 1}{a = n0 - 1} and \eqn{S^2_0}{S0^2} the within-cell mean square of the earlier trial (eq. 2.11; the values of the
#' earlier versions of the package).
#'
#' \strong{Conditional distributions} (paper, Section 2.3). With \eqn{\mathbf{Y}}{Y} the cell means of the trial,
#' \eqn{B}{B} its replicates and \eqn{n = BGE}{n = BGE}, the linear effects are normal with means
#' \eqn{(n\bar{y} + n_0GE\mu_0)/(n + n_0GE)}{(n ybar + n0 G E mu0)/(n + n0 G E)} and so on; the singular values are
#' truncated normal with mean \eqn{(B\,\boldsymbol{\alpha}_k^\top\mathbf{Y}\boldsymbol{\gamma}_k +
#' n_0\lambda_{0k})/(B + n_0)}{(B alpha_k' Y gamma_k + n0 lambda0_k)/(B + n0)} and variance \eqn{\sigma^2/(B +
#' n_0)}{sigma^2/(B + n0)} (eq. 2.15); the score matrices are von Mises-Fisher (eq. 2.13-2.14), drawn one column at a
#' time on the sphere orthogonal to \eqn{\mathbf{1}}{1} and the other columns (Hoff 2009) with parameters
#' \eqn{(B\lambda_k\mathbf{Y}\boldsymbol{\gamma}_k + n_0\lambda_{0k}\mathbf{Y}_0\boldsymbol{\gamma}_{0k})/\sigma^2}{(B
#' lambda_k Y gamma_k + n0 lambda0_k Y0 gamma0_k) / sigma^2} and the same for \eqn{\mathbf{V}}{V}; and
#' \deqn{\sigma^2 \sim \frac{aS^2_0 + \mathrm{SSW} + B\sum_{ij}(\bar{y}_{ij} - \hat{y}_{ij})^2}{\chi^2_{a +
#'   n}}.}{sigma^2 ~ (a S0^2 + SSW + B sum (ybar_ij - fit_ij)^2) / chi^2(a + n).}
#' As in the paper, the dependence of the priors on \eqn{\sigma^2}{sigma^2} is left out of its conditional
#' distribution (it would involve the normalising constants of the von Mises-Fisher priors). Without an earlier trial
#' the priors of \code{\link{bys_priors}()} apply, which do not depend on \eqn{\sigma^2}{sigma^2}, and all
#' conditionals are exact (Crossa et al. 2011).
#'
#' \strong{Starting values and signs.} With an earlier trial the first chain starts at its least-squares fit and every
#' draw keeps the orientation of its score vectors (\code{.align}); otherwise the first chain starts at the
#' least-squares fit of the trial. Other chains start at dispersed values. With all terms fitted and a prior trial that
#' disagrees with the data, the posterior of the scores can have several modes; compare chains (\code{.chains}) and
#' the potential scale reduction factors.
#'
#' \strong{Unbalanced trials.} Missing plots and cells are completed by data augmentation (see
#' \code{\link{bys_priors}}); an earlier trial with empty cells is completed by the EM algorithm before its fit.
#'
#' \strong{Changes from version 0.2.0.} The function no longer reads the earlier trial from objects in the workspace
#' (\code{n0}, \code{mu0}, \code{tau0}, ...; it stopped with "object 'n0' not found" without them): the earlier trial is
#' the argument \code{.prior_data}. The singular values were computed once and never sampled; \code{tau1} held a
#' precision drawn from a residual matrix that was never updated; the scores were drawn without their constraints
#' and orthogonalised afterwards; the columns of \code{delta1} were named for the genotypes. \code{tau1} now holds
#' the draws of the error variance \eqn{\sigma^2}{sigma^2}.
#'
#' @param .data a data frame with one row per plot.
#' @param .y,.gen,.env,.rep the columns of the response, genotypes, environments and replicates: bare names or
#'   strings. The same names are used in \code{.prior_data}.
#' @param .nIter number of iterations kept per chain before thinning (after the burn-in).
#' @param .prior_data \code{NULL} (default: vague priors) or a data frame of an earlier trial with the same columns,
#'   genotypes and environments, used as prior (Perez-Elizalde et al. 2012).
#' @param .n0 weight of the earlier trial: the number of its observations per cell in the prior (default its number of
#'   replicates).
#' @param .model \code{"AMMI"} (default) or \code{"SREG"} (sites regression, GGE; \code{"GGE"} is accepted).
#' @param .k number of bilinear terms \eqn{K}{K} (default \code{NULL}: all of them).
#' @param .burnin number of iterations discarded at the start of every chain (default 0).
#' @param .thin keep every \code{.thin}-th iteration (default 1).
#' @param .chains number of chains (default 1).
#' @param .priors the vague priors used without an earlier trial, from \code{\link{bys_priors}()}.
#' @param .level probability of the HPD intervals (default 0.95).
#' @param .seed integer seed for a reproducible run, or \code{NULL}; the random number state of the session is
#'   restored afterwards.
#' @param .align keep the sign of every bilinear term (default \code{TRUE}).
#' @param .waic compute WAIC during sampling (default \code{TRUE}).
#' @return An object of class \code{c("bys_bayes_ammi", "bys_result", "list")}; the components of version 0.2.0
#'   come first, the draws of all chains one after the other:
#' \describe{
#'   \item{\code{mu1}}{\code{data.table}, column \code{mu}: draws of the overall mean.}
#'   \item{\code{tau1}}{\code{data.table}, column \code{tau}: draws of the error variance \eqn{\sigma^2}{sigma^2}.}
#'   \item{\code{tao1}}{\code{data.table}, columns \code{tao1}, ..., \code{taoG}: draws of the genotype effects (no
#'     columns for SREG).}
#'   \item{\code{delta1}}{\code{data.table}, columns \code{delta1}, ..., \code{deltaE}: draws of the environment
#'     effects.}
#'   \item{\code{lambdas1}}{\code{data.table}, columns \code{lambdas1}, ...: draws of the singular values.}
#'   \item{\code{alphas1}, \code{gammas1}}{\code{data.table}s, one row per draw, columns \code{alphas1}, ...,
#'     \code{alphas<GK>} (genotype scores, genotypes within terms) and \code{gammas1}, ... (environment scores).}
#'   \item{\code{Mean_Alphas}, \code{Mean_Gammas}}{\code{data.table}s, columns \code{V1}, ..., \code{VK}: posterior
#'     means of the score matrices.}
#'   \item{\code{summary}}{\code{data.table}: posterior mean, standard deviation, Monte Carlo standard error, HPD
#'     interval, effective sample size and potential scale reduction factor of every parameter.}
#'   \item{\code{scores}}{\code{data.table}: posterior means of the biplot coordinates.}
#'   \item{\code{draws}}{\code{coda::mcmc.list} of all parameters.}
#'   \item{\code{diagnostics}}{Gelman-Rubin, Raftery-Lewis and Heidelberger-Welch diagnostics.}
#'   \item{\code{ls}}{least-squares fit of the trial; \code{prior}: the fit of the earlier trial used as prior.}
#'   \item{\code{fit}}{DIC and WAIC.}
#'   \item{\code{settings}}{iterations, burn-in, thinning, chains, terms, model, priors, level, seed, time used.}
#' }
#'   Methods: \code{print()}, \code{summary()}, \code{plot()}, \code{\link{bayes_ammi_biplot}()},
#'   \code{\link{bys_interpret}()}, \code{\link{bys_compare}()}, \code{\link{bys_cite}()}.
#' @references
#' Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear
#' models applied to plant breeding trials with genotype x environment interaction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, \strong{17}(1), 15--37. \doi{10.1007/s13253-011-0063-9}
#'
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with applications to
#' multivariate and relational data. \emph{Journal of Computational and Graphical Statistics}, \strong{18}(2),
#' 438--456. \doi{10.1198/jcgs.2009.07177}
#'
#' Wood, A. T. A. (1994). Simulation of the von Mises Fisher distribution. \emph{Communications in Statistics -
#' Simulation and Computation}, \strong{23}(1), 157--164.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com}), Diego Jarquin (\email{diego.jarquin@@gmail.com}), Sergio
#'   Perez-Elizalde (\email{sergiop@@colpos.mx}), Juan Burgueno, Jose Crossa (\email{j.crossa@@cgiar.org})
#' @seealso \code{\link{bayes_ammi_biplot}}, \code{\link{bys_priors}}, \code{\link{bys_interpret}}, \code{\link{bys_ls_ammi}}
#' @examples
#' data(cultivo2008)
#' data(cultivo2009)
#' fm3 <-
#'   bayes_ammi(
#'       .data       = cultivo2009
#'     , .y          = y
#'     , .gen        = entry
#'     , .env        = site
#'     , .rep        = rep
#'     , .nIter      = 300
#'     , .prior_data = cultivo2008
#'     , .n0         = NULL
#'     , .model      = "AMMI"
#'     , .k          = NULL
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
#' fm3
#' fm3$Mean_Alphas
#' plot(x = fm3, which = "biplot")
#' @export
bayes_ammi <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .prior_data = NULL, .n0 = NULL,
                       .model = c("AMMI", "SREG", "GGE"), .k = NULL, .burnin = 0L, .thin = 1L, .chains = 1L,
                       .priors = bys_priors(), .level = 0.95, .seed = NULL, .align = TRUE, .waic = TRUE) {
  UseMethod("bayes_ammi")
}

#' @rdname bayes_ammi
#' @export
bayes_ammi.default <- function(.data, .y, .gen, .env, .rep, .nIter = 1000L, .prior_data = NULL, .n0 = NULL,
                               .model = c("AMMI", "SREG", "GGE"), .k = NULL, .burnin = 0L, .thin = 1L,
                               .chains = 1L, .priors = bys_priors(), .level = 0.95, .seed = NULL, .align = TRUE,
                               .waic = TRUE) {
  fn  <- "bayes_ammi"
  t0  <- proc.time()[["elapsed"]]
  ex  <- list(y = substitute(.y), gen = substitute(.gen), env = substitute(.env), rep = substitute(.rep))
  cd  <- .bys_prep(fn, .data, ex$y, ex$gen, ex$env, ex$rep, parent.frame(), allow_missing = TRUE)
  G <- length(cd$G)
  E <- length(cd$E)
  .model  <- .bys_check_choice(.model[1L], c("AMMI", "SREG", "GGE"), ".model", fn)
  if (.model == "GGE") .model <- "SREG"
  gge     <- .model == "SREG"
  kmax    <- if (gge) min(G - 1L, E) else min(G, E) - 1L
  .k      <- if (is.null(.k)) kmax else .bys_check_number(.k, ".k", fn, min = 1, max = kmax, integer = TRUE)
  .nIter  <- .bys_check_number(.nIter, ".nIter", fn, min = 2, integer = TRUE)
  .burnin <- .bys_check_number(.burnin, ".burnin", fn, min = 0, integer = TRUE)
  .thin   <- .bys_check_number(.thin, ".thin", fn, min = 1, max = .nIter, integer = TRUE)
  .chains <- .bys_check_number(.chains, ".chains", fn, min = 1, integer = TRUE)
  .level  <- .bys_check_number(.level, ".level", fn, min = 0.5, max = 0.999)
  .align  <- .bys_check_flag(.align, ".align", fn)
  .waic   <- .bys_check_flag(.waic, ".waic", fn)
  .bys_check_class(.priors, "bys_priors", ".priors", fn)
  if (!is.null(.seed)) .seed <- .bys_check_number(.seed, ".seed", fn, integer = TRUE)
  K  <- .k
  ls <- .bys_ls(cd$Y, K, gge)
  ls$gge <- gge
  df_res <- cd$n - 1 - (E - 1) - (if (gge) 0 else G - 1) - K * (G + E - 1 - 2 * K + gge) - cd$missing
  resid  <- cd$ssw + sum((cd$N * (cd$Y - ls$fitted)^2)[cd$N > 0])
  # with all terms the cell means are fitted exactly: the within-cell mean square scales the prior
  s2_ls  <- if (df_res >= 1 && resid > 0) resid / df_res else cd$ssw / max(1, cd$n - G * E)
  tp <- if (is.null(.prior_data)) NULL else
    .bys_trial_prior(fn, .prior_data, ex, parent.frame(), cd, K, gge, .n0)
  pr <- .priors
  if (is.null(pr$s2)) pr$s2 <- s2_ls
  if (!is.null(tp)) {
    pr$nu <- tp$a
    pr$s2 <- tp$s2
  }
  if (!is.null(.seed)) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
      rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(.seed)
  }
  n_total <- .burnin + .nIter
  keep    <- seq.int(.burnin + .thin, n_total, by = .thin)
  flag    <- integer(n_total)
  flag[keep] <- 1L
  st0 <- if (is.null(tp)) ls else utils::modifyList(ls, tp$fit[c("mu", "tau", "delta", "lambda", "alpha", "gamma")])
  ref <- if (is.null(tp)) ls$alpha else tp$fit$alpha
  chains  <- lapply(seq_len(.chains), function(ch) {
    .bys_run_chain(cd, K, gge, FALSE, pr, .bys_start(st0, s2_ls, ch, FALSE, E), ref, n_total, flag, .align,
                   .waic, tp)
  })
  lay <- .bys_layout(cd$G, cd$E, K, gge, FALSE)
  out <- .bys_assemble(chains, cd, ls, K, keep, .level, lay)
  # the score draws in the form of version 0.2.0: one row per draw, genotypes (environments) within terms
  wide <- function(Arr, nm) {
    m <- matrix(Arr, dim(Arr)[1L])
    data.table::setnames(data.table::as.data.table(m), paste0(nm, seq_len(ncol(m))))
  }
  all <- .bys_all(out)
  A <- .bys_block(all, lay, "alpha")
  Gm <- .bys_block(all, lay, "gamma")
  out$alphas1 <- wide(A, "alphas")
  out$gammas1 <- wide(Gm, "gammas")
  vk <- function(Arr) data.table::setnames(data.table::as.data.table(apply(Arr, c(2L, 3L), mean)),
                                           paste0("V", seq_len(K)))
  out <- c(out[c("mu1", "tau1", "tao1", "delta1", "lambdas1", "alphas1", "gammas1")],
           list(Mean_Alphas = vk(A), Mean_Gammas = vk(Gm)),
           out[setdiff(names(out), c("mu1", "tau1", "tao1", "delta1", "lambdas1", "alphas1", "gammas1"))])
  out$ls <- list(mu = ls$mu, tau = stats::setNames(ls$tau, cd$G), delta = stats::setNames(ls$delta, cd$E),
                 lambda = ls$lambda, alpha = ls$alpha, gamma = ls$gamma, sigma2 = s2_ls)
  dimnames(out$ls$alpha) <- list(cd$G, paste0("axis", seq_len(K)))
  dimnames(out$ls$gamma) <- list(cd$E, paste0("axis", seq_len(K)))
  out$prior <- if (is.null(tp)) NULL else tp[c("n0", "a", "s2", "fit")]
  out$fit <- .bys_fit_stats(chains, out, cd, lay)
  out$settings <- list(nIter = .nIter, burnin = .burnin, thin = .thin, chains = .chains, kept = length(keep),
                       k = K, model = .model, het = FALSE, priors = pr, prior_trial = !is.null(tp),
                       n0 = if (is.null(tp)) NA_real_ else tp$n0, level = .level, seed = .seed, align = .align,
                       seconds = proc.time()[["elapsed"]] - t0)
  .bys_result(out, "bayes_ammi",
              label = sprintf("Bayesian %s with %d multiplicative term(s)%s: %d genotypes, %d environments, %s",
                              .model, K, if (is.null(tp)) "" else sprintf(", an earlier trial as prior (n0 = %g)", tp$n0),
                              G, E, if (cd$balanced) sprintf("%d replicates", cd$R) else
                                sprintf("%d plots (%d empty cells)", cd$n, cd$missing)),
              call = match.call(),
              meta = list(response = cd$names[["y"]], genotypes = G, environments = E, replicates = cd$B,
                          gen_labels = cd$G, env_labels = cd$E, k = K, level = .level, chains = .chains,
                          draws = length(keep) * .chains, model = if (gge) "GGE" else "AMMI", het = FALSE,
                          layout = lay, cells = list(Y = cd$Y, N = cd$N), balanced = cd$balanced,
                          missing = cd$missing, y_by_cell = cd$data$y[order(cd$data$cell)]))
}

# The earlier trial as prior (Perez-Elizalde et al. 2012, Section 2.2): its least-squares fit with the genotypes and
# environments in the order of the trial, the weight n0, the von Mises-Fisher parameters n0 Y0 V0 D0 and
# n0 Y0' U0 D0, and the prior of the error variance (a = n0 - 1, S0^2 its within-cell mean square)
.bys_trial_prior <- function(fn, data0, ex, env, cd, K, gge, n0) {
  c0 <- .bys_prep(fn, data0, ex$y, ex$gen, ex$env, ex$rep, env, allow_missing = TRUE)
  if (!setequal(c0$G, cd$G) || !setequal(c0$E, cd$E))
    .bys_stop(fn, "`.prior_data` must have the same genotypes and environments as `.data`.")
  Y0 <- c0$Y[cd$G, cd$E, drop = FALSE]
  f0 <- .bys_ls(Y0, K, gge)
  n0 <- if (is.null(n0)) c0$B else .bys_check_number(n0, ".n0", fn, min = 1e-6)
  # orient every prior term like the least-squares fit of the trial, so that the signs of the two agree
  ls1 <- .bys_ls(cd$Y, K, gge)
  for (k in seq_len(K)) if (sum(f0$alpha[, k] * ls1$alpha[, k]) < 0) {
    f0$alpha[, k] <- -f0$alpha[, k]
    f0$gamma[, k] <- -f0$gamma[, k]
  }
  Yf <- f0$filled
  s2 <- c0$ssw / max(1, c0$n - sum(c0$N > 0))
  list(n0 = n0, a = max(n0 - 1, 1), s2 = s2,
       P0 = n0 * Yf %*% f0$gamma %*% diag(f0$lambda, K),
       Q0 = n0 * crossprod(Yf, f0$alpha) %*% diag(f0$lambda, K),
       fit = list(mu = f0$mu, tau = f0$tau, delta = f0$delta, lambda = f0$lambda, alpha = f0$alpha,
                  gamma = f0$gamma, replicates = c0$B))
}

# the components of the result from the chains
.bys_assemble <- function(chains, cd, ls, K, keep, level, lay) {
  G  <- length(cd$G)
  E  <- length(cd$E)
  nk <- length(keep)
  mats <- lapply(chains, function(ch) {
    m <- cbind(ch$scalar, matrix(ch$alpha, nk), matrix(ch$gamma, nk))
    colnames(m) <- lay$parameter
    coda::mcmc(m, start = keep[1L], thin = if (nk > 1L) keep[2L] - keep[1L] else 1L)
  })
  draws <- coda::mcmc.list(mats)
  all   <- do.call(rbind, lapply(mats, unclass))
  N     <- nrow(all)
  long  <- function(Arr, n) {
    out <- matrix(0, N * n, K)
    for (k in seq_len(K)) out[, k] <- as.vector(t(Arr[, , k, drop = TRUE]))
    cbind(.nIter = rep(seq_len(N), each = n), rep(seq_len(n), N), out)
  }
  alphas1 <- long(.bys_block(all, lay, "alpha"), G)
  gammas1 <- long(.bys_block(all, lay, "gamma"), E)
  colnames(alphas1) <- c(".nIter", "Gen", paste0("alphas", seq_len(K)))
  colnames(gammas1) <- c(".nIter", "Env", paste0("gammas", seq_len(K)))
  dt <- function(M, nm) data.table::setnames(data.table::as.data.table(M), nm)
  sig <- .bys_block(all, lay, "sigma2")
  tau <- if (any(lay$block == "tau")) .bys_block(all, lay, "tau") else matrix(numeric(0), N, 0L)
  out <- list(
      mu1      = dt(all[, "mu", drop = FALSE], "mu")
    , tau1     = dt(sig, if (ncol(sig) == 1L) "tau" else paste0("tau", seq_len(ncol(sig))))
    , tao1     = if (ncol(tau)) dt(tau, paste0("tao", seq_len(G))) else data.table::data.table()
    , delta1   = dt(.bys_block(all, lay, "delta"), paste0("delta", seq_len(E)))
    , lambdas1 = dt(.bys_block(all, lay, "lambda"), paste0("lambdas", seq_len(K)))
    , alphas1  = alphas1
    , gammas1  = gammas1
    , alphas0  = dt(ls$alpha, paste0("alphas", seq_len(K)))
    , gammas0  = dt(ls$gamma, paste0("gammas", seq_len(K)))
    )
  out$summary     <- .bys_summary_table(draws, all, lay, level)
  out$scores      <- .bys_scores(all, lay, cd, K)
  out$draws       <- draws
  out$diagnostics <- .bys_diagnostics(draws, sum(!lay$block %in% c("alpha", "gamma")))
  out
}

# posterior summaries of every parameter: mean, sd, Monte Carlo standard error, HPD interval, effective sample size,
# potential scale reduction factor
.bys_summary_table <- function(draws, all, lay, level) {
  hpd <- .bys_hpd(all, level)
  ess <- coda::effectiveSize(draws)
  rhat <- if (length(draws) > 1L) {
    tryCatch(coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1L],
             error = function(e) rep(NA_real_, ncol(all)))
  } else rep(NA_real_, ncol(all))
  sdv <- apply(all, 2L, stats::sd)
  cbind(data.table::copy(lay),
        data.table::data.table(mean = colMeans(all), sd = sdv, mc_se = sdv / sqrt(pmax(ess, 1)), lower = hpd[, 1L],
                               upper = hpd[, 2L], ess = unname(ess), rhat = unname(rhat)))
}

# deviance information criterion and WAIC from the accumulators of the chains
.bys_fit_stats <- function(chains, out, cd, lay) {
  dev  <- unlist(lapply(chains, `[[`, "deviance"))
  all  <- .bys_all(out)
  pm   <- colMeans(all)
  # deviance at the posterior means of mu, tau, delta, sigma2 and of the cell means eta
  G <- length(cd$G)
  E <- length(cd$E)
  eta <- vapply(seq_len(E), function(j) colMeans(.bys_eta_env(all, lay, j, G)), numeric(G))
  s2  <- pm[lay$parameter[lay$block == "sigma2"]]
  s2  <- if (length(s2) == 1L) rep(s2, E) else s2
  N <- cd$N
  dhat <- sum(N * log(2 * pi * rep(s2, each = G)) + (cd$Q - 2 * eta * cd$S + N * eta^2) / rep(s2, each = G))
  pd   <- mean(dev) - dhat
  c(list(deviance = dev, Dbar = mean(dev), Dhat = dhat, pD = pd, DIC = dhat + 2 * pd), .bys_waic(chains))
}

# WAIC from the accumulators of the chains (log-sum-exp, running mean and sum of squares of the pointwise
# log-likelihood of every observed plot); an empty list when WAIC was not computed
.bys_waic <- function(chains) {
  res  <- list()
  lse  <- lapply(chains, `[[`, "lse")
  if (all(is.finite(lse[[1L]]))) {
    S  <- vapply(chains, function(ch) length(ch$deviance), numeric(1L))
    lt <- do.call(cbind, lse)
    mx <- apply(lt, 1L, max)
    lppd_i <- mx + log(rowSums(exp(lt - mx))) - log(sum(S))
    # pooled mean and variance of the pointwise log-likelihood over the chains
    mt <- do.call(cbind, lapply(chains, `[[`, "lmean"))
    m2 <- do.call(cbind, lapply(chains, `[[`, "lm2"))
    gm <- drop(mt %*% S) / sum(S)
    ss <- rowSums(m2) + drop(((mt - gm)^2) %*% S)
    p_i <- ss / (sum(S) - 1)
    res$lppd   <- sum(lppd_i)
    res$p_waic <- sum(p_i)
    res$WAIC   <- -2 * (res$lppd - res$p_waic)
    res$waic_i  <- -2 * (lppd_i - p_i)
    res$waic_se <- sqrt(length(p_i) * stats::var(res$waic_i))
  }
  res
}

# shortest intervals containing the proportion `level` of the draws (the HPD interval of a unimodal posterior;
# Chen and Shao 1999, the algorithm of Crossa et al. 2011)
.bys_hpd <- function(M, level) {
  M <- as.matrix(M)
  n <- nrow(M)
  s <- max(1L, min(n - 1L, ceiling(level * n) - 1L))
  t(apply(M, 2L, function(v) {
    v <- sort(v)
    w <- v[(s + 1L):n] - v[seq_len(n - s)]
    j <- which.min(w)
    c(v[j], v[j + s])
  }))
}

# biplot coordinates: posterior means of lambda_k^(1/2) alpha_ik and lambda_k^(1/2) gamma_jk
.bys_scores <- function(all, lay, cd, K) {
  lam <- .bys_block(all, lay, "lambda")
  A   <- .bys_block(all, lay, "alpha")
  Gm  <- .bys_block(all, lay, "gamma")
  ga <- vapply(seq_len(K), function(k) colMeans(A[, , k, drop = TRUE] * sqrt(lam[, k])), numeric(dim(A)[2L]))
  ea <- vapply(seq_len(K), function(k) colMeans(Gm[, , k, drop = TRUE] * sqrt(lam[, k])), numeric(dim(Gm)[2L]))
  ga <- matrix(ga, ncol = K)
  ea <- matrix(ea, ncol = K)
  colnames(ga) <- colnames(ea) <- paste0("PC", seq_len(K))
  data.table::data.table(type = rep(c("Genotype", "Environment"), c(nrow(ga), nrow(ea))), label = c(cd$G, cd$E),
                         mean = c(rowMeans(cd$Y, na.rm = TRUE), colMeans(cd$Y, na.rm = TRUE)), rbind(ga, ea))
}

# convergence diagnostics of the scalar parameters (Crossa et al. 2011: Raftery-Lewis, Heidelberger-Welch,
# Gelman-Rubin)
.bys_diagnostics <- function(draws, n_scalar) {
  sub <- coda::as.mcmc.list(lapply(draws, function(m) m[, seq_len(n_scalar), drop = FALSE]))
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  list(
      gelman  = if (length(sub) > 1L) safe(coda::gelman.diag(sub, autoburnin = FALSE, multivariate = FALSE))
    , raftery = safe(lapply(sub, coda::raftery.diag))
    , heidel  = safe(lapply(sub, coda::heidel.diag))
    )
}

.bys_print_bayes_ammi <- function(x, digits = 4L, ...) {
  s  <- x$summary
  st <- x$settings
  cat(attr(x, "label"), "\n")
  cat(sprintf("%d chain(s) x %d kept draws (iterations %d after a burn-in of %d, thinning %d); %.1f s\n\n",
              st$chains, st$kept, st$nIter, st$burnin, st$thin, st$seconds))
  cat(sprintf("Posterior means, standard deviations and %g%% HPD intervals\n", 100 * st$level))
  tab <- s[block %in% c("mu", "sigma2", "lambda", "tau", "delta"),
           list(parameter, mean, sd, lower, upper, ess = round(ess), rhat)]
  if (all(is.na(tab$rhat))) tab[, rhat := NULL]
  print(tab, digits = digits, nrows = 200L)
  cat("\nScores of the bilinear terms: x$summary (blocks alpha, gamma); biplot coordinates: x$scores\n")
  invisible(x)
}

.bys_plot_bayes_ammi <- function(x, which = c("trace", "density", "biplot", "effects"), pars = NULL, ...) {
  which <- match.arg(which)
  if (which == "biplot") return(.bys_plot_mean_biplot(x))
  if (which == "effects") return(.bys_plot_effects(x))
  if (is.null(pars)) pars <- c("mu", "sigma2", sprintf("lambda[%d]", seq_len(x$settings$k)))
  .bys_plot_draws(x, which, pars)
}

# traces or histograms of the draws of chosen parameters (a result with $draws and $summary)
.bys_plot_draws <- function(x, which, pars) {
  bad <- setdiff(pars, x$summary$parameter)
  if (length(bad)) .bys_stop("plot", "unknown parameter(s): ", paste(bad, collapse = ", "),
                             "; see x$summary$parameter.")
  d <- data.table::rbindlist(lapply(seq_along(x$draws), function(ch) {
    m  <- x$draws[[ch]][, pars, drop = FALSE]
    # traces show at most 1000 evenly spaced draws per chain (enough for the eye, small files)
    ix <- if (which == "trace" && nrow(m) > 1000L) unique(round(seq(1, nrow(m), length.out = 1000L))) else
      seq_len(nrow(m))
    m  <- m[ix, , drop = FALSE]
    data.table::data.table(chain = factor(ch), draw = rep(ix, length(pars)),
                           parameter = factor(rep(pars, each = nrow(m)), levels = pars), value = as.vector(m))
  }))
  if (which == "trace") {
    return(ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = draw, y = value, colour = chain)) +
             ggplot2::geom_line(linewidth = 0.25) +
             ggplot2::facet_wrap(ggplot2::vars(parameter), scales = "free_y", ncol = 1L) +
             ggplot2::scale_colour_manual(values = .bys_palette, name = "chain") +
             .bys_x_cont() +
             ggplot2::labs(x = "kept draw", y = NULL, title = "Traces of the Markov chains") +
             bys_theme())
  }
  ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 40L, fill = .bys_palette[1L], colour = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(ggplot2::vars(parameter), scales = "free", ncol = 2L) +
    .bys_x_cont() + .bys_y_cont(bars = TRUE) +
    ggplot2::labs(x = NULL, y = "draws", title = "Marginal posterior distributions") +
    bys_theme()
}

# the biplot of the posterior means of the coordinates (Crossa et al. 2011, Fig. 2)
.bys_plot_mean_biplot <- function(x) {
  s <- x$scores
  if (x$settings$k < 2L) .bys_stop("plot", "the biplot needs at least two bilinear terms (.k >= 2).")
  lim <- max(abs(c(s$PC1, s$PC2))) * 1.08
  sm  <- x$summary[block == "lambda"]
  e   <- s[type == "Environment"]
  ggplot2::ggplot(data = s, mapping = ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_segment(data = e, mapping = ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
                          colour = .bys_palette[2L], linewidth = 0.3,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))) +
    ggplot2::geom_point(mapping = ggplot2::aes(colour = type, shape = type), size = 2) +
    ggrepel::geom_text_repel(mapping = ggplot2::aes(label = label, colour = type), size = 3, family = "serif",
                             max.overlaps = 50L, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(Genotype = .bys_palette[1L], Environment = .bys_palette[2L]), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Genotype = 16L, Environment = 17L), name = NULL) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = expression(lambda[1]^{1/2} * alpha[i1] ~ ", " ~ lambda[1]^{1/2} * gamma[j1]),
                  y = expression(lambda[2]^{1/2} * alpha[i2] ~ ", " ~ lambda[2]^{1/2} * gamma[j2]),
                  title = "Bayesian AMMI biplot of the posterior means",
                  subtitle = sprintf("lambda_1 = %s, lambda_2 = %s (posterior means)",
                                     format(signif(sm$mean[1L], 4L)), format(signif(sm$mean[2L], 4L)))) +
    bys_theme()
}

# posterior means and HPD intervals of the main effects
.bys_plot_effects <- function(x) {
  s <- x$summary[block %in% c("tau", "delta")]
  s[, block := factor(block, levels = c("tau", "delta"), labels = c("genotype effects", "environment effects"))]
  s[, label := factor(label, levels = unique(label))]
  ggplot2::ggplot(data = s, mapping = ggplot2::aes(x = mean, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = lower, xmax = upper), height = 0.3,
                            colour = .bys_palette[1L]) +
    ggplot2::geom_point(size = 1.6, colour = .bys_palette[1L]) +
    ggplot2::facet_wrap(ggplot2::vars(block), scales = "free", ncol = 2L) +
    .bys_x_cont() +
    ggplot2::labs(x = sprintf("posterior mean and %g%% HPD interval", 100 * x$settings$level), y = NULL,
                  title = "Main effects") +
    bys_theme()
}

utils::globalVariables(c("sd", "block", "parameter", "lower", "upper", "ess", "rhat", "draw", "value", "chain", "type",
                         "PC1", "PC2", "label", ".nIter"))
