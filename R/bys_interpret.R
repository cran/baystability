#' Interpretation of a Bayesian AMMI analysis
#'
#' @description
#' \code{bys_interpret()} turns a result of \code{\link{bayes_ammi}()}, \code{\link{bayes_ammi_biplot}()} or
#' \code{\link{bys_ls_ammi}()} into a plain-language interpretation: whether the chains have converged, which genotypes
#' and environments differ from the average (main effects whose HPD interval excludes zero), how strong the bilinear
#' terms are, which genotypes and environments interact significantly and in which direction, which are stable, and
#' what to check before acting on the results. Every statement is computed from the numbers of the result.
#'
#' @details
#' The rules:
#' \itemize{
#'   \item convergence: potential scale reduction factor below 1.1 (two or more chains), effective sample size of
#'     at least 400, Raftery-Lewis dependence factor below 5 and the Heidelberger-Welch tests passed (Crossa et
#'     al. 2011);
#'   \item main effects: a genotype or environment differs from the average when the HPD interval of its effect
#'     excludes zero;
#'   \item bilinear terms: a term whose HPD interval of \eqn{\lambda_k}{lambda_k} reaches down to near zero (lower
#'     limit below a tenth of the posterior mean) is weakly supported;
#'   \item interaction: a genotype (environment) interacts significantly on the first term when the HPD interval of
#'     its first score excludes zero; genotypes and environments with significant scores of the same sign interact
#'     positively (the genotype does better there than its main effect predicts), of opposite sign negatively;
#'   \item stability (\code{\link{bayes_ammi_biplot}()}): a genotype is stable when its HPD region contains the origin.
#' }
#'
#' @param x a result of \code{\link{bayes_ammi}()}, \code{\link{bayes_ammi_biplot}()} or \code{\link{bys_ls_ammi}()}.
#' @param higher_better \code{TRUE} (default) when larger values of the trait are better.
#' @param top number of genotypes or environments named in lists (default 5).
#' @param ... unused.
#' @return An object of class \code{"bys_interpretation"}: a list with \code{analysis}, \code{title},
#'   \code{headline}, \code{findings}, \code{recommendations} and \code{cautions} (sentences), \code{facts} (the
#'   numbers and lists behind the sentences), \code{guide} (how the method is read) and \code{sources} (the full
#'   references of the method; \code{\link{bys_cite}()} gives them as BibTeX). \code{print()} shows the text;
#'   \code{format(x, style = "markdown")} returns it as Markdown.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{bayes_ammi_biplot}}
#' @examples
#' data(cultivo2009)
#' fm1 <-
#'   bayes_ammi(
#'       .data       = cultivo2009
#'     , .y          = y
#'     , .gen        = entry
#'     , .env        = site
#'     , .rep        = rep
#'     , .nIter      = 500
#'     , .prior_data = NULL
#'     , .n0         = NULL
#'     , .model      = c("AMMI", "GGE")
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
#' txt <-
#'   bys_interpret(
#'       x             = fm1
#'     , higher_better = TRUE
#'     , top           = 5
#'     )
#' txt
#' cat(format(x = txt, style = "markdown"))
#' @export
bys_interpret <- function(x, higher_better = TRUE, top = 5L, ...) UseMethod("bys_interpret")

#' @rdname bys_interpret
#' @export
bys_interpret.default <- function(x, higher_better = TRUE, top = 5L, ...) {
  fn  <- "bys_interpret"
  higher_better <- .bys_check_flag(higher_better, "higher_better", fn)
  top <- .bys_check_number(top, "top", fn, min = 1, integer = TRUE)
  if (!inherits(x, "bys_result")) .bys_stop(fn, "`x` must be a result of a function of the baystability package.")
  cls <- .bys_cls(x)
  f <- get0(paste0(".bys_int_", cls), envir = asNamespace("baystability"), inherits = FALSE)
  if (is.null(f)) .bys_stop(fn, "no interpretation is defined for ", cls, " results.")
  ctx <- list(hb = higher_better, top = top, good = if (higher_better) "above" else "below",
              bad = if (higher_better) "below" else "above")
  out <- f(x, ctx)
  out$analysis <- cls
  out$guide    <- .bys_guides[[cls]]
  out$sources  <- unname(vapply(.bys_refs_list()[unique(.bys_sources_of(x))], `[[`, "", "text"))
  structure(out[c("analysis", "title", "headline", "findings", "recommendations", "cautions", "facts", "guide",
                  "sources")],
            class = "bys_interpretation")
}

#' @rdname bys_interpret
#' @export
print.bys_interpretation <- function(x, ...) {
  cat(format(x, style = "text"), sep = "\n")
  invisible(x)
}

#' @rdname bys_interpret
#' @param style \code{"text"} (default) or \code{"markdown"}.
#' @param width line width of the text (default \code{getOption("width")}).
#' @export
format.bys_interpretation <- function(x, style = c("text", "markdown"), width = getOption("width", 80L), ...) {
  style <- match.arg(style)
  wrap <- function(s, prefix) {
    if (!length(s)) return(character(0))
    unlist(lapply(s, function(z) {
      w <- strwrap(z, width = max(40L, width - 2L) - nchar(prefix))
      c(paste0(prefix, w[1L]), if (length(w) > 1L) paste0(strrep(" ", nchar(prefix)), w[-1L]))
    }))
  }
  sec <- function(h, s) if (length(s)) c("", if (style == "markdown") paste("###", h) else paste0(h, ":"),
                                          wrap(s, "- ")) else character(0)
  c(if (style == "markdown") paste("##", x$title) else x$title,
    "", wrap(x$headline, ""),
    sec("Findings", x$findings), sec("Recommendations", x$recommendations), sec("Cautions", x$cautions),
    sec("Sources", x$sources))
}

.bys_list <- function(x, max = Inf) {
  x <- as.character(x)
  if (!length(x)) return("none")
  more <- length(x) - max
  if (more > 0) x <- c(x[seq_len(max)], sprintf("%d more", more))
  if (length(x) == 1L) return(x)
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}
.bys_n <- function(v, d = 4L) {
  vapply(v, function(z) format(signif(z, d), big.mark = ",", trim = TRUE, scientific = FALSE), "")
}

# convergence facts of a bayes_ammi result
.bys_conv_facts <- function(x) {
  s  <- x$summary
  sc <- s[block %in% c("mu", "sigma2", "tau", "delta", "lambda", "slope", "sigma2_d")]
  dg <- x$diagnostics
  rl <- if (!is.null(dg$raftery)) {
    unlist(lapply(dg$raftery, function(r) if (is.matrix(r$resmatrix)) r$resmatrix[, "I"] else NA_real_))
  }
  hw <- if (!is.null(dg$heidel)) {
    unlist(lapply(dg$heidel, function(h) c(h[, "stest"], h[, "htest"])))
  }
  list(max_rhat = if (all(is.na(s$rhat))) NA_real_ else max(s$rhat, na.rm = TRUE),
       min_ess = min(s$ess), min_ess_scalar = min(sc$ess), max_dependence = if (any(is.finite(rl))) max(rl, na.rm = TRUE) else NA_real_,
       heidel_failed = if (length(hw)) sum(hw == 0, na.rm = TRUE) else NA_integer_,
       heidel_tests = if (length(hw)) sum(!is.na(hw)) else NA_integer_, chains = x$settings$chains,
       draws = x$settings$kept * x$settings$chains)
}

.bys_int_bayes_ammi <- function(x, ctx) {
  s   <- x$summary
  K   <- x$settings$k
  lev <- 100 * x$settings$level
  cf  <- .bys_conv_facts(x)
  # Heidelberger-Welch: each test has level 0.05, so a few of the many tests fail by chance; more than 10% is a
  # sign of non-stationary chains
  hw_bad <- !is.na(cf$heidel_failed) && cf$heidel_tests > 0 && cf$heidel_failed / cf$heidel_tests > 0.1
  ok  <- (is.na(cf$max_rhat) || cf$max_rhat < 1.1) && cf$min_ess_scalar >= 400 &&
    (is.na(cf$max_dependence) || cf$max_dependence < 5) && !hw_bad
  lam <- s[block == "lambda"]
  weak <- lam$axis[lam$lower < 0.1 * lam$mean]
  sig_pos <- function(b) s[block == b & lower > 0, label]
  sig_neg <- function(b) s[block == b & upper < 0, label]
  gen_good <- if (ctx$hb) sig_pos("tau") else sig_neg("tau")
  gen_bad  <- if (ctx$hb) sig_neg("tau") else sig_pos("tau")
  env_hi   <- sig_pos("delta")
  env_lo   <- sig_neg("delta")
  a1 <- s[block == "alpha" & axis == 1L]
  g1 <- s[block == "gamma" & axis == 1L]
  ga_pos <- a1[lower > 0, label]
  ga_neg <- a1[upper < 0, label]
  ge_pos <- g1[lower > 0, label]
  ge_neg <- g1[upper < 0, label]
  stable_good <- setdiff(gen_good, c(ga_pos, ga_neg))
  mu  <- s[block == "mu"]
  s2  <- s[block == "sigma2"]
  cv  <- 100 * sqrt(s2$mean) / mu$mean
  f <- c(cf, list(converged = ok, lambda = stats::setNames(lam$mean, paste0("lambda", lam$axis)),
                  lambda_lower = lam$lower, weak_terms = weak, genotypes_above = gen_good,
                  genotypes_below = gen_bad, environments_high = env_hi, environments_low = env_lo,
                  axis1_genotypes_positive = ga_pos, axis1_genotypes_negative = ga_neg,
                  axis1_environments_positive = ge_pos, axis1_environments_negative = ge_neg,
                  stable_and_good = stable_good, cv = cv))
  conv <- if (ok) {
    sprintf("The chains have converged by the criteria of Crossa et al. (2011)%s; the smallest effective sample size of the main parameters is %s.",
            if (is.na(cf$max_rhat)) " that apply to one chain" else sprintf(" (largest potential scale reduction factor %.3f)", cf$max_rhat),
            .bys_n(cf$min_ess_scalar, 3L))
  } else {
    sprintf("Convergence is not established: %s.", paste(c(
        if (!is.na(cf$max_rhat) && cf$max_rhat >= 1.1) sprintf("potential scale reduction factor up to %.2f", cf$max_rhat)
      , if (cf$min_ess_scalar < 400) sprintf("effective sample size as small as %s", .bys_n(cf$min_ess_scalar, 3L))
      , if (!is.na(cf$max_dependence) && cf$max_dependence >= 5) sprintf("Raftery-Lewis dependence factor up to %.1f", cf$max_dependence)
      , if (hw_bad) sprintf("%d of %d Heidelberger-Welch tests failed", cf$heidel_failed, cf$heidel_tests)
      ), collapse = "; "))
  }
  inter <- c(
      sprintf("On the first bilinear term, genotypes with significantly positive scores: %s; negative: %s.",
              .bys_list(ga_pos, ctx$top), .bys_list(ga_neg, ctx$top))
    , sprintf("Environments with significantly positive first scores: %s; negative: %s.", .bys_list(ge_pos, ctx$top),
              .bys_list(ge_neg, ctx$top))
    , if (length(c(ga_pos, ga_neg)) && length(c(ge_pos, ge_neg)))
      "Genotypes and environments with significant first scores of the same sign interact positively: the genotype yields more there than its main effect predicts; of opposite sign, less.")
  list(title = "Bayesian AMMI analysis",
       headline = sprintf("%d bilinear term(s): lambda_1 = %s%s; %s", K, .bys_n(lam$mean[1L]),
                          if (K > 1L) sprintf(", lambda_2 = %s", .bys_n(lam$mean[2L])) else "",
                          if (ok) "the chains have converged." else "convergence is not yet established."),
       findings = c(
           conv
         , if (!is.na(cf$heidel_failed)) sprintf("%d of %d Heidelberger-Welch stationarity and half-width tests failed%s.",
                                                 cf$heidel_failed, cf$heidel_tests,
                                                 if (hw_bad) "" else " (about 5% fail by chance at level 0.05)")
         , sprintf("Overall mean %s (%g%% HPD %s to %s); error variance %s, a coefficient of variation of %.1f%%.",
                   .bys_n(mu$mean), lev, .bys_n(mu$lower), .bys_n(mu$upper), .bys_n(s2$mean), cv)
         , sprintf("Genotypes significantly %s average (HPD of the effect excludes zero): %s; %s average: %s.", ctx$good,
                   .bys_list(gen_good, ctx$top), ctx$bad, .bys_list(gen_bad, ctx$top))
         , sprintf("Environments with significantly high means: %s; low: %s.", .bys_list(env_hi, ctx$top),
                   .bys_list(env_lo, ctx$top))
         , sprintf("Singular values with %g%% HPD intervals: %s.", lev,
                   paste(sprintf("lambda_%d = %s (%s to %s)", lam$axis, .bys_n(lam$mean), .bys_n(lam$lower),
                                 .bys_n(lam$upper)), collapse = "; "))
         , inter),
       recommendations = c(
           if (length(stable_good)) sprintf("%s are %s average and do not interact significantly on the first term: candidates for wide recommendation.",
                                            .bys_list(stable_good, ctx$top), ctx$good)
         , if (length(c(ga_pos, ga_neg))) "Recommend the genotypes with significant first scores for the environments whose scores have the same sign (bayes_ammi_biplot() shows the regions)."
         , "Draw bayes_ammi_biplot() to see which genotypes and environments form homogeneous groups (overlapping HPD regions)."),
       cautions = c(
           if (!ok) "Run longer chains (and two or more chains) before relying on the estimates."
         , if (cf$chains < 2L) "With one chain the Gelman-Rubin diagnostic is not available; use .chains = 2 or more."
         , if (length(weak)) sprintf("The HPD interval of lambda_%s reaches down to near zero: the term is weakly supported and its scores are uncertain.",
                                     paste(weak, collapse = ", lambda_"))
         , if (!isTRUE(x$settings$align)) "Signs were not aligned (.align = FALSE): posterior means of scores of weak terms shrink towards zero."),
       facts = f)
}

.bys_int_biplots <- function(x, ctx) {
  b  <- x$biplot_data
  lv <- 100 * attr(x, "meta")$level
  gen <- b[byplot == "Genotype"]
  env <- b[byplot == "Environment"]
  f <- list(stable_genotypes = gen[stable == TRUE, name], unstable_genotypes = gen[stable == FALSE, name],
            stable_environments = env[stable == TRUE, name], unstable_environments = env[stable == FALSE, name],
            positive_genotypes = gen[group == "positive", name], negative_genotypes = gen[group == "negative", name],
            positive_environments = env[group == "positive", name], negative_environments = env[group == "negative", name])
  list(title = "Biplot with highest posterior density regions",
       headline = sprintf("%d of %d genotypes and %d of %d environments have a %g%% HPD region that excludes the origin (significant interaction).",
                          length(f$unstable_genotypes), nrow(gen), length(f$unstable_environments), nrow(env), lv),
       findings = c(
           if (nrow(gen)) sprintf("Genotypes whose region contains the origin (interaction not different from zero, stable): %s.",
                                  .bys_list(f$stable_genotypes, ctx$top))
         , if (nrow(gen)) sprintf("Groups on the first term (marginal HPD interval excludes zero): positive %s; negative %s.",
                                  .bys_list(f$positive_genotypes, ctx$top), .bys_list(f$negative_genotypes, ctx$top))
         , if (nrow(env)) sprintf("Environments on the first term: positive %s; negative %s; the others do not differ from zero.",
                                  .bys_list(f$positive_environments, ctx$top), .bys_list(f$negative_environments, ctx$top))
         , "Overlapping regions mark genotypes (environments) whose interaction parameters do not differ: a homogeneous group."),
       recommendations = c(
           if (length(f$positive_genotypes) && length(f$positive_environments))
             sprintf("Genotypes %s are adapted to environments %s (both positive on the first term).",
                     .bys_list(f$positive_genotypes, ctx$top), .bys_list(f$positive_environments, ctx$top))
         , if (length(f$negative_genotypes) && length(f$negative_environments))
             sprintf("Genotypes %s are adapted to environments %s (both negative on the first term).",
                     .bys_list(f$negative_genotypes, ctx$top), .bys_list(f$negative_environments, ctx$top))),
       cautions = c("The regions are estimated from the draws by a kernel density; with few draws their outline is rough.",
                    "The second term is usually estimated with much more uncertainty than the first: long regions along the second axis."),
       facts = f)
}

.bys_int_ge_ammi <- function(x, ctx) {
  pct <- 100 * x$lambdas^2 / sum(x$lambdas^2)
  f <- list(lambdas = x$lambdas, percent = pct, sigma2 = x$sigma2)
  list(title = "Least-squares AMMI",
       headline = sprintf("The first two of %d interaction axes hold %.1f%% and %.1f%% of the interaction sum of squares.",
                          x$k, pct[1L], pct[2L]),
       findings = sprintf("Singular values: %s; error variance %s.", paste(.bys_n(utils::head(x$lambdas, 4L)), collapse = ", "),
                          .bys_n(x$sigma2)),
       recommendations = "Use bayes_ammi() for the posterior distributions and bayes_ammi_biplot() for inference on the biplot.",
       cautions = "Least-squares singular values are biased upwards (Crossa et al. 2011); the Bayesian estimates correct this.",
       facts = f)
}







.bys_int_compare <- function(x, ctx) {
  t <- x$table
  sed   <- if (all(is.na(t$dWAIC_se))) t$se else t$dWAIC_se
  close <- t$model[t$dWAIC <= 2 * sed]
  f <- list(best = t$model[1L], close = close)
  list(title = "Model comparison",
       headline = sprintf("%s has the smallest WAIC%s.", t$model[1L],
                          if (length(close) > 1L) sprintf("; %s within two standard errors of the difference", .bys_list(setdiff(close, t$model[1L]))) else ""),
       findings = sprintf("%s: WAIC %s (%s), DIC %s, lower HPD limit of the last singular value %s.",
                          t$model, .bys_n(t$WAIC, 5L),
                          ifelse(t$dWAIC == 0, "the best", sprintf("%s above the best, standard error %s",
                                                                   .bys_n(t$dWAIC, 3L), .bys_n(sed, 3L))),
                          .bys_n(t$DIC, 5L), .bys_n(t$lambda_last_lower, 3L)),
       recommendations = "Among models within two standard errors of the best, prefer the simplest; a last singular value whose HPD interval reaches zero adds little.",
       cautions = "WAIC and DIC compare predictive accuracy for new plots of the same genotypes and environments.",
       facts = f)
}




.bys_guides <- list(
    bayes_ammi = "Bayesian AMMI: every parameter has a posterior distribution, summarised by its mean, standard deviation and highest posterior density (HPD) interval. Main effects whose HPD interval excludes zero differ from the average. Singular values measure the strength of each bilinear term. Genotype and environment scores with HPD intervals excluding zero interact significantly; same signs mean positive interaction. Check convergence first: potential scale reduction factor near 1, large effective sample sizes, Raftery-Lewis dependence factor below 5, Heidelberger-Welch tests passed.",
    biplots = "Biplot of lambda^(1/2) times the scores of the first two terms with HPD regions. A region containing the origin: no significant interaction (stable). Regions that do not overlap: genotypes (environments) that differ in interaction. The marginal HPD interval of the first coordinate decides significance on the first term, as in Crossa et al. (2011).",
    compare = "WAIC and DIC: smaller predicts new data better; a difference within two of its standard errors is not decisive.",
    ge_ammi = "Least-squares AMMI: singular value decomposition of the interaction effects. The share of each axis in the interaction sum of squares measures its importance; the estimates start the Bayesian sampler."
  )

utils::globalVariables(c("axis", "group"))
