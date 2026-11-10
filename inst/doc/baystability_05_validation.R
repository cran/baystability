## ----setup, include = FALSE---------------------------------------------------
# HTML figures are drawn as PNG and kept as WebP: lossless (the same pixels in about half the space of PNG), or
# lossy at quality 90 when that is less than half the size, as for traces and densities of thousands of draws. The
# PDF figures stay vector graphics. Any failure of magick (not installed, no WebP, or a broken installation whose
# system libraries are missing) keeps the PNG
webp <- function(path) {
  if (knitr::is_latex_output() || !grepl(pattern = "[.]png$", x = path)) return(path)
  tryCatch(
      {
        if (!requireNamespace("magick", quietly = TRUE) || !isTRUE(magick::magick_config()$webp)) return(path)
        img <- magick::image_read(path)
        out <- sub(pattern = "[.]png$", replacement = ".webp", x = path)
        tmp <- tempfile(fileext = ".webp")
        magick::image_write(
            image   = img
          , path    = out
          , format  = "webp"
          , quality = 100
          , defines = c(`webp:lossless` = "true")
          )
        magick::image_write(image = img, path = tmp, format = "webp", quality = 90)
        if (file.size(tmp) < file.size(out) / 2) file.copy(from = tmp, to = out, overwrite = TRUE)
        unlink(c(path, tmp))
        out
      }
    , error = function(e) path
    )
}
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 4.2,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("baystability", quietly = TRUE)) library(baystability) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----same---------------------------------------------------------------------
ns <- asNamespace(ns = "baystability")
ex <- list(y = quote(y), gen = quote(entry), env = quote(site), rep = quote(rep))
cd <-
  ns$.bys_prep(
      "check"
    , cultivo2009
    , ex$y
    , ex$gen
    , ex$env
    , ex$rep
    , environment()
    , allow_missing = TRUE
    )
K  <- 11
ls <- ns$.bys_ls(cd$Y, K)
tp <- ns$.bys_trial_prior("check", cultivo2008, ex, environment(), cd, K, FALSE, NULL)
pr <-
  bys_priors(
      mu_mean     = 0
    , mu_var      = 1e+8
    , tau_mean    = 0
    , tau_var     = 1e+8
    , delta_mean  = 0
    , delta_var   = 1e+8
    , lambda_mean = 0
    , lambda_var  = 1e+8
    , nu          = 1
    , s2          = NULL
    )
pr$nu <- tp$a
pr$s2 <- tp$s2
pl <- ns$.bys_prior_list(pr, 12, 25, K, tp)
st <- ns$.bys_start(ls, tp$s2, 1)
set.seed(seed = 9)
r  <- ns$.bys_chain(cd$Y, 3, cd$ssw, K, pl, st, list(alpha = ls$alpha), 5, 1:5)
set.seed(seed = 9)
cc <- ns$.bys_run_chain(cd, K, FALSE, FALSE, pr, st, ls$alpha, 5, rep(1L, 5), TRUE, FALSE, tp)
c(scalars = max(abs(cc$scalar - r$scalar) / abs(r$scalar)), alpha = max(abs(cc$alpha - r$alpha)),
  gamma = max(abs(cc$gamma - r$gamma)))


## ----pairs--------------------------------------------------------------------
fm <-
  bayes_ammi(
      .data       = cultivo2009
    , .y          = y
    , .gen        = entry
    , .env        = site
    , .rep        = rep
    , .nIter      = 1000
    , .prior_data = cultivo2008
    , .n0         = NULL
    , .model      = "AMMI"
    , .k          = NULL
    , .burnin     = 200
    , .thin       = 1
    , .chains     = 2
    , .priors     = bys_priors(
          mu_mean     = 0
        , mu_var      = 1e+8
        , tau_mean    = 0
        , tau_var     = 1e+8
        , delta_mean  = 0
        , delta_var   = 1e+8
        , lambda_mean = 0
        , lambda_var  = 1e+8
        , nu          = 1
        , s2          = NULL
        )
    , .level      = 0.95
    , .seed       = 3
    , .align      = TRUE
    , .waic       = TRUE
    )
fm$summary[block == "alpha", list(
    min_sd      = min(sd)
  , median_rhat = stats::median(rhat)
  , max_rhat    = max(rhat)
  )]


## ----sim----------------------------------------------------------------------
set.seed(seed = 8)
G <- 8
E <- 12
B <- 3
A  <- qr.Q(qr(cbind(1, matrix(data = rnorm(G * 2), nrow = G))))[, 2:3]
Gm <- qr.Q(qr(cbind(1, matrix(data = rnorm(E * 2), nrow = E))))[, 2:3]
tau <- as.vector(scale(x = rnorm(n = G, mean = 0, sd = 6), scale = FALSE))
del <- as.vector(scale(x = rnorm(n = E, mean = 0, sd = 10), scale = FALSE))
M <- 100 + outer(X = tau, Y = del, FUN = "+") + A %*% (c(80, 40) * t(Gm))
sim <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
sim$y <- M[cbind(sim$gen, sim$env)] + rnorm(n = nrow(sim), mean = 0, sd = 5)
fs <-
  bayes_ammi(
      .data       = sim
    , .y          = y
    , .gen        = gen
    , .env        = env
    , .rep        = rep
    , .nIter      = 3000
    , .prior_data = NULL
    , .n0         = NULL
    , .model      = "AMMI"
    , .k          = 2
    , .burnin     = 500
    , .thin       = 1
    , .chains     = 2
    , .priors     = bys_priors(
          mu_mean     = 0
        , mu_var      = 1e+8
        , tau_mean    = 0
        , tau_var     = 1e+8
        , delta_mean  = 0
        , delta_var   = 1e+8
        , lambda_mean = 0
        , lambda_var  = 1e+8
        , nu          = 1
        , s2          = NULL
        )
    , .level      = 0.95
    , .seed       = 9
    , .align      = TRUE
    , .waic       = TRUE
    )
s <- fs$summary
c(tau_covered   = mean(s[block == "tau", lower] <= tau & tau <= s[block == "tau", upper]),
  delta_covered = mean(s[block == "delta", lower] <= del & del <= s[block == "delta", upper]))
s[block %in% c("sigma2", "lambda"), list(parameter, truth = c(25, 80, 40), mean, lower, upper)]


## ----ref----------------------------------------------------------------------
bys_reference(
    source   = "Perez"
  , evaluate = TRUE
  , level    = 1
  )[, list(quantity, value, computed, agrees)]

