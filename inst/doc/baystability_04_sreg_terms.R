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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 5,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("baystability", quietly = TRUE)) library(baystability) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----sreg---------------------------------------------------------------------
fs <-
  bayes_ammi(
      .data       = cultivo2009
    , .y          = y
    , .gen        = entry
    , .env        = site
    , .rep        = rep
    , .nIter      = 3000
    , .prior_data = cultivo2008
    , .n0         = NULL
    , .model      = "SREG"
    , .k          = 2
    , .burnin     = 1000
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
    , .seed       = 2012
    , .align      = TRUE
    , .waic       = TRUE
    )
fs


## ----sreg-bip, fig.width = 9--------------------------------------------------
bs <-
  bayes_gge_biplot(
      fit           = fs
    , burnin        = 0
    , thin          = 1
    , pb            = 0.05
    , plot_stable   = FALSE
    , plot_unstable = TRUE
    , ncolors       = 5
    , .which        = "both"
    , .grid         = 100
    )
bs$biplot


## ----n0-----------------------------------------------------------------------
w <- c(0.5, 3, 12)
res <-
  data.table::rbindlist(
      l = lapply(
          X   = w
        , FUN = function(n0) {
            f <-
              bayes_ammi(
                  .data       = cultivo2009
                , .y          = y
                , .gen        = entry
                , .env        = site
                , .rep        = rep
                , .nIter      = 1500
                , .prior_data = cultivo2008
                , .n0         = n0
                , .model      = "AMMI"
                , .k          = 2
                , .burnin     = 500
                , .thin       = 1
                , .chains     = 1
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
                , .seed       = 4
                , .align      = TRUE
                , .waic       = TRUE
                )
            data.table::data.table(
                n0      = n0
              , mu      = f$summary[parameter == "mu", mean]
              , alpha1  = f$summary[block == "tau", mean][1]
              , lambda1 = f$summary[parameter == "lambda[1]", mean]
              )
          }
        )
    )
res


## ----n0f, include = FALSE-----------------------------------------------------
m08 <-
  bys_ls_means(
      .data = cultivo2008
    , .y    = y
    , .gen  = entry
    , .env  = site
    )$grand_mean
m09 <-
  bys_ls_means(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    )$grand_mean


## ----terms--------------------------------------------------------------------
fits <-
  lapply(
      X   = c(1, 2, 3, 5)
    , FUN = function(K) {
        bayes_ammi(
            .data       = cultivo2009
          , .y          = y
          , .gen        = entry
          , .env        = site
          , .rep        = rep
          , .nIter      = 1500
          , .prior_data = cultivo2008
          , .n0         = NULL
          , .model      = "AMMI"
          , .k          = K
          , .burnin     = 500
          , .thin       = 1
          , .chains     = 1
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
          , .seed       = 5
          , .align      = TRUE
          , .waic       = TRUE
          )
      }
    )
names(fits) <- paste0("K", c(1, 2, 3, 5))
cp <- do.call(what = bys_compare, args = fits)
cp


## ----termsf, include = FALSE--------------------------------------------------
t <- cp$table
close <- t$model[t$dWAIC <= 2 * t$dWAIC_se]

