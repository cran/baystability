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


## ----data---------------------------------------------------------------------
data(cultivo2008)
data(cultivo2009)
str(cultivo2009)
m08 <-
  bys_ls_means(
      .data = cultivo2008
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
m09 <-
  bys_ls_means(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
c(year_2008 = m08$grand_mean, year_2009 = m09$grand_mean)


## ----ls-----------------------------------------------------------------------
ls09 <-
  bys_ls_ammi(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
round(ls09$lambdas)


## ----fit----------------------------------------------------------------------
fm <-
  bayes_ammi(
      .data       = cultivo2009
    , .y          = y
    , .gen        = entry
    , .env        = site
    , .rep        = rep
    , .nIter      = 3000
    , .prior_data = cultivo2008
    , .n0         = NULL
    , .model      = "AMMI"
    , .k          = NULL
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
fm


## ----facts, include = FALSE---------------------------------------------------
s <- fm$summary
mu_post <- s[parameter == "mu", mean]


## ----trace, fig.height = 5.5--------------------------------------------------
plot(x = fm, which = "trace", pars = c("mu", "sigma2", "lambda[1]", "lambda[2]"))


## ----biplot, fig.width = 9, fig.height = 5.5----------------------------------
bp <-
  bayes_ammi_biplot(
      fit           = fm
    , burnin        = 0
    , thin          = 1
    , pb            = c(0.10, 0.05)
    , plot_stable   = FALSE
    , plot_unstable = TRUE
    , ncolors       = 5
    , .which        = "both"
    , .grid         = 100
    )
bp$biplot


## ----int----------------------------------------------------------------------
bys_interpret(
    x             = fm
  , higher_better = TRUE
  , top           = 5
  )

