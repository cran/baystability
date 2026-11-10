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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 6,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("baystability", quietly = TRUE)) library(baystability) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----fit----------------------------------------------------------------------
fm <-
  bayes_ammi(
      .data       = cultivo2009
    , .y          = y
    , .gen        = entry
    , .env        = site
    , .rep        = rep
    , .nIter      = 5000
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


## ----conv---------------------------------------------------------------------
s <- fm$summary
c(max_rhat = max(s$rhat), min_ess_linear = min(s[block %in% c("mu", "tau", "delta"), ess]))


## ----b1-----------------------------------------------------------------------
pub <- utils::read.csv(file = system.file("extdata", "perez2012.csv", package = "baystability"))
pos <- c("mu", sprintf("tau[%d]", 1:12), sprintf("delta[%d]", 1:25))
lin <- s[block %in% c("mu", "tau", "delta")]
b1 <- data.table::data.table(
    parameter = c("mu", sprintf("alpha%d", 1:12), sprintf("beta%d", 1:25))
  , published = pub$mean[match(x = pos, table = pub$parameter)]
  , here      = round(x = lin$mean, digits = 1)
  , sd_pub    = pub$sd[match(x = pos, table = pub$parameter)]
  , sd_here   = round(x = lin$sd, digits = 2)
  )
b1[, diff_in_sd := round(x = (here - published) / sd_pub, digits = 3)]
b1


## ----b1f, include = FALSE-----------------------------------------------------
mz <- max(abs(b1$diff_in_sd))


## ----avg----------------------------------------------------------------------
both <-
  bys_ls_means(
      .data = rbind(
          cultivo2008[c("site", "rep", "entry", "y")]
        , cultivo2009[c("site", "rep", "entry", "y")]
        )
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
c(average_mu = both$grand_mean, published = 3848)


## ----t2-----------------------------------------------------------------------
sig <- sqrt(unlist(fm$tau1))
lam <- s[block == "lambda"]
t2 <- data.table::data.table(
    parameter = c("sigma", sprintf("lambda%d", 1:11))
  , published = pub$mean[match(x = c("sigma", sprintf("lambda[%d]", 1:11)), table = pub$parameter)]
  , here      = round(c(mean(sig), lam$mean))
  , sd_pub    = pub$sd[match(x = c("sigma", sprintf("lambda[%d]", 1:11)), table = pub$parameter)]
  , sd_here   = round(c(sd(sig), lam$sd))
  )
t2[, ratio := round(x = here / published, digits = 2)]
t2


## ----t2f, include = FALSE-----------------------------------------------------
zs <- (mean(sig) - 1311) / 39.13
l01 <- fm$prior$fit$lambda[1]


## ----scores-------------------------------------------------------------------
sc <- s[block %in% c("alpha", "gamma") & axis %in% 1:2]
sc[(block == "alpha" & label == levels(factor(cultivo2009$entry))[12]) |
     (block == "gamma" & label %in% c("5", "9", "23", "25")), list(
         parameter
       , axis
       , mean
       , sd
       , lower
       , upper
       )]


## ----pubsc, include = FALSE---------------------------------------------------
g <- s[block == "gamma" & axis == 1 & label %in% c("5", "9", "23", "25"), mean]


## ----fig4, fig.width = 9------------------------------------------------------
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

