test_that("the application fits a model with the earlier trial and fills every page", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  app <- bys_app_object()
  expect_s3_class(app, "shiny.appobj")
  shiny::testServer(app, {
    session$setInputs(source = "cultivo", col_y = "y", col_gen = "entry", col_env = "site", col_rep = "rep",
                      hb = TRUE, prior_source = "cultivo", model = "AMMI", k = 2, n0 = NA, iter = 200, burnin = 50,
                      thin = 1, chains = 2, seed = 1, level = 0.95, lambda_var = 1e8, nu = 1, fit_plot = "trace",
                      pb = c("0.05"), which = "both", stable = TRUE, bp_thin = 1)
    session$setInputs(run = 1)
    expect_s3_class(fit(), "bys_bayes_ammi")
    expect_false(is.null(fit()$settings$prior_trial))
    expect_s3_class(bp(), "bys_biplots")
    expect_match(code(), "bayes_ammi\\(")
    expect_match(code(), ".prior_data = cultivo2008", fixed = TRUE)
    session$setInputs(cmp_k = c("1", "2"), cmp_prior = TRUE, compare = 1)
    expect_s3_class(cmp(), "bys_compare")
    expect_equal(nrow(cmp()$table), 4L)
    expect_match(paste(output$refs, collapse = "\n"), "Perez-Elizalde")
  })
})

test_that("the R code of the application runs without the earlier trial", {
  a <- list(.y = "y", .gen = "entry", .env = "site", .rep = "rep", .nIter = 100L, .n0 = NULL, .model = "SREG",
            .k = NULL, .burnin = 10L, .thin = 1L, .chains = 1L, .priors = bys_priors(), .level = 0.95, .seed = 1L)
  code <- baystability:::.bys_app_code(a = a, source = "cultivo", prior_source = "none")
  expect_match(code, ".prior_data = NULL", fixed = TRUE)
  expect_match(code, ".k          = NULL", fixed = TRUE)
  env <- new.env()
  expect_silent(eval(parse(text = sub("library(baystability)\n", "", code, fixed = TRUE)), envir = env))
  expect_s3_class(env$fit, "bys_bayes_ammi")
})

test_that("plots for the application are WebP when magick can write it", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  f <- baystability:::.bys_plot_file(p, 400, 300)
  expect_true(file.exists(f$src))
  if (isTRUE(tryCatch(requireNamespace("magick", quietly = TRUE) && isTRUE(magick::magick_config()$webp),
                      error = function(e) FALSE)))
    expect_equal(f$contentType, "image/webp")
})
