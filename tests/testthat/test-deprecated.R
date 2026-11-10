test_that("the names of baystability 0.2.0 still work, with a warning, and give the results of the new names", {
  data(cultivo2009)
  expect_warning(old <- ge_mean(.data = cultivo2009, .y = y, .gen = entry, .env = site), "bys_ls_means")
  expect_equal(old$ge_means, bys_ls_means(.data = cultivo2009, .y = y, .gen = entry, .env = site)$ge_means)
  expect_warning(old <- ge_ammi(cultivo2009, y, entry, site, rep), "bys_ls_ammi")
  expect_equal(old$lambdas,
               bys_ls_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep)$lambdas)
  expect_warning(k <- matrix_k(4), "bys_contrasts")
  expect_equal(k, bys_contrasts(4))
  for (f in c("ge_var", "g_eff", "e_eff", "ge_eff"))
    expect_warning(do.call(f, list(.data = cultivo2009, .y = "y", .gen = "entry", .env = "site")), "bys_ls_")
  expect_warning(orthnorm(u = cbind(1, 1:4)), "bys_orthonormal")
  expect_warning(ge_model(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep), "bys_ls_model")
})

test_that("bayes_gge() is bayes_ammi(.model = 'SREG'), with its own biplot", {
  data(cultivo2008)
  data(cultivo2009)
  fg <- bayes_gge(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 150, .k = 2,
                  .prior_data = cultivo2008, .seed = 5)
  fs <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 150, .k = 2,
                   .prior_data = cultivo2008, .model = "SREG", .seed = 5)
  expect_equal(fg$summary$mean, fs$summary$mean)
  expect_identical(deparse(attr(fg, "call")[[1L]]), "bayes_gge")
  expect_s3_class(bayes_gge_biplot(fit = fg, thin = 1), "bys_biplots")
  expect_error(bayes_ammi_biplot(fit = fg), "not an AMMI fit")
  for (f in c("bayes_ammi", "bayes_gge", "bayes_ammi_biplot", "bayes_gge_biplot", "bys_ls_means", "bys_compare",
              "bys_interpret", "bys_cite"))
    expect_true(is.function(getS3method(f, "default")), info = f)
})
