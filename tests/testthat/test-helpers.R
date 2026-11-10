test_that("the tables and effects agree with the released version", {
  f <- file.path(legacy_dir(), "helpers.rds")
  skip_if(!nzchar(legacy_dir()) || !file.exists(f), "no legacy snapshot (source tree only)")
  old <- readRDS(f)
  expect_equal(bys_ls_means(.data = cultivo2009, .y = y, .gen = entry, .env = site)$ge_means, old$ge_mean$ge_means,
               ignore_attr = TRUE)
  expect_equal(bys_ls_means(.data = cultivo2009, .y = y, .gen = entry, .env = site)$grand_mean, old$ge_mean$grand_mean)
  expect_equal(bys_ls_var(.data = cultivo2009, .y = y, .gen = entry, .env = site)$ge_variances, old$ge_var$ge_variances,
               ignore_attr = TRUE)
  expect_equal(bys_ls_gen_effects(.data = cultivo2009, .y = y, .gen = entry, .env = site), old$g_eff, ignore_attr = TRUE)
  expect_equal(bys_ls_env_effects(.data = cultivo2009, .y = y, .gen = entry, .env = site), old$e_eff, ignore_attr = TRUE)
  ge <- bys_ls_ge_effects(.data = cultivo2009, .y = y, .gen = entry, .env = site)
  expect_equal(ge$ge_effects, old$ge_eff$ge_effects, ignore_attr = TRUE)
  expect_equal(ge$ge_svd$d, old$ge_eff$ge_svd$d)
  expect_equal(bys_contrasts(5), old$matrix_k)
  # the least-squares AMMI of version 0.2.0: the same decomposition (its error variance came from lme4)
  ga <- bys_ls_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_equal(ga$lambdas, old$ge_ammi$lambdas, ignore_attr = TRUE)
  expect_equal(abs(ga$alphas), abs(old$ge_ammi$alphas), ignore_attr = TRUE)
  expect_equal(ga$sigma2, old$ge_ammi$sigma2, tolerance = 1e-4)
  expect_equal(bys_orthonormal(cbind(1, 1:5, (1:5)^2)), old$orthnorm)
})

test_that("tables and effects from their definitions; columns by name, string or vector", {
  M <- tapply(cultivo2009$y, list(cultivo2009$entry, cultivo2009$site), mean)
  m <- bys_ls_means(.data = cultivo2009, .y = "y", .gen = "entry", .env = "site")
  expect_s3_class(m, "bys_ge_mean")
  expect_equal(m$ge_means, M, ignore_attr = TRUE)
  expect_equal(unname(drop(bys_ls_gen_effects(.data = cultivo2009, .y = y, .gen = entry, .env = site))), unname(rowMeans(M) - mean(M)))
  expect_equal(unname(drop(bys_ls_env_effects(.data = cultivo2009, .y = cultivo2009$y, .gen = cultivo2009$entry, .env = cultivo2009$site))),
               unname(colMeans(M) - mean(M)))
  Z <- bys_ls_ge_effects(.data = cultivo2009, .y = y, .gen = entry, .env = site)$ge_effects
  expect_equal(unname(rowSums(Z)), rep(0, 12))
  expect_equal(unname(colSums(Z)), rep(0, 25))
  v <- bys_ls_var(.data = cultivo2009, .y = y, .gen = entry, .env = site)$ge_variances
  expect_equal(unname(v), unname(tapply(cultivo2009$y, list(cultivo2009$entry, cultivo2009$site), var)))
})

test_that("ge_ammi: closed-form error variance equals the REML fit of bys_ls_model()", {
  a <- bys_ls_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_s3_class(a, "bys_ge_ammi")
  expect_equal(a$k, 11L)
  expect_equal(round(a$lambdas[1:2]), c(5204, 3886))
  expect_equal(a$tau, 1 / a$sigma2)
  skip_if_not_installed("lme4")
  fm <- bys_ls_model(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_equal(a$sigma2, stats::sigma(fm)^2, tolerance = 1e-6)
  # the boundary case: replicate effects removed, so the REML replicate variance is zero
  d <- cultivo2009
  d$y <- d$y - ave(d$y, d$site, d$rep) + ave(d$y, d$site)
  a2 <- bys_ls_ammi(.data = d, .y = y, .gen = entry, .env = site, .rep = rep)
  expect_equal(a2$sigma2, stats::sigma(bys_ls_model(.data = d, .y = y, .gen = entry, .env = site, .rep = rep))^2,
               tolerance = 1e-4)
})

test_that("matrix_k and orthnorm", {
  for (n in c(2, 3, 9, 20)) {
    K <- bys_contrasts(n)
    expect_equal(dim(K), c(n - 1, n))
    expect_equal(K %*% t(K), diag(n - 1))
    expect_equal(drop(K %*% rep(1, n)), rep(0, n - 1))
    expect_equal(t(K) %*% K, diag(n) - 1 / n)
  }
  expect_error(bys_contrasts(1))
  V <- bys_orthonormal(cbind(1, 1:6), basis = TRUE, norm = TRUE)
  expect_equal(crossprod(V), diag(6))
  V2 <- bys_orthonormal(cbind(1, 1:6), basis = FALSE, norm = FALSE)
  expect_equal(sum(V2[, 1] * V2[, 2]), 0)
  expect_null(bys_orthonormal(NULL))
  expect_warning(bys_orthonormal(matrix(1:6, 2)))
})

test_that("input errors are reported by name", {
  expect_error(bys_ls_means(.data = cultivo2009, .y = yield, .gen = entry, .env = site), "\\.y")
  d <- cultivo2009[!(cultivo2009$entry == 2 & cultivo2009$site == 1), ]
  expect_error(bys_ls_means(.data = d, .y = y, .gen = entry, .env = site), "genotype 2 in environment 1")
  d <- cultivo2009[-1, ]
  # unequal replicates are completed by data augmentation since version 1.0.0
  expect_s3_class(bayes_ammi(.data = d, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10), "bys_bayes_ammi")
  expect_error(bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10,
                          .model = "nope"))
  expect_error(bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10, .k = 12))
  other <- cultivo2008[cultivo2008$entry != 2, ]
  expect_error(bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 10,
                          .prior_data = other), "same genotypes and environments")
})
