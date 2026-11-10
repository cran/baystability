test_that("von Mises-Fisher and von Mises draws have the expected moments", {
  set.seed(1)
  for (p in c(2L, 5L, 18L)) {
    mu <- rep(1, p) / sqrt(p)
    for (kap in c(0.5, 5, 50)) {
      X <- t(replicate(3000, baystability:::.bys_rvmf(kap * mu)))
      expect_equal(rowSums(X^2), rep(1, 3000))
      r <- drop(X %*% mu)
      A <- besselI(kap, p / 2, expon.scaled = TRUE) / besselI(kap, p / 2 - 1, expon.scaled = TRUE)
      expect_equal(mean(r), A, tolerance = 4 * sqrt(var(r) / 3000) / A + 1e-3, label = paste(p, kap))
    }
  }
  for (kap in c(0, 0.5, 4, 40)) {
    th <- replicate(4000, baystability:::.bys_rvonmises(1, kap))
    A <- if (kap == 0) 0 else besselI(kap, 1, TRUE) / besselI(kap, 0, TRUE)
    expect_lt(abs(mean(cos(th - 1)) - A), 4 * stats::sd(cos(th - 1)) / sqrt(4000) + 1e-3)
  }
})

test_that("truncated normal draws stay in their interval", {
  set.seed(2)
  for (b in list(c(0, Inf), c(-1, 1), c(3, 4), c(8, Inf))) {
    x <- replicate(3000, baystability:::.bys_rtnorm(0, 1, b[1], b[2]))
    expect_true(all(x >= b[1] & x <= b[2]))
  }
})

test_that("the compiled sampler draws the same numbers as the R reference, with and without an earlier trial", {
  ns <- asNamespace("baystability")
  ex <- list(y = quote(y), gen = quote(entry), env = quote(site), rep = quote(rep))
  cd <- ns$.bys_prep("t", cultivo2009, ex$y, ex$gen, ex$env, ex$rep, environment(), allow_missing = TRUE)
  for (prior in c(FALSE, TRUE)) for (K in c(2L, 11L)) {
    ls <- ns$.bys_ls(cd$Y, K)
    pr <- bys_priors()
    pr$s2 <- 7e5
    tp <- if (prior) ns$.bys_trial_prior("t", cultivo2008, ex, environment(), cd, K, FALSE, NULL)
    if (prior) {
      pr$nu <- tp$a
      pr$s2 <- tp$s2
    }
    pl <- ns$.bys_prior_list(pr, 12, 25, K, tp)
    st <- ns$.bys_start(ls, 7e5, 1)
    set.seed(4)
    r <- ns$.bys_chain(cd$Y, 3, cd$ssw, K, pl, st, list(alpha = ls$alpha), 4, 1:4)
    set.seed(4)
    cc <- ns$.bys_run_chain(cd, K, FALSE, FALSE, pr, st, ls$alpha, 4, rep(1L, 4), TRUE, FALSE, tp)
    expect_equal(cc$scalar, r$scalar, tolerance = 1e-10, label = paste("scalars", prior, K))
    expect_equal(cc$alpha, r$alpha, tolerance = 1e-10, label = paste("alpha", prior, K))
    expect_equal(cc$gamma, r$gamma, tolerance = 1e-10, label = paste("gamma", prior, K))
  }
})

test_that("with all terms the genotype scores move (pair moves) and the chains agree", {
  fm <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 400,
                   .prior_data = cultivo2008, .burnin = 100, .chains = 2, .seed = 3)
  s <- fm$summary
  expect_equal(fm$settings$k, 11L)
  expect_gt(min(s[block == "alpha", sd]), 0.01)
  expect_lt(stats::median(s[block == "alpha", rhat]), 1.05)
  # constraints of every draw
  A <- baystability:::.bys_block(baystability:::.bys_all(fm), attr(fm, "meta")$layout, "alpha")
  a1 <- A[1L, , ]
  expect_equal(crossprod(a1), diag(11), tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(colSums(a1), rep(0, 11), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("the sampler recovers the parameters of simulated data", {
  set.seed(3)
  G <- 8; E <- 10; B <- 3
  a <- baystability:::.bys_orth1(matrix(rnorm(G * 2), G)); g <- baystability:::.bys_orth1(matrix(rnorm(E * 2), E))
  lam <- c(60, 30); tau <- rnorm(G, 0, 5); tau <- tau - mean(tau); del <- rnorm(E, 0, 10); del <- del - mean(del)
  M <- 100 + outer(tau, del, "+") + a %*% (lam * t(g))
  d <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
  d$y <- M[cbind(d$gen, d$env)] + rnorm(nrow(d), 0, 4)
  fit <- bayes_ammi(.data = d, .y = y, .gen = gen, .env = env, .rep = rep, .nIter = 2500, .k = 2, .burnin = 500,
                    .chains = 2, .seed = 4)
  s <- fit$summary
  expect_true(all(s[block %in% c("mu", "sigma2", "lambda", "tau", "delta"), rhat] < 1.05))
  covered <- function(b, truth) mean(s[block == b, lower] <= truth & truth <= s[block == b, upper])
  expect_gte(covered("tau", tau), 0.75)
  expect_gte(covered("delta", del), 0.75)
  expect_true(s[parameter == "sigma2", lower] < 16 & 16 < s[parameter == "sigma2", upper])
  expect_true(all(abs(s[block == "lambda", mean] - lam) < 4 * s[block == "lambda", sd]))
})

test_that("results keep the components of version 0.2.0; seeds reproduce runs", {
  set.seed(10)
  before <- .Random.seed
  f1 <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 30, .k = 2,
                   .seed = 7)
  expect_identical(.Random.seed, before)
  f2 <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 30, .k = 2,
                   .seed = 7)
  expect_equal(f1$mu1, f2$mu1)
  expect_equal(names(f1)[1:9], c("mu1", "tau1", "tao1", "delta1", "lambdas1", "alphas1", "gammas1", "Mean_Alphas",
                                 "Mean_Gammas"))
  expect_equal(names(f1$delta1), paste0("delta", 1:25))
  expect_equal(dim(f1$alphas1), c(30L, 24L))
  expect_equal(dim(f1$Mean_Gammas), c(25L, 2L))
  expect_true(all(f1$tau1$tau > 1e5))
  expect_output(print(f1), "Sources:")
})

test_that("SREG: no genotype effects, genotype scores sum to zero, environment scores do not", {
  fs <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 200, .k = 2,
                   .prior_data = cultivo2008, .model = "SREG", .seed = 5)
  expect_equal(ncol(fs$tao1), 0L)
  expect_equal(fs$settings$model, "SREG")
  A <- baystability:::.bys_block(baystability:::.bys_all(fs), attr(fs, "meta")$layout, "alpha")
  G <- baystability:::.bys_block(baystability:::.bys_all(fs), attr(fs, "meta")$layout, "gamma")
  expect_equal(colSums(A[5L, , ]), c(0, 0), tolerance = 1e-8)
  expect_gt(max(abs(colSums(G[5L, , ]))), 1e-4)
  expect_s3_class(bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 20,
                             .model = "GGE", .k = 2, .seed = 1), "bys_bayes_ammi")
  bp <- bayes_gge_biplot(fit = fs, burnin = 0, thin = 1, pb = 0.05)
  expect_s3_class(bp$biplot, "ggplot")
})
