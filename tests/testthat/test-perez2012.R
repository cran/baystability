test_that("Table B.1 and sigma of Perez-Elizalde et al. (2012) are reproduced", {
  skip_on_cran()
  ref <- utils::read.csv(system.file("extdata", "perez2012.csv", package = "baystability"))
  fm <- bayes_ammi(.data = cultivo2009, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 2500,
                   .prior_data = cultivo2008, .burnin = 500, .chains = 2, .seed = 2012)
  s <- fm$summary
  lin <- s[block %in% c("mu", "tau", "delta")]
  pub <- ref[ref$table == "B.1", ]
  pos <- c("mu", sprintf("tau[%d]", 1:12), sprintf("delta[%d]", 1:25))
  z <- (lin$mean - pub$mean[match(pos, pub$parameter)]) / pub$sd[match(pos, pub$parameter)]
  expect_lt(max(abs(z)), 0.15)
  expect_equal(lin$sd, pub$sd[match(pos, pub$parameter)], tolerance = 0.1)
  sig <- sqrt(unlist(fm$tau1))
  expect_lt(abs(mean(sig) - 1311) / 39.13, 0.5)
  expect_equal(sd(sig), 39.13, tolerance = 0.15)
})
