test_that("glm.fit BIC Bayes factor matches its deviance expression", {
  data <- glm_test_data()
  X <- model.matrix(~ x1, data = data)
  null_fit <- glm(y ~ 1, data = data, family = binomial())
  model_fit <- glm.fit(x = X, y = data$y, family = binomial())
  n <- nrow(data)
  k <- 1

  result <- MissingBVS:::BF.BIC.glm.fit(
    y = data$y,
    X = X,
    family = binomial(),
    devnull = null_fit$deviance,
    n = n,
    k = k
  )
  expected <- (null_fit$deviance - model_fit$deviance - k * log(n)) / 2

  expect_equal(as.numeric(result), as.numeric(expected))
})

test_that("glm.fit TBF Bayes factor passes k and deviance to its method", {
  data <- glm_test_data()
  X <- model.matrix(~ x1, data = data)
  received <- list()
  lTBF.method <- function(k, dev) {
    received <<- list(k = k, dev = dev)
    as.numeric(k + dev)
  }

  result <- MissingBVS:::BF.TBF.glm.fit(
    y = data$y,
    X = X,
    family = binomial(),
    n = nrow(data),
    k = 1,
    lTBF.method = lTBF.method
  )
  expected_dev <- glm.fit(x = X, y = data$y, family = binomial())$deviance

  expect_identical(received$k, 1)
  expect_equal(as.numeric(received$dev), as.numeric(expected_dev))
  expect_equal(as.numeric(result), 1 + as.numeric(expected_dev))
})

test_that("BAS glm Bayes factors are finite and respect log-marginal shifts", {
  data <- glm_test_data()
  X <- model.matrix(~ x1, data = data)
  null_fit <- glm(y ~ 1, data = data, family = binomial())
  n <- nrow(data)

  bic <- function(logmargnull) {
    MissingBVS:::BF.BIC.glm(
      y = data$y,
      X = X,
      family = binomial(),
      logmargnull = logmargnull,
      n = n,
      k = 1,
      p0 = 1
    )
  }
  tbf_prior <- BAS::testBF.prior(g = n)
  tbf_prior$hyper.parameters$loglik_null <- as.numeric(-0.5 * null_fit$deviance)
  tbf <- function(logmargnull) {
    MissingBVS:::BF.TBF.glm(
      y = data$y,
      X = X,
      family = binomial(),
      prior.betas = tbf_prior,
      logmargnull = logmargnull,
      k = 1,
      p0 = 1
    )
  }
  gprior <- function(logmargnull) {
    MissingBVS:::BF.gprior.glm(
      y = data$y,
      X = X,
      family = binomial(),
      prior.betas = BAS::g.prior(g = n),
      logmargnull = logmargnull,
      k = 1,
      p0 = 1
    )
  }

  for (method in list(bic, tbf, gprior)) {
    at_zero <- method(0)
    shifted <- method(2)
    expect_length(at_zero, 1)
    expect_true(is.finite(at_zero))
    expect_equal(as.numeric(at_zero - shifted), 2)
  }
})

test_that("checkforprior.betas.glm selects BAS methods", {
  data <- glm_test_data()
  null_fit <- glm(y ~ 1, data = data, family = binomial(), x = TRUE, y = TRUE)
  X <- model.matrix(~ x1, data = data)
  y <- as.numeric(null_fit$y)
  n <- nrow(data)

  bic <- MissingBVS:::checkforprior.betas.glm(
    "BIC", NULL, TRUE, n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
  )
  tbf <- MissingBVS:::checkforprior.betas.glm(
    "TBF", "gZellner", TRUE, n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
  )
  gprior <- MissingBVS:::checkforprior.betas.glm(
    "gprior", "gZellner", TRUE, n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
  )

  expect_true(is.finite(bic(k = 1, X = X)))
  expect_true(is.finite(tbf(k = 1, X = X)))
  expect_true(is.finite(gprior(k = 1, X = X)))
})

test_that("checkforprior.betas.glm rejects unsupported choices", {
  data <- glm_test_data()
  null_fit <- glm(y ~ 1, data = data, family = binomial())

  expect_test_error(
    MissingBVS:::checkforprior.betas.glm(
      "invalid", NULL, TRUE, n = nrow(data), p = 1, p0 = 1,
      y = data$y, glmnull = null_fit, laplace = 0
    ),
    "Only BF approximations"
  )
  expect_test_error(
    MissingBVS:::checkforprior.betas.glm(
      "TBF", "invalid", FALSE, n = nrow(data), p = 1, p0 = 1,
      y = data$y, glmnull = null_fit, laplace = 0
    ),
    "Prior.betas must be one of"
  )
})
