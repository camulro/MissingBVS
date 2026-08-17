test_that("BF.BIC.lm matches the BIC expression", {
  data <- lm_test_data()
  X <- lm_test_complete_design(data)[, c("(Intercept)", "x1")]
  null_fit <- lm(y ~ 1, data = data)
  model_fit <- lm.fit(x = X, y = data$y)
  n <- nrow(data)
  k <- 1
  SS0 <- crossprod(null_fit$residuals)
  SSE <- crossprod(model_fit$residuals)

  result <- MissingBVS:::BF.BIC.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    n = n,
    k = k,
    p0 = 1
  )
  expected <- n / 2 * log(SS0 / SSE) - k / 2 * log(n)

  expect_equal(as.numeric(result), as.numeric(expected))
})

test_that("BF.TBF.lm passes the model dimension and deviance to its method", {
  data <- lm_test_data()
  X <- lm_test_complete_design(data)[, c("(Intercept)", "x1")]
  null_fit <- lm(y ~ 1, data = data)
  SSE <- crossprod(lm.fit(x = X, y = data$y)$residuals)
  SS0 <- crossprod(null_fit$residuals)
  received <- list()
  lTBF.method <- function(k, dev) {
    received <<- list(k = k, dev = dev)
    as.numeric(k + dev)
  }

  result <- MissingBVS:::BF.TBF.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    lTBF.method = lTBF.method,
    n = nrow(data),
    k = 1,
    p0 = 1
  )

  expected_dev <- nrow(data) * log(as.numeric(SSE / SS0))
  expect_identical(received$k, 1)
  expect_equal(as.numeric(received$dev), expected_dev)
  expect_equal(as.numeric(result), 1 + expected_dev)
})

test_that("g-prior and FLS Bayes factors return zero for the null model", {
  data <- lm_test_data()
  X <- lm_test_complete_design(data)[, "(Intercept)", drop = FALSE]
  X_model <- lm_test_complete_design(data)[, c("(Intercept)", "x1")]
  SS0 <- crossprod(lm(y ~ 1, data = data)$residuals)

  gprior <- MissingBVS:::BF.gprior.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    prior.betas = "gBF",
    n = nrow(data),
    k = 0,
    p0 = 1
  )
  fls <- MissingBVS:::BF.FLS.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    dmax = 2,
    n = nrow(data),
    k = 0,
    p0 = 1
  )

  expect_equal(as.numeric(gprior), 0)
  expect_equal(as.numeric(fls), 0)

  expect_true(is.finite(MissingBVS:::BF.gprior.lm(
    y = data$y,
    X = X_model,
    SS0 = SS0,
    prior.betas = "gBF",
    n = nrow(data),
    k = 1,
    p0 = 1
  )))
  expect_true(is.finite(MissingBVS:::BF.FLS.lm(
    y = data$y,
    X = X_model,
    SS0 = SS0,
    dmax = 2,
    n = nrow(data),
    k = 1,
    p0 = 1
  )))
})

test_that("checkforprior.betas.lm selects the requested methods", {
  data <- lm_test_data()
  X <- lm_test_complete_design(data)[, c("(Intercept)", "x1")]
  SS0 <- crossprod(lm(y ~ 1, data = data)$residuals)
  n <- nrow(data)

  bic <- MissingBVS:::checkforprior.betas.lm(
    "BIC", NULL, n, p = 1, p0 = 1, y = data$y, SS0 = SS0
  )
  direct_bic <- MissingBVS:::BF.BIC.lm(
    y = data$y, X = X, SS0 = SS0, n = n, k = 1, p0 = 1
  )
  expect_equal(bic(k = 1, X = X), direct_bic)

  gprior <- MissingBVS:::checkforprior.betas.lm(
    "gprior", "gZellner", n, p = 1, p0 = 1, y = data$y, SS0 = SS0
  )
  expect_length(gprior(k = 1, X = X), 1)
  expect_true(is.finite(gprior(k = 1, X = X)))
})

test_that("checkforprior.betas.lm rejects unsupported choices", {
  expect_test_error(
    MissingBVS:::checkforprior.betas.lm(
      "invalid", NULL, n = 10, p = 1, p0 = 1, y = 1:10, SS0 = 1
    ),
    "Only BF approximations"
  )
  expect_test_error(
    MissingBVS:::checkforprior.betas.lm(
      "TBF", "invalid", n = 10, p = 1, p0 = 1, y = 1:10, SS0 = 1
    ),
    "Prior.betas must be one of"
  )
})
