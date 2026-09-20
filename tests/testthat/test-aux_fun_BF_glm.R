test_that("lTBF.gfixed is zero for the null model", {
  g <- 30
  devnull <- rnorm(1)
  expect_equal(as.numeric(MissingBVS:::lTBF.gfixed(g, k = 0, dev = devnull, devnull = devnull)), 0)
})

for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  test_that(sprintf(
    "glm.fit BIC and TBF Bayes factors match their expressions for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)
    n <- nrow(data)
    X <- model.matrix(~ x1, data = data)
    null_fit <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj))
    model_fit <- suppressWarnings(glm.fit(x = X, y = data$y, family = family_obj))
    k <- 1

    bic <- suppressWarnings(MissingBVS:::BF.BIC.glm.fit(
      y = data$y,
      X = X,
      family = family_obj,
      devnull = null_fit$deviance,
      n = n,
      k = k
    ))
    expected_bic <- (null_fit$deviance - model_fit$deviance - k * log(n)) / 2

    expect_equal(as.numeric(bic), as.numeric(expected_bic), tolerance = 1e-6)

    lTBF.method <- function(k, dev) k + dev
    tbf <- suppressWarnings(MissingBVS:::BF.TBF.glm.fit(
      y = data$y,
      X = X,
      family = family_obj,
      n = n,
      k = k,
      lTBF.method = lTBF.method
    ))

    expect_equal(as.numeric(tbf), k + as.numeric(model_fit$deviance), tolerance = 1e-6)
  })
}

for (fam_string in names(glm_bas_families())) {
  family_obj <- glm_bas_family(fam_string)

  test_that(sprintf(
    "BAS glm Bayes factors are finite and respect log-marginal shifts for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)
    n <- nrow(data)
    X <- model.matrix(~ x1, data = data)
    null_fit <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj))

    bic <- function(logmargnull) {
      MissingBVS:::BF.BIC.glm(
        y = data$y,
        X = X,
        family = family_obj,
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
        family = family_obj,
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
        family = family_obj,
        prior.betas = BAS::g.prior(g = n),
        logmargnull = logmargnull,
        k = 1,
        p0 = 1
      )
    }

    for (method in list(bic, tbf, gprior)) {
      at_zero <- suppressWarnings(method(0))
      shifted <- suppressWarnings(method(2))
      expect_length(at_zero, 1)
      expect_true(is.finite(at_zero))
      expect_equal(as.numeric(at_zero - shifted), 2)
    }
  })
}

test_that("checkforfamily accepts the BAS families for gprior", {
  for (fam_string in names(glm_bas_families())) {
    expect_error(MissingBVS:::checkforfamily(glm_bas_family(fam_string), "gprior"), NA)
  }
})

test_that("checkforfamily rejects gprior for non-BAS families", {
  for (fam_string in setdiff(names(glm_families()), names(glm_bas_families()))) {
    family_obj <- glm_families()[[fam_string]]
    expect_test_error(
      MissingBVS:::checkforfamily(family_obj, "gprior"),
      "family not implemented in BAS"
    )
  }
})

test_that("checkforfamily rejects gprior for non-default links", {
  bad_families <- list(
    binomial(link = "probit"),
    binomial(link = "cloglog"),
    poisson(link = "identity"),
    poisson(link = "sqrt"),
    Gamma(link = "inverse"),
    Gamma(link = "identity")
  )
  for (family_obj in bad_families) {
    expect_test_error(
      MissingBVS:::checkforfamily(family_obj, "gprior"),
      "family not implemented in BAS"
    )
  }
})

test_that("checkforfamily allows BIC and TBF for every family", {
  for (family_obj in glm_families()) {
    expect_error(MissingBVS:::checkforfamily(family_obj, "BIC"), NA)
    expect_error(MissingBVS:::checkforfamily(family_obj, "TBF"), NA)
  }
})

for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  test_that(sprintf(
    "checkforprior.betas.glm selects working BIC and TBF methods for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)
    n <- nrow(data)
    X <- model.matrix(~ x1, data = data)
    null_fit <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj, y = TRUE, x = TRUE))
    y <- as.numeric(null_fit$y)

    bic <- MissingBVS:::checkforprior.betas.glm(
      "BIC", NULL, n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    )
    tbf <- MissingBVS:::checkforprior.betas.glm(
      "TBF", "gZellner", n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    )

    expect_true(is.finite(suppressWarnings(bic(k = 1, X = X, fitstart = NULL))))
    expect_true(is.finite(suppressWarnings(tbf(k = 1, X = X, fitstart = NULL))))
  })
}

for (fam_string in names(glm_bas_families())) {
  family_obj <- glm_bas_family(fam_string)

  test_that(sprintf(
    "checkforprior.betas.glm supports the gprior method for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)
    n <- nrow(data)
    X <- model.matrix(~ x1, data = data)
    null_fit <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj, y = TRUE, x = TRUE))
    y <- as.numeric(null_fit$y)

    gprior <- MissingBVS:::checkforprior.betas.glm(
      "gprior", "gZellner", n, p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    )

    expect_true(is.finite(suppressWarnings(gprior(k = 1, X = X))))
  })
}

test_that("checkforprior.betas.glm rejects unsupported choices", {
  data <- glm_family_data("poisson")
  null_fit <- suppressWarnings(glm(y ~ 1, data = data, family = poisson()))
  y <- as.numeric(null_fit$y)

  expect_test_error(
    MissingBVS:::checkforprior.betas.glm(
      "invalid", NULL, n = nrow(data), p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    ),
    "Only BF approximations"
  )
  expect_test_error(
    MissingBVS:::checkforprior.betas.glm(
      "TBF", "invalid", n = nrow(data), p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    ),
    "Prior.betas must be one of"
  )
  expect_test_error(
    MissingBVS:::checkforprior.betas.glm(
      "gprior", "invalid", n = nrow(data), p = 1, p0 = 1, y = y, glmnull = null_fit, laplace = 0
    ),
    "Prior.betas must be one of"
  )
})

for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  test_that(sprintf(
    "lBF.av.glm.fit equals the single-imputation BIC Bayes factor for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string, missing = TRUE)
    x1imp <- data$x1
    x1imp[is.na(x1imp)] <- mean(x1imp, na.rm = TRUE)
    design <- cbind(`(Intercept)` = 1, x1 = x1imp, x2 = data$x2)
    arr <- array(design, dim = c(nrow(design), 3, 1),
                 dimnames = list(NULL, colnames(design), "1"))
    y <- as.numeric(data$y)
    n <- nrow(data)
    glmnull <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj, y = TRUE, x = TRUE))
    lBF_fun <- MissingBVS:::checkforprior.betas.glm(
      "BIC", NULL, n, p = 2, p0 = 1, y = y, glmnull = glmnull, laplace = 0
    )

    av <- suppressWarnings(MissingBVS:::lBF.av.glm.fit(
      model = 1,
      imputation.array = arr,
      lBF = lBF_fun,
      p0 = 1,
      n.imp = 1,
      y = y,
      glmnull = glmnull
    ))

    fit <- suppressWarnings(glm.fit(
      y = y,
      x = arr[, 1:2, 1],
      family = glmnull$family,
      weights = glmnull$prior.weights,
      offset = glmnull$offset,
      control = glmnull$control
    ))
    expected <- suppressWarnings(lBF_fun(k = 1, X = arr[, 1:2, 1], fitstart = fit$coefficients))

    expect_equal(as.numeric(av), as.numeric(expected), tolerance = 1e-6)
  })
}

for (fam_string in names(glm_bas_families())) {
  family_obj <- glm_bas_family(fam_string)

  test_that(sprintf(
    "lBF.av equals the single-imputation gprior Bayes factor for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string, missing = TRUE)
    x1imp <- data$x1
    x1imp[is.na(x1imp)] <- mean(x1imp, na.rm = TRUE)
    design <- cbind(`(Intercept)` = 1, x1 = x1imp, x2 = data$x2)
    arr <- array(design, dim = c(nrow(design), 3, 1),
                 dimnames = list(NULL, colnames(design), "1"))
    y <- as.numeric(data$y)
    n <- nrow(data)
    glmnull <- suppressWarnings(glm(y ~ 1, data = data, family = family_obj, y = TRUE, x = TRUE))
    lBF_fun <- MissingBVS:::checkforprior.betas.glm(
      "gprior", "gZellner", n, p = 2, p0 = 1, y = y, glmnull = glmnull, laplace = 0
    )

    av <- suppressWarnings(MissingBVS:::lBF.av(
      model = 1,
      imputation.array = arr,
      lBF = lBF_fun,
      p0 = 1,
      n.imp = 1
    ))
    expected <- suppressWarnings(lBF_fun(k = 1, X = arr[, 1:2, 1]))

    expect_equal(as.numeric(av), as.numeric(expected), tolerance = 1e-6)
  })
}