for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  for (BF.method in c("BIC", "TBF")) {
    test_that(sprintf(
      "missingBtest.glm computes posterior probabilities for %s with %s",
      encode_family(family_obj), BF.method
    ), {
      data <- glm_family_data(fam_string)
      models <- list(
        null = y ~ 1,
        x1 = y ~ x1,
        full = y ~ x1 + x2
      )

      capture.output({
        result <- suppressWarnings(missingBtest.glm(
          data = data,
          models = models,
          family = family_obj,
          BF.method = BF.method
        ))
      })

      expect_s3_class(result, "MissingBtest")
      expect_identical(result$nullmodel, "null")
      expect_identical(names(result$lBFi0), c("null.to.null", "x1.to.null", "full.to.null"))
      expect_equal(result$lBFi0[["null.to.null"]], 0, tolerance = 1e-6)
      expect_true(all(is.finite(result$lBFi0)))
      expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
      expect_length(result$modelspool, 3)
      expect_true(all(vapply(result$modelspool, inherits, logical(1), "glm")))
      expect_identical(result$family$family, family_obj$family)
      expect_identical(result$BF.method, BF.method)
      if (BF.method == "TBF") {
        expect_identical(result$prior.betas, "gZellner")
      }
    })
  }
}

for (fam_string in names(glm_bas_families())) {
  family_obj <- glm_bas_family(fam_string)

  test_that(sprintf(
    "missingBtest.glm supports the gprior method for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)
    models <- list(
      null = y ~ 1,
      x1 = y ~ x1,
      full = y ~ x1 + x2
    )

    capture.output({
      result <- suppressWarnings(missingBtest.glm(
        data = data,
        models = models,
        family = family_obj,
        BF.method = "gprior"
      ))
    })

    expect_s3_class(result, "MissingBtest")
    expect_identical(result$nullmodel, "null")
    expect_equal(result$lBFi0[["null.to.null"]], 0, tolerance = 1e-6)
    expect_true(all(is.finite(result$lBFi0)))
    expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
    expect_identical(result$family$family, family_obj$family)
    expect_identical(result$BF.method, "gprior")
    expect_identical(result$prior.betas, "gZellner")
  })
}

for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  test_that(sprintf(
    "missingBtest.glm accepts external imputations for %s with BIC",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string, missing = TRUE)
    imputations <- glm_test_imputations(data)
    models <- list(
      null = y ~ 1,
      x1 = y ~ x1,
      full = y ~ x1 + x2
    )

    capture.output({
      result <- suppressWarnings(missingBtest.glm(
        data = data,
        models = models,
        family = family_obj,
        BF.method = "BIC",
        imp.datasets = imputations
      ))
    })

    expect_s3_class(result, "MissingBtest")
    expect_true("compress.imp.array" %in% names(result))
    expect_equal(result$imp.info$n.imp, 1)
    expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
    expect_length(result$modelspool, 3)
    expect_identical(result$family$family, family_obj$family)
  })
}

test_that("missingBtest.glm rejects an unknown null model", {
  data <- glm_test_data()
  models <- list(null = y ~ 1, x1 = y ~ x1)

  expect_test_error(
    missingBtest.glm(
      data = data,
      models = models,
      family = binomial(),
      null.model = "unknown"
    ),
    "null model provided is not in the list"
  )
})