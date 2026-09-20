for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  for (BF.method in c("BIC", "TBF")) {
    test_that(sprintf(
      "missingBVS.glm returns the full posterior for %s with %s",
      encode_family(family_obj), BF.method
    ), {
      data <- glm_family_data(fam_string)

      capture.output({
        result <- suppressWarnings(missingBVS.glm(
          y ~ x1 + x2,
          data = data,
          family = family_obj,
          BF.method = BF.method,
          n.keep = "all"
        ))
      })

      expect_s3_class(result, "MissingBvs")
      expect_identical(result$method, "Full")
      expect_identical(result$variables, c("x1", "x2"))
      expect_equal(result$p, 2)
      expect_equal(nrow(result$modelsprob), 4)
      expect_equal(sum(result$modelsprob[, "Post"]), 1, tolerance = 1e-6)
      expect_s3_class(result$glmnull, "glm")
      expect_s3_class(result$glmfull, "glm")
      expect_identical(result$family$family, family_obj$family)
      expect_identical(result$BF.method, BF.method)
      expect_length(result$inclprob, 2)
      expect_length(result$postprobdim, 3)
      expect_true(is.finite(result$C))
      if (BF.method == "TBF") {
        expect_identical(result$prior.betas, "gZellner")
      }
    })
  }
}

for (fam_string in names(glm_bas_families())) {
  family_obj <- glm_bas_family(fam_string)

  test_that(sprintf(
    "missingBVS.glm supports the gprior method for %s",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string)

    capture.output({
      result <- suppressWarnings(missingBVS.glm(
        y ~ x1 + x2,
        data = data,
        family = family_obj,
        BF.method = "gprior",
        n.keep = "all"
      ))
    })

    expect_s3_class(result, "MissingBvs")
    expect_identical(result$method, "Full")
    expect_identical(result$family$family, family_obj$family)
    expect_identical(result$BF.method, "gprior")
    expect_identical(result$prior.betas, "gZellner")
    expect_equal(sum(result$modelsprob[, "Post"]), 1, tolerance = 1e-6)
    expect_true(is.finite(result$C))
  })
}

for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  test_that(sprintf(
    "missingBVS.glm uses supplied imputations for %s with BIC",
    encode_family(family_obj)
  ), {
    data <- glm_family_data(fam_string, missing = TRUE)
    imputations <- glm_test_imputations(data)

    capture.output({
      result <- suppressWarnings(missingBVS.glm(
        y ~ x1 + x2,
        data = data,
        family = family_obj,
        BF.method = "BIC",
        n.keep = "all",
        imp.datasets = imputations
      ))
    })

    expect_s3_class(result, "MissingBvs")
    expect_identical(result$method, "Full")
    expect_true("compress.imp.array" %in% names(result))
    expect_equal(result$imp.info$n.imp, 1)
    expect_equal(sum(result$modelsprob[, "Post"]), 1, tolerance = 1e-6)
    expect_identical(result$family$family, family_obj$family)
  })
}

test_that("missingBVS.glm rejects gprior for non-BAS families", {
  data <- glm_test_data()

  expect_test_error(
    missingBVS.glm(y ~ x1 + x2, data = data, family = gaussian(), BF.method = "gprior"),
    "family not implemented in BAS"
  )
})

test_that("missingBVS.glm reports a mismatched response", {
  data <- glm_test_data()

  expect_test_error(
    missingBVS.glm(y ~ x1, data = data, family = binomial(), null.model = z ~ 1),
    "response in the full and null model does not coincide"
  )
})