for (i in seq_along(glm_families())) {
  fam_string <- names(glm_families())[i]
  family_obj <- glm_families()[[i]]

  for (BF.method in c("BIC", "TBF")) {
    test_that(sprintf(
      "missingGibbsBVS.glm returns Gibbs summaries for %s with %s",
      encode_family(family_obj), BF.method
    ), {
      data <- glm_family_data(fam_string)

      capture.output({
        result <- suppressWarnings(missingGibbsBVS.glm(
          y ~ x1 + x2,
          data = data,
          family = family_obj,
          BF.method = BF.method,
          init.model = "Null",
          n.iter = 12,
          n.burnin = 3,
          n.thin = 2,
          Gibbs.seed = 123
        ))
      })

      expect_s3_class(result, "MissingBvs")
      expect_identical(result$method, "Gibbs")
      expect_equal(result$p, 2)
      expect_equal(nrow(result$modelslogBF), 6)
      expect_equal(ncol(result$modelslogBF), 3)
      expect_length(result$inclprob, 2)
      expect_length(result$inclprobRB, 2)
      expect_length(result$postprobdim, 3)
      expect_true(all(is.finite(result$modelslogBF)))
      expect_true(is.finite(result$C))
      expect_s3_class(result$glmnull, "glm")
      expect_s3_class(result$glmfull, "glm")
      expect_identical(result$family$family, family_obj$family)
      expect_identical(result$BF.method, BF.method)
      if (BF.method == "TBF") {
        expect_identical(result$prior.betas, "gZellner")
      }
    })
  }
}

test_that("missingGibbsBVS.glm validates the initial model", {
  data <- glm_test_data()

  expect_test_error(
    suppressWarnings(missingGibbsBVS.glm(
      y ~ x1 + x2,
      data = data,
      family = binomial(),
      BF.method = "BIC",
      init.model = c(1, 0, 1),
      n.iter = 2,
      n.burnin = 0,
      Gibbs.seed = 123
    )),
    "Initial model with incorrect length"
  )
})