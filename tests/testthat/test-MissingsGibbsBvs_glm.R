test_that("missingGibbsBVS.glm returns Gibbs summaries", {
  data <- glm_test_data()

  capture.output({
    result <- suppressWarnings(missingGibbsBVS.glm(
      y ~ x1 + x2,
      data = data,
      family = binomial(),
      BF.method = "BIC",
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
  expect_equal(length(result$postprobdim), 3)
  expect_true(all(is.finite(result$modelslogBF)))
  expect_s3_class(result$glmnull, "glm")
  expect_s3_class(result$glmfull, "glm")
})

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
