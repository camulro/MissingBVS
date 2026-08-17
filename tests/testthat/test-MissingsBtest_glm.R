test_that("missingBtest.glm computes posterior probabilities", {
  data <- glm_test_data()
  models <- list(
    null = y ~ 1,
    x1 = y ~ x1,
    full = y ~ x1 + x2
  )

  capture.output({
    result <- missingBtest.glm(
      data = data,
      models = models,
      family = binomial(),
      BF.method = "BIC"
    )
  })

  expect_s3_class(result, "MissingBtest")
  expect_identical(result$nullmodel, "null")
  expect_identical(names(result$lBFi0), c("null.to.null", "x1.to.null", "full.to.null"))
  expect_identical(result$lBFi0[["null.to.null"]], 0)
  expect_equal(sum(result$PostProbi), 1)
  expect_length(result$modelspool, 3)
  expect_true(all(vapply(result$modelspool, inherits, logical(1), "glm")))
  expect_identical(result$family$family, "binomial")
})

test_that("missingBtest.glm accepts external imputations", {
  data <- glm_test_data(missing = TRUE)
  imputations <- glm_test_imputations(data)
  models <- list(
    null = y ~ 1,
    x1 = y ~ x1,
    full = y ~ x1 + x2
  )

  capture.output({
    result <- missingBtest.glm(
      data = data,
      models = models,
      family = binomial(),
      BF.method = "BIC",
      imp.datasets = imputations
    )
  })

  expect_s3_class(result, "MissingBtest")
  expect_true("compress.imp.array" %in% names(result))
  expect_equal(result$imp.info$n.imp, 1)
  expect_equal(sum(result$PostProbi), 1)
})

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
