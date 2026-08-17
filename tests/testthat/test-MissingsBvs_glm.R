test_that("missingBVS.glm returns the full posterior for complete data", {
  data <- glm_test_data()

  capture.output({
    result <- missingBVS.glm(
      y ~ x1 + x2,
      data = data,
      family = binomial(),
      BF.method = "BIC",
      n.keep = "all"
    )
  })

  expect_s3_class(result, "MissingBvs")
  expect_identical(result$method, "Full")
  expect_identical(result$variables, c("x1", "x2"))
  expect_equal(result$p, 2)
  expect_equal(nrow(result$modelsprob), 4)
  expect_equal(sum(result$modelsprob[, "Post"]), 1)
  expect_s3_class(result$glmnull, "glm")
  expect_s3_class(result$glmfull, "glm")
  expect_identical(result$family$family, "binomial")
})

test_that("missingBVS.glm uses supplied imputations", {
  data <- glm_test_data(missing = TRUE)
  imputations <- glm_test_imputations(data)

  capture.output({
    result <- missingBVS.glm(
      y ~ x1 + x2,
      data = data,
      family = binomial(),
      BF.method = "BIC",
      n.keep = "all",
      imp.datasets = imputations
    )
  })

  expect_s3_class(result, "MissingBvs")
  expect_identical(result$method, "Full")
  expect_true("compress.imp.array" %in% names(result))
  expect_equal(result$imp.info$n.imp, 1)
  expect_equal(sum(result$modelsprob[, "Post"]), 1)
})

test_that("missingBVS.glm reports a mismatched response", {
  data <- glm_test_data()

  expect_test_error(
    missingBVS.glm(y ~ x1, data = data, family = binomial(), null.model = z ~ 1),
    "response in the full and null model does not coincide"
  )
})
