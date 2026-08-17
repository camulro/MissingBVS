test_that("missingBVS.lm returns the full posterior for complete data", {
  data <- lm_test_data()

  output <- capture.output({
    result <- missingBVS.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "BIC",
      n.keep = "all"
    )
  })

  expect_gt(length(output), 0)
  expect_s3_class(result, "MissingBvs")
  expect_identical(result$method, "Full")
  expect_identical(result$variables, c("x1", "x2"))
  expect_equal(result$p, 2)
  expect_equal(nrow(result$modelsprob), 4)
  expect_equal(sum(result$modelsprob[, "Post"]), 1)
  expected_inclprob <- colSums(result$modelsprob[, c("x1", "x2")] *
    result$modelsprob[, "Post"])
  expect_equal(unname(result$inclprob), unname(expected_inclprob))
  expect_s3_class(result$lmnull, "lm")
  expect_s3_class(result$lmfull, "lm")
})

test_that("missingBVS.lm uses supplied imputations", {
  data <- lm_test_data(missing = TRUE)
  imputations <- lm_test_imputations(data)

  capture.output({
    result <- missingBVS.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "BIC",
      n.keep = "all",
      imp.datasets = imputations
    )
  })

  expect_s3_class(result, "MissingBvs")
  expect_identical(result$method, "Full")
  expect_true("compress.imp.array" %in% names(result))
  expect_equal(result$imp.info$n.imp, 1)
  expect_equal(nrow(result$modelsprob), 4)
  expect_equal(sum(result$modelsprob[, "Post"]), 1)
})

test_that("missingBVS.lm reports a mismatched response", {
  data <- lm_test_data()

  expect_test_error(
    missingBVS.lm(y ~ x1, data = data, null.model = z ~ 1),
    "response in the full and null model does not coincide"
  )
})
