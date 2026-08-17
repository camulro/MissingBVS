lm_test_data <- function(missing = FALSE) {
  n <- 24
  angle <- seq(0, 2 * pi, length.out = n)
  x1 <- seq(-1, 1, length.out = n)
  x2 <- cos(angle)
  y <- 1 + 0.8 * x1 - 0.4 * x2 + sin(angle) / 10

  data <- data.frame(y = y, x1 = x1, x2 = x2)
  if (missing) data$x1[c(4, 11)] <- NA_real_
  data
}

lm_test_imputations <- function(data) {
  imputed <- data
  imputed$x1[is.na(imputed$x1)] <- mean(data$x1, na.rm = TRUE)

  result <- array(
    rep(as.matrix(imputed), 1),
    dim = c(nrow(imputed), ncol(imputed), 1),
    dimnames = list(NULL, names(imputed), "1")
  )
  result
}

lm_test_complete_design <- function(data = lm_test_data()) {
  cbind(`(Intercept)` = 1, x1 = data$x1, x2 = data$x2)
}

expect_test_error <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  testthat::expect_s3_class(error, "error")
  testthat::expect_match(conditionMessage(error), pattern)
}

glm_test_data <- function(missing = FALSE) {
  n <- 30
  angle <- seq(0, 2 * pi, length.out = n)
  x1 <- seq(-1.5, 1.5, length.out = n)
  x2 <- sin(angle)
  y <- as.integer(seq_len(n) %% 4 %in% c(1, 2))

  data <- data.frame(y = y, x1 = x1, x2 = x2)
  if (missing) data$x1[c(5, 16)] <- NA_real_
  data
}

glm_test_imputations <- function(data) {
  imputed <- data
  imputed$x1[is.na(imputed$x1)] <- mean(data$x1, na.rm = TRUE)

  array(
    rep(as.matrix(imputed), 1),
    dim = c(nrow(imputed), ncol(imputed), 1),
    dimnames = list(NULL, names(imputed), "1")
  )
}
