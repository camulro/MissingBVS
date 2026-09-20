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

glm_families <- function() {
  list(
    "binomial" = binomial(),
    "gaussian" = gaussian(),
    "Gamma" = Gamma(),
    "inverse.gaussian" = inverse.gaussian(),
    "poisson" = poisson(),
    "quasi" = quasi(),
    "quasibinomial" = quasibinomial(),
    "quasipoisson" = quasipoisson()
  )
}

glm_bas_families <- function() {
  c("binomial" = NA, "poisson" = NA, "Gamma" = NA)
}

glm_bas_family <- function(family) {
  switch(family,
    "binomial" = binomial(link = "logit"),
    "poisson" = poisson(link = "log"),
    "Gamma" = Gamma(link = "log"),
    stop("unknown BAS family '", family, "'", call. = FALSE)
  )
}

glm_family_data <- function(family, missing = FALSE) {
  n <- 30
  angle <- seq(0, 2 * pi, length.out = n)
  x1p <- seq(-1.5, 1.5, length.out = n)
  x2 <- sin(angle)
  family_name <- switch(
    family,
    "binomial" = "binomial",
    "quasibinomial" = "binomial",
    "quasipoisson" = "poisson",
    "Gamma" = "inverse",
    "inverse.gaussian" = "inverse",
    "quasi" = "gaussian",
    family
  )

  y <- switch(
    family_name,
    "binomial" = as.integer(seq_len(n) %% 4 %in% c(1, 2)),
    "poisson" = as.integer(round(2 + 0.8 * x1p - 0.2 * x2)),
    "inverse" = 2 + 0.5 * x1p^2 + 0.3 * abs(x2),
    "gaussian" = 0.5 + 0.8 * x1p - 0.4 * x2 + 0.25 * sin(angle),
    stop("unknown family '", family, "'", call. = FALSE)
  )

  data <- data.frame(y = y, x1 = x1p, x2 = x2)
  if (missing) data$x1[c(5, 16)] <- NA_real_
  data
}

encode_family <- function(family) {
  paste(family$family, family$link, sep = "/")
}
