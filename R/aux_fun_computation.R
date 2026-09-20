
#' Build design matrix (fixed and non-fixed), matrix to use -if needed- for
#' imputation step and factor objects -if needed- to use for computation.
#'
#' @keywords internal
buildmatrices <- function(formula, null.model, data, marginal.factors) {

  ##Null model (fixed vars, including the intercept term)
  #for imputation:
  framenull <- model.frame(null.model, data = data, na.action = NULL)
  #for computation when no missing data present
  X0rdf <- X0 <- model.matrix(framenull, data)
  namesnull <- colnames(X0)

  q0 <- ncol(framenull) #covariate-factor level
  p0 <- ncol(X0) #covariate-dummy level

  ##Full model
  framefull <- model.frame(formula, data = data, na.action = NULL)

  ##Factors' treatment: rank deficient decomposition
  if (marginal.factors) {

    X0rdf <- model.matrix.rankdef(framenull) #to remove vars from full
    X.full <- model.matrix.rankdef(framefull) #rank defficient full design

    namesx <- colnames(X.full)
    namesxnotnull <- setdiff(namesx, colnames(X0rdf)) #non-fixed var names

  } else {
    X.full <- model.matrix(formula, framefull)

    namesx <- colnames(X.full)
    namesxnotnull <- setdiff(namesx, namesnull)
  }
  #Full design matrix without fixed vars
  X.full <- X.full[, namesxnotnull, drop = FALSE]
  p <- length(namesxnotnull) #covariate-dummy level

  ordvars <- c(namesnull, namesxnotnull) #[X0, X.full]

  if (!marginal.factors) { #If no factors' treatment

    return(list(
      q0 = p0,
      p0 = p0,
      X0 = X0,
      namesnull = namesnull,
      framenull = framenull,
      q = p,
      p = p,
      X.full = X.full,
      namesxnotnull = namesxnotnull,
      namesx = namesx,
      framefull = framefull,
      ordvars = ordvars,
      depvars = namesxnotnull,
      L = 0
    ))
  }

  ##Get covariates and factors to select from if needed
  #Info from framefull to get factors and dummy names:
  depvars <- setdiff(attr(terms(framefull), "term.labels"),
                     attr(terms(framenull), "term.labels"))

  factorsfull <- attr(terms(framefull), "factors")
  classesfull <- attr(terms(framefull), "dataClasses")

  #positions has number of rows equal to the number of regressors and p columns.
  #A 1 in a row denotes the position in X of a regressor (several positions for
  #the dummies of a factor).
  positions <- t(sapply(depvars,
    function(var) {

      ##OLD: (does not select interactions of factors)
      # if(is.factor(data[[var]])) {
      #   levs <- levels(data[[var]])
      #   ind <- which(namesxnotnull %in% paste0(var,levs)) #1 if the namelevel matches
      # } else ind <- which(namesxnotnull == var) #1 if the name matches

      if (classesfull[var] %in% c("factor", "ordered")) {

        #get vars corresponding to that factor (to include interactions)
        terms_var <- names(which(factorsfull[, var] > 0))
        levels_var <- lapply(terms_var,
          function(x) paste0(x, levels(data[[x]]))
        )
        #get all possible combinations for factor and levels
        gridlevs <- expand.grid(levels_var, stringsAsFactors = FALSE)

        #for interaction terms
        levelnames <- apply(gridlevs, 1, paste0, collapse = ":")

        #select dummies corresponding to each factor
        ind <- which(namesxnotnull %in% levelnames)

      } else ind <- which(namesxnotnull == var)

      #build position row for these ind
      pos <- numeric(p)
      pos[ind] <- 1
      pos
    }
  ))

  colnames(positions) <- namesxnotnull

  #Identify factors
  tmp <- colSums(positions %*% t(positions))

  positionsx <- tmp == 1 #TRUE if numeric covariate
  L <- sum(!positionsx) #Number of factors to select from

  ##Factor information needed for model prior and posterior computation
  if (L > 0) {

    #1 if the dummy column variable belongs to the row factor
    positionsfac <- positions[!positionsx, , drop = FALSE]

    l <- tmp[tmp > 1] #number of levels for each factor

    #save position of the first dummy in each factor
    #used on pool estimation step to get representant over repeated models
    indf <- apply(positionsfac, 1,
      function(x) which(x == 1)[1]
    )

    #save hash for representant over saturated models
    satmodels.repr <- sapply(seq_len(L),
      function(j) {
        fac <- positionsfac[j, ]
        fac[indf[j]] <- 0
        digest::digest(fac)
      }
    )

    q <- p - sum(l) + L #Total number of factors and covariates

  } else { #if there are no factors

    positionsfac <- NULL
    l <- integer(0)
    indf <- integer(0)
    satmodels.repr <- character(0)

    q <- p
  }

  #return:
  list(
    q0 = q0,
    p0 = p0,
    X0 = X0,
    namesnull = namesnull,
    framenull = framenull,

    q = q,
    p = p,
    X.full = X.full,
    namesxnotnull = namesxnotnull,
    namesx = namesx,

    framefull = framefull,
    ordvars = ordvars,
    depvars = depvars,

    positions = positions,
    positionsx = positionsx,
    positionsfac = positionsfac,
    L = L,
    l = l,
    indf = indf,
    satmodels.repr = satmodels.repr
  )
}

#' Builds the matrix for the whole model space and adds an aditional column
#' to associated probabilites
#'
#' @keywords internal
buildmodelsmatrix <- function(q) {
  models.ord <- c(seq_len(2^q-1),0)
  zeros <- numeric(2^q)

  tmodels.mat <- sapply(models.ord,
                        function(j) num2bin.model(j, q, NULL)["bin",])

  cbind(t(tmodels.mat), zeros)
}

#' Original code from package 2.6.0 \pkg{BayesVarSel} (distributed under GPL-2),
#' by Gonzalo García-Donato and Anabel Forte.
#'
#' Adapted by Carolina Mulet to fit \pkg{MissingBVS}'s code and return a matrix
#' with the binary expression of a model and the active variables with NA
#' to reduce computational burden
#'
#' @keywords internal
num2bin.model <- function(x, p, NAvars) {
  #x is the number to get its binary expression, p is the number of variables,
  #namesxnotnull is the name of the p competing vars and NAvars, the ones with NAs.
  #If NAvars = NULL, the row "bin" equals integer.base.b_C

  if (x == 0) {
    res <- numeric(p)
  } else {
    ndigits <- (floor(logb(x, base = 2)) + 1)
    res <- numeric(ndigits)
    for (i in 1:ndigits) {
      res[i] <- (x %% 2)
      x <- (x %/% 2)
    }

    res <- c(res, numeric(p - ndigits)) #variables active
  }

  resNA <- res * NAvars #variables in model with NA

  matrix(c(res, resNA), byrow = T, nrow = 2, ncol = p,
         dimnames = list(c("bin","NA"), names(NAvars)))
}

#' Original code from package 3.1-0 \pkg{lmerTest} (distributed under GPL-2, GPL-3),
#' by Alexandra Kuznetsova, Per Bruun Brockhoff and Rune Haubo Bojesen Christensen.
#'
#' Adapted by Carolina Mulet to fit \pkg{MissingBVS}'s code and obtain just
#' rank defficient matrices.
#'
#' @keywords internal
model.matrix.rankdef <- function (model.frame.aux) {
  #internal function to create rank defficient matrices from a given dataframe
  #created from a model.frame call

  if (ncol(model.frame.aux) == 1) { #just the response
    Xnull.def <- cbind(`(Intercept)` = rep.int(1, nrow(model.frame.aux)))
    return(Xnull.def)
  }

  terms <- attr(terms(model.frame.aux), "term.labels")

  Xi.rdef <- sapply(terms, function(var) {

    f <- as.formula(paste0("~ 0 + ", var))
    #without intercept produces one columns per level
    model.matrix(f, data = model.frame.aux)

  }, simplify = FALSE)

  Xfull.def <- do.call(cbind, Xi.rdef)

  #Always include intercept termn:
  Xfull.def <- cbind(`(Intercept)` = rep.int(1, nrow(Xfull.def)), Xfull.def)
  return(Xfull.def)
}

#' Stable version of log(sum(exp(vector of logarithms))).
#' To avoid infite Bayes factors at the average step.
#'
#' If maximum is not infinite, it returns:
#'  log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x)))),
#' a more stable computational expression for the sum in logarithmic scale.
#'
#' @keywords internal
logsumexp.stable <- function (x) {
  max.x <- max(x)
  if (!is.finite(max.x)) {

    warning("A Bayes factor in infinite.\n", immediate. = TRUE)
    return(max.x)

  } else logsum <- max.x + log(sum(exp(x - max.x)))

  return(logsum)
}

#' Get computational time from an starting event
#'
#' @keywords internal
get.time <- function (start.time) {

  t <- difftime(Sys.time(), start.time, units = "secs")

  as.numeric(t)
}

#' Show computational time
#'
#' @keywords internal
format_time <- function(x) {
  if (is.na(x)) return("--:--:--")

  x <- round(x)
  h <- x %/% 3600
  m <- (x %% 3600) %/% 60
  s <- x %% 60

  sprintf("%02d:%02d:%02d", h, m, s)
}

#' Personalized progress bar. If update == FALSE, just comutes models rate.
#'
#' @keywords internal
update.progress <- function(i, total, start.time, width = 20) {

  #current computation time
  elapsed <- get.time(start.time)
  #models rate computed
  mrate <- if (elapsed > 0) i / elapsed else 0
  #update each ceiling(mrate / 4) models (approx 4 times per second)
  update <- i %% max(1, ceiling(mrate / 4)) == 0 || i == total

  if (update) { #check whether or not to update progress bar
    #estimated remaining time
    eta <- if (mrate > 0) (total - i) / mrate else NA
    #progress percent
    percent <- i / total * 100

    #build bar
    numcomplete <- floor(width * i / total)
    bar <- paste0(
      strrep("=", max(0, numcomplete - 1)),
      if (i < total) ">" else "=",
      strrep("-", max(0, width - numcomplete))
    )

    #write bar and info
    cat(
      sprintf(
        "\r\033[K[%s] %3.0f%% [Elapsed time: %s | Estim. time: %s | ≈ %.0f models/s]",
        bar,
        percent,
        format_time(elapsed),
        format_time(eta),
        mrate
      )
    )

    flush.console()
  }

}

#' !%in% function
#'
#' @keywords internal
"%notin%" <- function(x, table) {

  match(x, table, nomatch = 0) == 0

}
