#' Logarithm of the Scott and Berger model prior
#'
#' Computes the logarithm of the prior that assigns an uniform distribution to
#' model size (Scott and Berger, 2010)
#'
#' It is derived by assigning an independent Bernouilli distribution with
#' probability Beta(1,1) to the a priori inclusion probability of each
#' competing variable.
#'
#' @export
#' @param p Number of covariates to select from.
#' @param model Binary vector of length \code{p} specifying a model.
#'
#' @return \code{logScottBerger} returns the logarithm of the Scott and Berger
#' prior for a model of size of the given \code{model}.
#'
#' @author Carolina Mulet
#'
#' @seealso \code{\link[MissingBVS]{MissingGibbsBvs.lm}},
#' \code{\link[MissingBVS]{logConstant}}, \code{\link[MissingBVS]{logUser}}
#'
#' @references Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes
#' multiplicity adjustment in the variable-selection problem. The Annals of
#' Statistics. 38: 2587–2619.
#'
#' @examples
#' logScottBerger(7, c(rep(1,4),rep(0,3))) # log((p + 1)^(-1) / choose(p, sum(model)))
#'
logScottBerger <- function(p, model) {
  # log((p + 1)^(-1) / choose(p, sum(model)))
  -log(p + 1) - lchoose(p, sum(model)) #logaritmic scale
}

#' Logarithm of the Constant model prior
#'
#' Computes the logarithm of the prior that assigns an uniform distribution to
#' the model space
#'
#' It is derived by assigning an independent Bernouilli distribution with fixed
#' probability 0.5 to the a priori inclusion probability of each
#' competing variable.
#'
#' @export
#' @param p Number of covariates to select from.
#'
#' @return \code{logConstant} returns the logarithm of the Constant prior for a
#' problem with \code{p} competing variables.
#'
#' @author Carolina Mulet
#'
#' @seealso \code{\link[MissingBVS]{MissingGibbsBvs.lm}},
#' \code{\link[MissingBVS]{logScottBerger}}, \code{\link[MissingBVS]{logUser}}
#'
#' @examples
#' logConstant(7) # log(1 / (2^p))
#'
logConstant <- function(p) {
  # 1 / (2^p)
  -p * log(2) #logaritmic scale
}

#' Logarithm of a User model prior
#'
#' Computes the logarithm of the prior that assigns a user given prior for model
#' size
#'
#' It is derived by assigning a fixed prior for model sizes given by user
#' through the numeric vector \code{priorprobs}. Its first entry correspond to
#' the null model, the second to the ones with exactly one active variable, and
#' so on. Models of a given size receive the same amount of probability, such
#' as Scott and Berger (2010).
#'
#' @export
#' @param p Number of covariates to select from.
#' @param model Binary vector of length \code{p} specifying a model.
#' @param priorprobs Vector of length \code{p}+1 made up by the prior
#' probabilities for each model size: 0,...,\code{p}.
#'
#' @return \code{logUser} returns the logarithm of a given User prior for a
#' problem with \code{p} competing variables and prior model sizes given by
#' \code{priorprobs}.
#'
#' @author Carolina Mulet
#'
#' @seealso \code{\link[MissingBVS]{MissingGibbsBvs.lm}},
#' \code{\link[MissingBVS]{logScottBerger}}, \code{\link[MissingBVS]{logConstant}}
#'
#' @references Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes
#' multiplicity adjustment in the variable-selection problem. The Annals of
#' Statistics. 38: 2587–2619.
#'
#' @examples
#' logUser(7, c(rep(1,4),rep(0,3)), c(0.3,rep(0.2,3),rep(0.1,4)))
#' # log(priorprobs[sum(model)]/sum(priorprobs) / choose(p, sum(model)))
#'
logUser <- function(p, model, priorprobs) {
  # priorprobs[sum(model)]/sum(priorprobs) / choose(p, sum(model))
  log(priorprobs[sum(model)]) - log(sum(priorprobs)) - lchoose(p, sum(model)) #logaritmic scale
  #prior prob for each model size divided by the number of models with that size
}

#auxiliar function for lprior.model.dummies computation
#returns lchoose(n, k) if k < n -1 and 1 otherwise (k <= n)
mylchoose <- function(n, k) {
  ifelse(k < n - 1, lchoose(n, k), 1)
}

#auxiliar function to calculate the number of models with r dummies active
lG <- function(r, ltau) {
  m2 <- length(ltau) #number of factors active

  mat.ind <- matrix(1:ltau[1], nc = 1)
  if(m2 > 1) {
    for (i in 2:m2) {
      mat.ind <- merge(mat.ind, 1:ltau[i], by = NULL)
    }
  }

  ind.r <- which(rowSums(mat.ind) == r)
  Gtau.r <- mat.ind[ind.r,]

  if (length(ind.r) > 1) {
    res <- 0
    for (j in 1:dim(Gtau.r)[1]) {
      res <- res + exp(sum(mylchoose(ltau + 1, as.numeric(Gtau.r[j,]))))
    }
  } else res <- exp(sum(mylchoose(ltau + 1, as.numeric(Gtau.r))))
  log(res)
}

#l is a vector of length L with the number of levels - 1 for each factor,
#where L is the number of factors
#delta is a vector of length L with the number of active dummies for each factor
#tau is a vector of length L with 1 if the factor is active
logScottBerger.d <- function(delta, tau, l) {
  if (sum(tau) == 0) {
    return(0)
  } else {
    ltau <- l[tau] #levels - 1 of active factors
    -lG(sum(delta), ltau) - log(sum(ltau) - sum(tau) + 1) #logarithmic scale
  }
}

logConstant.d <- function(tau, l) {
  if (sum(tau) == 0) {
    return(0)
  } else -sum(l) * log(2) #logaritmic scale of 1 / (2^sum(l))
}
