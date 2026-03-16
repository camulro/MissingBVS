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
