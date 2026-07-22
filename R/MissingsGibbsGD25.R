#' Bayesian Variable Selection with Missing data for linear regression models
#' using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior model distribution using a
#' Gibbs sampling algorithm to explore the model space. Each posterior model probability
#' is computed under the g'-imputation prior of García-Donato et al. (2025) in
#' linear models with normally distributed covariates.
#'
#' Gibbs sampling search algorithm to avoid exhaustive enumeration of model space
#' when it is unfeasible. It draws from the model posterior distribution and uses
#' frequency of "visits" to construct the estimates. The algorithm was originally
#' proposed by  George and McCulloch (1997). Later, Garcia-Donato and Martinez-Beneito (2013)
#' shown that the sampling strategy in combination with estimates based on frequency of
#' visits provides very reliable results.
#'
#' \code{\link[MissingBVS]{missingGibbsGD25}} is a heuristic approximation of
#' \code{\link[MissingBVS]{missingGD25}}. See the later for common details.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param prior.models Model prior distribution over the covariates and/or factors
#' model space (to be literally specified). Possible choices are "Constant",
#' "ScottBerger" and "User" (see details).
#' @param priorprobs A p+1 (being p the number of competing variables) dimensional
#' vector defining the prior model size probabilities (if \code{prior.models}= "User";
#' see details).
#' @param init.model The model at which the simulation process starts.
#' It can be either a string: "Null" for \code{null.model}, "Full" for \code{formula}
#' and "Random" for a randomly selected model; or a vector with p (the number of factors
#' and/or covariates to select from) zeros and ones defining a model.
#' @param n.iter The total number of iterations performed after burn in.
#' @param n.burnin Number of iterations to discard at the beginning.
#' @param n.thin Positive integer that states the number of models to discard before one
#' is saved. Default is 1, larger values are suggested if needed less memory and computation
#' but they can reduce accuracy because estimates are based on fewer simulations.
#' @param n.imp Number of imputed datasets for model posterior computation.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGibbsGD25} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{Time lasted solving the problem}
#' \item{lmfull}{If missings on the \code{formula} competing variables, combination
#' of the estimates of fitted full model over the \code{n.imp} imputed datasets.
#' Otherwise, it is the \code{\link[stats]{lm}} object}
#' \item{lmnull}{The \code{lm} class object that results when \code{null.model}
#' is fitted by \code{\link[stats]{lm}}}
#' \item{variables}{Names of all the competing variables given by \code{formula}}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory covariates to select from}
#' \item{k}{Number of fixed variables, which is fixed to 1}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model}
#' \item{modelslogBF}{A floor(\code{n.iter}/\code{n.thin}) x (p+1) matrix which summarizes
#' the floor(\code{n.iter}/\code{n.thin}) visited models and their associated log-Bayes factors}
#' \item{inclprob}{Named vector with the inclusion probabilities of p competing variables}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Estimated posterior probabilities over the true model size}
#' \item{C}{The value of the estimated normalizing constant}
#' \item{postprobs}{Estimated posterior probability}
#' \item{call}{The \code{call} to the function}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{imp.info}{List of arguments used for the imputation step and other information}
#' \item{compress.imp.array}{Compressed array of imputed datasets}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{prior.models}{Argument chosen for \code{prior.models}}
#' \item{method}{String "Gibbs" denoting Gibbs sampling model search}
#'
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGD25}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' George E. and McCulloch R. (1997) Approaches for Bayesian variable
#' selection. Statistica Sinica, 7, 339:372.
#'
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238> Optimal
#' Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
#' @examples
#' \dontrun{
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' #Here we keep the 8 competing models:
#' f <- Ozone ~ 1 + Wind + Temp + Solar.R
#' airq.mGBVS <- missingGibbsGD25(formula = f, data = airquality)
#'
#' #Show the results:
#' airq.mGBVS
#'
#' #Summ up the results:
#' summary(airq.mGBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(airq.mGBVS)
#' }
#'
missingGibbsGD25 <- function (formula,
                              data,
                              prior.models = "ScottBerger",
                              priorprobs = NULL,
                              init.model = "Full",
                              n.iter = 10000,
                              n.burnin = 500,
                              n.thin = 1,
                              n.imp = 039E1,
                              Gibbs.seed = runif(1,0,26061970),
                              imp.seed = runif(1,0,09011975)) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(paste(formula[[2]], " ~ 1", sep=""))

  #Check for numeric covariates
  aux <- model.frame(formula, data)
  isnum <- sapply(aux, is.numeric)
  # isint <- sapply(aux, is.integer)
  # if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
  if (sum(isnum) < dim(aux)[2]) {
    stop("This method is only for continuous covariates.\nTry missingGibbsBVS.lm instead.\n")
  }

  #Full design matrix
  framefull <- model.frame(formula, data, na.action = NULL)
  X.full <- framefull[,-1] #remove intercept
  namesx <- dimnames(X.full)[[2]]
  p <- length(namesx) #Number of covariates to select from

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, p)

  #Check arguments and compute init.model
  init.model <- checkGibbsarguments(p, 1, "(Intercept)", c("(Intercept)", namesx),
                                    init.model, FALSE, NULL, NULL, NULL, NULL)

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #check for missings
  NAvars <- checkformissings(y = framefull[,1], X.full = X.full[obsnotNA,])

  #BF function
  BF.miss.aux <- function (X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
                                                           y = y, SS0 = SS0, n = n, k)

  #Imputation of missing data
  if (p*n > 10000 | n.imp > 039E1) cat("Imputation step could take a while.\n",
    "Consider reducing the number of imputed datasets if that is the case.\n")

  cat("Performing imputation of missing data with Garcia-Donato's 2025 method.\n",
      "Please wait . . . \n")
  imputation.list <- MC.imputation(X = X.full, nMC = n.imp, seed = imp.seed)

  #remove observations with missings on the response
  imputation.list$rX.imput <- imputation.list$rX.imput[obsnotNA, , , drop = FALSE]
  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model with García-Donato's 2025 method
    lBF.method <- function (model) lBF.miss(model,
                                            imputation.list = imputation.list,
                                            BF.miss.aux = BF.miss.aux,
                                            n = n, nMC = n.imp)
  } else lBF.method <- function (model) BF.miss.aux(X.center = imputation.list$rX.imput[,model,],
                                                    Sigma11 = imputation.list$rSigma[model, model,],
                                                    k = length(model))

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", p + 1, "single covariates.\n")
  cat("From those 1 is fixed (the intercept) and we should select from the remaining",
      p, ":\n")
  cat(paste(namesx, collapse = ", ", sep = ""))

  cat("\nThe problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  set.seed(Gibbs.seed)
  gibbs.list <- GM97.Gibbs(1, X.full, p, namesx, namesx, lprior.models, lBF.method, BF.miss.aux,
                           FALSE, NULL, init.model, n.iter, n.burnin, n.thin)
  list2env(gibbs.list, envir = environment())

  summ.Gibbs.list <- summ.Gibbs(cf.models.lBF, all.lBF.PM, inclprobRB, p, n.iter)
  list2env(summ.Gibbs.list, envir = environment())

  if (!is.null(NAvars)) {#Pool results for imputed datasets
    #Evaluate lm of full model with missings using Rubin's rule
    fit <- list(); mt <- attr(framefull, "terms")
    for (i in 1:n.imp) {
      z <- lm.fit(x = cbind(1, imputation.list$rX.imput[,,i]), y = y)
      z$terms <- mt; class(z) <- "lm"; fit[[i]] <- z
    }
    lmfull <- mice::pool(fit)
    lmfull$call <- NULL #otherwise, Rstudio returns a warning trying to read lmfull$call
  } else lmfull <- lm(formula, data, x = TRUE, y = TRUE)

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the programm to finish
  result$lmfull <- lmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, lmfull is the lm object for the full model
  result$lmnull <- lmnull #The lm object for the null model (omits NAs)

  result$variables <- namesx #The name of the competing variables
  result$n <- n #number of observations
  result$p <- p #number of competing variables
  result$k <- 1 #number of fixed covariates
  result$HPMbin <- hpm #The binary code for the HPM model and its BF.PM
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- namesx

  #The binary code for all the visited models (after n.thin is applied) and the correspondent post
  result$modelslogBF <- cf.models.lBF

  result$inclprob <- inclprob #inclusion probability for each variable
  result$inclprobRB <- inclprobRB[n.iter, ] #Rao-Blackwellized inclusion probability
  names(result$inclprobRB) <- depvars

  result$postprobdim <- probdim #vector with the estimated posterior dimension probability
  names(result$postprobdim) <- 0:p + 1 #dimension of the true model
  result$C <- C #estimated normilizing constant
  #Estimation of posterior probabilities based on C
  result$postprobs <- post

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- numeric(p+1)
    priorprobs[1] <- exp(lprior.models(numeric(p))) #prior inclusion prob for dimension 0
    for (i in seq_len(p)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep.int(1, i), rep.int(0, p - i))) + lchoose(p, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:p + 1 #dimension prior probability

  #arguments used for imputation
  result$imp.info <- list(n.imp = n.imp, imp.seed = imp.seed)
  #save the imputed datasets for sensitivity analysis
  raw.imp.list <- serialize(imputation.list, NULL)
  result$compress.imp.list <- memCompress(raw.imp.list, type = "xz")

  result$logprior.models <- lprior.models #function used for model prior
  result$prior.models <- prior.models

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}
