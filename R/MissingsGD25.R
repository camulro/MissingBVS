#' Bayesian Variable Selection with Missing data for linear regression models
#'
#' Computation and summaries of posterior distribution over the model space
#' in problems of small to moderate size, when missingness occurs in linear
#' models with normally distributed covariates under the g'-imputation prior of
#' García-Donato et al. (2025).
#'
#' The set of competing models is made up by all the possible subsets of regressors
#' specified by \code{formula}: Mi for i in 1,...,2^p, being p the number of potential
#' regressors in the variable selection problem. It is assumed that the intercept term
#' is present in all models, and the simplest one, the null, contains only this term.
#' \code{\link[MissingBVS]{missingGD25}} performs for the regressors \code{n.imp}
#' imputations from an objective Bayesian multivariate normal model (Sun and Berger, 2006),
#' along with their variance-covariance matrices. Hence, the posterior distribution
#' over the model space is given through Bayes' theorem:
#'
#' Pr(Mi | \code{data}) = Pr(Mi) * Bi / C,
#'
#' where Bi is the Bayes factor (BF) of Mi to M0 derived a the expected value of g'BF of
#' García-Donato et al. (2025) for normal random regressors, Pr(Mi) is the prior
#' probability of Mi and C is the normalizing constant. The g'BFs are computed
#' with \code{\link[MissingBVS]{BF.miss.X}}, which resemble g-prior BFs (Zellner, 1986)
#' for random covariates, for each pair of imputed dataset and  variance-covariance
#' matrix simulated. The integral of g'BF is approximated with a Monte Carlo scheme.
#'
#' The prior over the model space Pr(Mi) offers three options throuh \code{prior.models}:
#' -"Constant" assigns the same prior probability to every model.
#' -"ScottBerger" is the default choice. It assigns the same prior probability to
#' every possible model size and, therefore, accounts for multiplicity issues
#' (Scott and Berger, 2010).
#' -"User": if chosen, user has to provide a p+1 dimensional vector with the model size
#' prior probabilities through \code{priorprobs}. The first component must contain the
#' probability of the null model M0 and next p components correspond to the p prior
#' probabilities of model sizes 1,...,p.
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
#' @param n.keep It can be either the character "all" to return the whole model
#' space or a numeric for the exact number of the most probable models to keep.
#' By default it is set to 10 and automatically adjusted if 10 is greater than
#' the total number of models.
#' @param n.imp Number of imputed datasets for model posterior computation.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGD25} returns an object of class \code{missingBVS}
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
#' \item{modelsprob}{A \code{n.keep} x (p+1) matrix which summarizes the \code{n.keep}
#' most probable a posteriori models and their associated probability}
#' \item{inclprob}{Named vector with the inclusion probabilities of p competing variables}
#' \item{postprobdim}{Posterior probabilities over the true model size}
#' \item{C}{The value of the normalizing constant C=sum_i AvBi * Pr(Mi)}
#' \item{call}{The \code{call} to the function}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{imp.info}{List of arguments used for the imputation step and other information}
#' \item{compress.imp.list}{Compressed list of imputed datasets and covariance
#' matrices}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{prior.models}{Argument chosen for \code{prior.models}}
#' \item{method}{String "Full" denoting exhaustive model search}
#'
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGibbsGD25}} for a heuristic
#' approximation based on Gibbs sampling (recommended when p>20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238> Optimal
#' Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' Sun, D. and Berger, J. O. (2006) Objective Bayesian analysis for the multivariate
#' normal model. In Bernardo, J. M., Bayarri, M. J., Berger, J. O., Dawid, A. P.,
#' Heckerman, D., Smith, A. F. M., and West, M. (eds.), Proc. Valencia / ISBA 8th
#' World Meeting on Bayesian statistics. Oxford university Press. MR2433206. 16
#'
#' @examples
#' \dontrun{
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' #Here we keep the 8 competing models:
#' f <- Ozone ~ 1 + Wind + Temp + Solar.R
#' airq.mBVS <- missingGD25(formula = f, data = airquality, n.keep = 8)
#'
#' #Show the results:
#' airq.mBVS
#'
#' #Summ up the results:
#' summary(airq.mBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(airq.mBVS)
#' }
#'
missingGD25 <- function (formula,
                         data,
                         prior.models = "ScottBerger",
                         priorprobs = NULL,
                         n.keep = 10,
                         n.imp = 039E1,
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
    stop("This method is only for continuous covariates.\nTry missingBvs.lm instead.\n")
  }

  #Full design matrix
  framefull <- model.frame(formula, data, na.action = NULL)
  X.full <- framefull[,-1] #remove intercept
  namesx <- dimnames(X.full)[[2]]
  p <- length(namesx) #Number of covariates to select from

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(p, 1, "(Intercept)", c("(Intercept)", namesx), n.keep, p)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, p)

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
  #function to compute log(BFa0) for a given model with García-Donato's 2025 method
  lBF.method <- function (model) lBF.miss(model,
                                          imputation.list = imputation.list,
                                          BF.miss.aux = BF.miss.aux,
                                          n = n, nMC = n.imp)

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
  cat("Most complex model has a total of", p + 1, "covariates.\n")
  cat("From those 1 is fixed (the intercept) and we should select from the remaining",
      p, ":\n")
  cat(paste(namesx, collapse = ", ", sep = ""))

  cat("\nThe problem has a total of", 2^p, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #progress bar for loop
  pb <- txtProgressBar(min = 0, max = 2^p, style = 3, width = 50, char = "=")

  #Posterior computation
  all.models.lPM <- matrix(0, nr = 2^p, nc = p+1) #last column contains log(BF_a0*Pr(M))
  for (i in seq_len(2^p-1)){
    setTxtProgressBar(pb, i)

    #transform the number of the model into a binary number
    current.model <- BayesVarSel:::integer.base.b_C(i, p)
    all.models.lPM[i, seq_len(p)] <- current.model

    all.models.lPM[i, p+1] <- lBF.method(model = which(current.model == 1)) +
      lprior.models(current.model) #log(BF_a0*Pr(M))
  }
  setTxtProgressBar(pb, 2^p)
  #null model
  all.models.lPM[2^p, seq_len(p)] <- numeric(p)
  all.models.lPM[2^p, p+1] <- lprior.models(numeric(p)) #BF = 1 for null model

  #renormalize
  C <- sum(exp(all.models.lPM[, p+1]))
  all.models.PM <- all.models.lPM
  all.models.PM[, p+1] <- exp(all.models.lPM[, p+1] - log(C))
  colnames(all.models.PM) <- c(namesx, "Post")

  #Summ up the posterior distribution
  summ.posterior.list <- summ.posterior(all.models.PM, p, p, FALSE, NULL)
  list2env(summ.posterior.list, envir = environment())

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

  ##result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$lmfull <- lmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, lmfull is the lm object for the full model
  result$lmnull <- lmnull # The lm object for the null model (omits NAs)

  result$variables <- namesx #The name of the competing variables
  result$n <- n #number of observations
  result$p <- p #number of competing variables
  result$k <- 1 # intercept #number of fixed covariates
  result$HPMbin <- hpm #The binary code for the HPM model
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- namesx

  #The binary code for the n.keep best models (after n.thin is applied) and the correspondent post
  result$modelsprob <- all.models.PM[order(all.models.PM[,p+1],
                                           decreasing = TRUE)[seq_len(n.keep)],]

  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- namesx

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:p + 1 #dimension of the true model
  result$C <- C #normalizing constant

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
  names(result$priorprobs) <- 0:p + 1 #prior dimension probability

  #arguments used for imputation
  result$imp.info <- list(n.imp = n.imp, imp.seed = imp.seed)
  #save the imputed datasets for BMA or sensitivity analysis
  raw.imp.list <- serialize(imputation.list, NULL)
  result$compress.imp.list <- memCompress(raw.imp.list, type = "xz")

  result$logprior.models <- lprior.models #function used for model prior
  result$prior.models <- prior.models

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}
