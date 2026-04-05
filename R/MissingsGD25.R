#' Bayesian Variable Selection with Missing data for linear regression models
#'
#' Computation and summaries of posterior distribution over the model space
#' in problems of small to moderate size when missingness occurs in linear
#' models with normally distributed covariates.
#'
#' The set of competing models is made up by all the possible subsets of
#' regressors specified by \code{formula}: Mi for i in 1,...,2^p, being p the number
#' of potential (non-fixed) covariates in the variable selection problem. The simplest,
#' nested in all of them, contains only the intercept. \code{MissingGD25} performs
#' \code{n.imp} imputations given by the result of \code{\link[MissingBVS]{MC.imputation}}
#' and computes the posterior distribution over this model space through Bayes' theorem:
#'
#' Pr(Mi | \code{data})=Pr(Mi)*Bi/C,
#'
#' where Bi is the Bayes factor(BF) of Mi to M0 under missing data proposed by
#' García-Donato et al (2025), Pr(Mi) is the prior probability of Mi and C is
#' the normalizing constant.
#'
#' Bi is computed as the MonteCarlo approximation of the integral defining the
#' BF, for Mi to M0 calculated for the jth imputed data set.
#'
#' The prior over the model space Pr(Mi) offers three possibilities:
#' "Constant" assigns the same prior probability to every model.
#' "ScottBerger" is the default choice. It assigns the same prior probability to
#' every possible model dimension and, therefore, accounts for multiplicity issues
#' (Scott and Berger 2010).
#' "User" (see below).
#'
#' If \code{prior.models}="User" is chosen, user has to provide a p+1 dimensional
#' parameter vector with the model dimension prior probabilities through \code{priorprobs}.
#' The first component of \code{priorprobs} must contain the probability of the
#' model with fixed covariates; next p components correspond to the p prior probabilities
#' of the possible model dimensions.
#'
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param prior.models Prior distribution over the model space (to be literally specified).
#' Possible choices are "Constant", "ScottBerger" and "User" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed covariates)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details).
#' @param n.keep Number of the most probable models kept. By default it is set to
#' 10 and automatically adjusted if 10 is greater than the total number of models.
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 30} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param initialimp.mice.method Method for mice's imputation.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGD25} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{lmfull}{Object of class \code{\link[mice]{mipo}} that combines the estimates
#' for the model defined by \code{formula} fitted by \code{\link[stats]{lm}} over
#' the \code{n.imp} imputed datasets. See \code{\link[mice]{pool}} for details}
#' \item{lmnull}{The \code{lm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{lm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model}
#' \item{modelsprob}{A (n.keep)x(p+1) \code{matrix} which summaries the \code{n.keep}
#' most probable a posteriori models and their associated probability}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables.}
#' \item{postprobdim}{Posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the
#' model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{method}{\code{Full}}
#'
#' @author Carolina Mulet and Gonzalo Garcia-Donato
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
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @keywords package
#'
#' @examples #To be completed
#'
missingGD25 <- function (formula,
                         data,
                         prior.models = "ScottBerger",
                         priorprobs = NULL, #needed if prior.models = User
                         n.keep = 10,
                         imp.time.test = TRUE,
                         initialimp.mice.method = "pmm", #mice's default
                         n.imp = 039E1, #number of imputed datasets for BF
                         imp.seed = runif(1,0,09011975)) { #seed for the imputation

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(paste(formula[[2]], " ~ 1", sep=""))

  #Check for numeric covariates
  aux <- model.frame(formula, data)
  isnum <- sapply(aux, is.numeric)
  isint <- sapply(aux, is.integer)
  if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
    stop("This method is only for continuous covariates.\n")
  }
  cat("Be careful, this method is only for normally distributed covariates.\n",
      "Do you want to continue? (y/n)\n")
  if (tolower(readline()) != "y") {
    stop("Try the missingBVS.lm function instead.\n")
  }

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #Full design matrix
  framefull <- model.frame(formula, data, na.action = NULL)
  X.full <- framefull[,-1] #remove intercept
  namesx <- dimnames(X.full)[[2]]
  p <- length(namesx) #Number of covariates to select from

  #Is there any variable to select from?
  if (p == 0) { #only the intercept can be fixed
    stop(paste0("The number of fixed covariates is equal to the number of\n",
                "covariates in the full model. No model selection can be done.\n"))
  }

  #check if the number of regressors is too big.
  if (p > 20) {
    warning("Number of covariates too big. . . consider using missingGibbsBvs.lm.\n",
            immediate. = TRUE)
  }

  #n.keep > 2^p, the number of models?
  if (n.keep > 2^p) {
    cat(paste0("The number of models to keep (", n.keep,
               ") is larger than the total number of models (",
               2^p, ") and it has been set to ", 2^p))
    n.keep <- 2^p
  }

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings
  checkformissings.lm(y = framefull[,1], X.full = X.full)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, p)

  #Check methods and options
  BF.miss.aux <- function (X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
                                                           y = y, SS0 = SS0,
                                                           n = n, k)

  #Imputation of missing data
  if (imp.time.test & (p*n > 10000 | n.imp > 039E1)) {
    #test imputation time
    cat("Time test . . . \n")
    time.test <- MC.imputation(X = X.full,
                               time.test = TRUE)

    estim.time <- time.test * n.imp / (60 * 30) #30 imputed datasets used to test
    cat("The whole imputation would take ", estim.time,
        "minutes (approx.) to run.\n Do you want to continue? (y/n)\n")
    if (tolower(readline()) != "y") {
       stop("Reduce the number of imputed datasets.\n")
    }
  }
  cat("Performing imputation of missing data with Garcia-Donato's 2025 method.\n",
      "Please wait . . . \n")
  imputation.list <- MC.imputation(X = X.full,
                                   nMC = n.imp,
                                   seed = imp.seed,
                                   initialimp.mice.method = initialimp.mice.method)

  #remove observations with missings on the response
  imputation.list$rX.imput <- imputation.list$rX.imput[obsnotNA,,]
  #function to compute log(BFa0) for a given model with García-Donato's 2025 method
  lBF.method <- function (model) lBF.miss(model,
                                          imputation.list = imputation.list,
                                          BF.miss.aux = BF.miss.aux,
                                          n = n, nMC = n.imp)

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", p + 1, "single covariates.\n")
  cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
             p, ":\n"))
  cat(paste(paste(namesx, collapse = ", ", sep = ""), "\n", sep = ""))

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #progress bar for loop
  pb <- txtProgressBar(min = 0,
                       max = 2^p,
                       style = 3,
                       width = 50,
                       char = "=")

  #Posterior computation
  all.models.lPM <- matrix(0, nr = 2^p, nc = p+1) #last column contains log(BF_a0*Pr(M))
  for (i in seq_len(2^p-1)){
    setTxtProgressBar(pb, i)

    #transform the number of the model into a binary number
    current.model <- BayesVarSel:::integer.base.b_C(i, p)
    all.models.lPM[i, seq_len(p)] <- current.model

    lBF.PM <- lBF.method(model = which(current.model == 1)) +
      lprior.models(current.model) #log(BF_a0*Pr(M))

    all.models.lPM[i, p+1] <- lBF.PM
  }
  setTxtProgressBar(pb, 2^p)
  #null model
  all.models.lPM[2^p, seq_len(p)] <- rep(0, p)
  all.models.lPM[2^p, p+1] <- lprior.models(rep(0, p)) #BF = 1 for null model

  #renormalize
  C <- sum(exp(all.models.lPM[, p+1]))
  all.models.PM <- all.models.lPM
  all.models.PM[, p+1] <- exp(all.models.lPM[, p+1] - log(C))
  colnames(all.models.PM) <- c(namesx, "Post")

  inclprob <- rep(0, p)
  probdim <- rep(0, p + 1)
  #compute inclusion probabilities (except for fixed variables) and
  #posterior probability of the dimension of the true model
  for (i in seq_len(2^p)) {
    inclprob[which(all.models.PM[i, seq_len(p)] == 1)] <-
      inclprob[which(all.models.PM[i, seq_len(p)] == 1)] + all.models.PM[i, p + 1]
    probdim[sum(all.models.PM[i, seq_len(p)]) + 1] <-
      probdim[sum(all.models.PM[i, seq_len(p)]) + 1] + all.models.PM[i, p + 1]
  }

  #HPM
  nPmax <- which.max(all.models.PM[, p+1])
  hpm <- all.models.PM[nPmax, ]

  #MPM
  mpm <- rep(0,p)
  mpm[which(inclprob >= 0.5)] <- 1

  #Evaluate lm of full model with missings using Rubin's rule
  fit <- list()
  mt <- attr(framefull, "terms")
  for (i in 1:n.imp) {
    z <- lm.fit(x = cbind(1, imputation.list$rX.imput[,,i]), y = y)
    z$terms <- mt
    class(z) <- "lm"

    fit[[i]] <- z
  }
  lmfull <- mice::pool(fit)

  ##result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$lmfull <- lmfull # Object of class mipo combining the estimates for the
  # n.imp imputed datasets for the fitted full model
  result$lmnull <- lmnull # The lm object for the null model (omits NAs)

  result$variables <- namesx #The name of the competing variables
  result$n <- n #number of observations
  result$p <- p #number of competing variables
  result$k <- 1 # intercept #number of fixed covariates
  result$HPMbin <- hpm #The binary code for the HPM model
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- namesx

  result$modelsprob <- all.models.PM[order(all.models.PM[,p+1],
                                           decreasing = TRUE)[seq_len(n.keep)],]
  #The binary code for the n.keep best models (after n.thin is applied) and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- namesx

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:p + 1 #dimension of the true model

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, p + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, p))) #prior inclusion prob for dimension 0
    for (i in seq_len(p)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, p - i))) + lchoose(p, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:p + 1 #prior dimension probability

  result$C <- C #normilizing constant

  #arguments used for imputation
  result$imp.args <- list(initialimp.mice.method = initialimp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  #save the imputed datasets for sensitivity analysis
  # raw.imp.array <- serialize(imputation.array, NULL)
  # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$logprior.models <- lprior.models #function used for model prior

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}
