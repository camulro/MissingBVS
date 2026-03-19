#' Bayesian Variable Selection with Missing data for linear regression models
#' using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior distribution using a
#' Gibbs sampling algorithm to explore the model space and frequency of
#' "visits" to construct the estimates. It was originally proposed by George
#' and McCulloch (1997) and further studied by Garcia-Donato and
#' Martinez-Beneito (2013). They shown that the sampling strategy in
#' combination with estimates based on frequency of visits provides very
#' reliable results.
#'
#' \code{\link[MissingBVS]{MissingGibbsBvs.lm}} is a heuristic approximation of
#' \code{\link[MissingBVS]{MissingBvs.lm}}. See the later for common details.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param null.model Formula defining which is the simplest (null) model, which.
#' should be nested in the full one. By default, it is defined to be the one
#' with just the intercept.
#' @param data Data frame containing the data.
#' @param BF.approx.method Method used to approximate Bayes factors with missing
#' data (to be literally specified). Possible choices include "BIC", "TBF" and
#' "gprior" (see details).
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices are the ones implemented
#' in \pkg{BayesVarSel}: "Robust", "Liangetal", "gZellner", "ZellnerSiow",
#' "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified). Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed covariates)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details.)
#' @param init.model The model at which the simulation process starts. Options
#' include "Null" (the model only with the covariates specified in
#' \code{null.model}), "Full" (the model defined by \code{formula}), "Random" (a
#' randomly selected model) and a vector with p (the number of covariates to
#' select from) zeros and ones defining a model.
#' @param n.iter The total number of iterations performed after the burn in
#' process.
#' @param n.burnin Length of burn in, i.e. number of iterations to discard at
#' the beginning.
#' @param n.thin Positive integer defining the thinning rate. Default is 1.
#' Set 'n.thin' > 1 to save memory and computation time if 'n.iter' is large.
#' A large \code{n.thin} can reduce the accuracy because it, along with \code{n.iter},
#' sets the number of simulations used to construct the estimates.
#' @param parallelmice Logital to indicate whether or not to use parallel
#' \code{\link[mice]{mice}} imputation. If \code{NULL}, automatically performs
#' the parallel mice imputation if the number of imputations or the dataset are
#' big enough.
#' @param n.core See \code{\link[mice]{futuremice}} for details.
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 10} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>300}).
#' @param imp.mice.method Method for mice's imputation.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGibbsBVS.lm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{lmfull}{The \code{lm} class object that results when the model
#' defined by \code{formula} is fitted by \code{\link[stats]{lm}}}
#' \item{lmnull}{The \code{lm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{lm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model using
#' \code{inclprobRB}}
#' \item{modelsprob}{\code{data.frame} which summaries the \code{n.keep}
#' most probable a posteriori models and their associated Bayes factor in
#' logaritmic scale}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{method}{\code{Gibbs}}
#'
#' @author Carolina Mulet and Gonzalo Garcia-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
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
#' García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
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
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6(2): 461–464.
#'
#' Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
#' @keywords package
#'
#' @examples #To be completed
#'
missingGibbsBVS.lm <- function (formula,
                                null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""), #only the intercept for now
                                data,
                                BF.approx.method = "gprior",
                                prior.betas = "Robust", #if BF.approx.method = "gprior"
                                prior.models = "ScottBerger",
                                priorprobs = NULL, #needed if prior.models = User
                                init.model = "Full",
                                n.iter = 10000, #number of iterations for Gibbs Sampling algorithm
                                n.burnin = 500,
                                n.thin = 1,
                                parallelmice = NULL,
                                n.core = NULL,
                                imp.time.test = TRUE,
                                imp.mice.method = "pmm", #mice's default
                                n.imp = 3E2, #number of imputed datasets for BF
                                Gibbs.seed = runif(1,0,26061970), #seed for the Gibbs sampling
                                imp.seed = runif(1,0,09011975)) { #seed for the imputation

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #The response in the null model and in the full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #Response and fixed vars for imputation
  auxnull <- model.frame(null.model, data, na.action = NULL)
  namesnull.toimp <- dimnames(auxnull)[[2]][-1] #name of fixed variables to imputation

  #Missing model matrix of fixed covariates
  X0 <- model.matrix(null.model, auxnull)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of covariates to select from

  #Eval the full model
  lmfull <- lm(formula, data = data, y = TRUE, x = TRUE) #omits NA observations

  #Full design matrix for imputation
  auxfull <- model.frame(formula, data, na.action = NULL)
  namesx.toimp <- dimnames(auxfull)[[2]][-1] #name of variables to imputation
  namesxnotnull.toimp <- namesx.toimp[namesx.toimp %notin% namesnull.toimp]
  X.toimp <- data[,c(namesnull.toimp, namesxnotnull.toimp)] #design matrix with missing data with fixed cov

  #Model matrix data with missings
  X.full <- model.matrix(formula, auxfull)
  namesx <- dimnames(X.full)[[2]]
  namesxnotnull <- namesx[namesx %notin% namesnull]
  X.full <- X.full[, namesxnotnull]
  p <- dim(X.full)[2] #Number of covariates to select from

  #Is there any variable to select from?
  if (p == p0) {
    stop(paste0("The number of fixed covariates is equal to the number of\n",
                "covariates in the full model. No model selection can be done.\n"))
  }

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #check for missings and define covariates with NAs
  NAvars <- checkformissings.lm(y = auxnull[,1], X0, X.full, obsnotNA)

  if (p <= 20) {
    warning(paste0("The number of variables is small enough to visit every model.\n",
                   "Consider using MissingBvs.lm.\n"), immediate. = TRUE)
  }

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, p)

  #Check the initial model:
  if (is.character(init.model) == TRUE) {
    im <- substr(tolower(init.model), 1, 1)
    if (im %notin% c("n", "f", "r")) stop("Initial model not valid.\n")
    if (im == "n") init.model <- rep(0, p) #null model
    if (im == "f") init.model <- rep(1, p) #full model
    if (im == "r") init.model <- rbinom(n = p,
                                        size = 1,
                                        prob = .5)
  } else {
    init.model <- as.numeric(init.model > 0)
    if (length(init.model) != p) stop("Initial model with incorrect length.\n")
  }

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.lm(BF.approx.method,
                                             prior.betas,
                                             n, p, p0, y, SS0)

  #Imputation of missing data
  if (is.null(parallelmice)) {
    if (n.imp > 120 | n*p > 50000) {
      parallelmice <- TRUE #faster
    } else parallelmice <- FALSE
  }

  if (imp.time.test & (n*p > 10000 | n.imp > 3E2)) {
    #test imputation time
    cat("Time test . . . \n")
    time.test <- mice.imputation(X = X.toimp,
                                 formula,
                                 imp.mice.method = imp.mice.method,
                                 parallel = parallelmice,
                                 n.core = n.core,
                                 time.test = TRUE)

    estim.time <- time.test * n.imp / (60 * 30) #30 imputed datasets used to time
    cat("The whole imputation can take ", estim.time,
        "minutes (approx.) to run.\n Do you want to continue? (y/n)\n")
    if (tolower(readline()) != "y") {
      if (!parallelmice) {
        cat("Do you want to faster imputation running a parallel version of mice? (y/n)\n")
        if (tolower(readline()) == "y") {
          parallelmice <- TRUE
        } else stop("Reduce the number of imputed datasets.\n")
      } else stop("Reduce the number of imputed datasets.\n")
    }
  }

  cat("Performing imputation of missing data with mice's", imp.mice.method)
  if (parallelmice) cat(" parallel")
  cat(" method.\n", "Please wait . . . \n")

  imputation.array <-  mice.imputation(X = X.toimp,
                                       formula,
                                       n.imp = n.imp,
                                       imp.mice.method = imp.mice.method,
                                       seed = imp.seed,
                                       parallel = parallelmice,
                                       n.core = n.core)

  #remove observations with missings on the response
  imputation.array <- imputation.array[obsnotNA,,]
  #function to compute log(BFa0) for a given model as an average of BF computed
  #by BF.approx.method over the imputed datasets
  lBF.method <- function (model) lBF.approx(model,
                                            imputation.array = imputation.array,
                                            BF.approx.method = BF.approx.method,
                                            p0 = p0, n.imp = n.imp)

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", p + p0, "single covariates.\n")
  if (p0 == 1) {
    cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
               p, ":\n"))
  } else {
    cat(paste0("From those ", p0, " are fixed and we should select from the remaining ",
               p, ":\n"))
  }
  cat(paste(paste(namesxnotnull, collapse = ", ", sep = ""), "\n", sep = ""))

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  gibbs.list <- GM97.Gibbs(y, X0, X.full, p, namesxnotnull, NAvars,
                           lprior.models, lBF.method,
                           init.model, n.iter, n.burnin, n.thin, Gibbs.seed)

  all.models.lPM <- gibbs.list$all.models.lPM

  inclprob <- colMeans(all.models.lPM[,-(p+1)]) #inclusion probabilities except for fixed variables

  all.models.PM <- all.models.lPM[,seq_len(p)]
  all.models.PM <- cbind(all.models.PM, exp(all.models.lPM[, p + 1]))
  colnames(all.models.PM) <- c(namesxnotnull, "BF.PM")

  #Estimation of the normalizing constant:
  K <- round(dim(all.models.PM)[1]/2)
  Aset <- sample(x = 1:dim(all.models.PM)[1], size = K, replace = F)
  Bset <- (1:dim(all.models.PM)[1])[-Aset]
  #Bayes factors multiplied by prior probs of the models in A
  BF.PMAset <- all.models.PM[Aset, "BF.PM"]
  #Remove duplicates
  BF.PMAset <- unique(BF.PMAset)
  gAset <- sum(BF.PMAset)
  #How many of the models in Bset are in A?
  sumIA <- sum(all.models.PM[Bset,"BF.PM"] %in% all.models.PM[Aset,"BF.PM"])
  #Normalizing constant estimation
  C <- gAset*K / sumIA

  #compute estimated posterior probabilities
  all.models.PM[, p + 1] <- exp(all.models.lPM[, p + 1] - log(C))
  colnames(all.models.PM) <- c(namesxnotnull, "Post")

  probdim <- rep(0, p + 1)
  all.models.lBF <- all.models.lPM
  #compute posterior probability of the dimension of the true model and
  #save logBF for each model
  for (i in seq_len(floor(n.iter / n.thin))) {
    probdim[sum(all.models.PM[i, seq_len(p)])] <-
      probdim[sum(all.models.PM[i, seq_len(p)])] + all.models.PM[i, p + 1]

    all.models.lBF[i, p + 1] <- all.models.lPM[i, p + 1] -
      lprior.models(all.models.lPM[i, seq_len(p)]) # lBF
  }
  colnames(all.models.lBF) <- c(namesxnotnull, "logBF")

  #HPM
  nPmax <- which.max(all.models.lBF[, p+1])
  hpm <- all.models.lBF[nPmax, ]

  #MPM
  mpm <- rep(0,p)
  mpm[which(gibbs.list$inclprobRB[n.iter, ] >= 0.5)] <- 1

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the programm to finish
  result$lmfull <- lmfull #The lm object for the full model (without NAs)
  result$lmnull <- lmnull #The lm object for the null model

  result$variables <- namesxnotnull #The name of the competing variables
  result$n <- n #number of observations
  result$p <- p #number of competing variables
  result$k <- p0 #number of fixed covariates
  result$HPMbin <- hpm #The binary code for the HPM model and its BF.PM
  names(result$HPMbin) <- c(namesxnotnull, "logBF")
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- namesxnotnull

  result$modelsprob <- all.models.lBF

  #The binary code for all the visited models (after n.thin is applied) and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- namesxnotnull

  result$inclprobRB <- gibbs.list$inclprobRB #Rao-Blackwellized inclusion probability

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:p + p0 #dimension of the true model

  result$call <- match.call()
  result$C <- C #estimated normilizing constant

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, p + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, p))) #prior inclusion probability for dimension 0
    for (i in seq_len(p)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, p - i))) + lchoose(p, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:p + p0 #dimension prior probability

  #Estimation of posterior probabilities based on C
  result$postprobs <- all.models.PM[,"Post"]

  #arguments used for imputation
  result$imp.args <- list(parallelmice = parallelmice,
                          imp.mice.method = imp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lprior.models #function used for model prior

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}

GM97.Gibbs <- function (y, X0, X.full, p, namesxnotnull, NAvars,
                        lprior.models, lBF.method,
                        init.model, n.iter, n.burnin, n.thin, Gibbs.seed) {
  #Gibbs sampling algorithm, originally proposed by George and McCulloch (1997)
  #and further studied by Garcia-Donato and Martinez-Beneito (2013), to explore
  #the model space and approximate the model posterior distribution
  #progress bar for loop
  pb <- txtProgressBar(min = 0,
                       max = n.iter + n.burnin,
                       style = 3,
                       width = 50,
                       char = "=")

  set.seed(Gibbs.seed)
  all.models.lPM <- matrix(0, nr = n.iter + n.burnin, nc = p+1) #last column contains log(BF_a0*Pr(M))
  #Rao-Blackwellized inclusion probabilities:
  inclprobRB <- matrix(0, nr = n.iter + n.burnin, nc = p)

  current.model <- init.model
  if (sum(current.model) == 0) { #null
    lBF.PMcurrent <- lprior.models(current.model)
  } else {
    lBF.PMcurrent <- lBF.method(model = which(current.model == 1)) +
                     lprior.models(current.model) #log(BF_a0*Pr(M))
  }
  #visited models in decimal notation and the corresponding log(BF_a0*Pr(M))
  visited.models.lBF.PM <- matrix(c(sum(current.model * 2^(0:(p-1))),
                                    lBF.PMcurrent), nc = 2)
  for (i in seq_len(n.iter + n.burnin)){
    setTxtProgressBar(pb, i)
    for (j in seq_len(p)){
      proposal.model <- current.model; proposal.model[j] <- 1 - current.model[j]

      #avoiding recomputing the BF for models already visited
      already.visited <- which(visited.models.lBF.PM[,1] ==
                                 sum(proposal.model * 2^(0:(p-1))))
      if (length(already.visited) > 0) {
        lBF.PMproposal <- visited.models.lBF.PM[already.visited, 2]
      } else {
        #Check if proposal.model is the null model
        if(sum(proposal.model) > 0){
          #check if there are NAs in the model considered to save computation time
          if (any(namesxnotnull[which(proposal.model == 1)] %in% NAvars)) {
            lBF.PMproposal <- lBF.method(model = which(proposal.model == 1)) +
                              lprior.models(proposal.model) #log(BF_a0*Pr(M))
          } else { #if there are no missings, compute the BF by the method selected
            X.i <- cbind(X0[names(y),], X.full[names(y), which(proposal.model == 1)])
            lBF.PMproposal <- BF.approx.method(k = sum(proposal.model),
                                               X = as.matrix(X.i)) + lprior.models(proposal.model)
          }
        } else {
          lBF.PMproposal <- lprior.models(proposal.model) #BF_a0 = 1
        }
        visited.models.lBF.PM <- rbind(visited.models.lBF.PM,
                                       c(sum(proposal.model * 2^(0:(p-1))), lBF.PMproposal))
      }

      ratio <- exp(lBF.PMproposal - log(exp(lBF.PMproposal) + exp(lBF.PMcurrent)))
      if (runif(1) < ratio) {
        current.model[j] <- proposal.model[j]; lBF.PMcurrent <- lBF.PMproposal
      }

      if(i > 1) {
        inclprobRB[i,j] <- inclprobRB[i-1,j] + proposal.model[j]*ratio +
                           (1 - proposal.model[j])*(1 - ratio)
      }
    }

    all.models.lPM[i,] <-  c(current.model, lBF.PMcurrent)
  }
  colnames(all.models.lPM) <- c(namesxnotnull, "logBF.PM")
  cat("\n")
  for(j in seq_len(p)) inclprobRB[,j] <- inclprobRB[,j] / seq(1,(n.iter + n.burnin))

  if (n.burnin > 0) all.models.lPM <- all.models.lPM[-seq_len(n.burnin),] #remove burnin
  all.models.lPM <- all.models.lPM[seq(1, n.iter, by = n.thin), ] #keep 1 each n.thin iterations

  return(list(all.models.lPM = all.models.lPM, inclprobRB = inclprobRB))
}
