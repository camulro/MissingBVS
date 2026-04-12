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
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 30} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param initialimp.mice.method Method for mice's imputation.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGibbsGD25} returns an object of class \code{missingBVS}
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
#' \item{MPMbin}{Binary expression of the Median Probability model using
#' \code{inclprobRB}}
#' \item{modelsprob}{A (n.keep)x(p+1) \code{matrix} which summaries the \code{n.keep}
#' most probable a posteriori models and their associated probability}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables.}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Estimated posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the
#' model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{method}{\code{Gibbs}}
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
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @keywords package
#'
#' @examples
#' \dontrun{
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' #Here we keep the 8 competing models:
#' f <- Solar.R ~ 1 + Ozone + Wind + Temp
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
missingGibbsGD25 <- function (formula,
                              data,
                              prior.models = "ScottBerger",
                              priorprobs = NULL, #needed if prior.models = User
                              init.model = "Full",
                              n.iter = 10000, #number of iterations for Gibbs Sampling algorithm
                              n.burnin = 500,
                              n.thin = 1,
                              imp.time.test = TRUE,
                              initialimp.mice.method = "pmm", #mice's default
                              n.imp = 039E1, #number of imputed datasets for BF
                              Gibbs.seed = runif(1,0,26061970), #seed for the Gibbs sampling
                              imp.seed = runif(1,0,09011975)) { #seed for the imputation

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(paste(formula[[2]], " ~ 1", sep=""))

  #Check for numeric covariates
  aux <- model.frame(formula, data)
  isnum <- sapply(aux, is.numeric)
  # isint <- sapply(aux, is.integer)
  # if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
  if (sum(isnum) < dim(aux)[2]) {
    stop("This method is only for continuous covariates.\n")
  }
  # cat("Be careful, this method is only for normally distributed covariates.\n",
  #     "Do you want to continue? (y/n)\n")
  # if (tolower(readline()) != "y") {
  #   stop("Try the missingGibbsBVS.lm function instead.\n")
  # }

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

  if (p <= 20) {
    warning(paste0("The number of variables is small enough to visit every model.\n",
                   "Consider using MissingBvs.lm.\n"), immediate. = TRUE)
  }

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #check for missings
  NAvars <- checkformissings(y = framefull[,1], X.full = X.full[obsnotNA,])

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
  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model with García-Donato's 2025 method
    lBF.method <- function (model) lBF.miss(model,
                                            imputation.list = imputation.list,
                                            BF.miss.aux = BF.miss.aux,
                                            n = n, nMC = n.imp)
  } else lBF.method <- function (model) BF.miss.aux(X.center = imputation.list$rX.imput[,model],
                                                    Sigma11 = imputation.list$rSigma[model, model,],
                                                    k = length(model))

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", p + 1, "single covariates.\n")
  cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
      p, ":\n"))
  cat(paste(paste(namesx, collapse = ", ", sep = ""), "\n", sep = ""))

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #progress bar for loop
  pb <- txtProgressBar(min = 0,
                       max = n.iter + n.burnin,
                       style = 3,
                       width = 50,
                       char = "=")

  #George and McCulloch's Gibbs exploration
  set.seed(Gibbs.seed)
  all.models.lPM <- matrix(0, nr = n.iter + n.burnin, nc = p+1) #last column is log(BF_a0*Pr(M))
  #Rao-Blackwellized inclusion probabilities:
  inclprobRB <- matrix(0, nr = n.iter + n.burnin, nc = p)

  current.model <- init.model
  #If init.model is the null
  if(sum(current.model) > 0){
    lBF.PMcurrent <- lBF.method(model = which(current.model == 1)) +
      lprior.models(current.model) #log(BF_a0*Pr(M))
  } else { #null
    lBF.PMcurrent <- lprior.models(current.model) #BF_a0 = 1
  }

  #visited models with hash and the corresponding log(BF_a0*Pr(M))
  visited.models.lBF.PM <- list()
  visited.models.lBF.PM$models <- digest::digest(current.model)
  visited.models.lBF.PM$lBF.PM <- lBF.PMcurrent
  for (i in seq_len(n.iter + n.burnin)){
    setTxtProgressBar(pb, i)
    for (j in seq_len(p)){
      proposal.model <- current.model; proposal.model[j] <- 1 - current.model[j]
      hash.proposal <-  digest::digest(proposal.model)

      #avoiding recomputing the BF for models already visited
      already.visited <- which(visited.models.lBF.PM$models == hash.proposal)
      if (length(already.visited) > 0) {
        lBF.PMproposal <- visited.models.lBF.PM$lBF.PM[already.visited]
      } else {
        #Check if proposal.model is the null model
        if(sum(proposal.model) > 0){
          lBF.PMproposal <- lBF.method(model = which(proposal.model == 1)) +
                            lprior.models(proposal.model) #log(BF_a0*Pr(M))
        } else { #null
          lBF.PMproposal <- lprior.models(proposal.model) #BF_a0 = 1
        }
        visited.models.lBF.PM$models <- c(visited.models.lBF.PM$models, hash.proposal)
        visited.models.lBF.PM$lBF.PM <- c(visited.models.lBF.PM$lBF.PM, lBF.PMproposal)
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
  colnames(all.models.lPM) <- c(namesx, "logBF.PM")

  for(j in seq_len(p)) inclprobRB[,j] <- inclprobRB[,j] / seq(1,(n.iter + n.burnin))
  colnames(inclprobRB) <- namesx

  if (n.burnin > 0) all.models.lPM <- all.models.lPM[-seq_len(n.burnin),] #remove burnin
  all.models.lPM <- all.models.lPM[seq(1, n.iter, by = n.thin), ] #keep 1 each n.thin iterations

  inclprob <- colMeans(all.models.lPM[,-(p+1)]) #inclusion probabilities except for fixed variables

  all.models.PM <- all.models.lPM[,seq_len(p)]
  all.models.PM <- cbind(all.models.PM, exp(all.models.lPM[, p + 1]))
  colnames(all.models.PM)[p+1] <- "BF.PM"

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
  colnames(all.models.PM)[p+1] <- "Post"

  probdim <- rep(0, p + 1)
  all.models.lBF <- all.models.lPM
  #compute posterior probability of the dimension of the true model and
  #save logBF for each model
  for (i in seq_len(floor(n.iter / n.thin))) {
      probdim[sum(all.models.PM[i, seq_len(p)]) + 1] <-
        probdim[sum(all.models.PM[i, seq_len(p)]) + 1] + all.models.PM[i, p + 1]

    all.models.lBF[i, p + 1] <- all.models.lPM[i, p + 1] -
      lprior.models(all.models.lPM[i, seq_len(p)]) # lBF
  }
  colnames(all.models.lBF)[p+1] <- "logBF"

  #HPM
  nPmax <- which.max(all.models.lBF[, p+1])
  hpm <- all.models.lBF[nPmax, ]

  #MPM
  mpm <- rep(0,p)
  mpm[which(inclprobRB[n.iter, ] >= 0.5)] <- 1

  if (!is.null(NAvars)) {
    if (n.imp > 1) {
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
    } else {
      #compute lm.fit for the unique imputation
      lmfull <- lm.fit(x = cbind(1, imputation.list$rX.imput), y = y)
    }
  } else lmfull <- lm(formula, data)

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

  result$modelsprob <- all.models.lBF

  #The binary code for all the visited models (after n.thin is applied) and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  result$inclprobRB <- inclprobRB #Rao-Blackwellized inclusion probability

  result$postprobdim <- probdim/sum(probdim) #vector with the estimated posterior dimension probability
  names(result$postprobdim) <- 0:p + 1 #dimension of the true model

  result$call <- match.call()
  result$C <- C #estimated normilizing constant

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, p + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, p))) #prior inclusion prob for dimension 0
    for (i in seq_len(p)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, p - i))) + lchoose(p, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:p + 1 #dimension prior probability

  #Estimation of posterior probabilities based on C
  result$postprobs <- all.models.PM[,"Post"]

  #arguments used for imputation
  result$imp.args <- list(initialimp.mice.method = initialimp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  #save the imputed datasets for sensitivity analysis
  # raw.imp.array <- serialize(imputation.array, NULL)
  # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$logprior.models <- lprior.models #function used for model prior

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}
