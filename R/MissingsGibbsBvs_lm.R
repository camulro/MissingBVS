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
#' @param data Data frame containing the data.
#' @param null.model Formula defining which is the simplest (null) model, which.
#' should be nested in the full one. By default, it is defined to be the one
#' with just the intercept.
#' @param BF.approx.method Method used to approximate Bayes factors with missing
#' data (to be literally specified). Possible choices include "BIC", "TBF" and
#' "gprior" (see details).
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices are the ones implemented
#' in \pkg{BayesVarSel}: "Robust", "Liangetal", "gZellner", "ZellnerSiow",
#' "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified). Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the model space of the
#' factor levels (to be literally specified). Possible choices are "Constant" and
#' "ScottBerger" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed variables)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details.)
#' @param init.model The model at which the simulation process starts. Options
#' include "Null" (the model only with the variables specified in
#' \code{null.model}), "Full" (the model defined by \code{formula}), "Random" (a
#' randomly selected model) and a vector with p (the number of factors and/or
#' covariates to select from) zeros and ones defining a model.
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
#' of the imputation process with \code{n.imp = 30} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param imp.mice.method Method for mice's imputation. Can be a string or a
#' vector of strings of length p1, where p1 is the number of independent
#' variables given by \code{formula} with \code{NAs}.
#' @param imp.predict.mat \code{matrix} with p1 rows and p2 columns, where p2
#' is the number of variables used to impute. Each entry equals 1 if the
#' column variable is used as a predictor for the corresponding row variable in the
#' imputation step. By default, the \code{\link[mice]{quickpred}} function
#' is used for the variables with NAs in formula and 0s for the rest of them.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGibbsBVS.lm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{lmfull}{If missings on the competing variables, object of class
#' \code{\link[mice]{mipo}} that combines the estimates
#' for the model defined by \code{formula} fitted by \code{\link[stats]{lm}} over
#' the \code{n.imp} (if \code{> 1}) imputed datasets; see \code{\link[mice]{pool}}
#' for details. Otherwise, it is the \code{lm} object for the \code{formula} model}
#' \item{lmnull}{The \code{lm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{lm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model using
#' \code{inclprobRB}}
#' \item{positions}{\code{matrix} with L rows and p plus the number of dummies
#' resulting from factors columns - L, where L is the number of factors, with 1
#' if the column dummy corresponds to the row factor and 0 otherwise}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate}
#' \item{modelsprob}{A \code{floor(n.iter/n.thin)}x(p+1) \code{matrix} which
#' summaries the keeped models and their associated Bayes factor in logaritmic scale}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Estimated posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{logprior.models.dumm}{Function used to compute the log-prior over the
#' model space of the factor levels}
#' \item{method}{\code{Gibbs}}
#'
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
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
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Default choices are: robust and ScottBerger priors, 10000 iterations with 500
#' #of burn in period. Here, 100 imputations with mice's pmm method.
#' dataS97.mGBVS <- missingGibbsBVS.lm(formula = gr56092 ~ ., data = dataS97,
#'                                     n.iter = 1000, n.burnin = 50, n.imp = 100)
#'
#' #Show the results:
#' dataS97.mGBVS
#'
#' #Summ up the results:
#' summary(dataS97.mGBVS)
#'
#' #A plot with the estimated posterior inclusion probabilities for each
#' #competing variable and the dimension probability of the true model:
#' plot(dataS97.mGBVS)
#' }
#'
missingGibbsBVS.lm <- function (formula,
                                data,
                                null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                                BF.approx.method = "gprior",
                                prior.betas = "Robust",
                                prior.models = "ScottBerger",
                                prior.models.dummies = "ScottBerger",
                                priorprobs = NULL,
                                init.model = "Full",
                                n.iter = 10000,
                                n.burnin = 500,
                                n.thin = 1,
                                parallelmice = NULL,
                                n.core = NULL,
                                imp.time.test = TRUE,
                                imp.mice.method = "pmm",
                                imp.predict.mat = NULL,
                                n.imp = 039E1,
                                Gibbs.seed = runif(1,0,26061970),
                                imp.seed = runif(1,0,09011975)) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #Response in the null model and full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(formula, null.model, data)
  list2env(buildmatrices.list, envir = environment())

  #Check arguments and compute init.model
  init.model <- checkGibbsarguments(p, p0, namesnull, namesx, init.model)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  if (L > 0) {
    lprior.models.dummies <- checkforprior.models.dummies(prior.models.dummies, l)
  } else lprior.models.dummies <- function(delta, tau) 0

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.lm(BF.approx.method,
                                             prior.betas,
                                             n, p, p0, y, SS0)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define competing variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)

  #Imputation step
  if (!is.null(NAvars)) {
    buildimputation.list <- buildimputation(NAvars, formula, data, imp.predict.mat,
                                            n.imp, n, q, p0, imp.time.test, imp.mice.method, imp.seed,
                                            parallelmice, n.core, obsnotNA, ordvars, BF.approx.method)
    list2env(buildimputation.list, envir = environment())
  }

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", q + q0, "covariates and/or factors.\n")
  if (q0 == 1) {
    cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
               q, ".\n"))
  } else  cat(paste0("From those ", q0, " are fixed and we should select from the remaining ",
               q, ".\n"))

  cat("  Numerical covariates:", depvars[positionsx], "\n")
  if (L > 0) cat(" Factors:", depvars[!positionsx], "\n")

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  gibbs.list <- GM97.Gibbs(y, X0, X.full, p, namesxnotnull, NAvars,
                           lprior.models, lprior.models.dummies, lBF.method, BF.approx.method,
                           positions, positionsfac, l, L, lastd,
                           init.model, n.iter, n.burnin, n.thin, Gibbs.seed)
  list2env(gibbs.list, envir = environment())

  #Summ up Gibbs sampling results
  summ.Gibbs.list <- summ.Gibbs(cf.models.lPM, inclprobRB, q, n.iter, n.thin,
                                lprior.models, lprior.models.dummies, positionsx)
  list2env(summ.Gibbs.list, envir = environment())

  if (!is.null(NAvars)) {
    #Pool results for imputed datasets
    if (n.imp > 1) {
      #remove last dummy for each factor, first p0 vars are the fixed ones
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0),]
      #Evaluate lm of full model with missings using Rubin's rule
      fit <- list()
      mt <- attr(framefull, "terms")
      for (i in 1:n.imp) {
        z <- lm.fit(x = imputation.array[,,i], y = y)
        z$terms <- mt
        class(z) <- "lm"

        fit[[i]] <- z
      }
      lmfull <- mice::pool(fit)
    } else {
      #compute lm.fit for the unique imputation
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0)]
      lmfull <- lm.fit(x = imputation.array, y = y)
    }
  } else lmfull <- lm(formula, data)

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the programm to finish
  result$lmfull <- lmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, lmfull is the lm object for the full model
  result$lmnull <- lmnull #The lm object for the null model (omits NAs)

  result$variables <- depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- q #number of competing vars
  result$k <- q0 #number of fixed vars
  result$HPMbin <- hpm #The binary code for the HPM model
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- depvars

  if (L > 0) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  result$modelsprob <- cf.models.lBF

  #The binary code for all the visited models (after n.thin is applied) and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  result$inclprobRB <- inclprobRB #Rao-Blackwellized inclusion probability

  result$postprobdim <- probdim/sum(probdim) #vector with the estimated posterior dimension probability
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model

  result$call <- match.call()
  result$C <- C #estimated normilizing constant

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, q + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, q))) #prior inclusion prob for dimension 0
    for (i in seq_len(q)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, q - i))) + lchoose(q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:q + q0 #dimension prior probability

  #Estimation of posterior probabilities based on C
  result$postprobs <- post

  if (!is.null(NAvars)) {
    #arguments used for imputation
    result$imp.args <- imp.args

    #save the imputed datasets for sensitivity analysis
    # raw.imp.array <- serialize(imputation.array, NULL)
    # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lprior.models #function used for model prior
  if (L > 0) result$logprior.models.dumm <- lprior.models.dummies

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}

#' @keywords internal
GM97.Gibbs <- function (y, X0, X.full, p, namesxnotnull, NAvars,
                        lprior.models, lprior.models.dummies, lBF.method, BF.approx.method,
                        positions, positionsfac, l, L, lastd,
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
  all.models.lPM <- matrix(0, nr = n.iter + n.burnin, nc = p+1) #last column is log(BF_a0*Pr(M))
  #Rao-Blackwellized inclusion probabilities:
  inclprobRB <- matrix(0, nr = n.iter + n.burnin, nc = p)

  current.model <- init.model
  gt <- positions %*% current.model > 0 #covariates and/or factors active
  d <- positionsfac %*% current.model #levels of factors
  t <- d > 0 #active factors

  #change saturated or oversaturated model for c(1,...,1,0)
  if (sum(t) > 0) {
    M <- t(apply(positionsfac, 1, function(x) x * current.model))
    f.check <- (d == l) | ((d == l - 1) & lastd(M))
    if(any(f.check)) {
      for (f in which(f.check)) {
        current.model[as.logical(positionsfac[f,])] <- c(rep(1,l[f]-1), 0)
      }
    }
  }

  if (sum(current.model) == 0) { #null
    lBF.PMcurrent <- lprior.models(gt)
  } else {
    #check if there are NAs in the model considered to save computation time
    if (any(namesxnotnull[which(current.model == 1)] %in% NAvars)) {
      lBF.PMcurrent <- lBF.method(model = which(current.model == 1)) +
        lprior.models(gt) + lprior.models.dummies(d, t) #log(BF_a0*Pr(M))
    } else { #if there are no missings, compute the BF by the method selected
      X.i <- cbind(X0, X.full[, which(current.model == 1)])
      lBF.PMcurrent <- BF.approx.method(k = sum(current.model == 1), X = X.i) +
        lprior.models(gt) + lprior.models.dummies(d, t) #log(BF_a0*Pr(M))
    }
  }

  #visited models with hash and the corresponding log(BF_a0*Pr(M))
  visited.models.lBF.PM <- list()
  visited.models.lBF.PM$models <- digest::digest(current.model)
  visited.models.lBF.PM$lBF.PM <- lBF.PMcurrent
  for (i in seq_len(n.iter + n.burnin)){
    setTxtProgressBar(pb, i)
    for (j in seq_len(p)){
      proposal.model <- current.model; proposal.model[j] <- 1 - current.model[j]

      gt <- positions %*% proposal.model > 0 #covariates and/or factors active
      d <- positionsfac %*% proposal.model #levels of factors
      t <- d > 0 #active factors

      #change saturated or oversaturated model for c(1,...,1,0)
      if (sum(t) > 0) {
        M <- t(apply(positionsfac, 1, function(x) x * proposal.model))
        f.check <- (d == l) | ((d == l - 1) & lastd(M))
        if(any(f.check)) {
          for (f in which(f.check)) {
            proposal.model[as.logical(positionsfac[f,])] <- c(rep(1,l[f]-1), 0)
          }
        }
      }
      hash.proposal <- digest::digest(proposal.model)

      #avoiding recomputing the BF for models already visited
      already.visited <- which(visited.models.lBF.PM$models == hash.proposal)
      if (length(already.visited) > 0) {
        lBF.PMproposal <- visited.models.lBF.PM$lBF.PM[already.visited]
      } else {
        #Check if proposal.model is the null model
        if(sum(proposal.model) > 0){
          #check if there are NAs in the model considered to save computation time
          if (any(namesxnotnull[which(proposal.model == 1)] %in% NAvars)) {
            lBF.PMproposal <- lBF.method(model = which(proposal.model == 1)) +
              lprior.models(gt) + lprior.models.dummies(d, t) #log(BF_a0*Pr(M))
          } else { #if there are no missings, compute the BF by the method selected
            X.i <- cbind(X0, X.full[, which(proposal.model == 1)])
            lBF.PMproposal <- BF.approx.method(k = sum(proposal.model == 1), X = X.i) +
              lprior.models(gt) + lprior.models.dummies(d, t) #log(BF_a0*Pr(M))
          }
        } else { #null
          lBF.PMproposal <- lprior.models(gt) + lprior.models.dummies(d, t) #BF_a0 = 1
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
  cat("\n")
  for(j in seq_len(p)) inclprobRB[,j] <- inclprobRB[,j] / seq(1,(n.iter + n.burnin))

  if (n.burnin > 0) all.models.lPM <- all.models.lPM[-seq_len(n.burnin),] #remove burnin
  all.models.lPM <- all.models.lPM[seq(1, n.iter, by = n.thin), ] #keep 1 each n.thin iterations

  #models matrix at the covariate-factor level
  cf.models.lPM <- all.models.lPM[,seq_len(p)] %*% t(positions)
  cf.models.lPM <- cbind(cf.models.lPM, all.models.lPM[, p+1])
  colnames(cf.models.lPM)[ncol(cf.models.lPM)] <- "logBF.PM"
  #cf.models.lPM is exactly all.models.lPM if there are no factors

  return(list(cf.models.lPM = cf.models.lPM, inclprobRB = inclprobRB %*% t(positions)))
}

#' @keywords internal
summ.Gibbs <- function (cf.models.lPM, inclprobRB, q, n.iter, n.thin,
                        lprior.models, lprior.models.dummies, positionsx) {
  #Summ up Gibbs sampling results

  inclprob <- colMeans(cf.models.lPM[,-(q+1)]) #inclusion probabilities except for fixed variables

  cf.models.PM <- cf.models.lPM[,seq_len(q)]
  cf.models.PM <- cbind(cf.models.PM, exp(cf.models.lPM[, q + 1]))
  colnames(cf.models.PM)[q+1] <- "BF.PM"

  #Estimation of the normalizing constant:
  K <- round(dim(cf.models.PM)[1]/2)
  Aset <- sample(x = 1:dim(cf.models.PM)[1], size = K, replace = F)
  Bset <- (1:dim(cf.models.PM)[1])[-Aset]
  #Bayes factors multiplied by prior probs of the models in A
  BF.PMAset <- cf.models.PM[Aset, "BF.PM"]
  #Remove duplicates
  BF.PMAset <- unique(BF.PMAset)
  gAset <- sum(BF.PMAset)
  #How many of the models in Bset are in A?
  sumIA <- sum(cf.models.PM[Bset,"BF.PM"] %in% cf.models.PM[Aset,"BF.PM"])
  #Normalizing constant estimation
  C <- gAset*K / sumIA

  #compute estimated posterior probabilities
  cf.models.PM[, q + 1] <- exp(cf.models.lPM[, q + 1] - log(C))
  post <- cf.models.PM[, q + 1]

  probdim <- rep(0, q + 1)
  cf.models.lBF <- cf.models.lPM
  #compute posterior probability of the dimension of the true model and
  #save logBF for each model
  for (i in seq_len(floor(n.iter / n.thin))) {
    probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] <-
      probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] + cf.models.PM[i, q + 1]

    gt <- cf.models.lPM[i, seq_len(q)] > 0
    d <- cf.models.lPM[i, !positionsx]; t <- d > 0
    cf.models.lBF[i, q + 1] <- cf.models.lPM[i, q + 1] -
      lprior.models(gt) - lprior.models.dummies(d, t) # lBF
  }
  colnames(cf.models.lBF)[q+1] <- "logBF"

  #HPM
  nPmax <- which.max(cf.models.lBF[, q+1])
  hpm <- cf.models.lBF[nPmax, ]

  #MPM
  mpm <- rep(0,q)
  mpm[which(inclprobRB[n.iter, ] >= 0.5)] <- 1

  return(list(cf.models.lBF = cf.models.lBF, post = post, C = C,
              inclprob = inclprob, probdim = probdim, hpm = hpm, mpm = mpm))
}

#' @keywords internal
checkGibbsarguments <- function (p, p0, namesnull, namesx, init.model) {
  #check arguments
  #Is there any variable to select from?
  if (p == 0) {
    stop(paste0("The number of fixed variables is equal to the number of\n",
                "regressors in the full model. No model selection can be done.\n"))
  }

  if (p <= 20) {
    warning(paste0("The number of variables is small enough to visit every model.\n",
                   "Consider using the Bvs version.\n"), immediate. = TRUE)
  }

  #check if null model is contained in the full one:
  for (i in 1:p0){
    if (namesnull[i] %notin% namesx) {
      stop(paste0("Error in var: ", namesnull[i],"; null model not nested in full model.\n"))
    }
  }

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

  return(init.model)
}
