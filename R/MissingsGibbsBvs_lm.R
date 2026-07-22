#' Bayesian Imputation Averaging for Variable Selection with Missing data in
#' linear regression models using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior model distribution using a
#' Gibbs sampling algorithm to explore the model space. Each posterior model probability
#' is computed following the Bayesian Imputation Averaging (BIA) framework, using
#' standard priors for model coefficients and the hierarchical approach of
#' García-Donato and Paulo (2022) with factors.
#'
#' Gibbs sampling search algorithm to avoid exhaustive enumeration of model space
#' when it is unfeasible. It draws from the model posterior distribution and uses
#' frequency of "visits" to construct the estimates. The algorithm was originally
#' proposed by  George and McCulloch (1997). Later, Garcia-Donato and Martinez-Beneito (2013)
#' shown that the sampling strategy in combination with estimates based on frequency of
#' visits provides very reliable results.
#'
#' \code{\link[MissingBVS]{MissingGibbsBvs.lm}} is a heuristic approximation of
#' \code{\link[MissingBVS]{MissingBvs.lm}}. See the later for common details.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param null.model Formula defining which is the simplest (null) model, nested in
#' the full one with possible fixed variables. By default, it is defined to be the one
#' with just the intercept.
#' @param BF.approx.method Method used to compute or approximate data-driven Bayes factors
#' (to be literally specified). Possible choices include "BIC", "TBF" and "gprior"
#' (see details).
#' @param prior.betas Prior distribution for model coefficients if "gprior" method is
#' chosen (to be literally specified). Possible choices are: "Robust", "Liangetal",
#' "gZellner", "ZellnerSiow", "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Model prior distribution over the covariates and/or factors
#' model space (to be literally specified). Possible choices are "Constant",
#' "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the dummies submodel space
#' given by the active factors (to be literally specified). Possible choices are
#' "Constant" and "ScottBerger" (see details).
#' @param marginal.factors Logical to indicate whether or not to marginalize factors'
#' probabilities such as García-Donato and Paulo (2022). By default, it is set to TRUE.
#' @param priorprobs A p+1 (being p the number of non-fixed variables) dimensional
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
#' @param imp.mice.method Method for \pck{mice}'s imputation. Can be either a string
#' or a vector of strings of length the number of variables in data, except the response.
#' @param imp.predict.mat Matrix with \code{formula}'s competing variables in rows
#' and some \code{data}'s variables in columns. Each entry equals 1 if the column variable
#' is used as a predictor for the corresponding row variable in the imputation step. Order
#' in columns defines the imputation visit sequence. By default, a shortcut is used to
#' define the most important predictors for each variable based on correlations.
#' @param n.imp Number of imputed datasets for model posterior computation.
#' @param maxit Number of iterations for \pck{mice}'s imputation. By default, it is 5.
#' @param parallelmice Logical to indicate whether or not to use parallelization on
#' \code{\link[mice]{mice}}'s imputation. By default, automatically performs it if the
#' number of imputations or competing variables given by \code{formula} are big enough.
#' @param n.core Number of cores for parallel imputation.
#' @param imp.datasets Array or list for imputed datasets if given by user. By default
#' it is set to NULL and imputation is performed following other imputation arguments.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingGibbsBVS.lm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{Time lasted solving the problem}
#' \item{lmfull}{If missings on the \code{formula} competing variables, combination
#' of the estimates of fitted full model over the \code{n.imp} imputed datasets.
#' Otherwise, it is the \code{\link[stats]{lm}} object}
#' \item{lmnull}{The \code{lm} class object that results when \code{null.model}
#' is fitted by \code{\link[stats]{lm}}}
#' \item{variables}{Names of all the competing variables given by \code{formula}}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables (covariates and/or factors) to select from}
#' \item{k}{Number of fixed variables given by \code{null.model}}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model using \code{inclprobRB}}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{modelsrankdefprob}{A floor(\code{n.iter}/\code{n.thin}) x (p1 * (sum_j l_j - L)+1)
#' matrix which summarizes the floor(\code{n.iter}/\code{n.thin}) visited submodels
#' and their log-Bayes factors (when relevant)}
#' \item{modelslogBF}{A floor(\code{n.iter}/\code{n.thin}) x (p+1) matrix which summarizes
#' the floor(\code{n.iter}/\code{n.thin}) visited models and their associated log-Bayes factors}
#' \item{inclprob}{Named vector with the inclusion probabilities of p competing variables}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Estimated posterior probabilities over the true model size}
#' \item{C}{The value of the estimated normalizing constant}
#' \item{postprobs}{Estimated posterior probability}
#' \item{call}{The \code{call} to the function}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information (when relevant)}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{BF.approx.method}{Function used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
#' \item{logprior.models}{Function used to compute the log-prior over the model space
#' defined by covariates and/or factors}
#' \item{prior.models}{Two-dimensional vector with \code{prior.models} and
#' \code{prior.models.dummies} chosen. If there are no factors or \code{marginal.factors}
#' is set to FALSE, it saves the only argument used, \code{prior.models}}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' model space probabilities such as García-Donato and Paulo (2022)}
#' \item{method}{String "Gibbs" denoting Gibbs sampling model search}
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
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Default choices are: robust and ScottBerger priors, 10000 iterations with 500
#' #of burn in period. Here, 100 imputations with mice's pmm method.
#' dataS97.mGBVS <- missingGibbsBVS.lm(formula = gr56092 ~ ., data = dataS97,
#'                                     n.iter = 1000, n.burnin = 50, n.imp = 10)
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
                                marginal.factors = TRUE,
                                priorprobs = NULL,
                                init.model = "Full",
                                n.iter = 10000,
                                n.burnin = 500,
                                n.thin = 1,
                                imp.mice.method = "pmm",
                                imp.predict.mat = NULL,
                                n.imp = 039E1,
                                maxit = 5,
                                parallelmice = NULL,
                                n.core = NULL,
                                imp.datasets = NULL,
                                Gibbs.seed = runif(1,0,26061970),
                                imp.seed = runif(1,0,09011975)) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #Response in the null model and full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  env <- environment()

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(formula, null.model, data, marginal.factors)
  list2env(buildmatrices.list, envir = env)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  mF <- L > 0 & marginal.factors
  #Check arguments and compute init.model
  init.model <- checkGibbsarguments(p, p0, namesnull, namesx, init.model, mF,
                                    positions, positionsfac, l, firstd)

  #Check if factors present and if marginalization of their probabilities. Define model prior
  lp.model <- checkmarg.factorsprior(mF, prior.models.dummies, l, positions, positionsfac,
                                     firstd, lprior.models)

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.lm(BF.approx.method, prior.betas, n, p, p0, y, SS0)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define competing variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)
  #Imputation step
  if (anyNAvar <- !is.null(NAvars)) {
    if (is.null(imp.datasets)) { #if there are no given imputations, build them
      imputation.list <- buildimputation(NAvars, formula, data, imp.predict.mat, n.imp, maxit,
                                         n, q, p0, imp.mice.method, imp.seed,
                                         parallelmice, n.core, obsnotNA, ordvars, BF.approx.method)
    } else imputation.list <- extimputation(formula, imp.datasets, n0 = dim(data)[1], framefull,
                                            ordvars, obsnotNA, p0, BF.approx.method, NAvars)
    list2env(imputation.list, envir = env)
  }

  #Info:
  cat("Info. . .\n")
  if (mF) {
    cat("Most complex model has a total of", q + q0, "covariates and/or factors.\n")
  } else cat("Most complex model has a total of", q + q0, "competing variables.\n")
  if (q0 == 1) {
    cat("From those 1 is fixed (the intercept) and we should select from the remaining",
        q, ".\n")
  } else cat("From those", q0, "are fixed and we should select from the remaining",
             q, ".\n")
  if (mF) {
    cat("  Numerical covariates:", depvars[positionsx], "\n")
    cat(" Factors:", depvars[!positionsx], "\n")
  } else  cat("  Competing variables:", depvars, "\n")

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  set.seed(Gibbs.seed)
  gibbs.list <- GM97.Gibbs(X0, X.full, p, namesxnotnull, NAvars, lp.model, lBF.method,
                           BF.approx.method, mF, positions, init.model, n.iter, n.burnin, n.thin)
  list2env(gibbs.list, envir = env)

  #Summ up Gibbs sampling results
  summ.Gibbs.list <- summ.Gibbs(cf.models.lBF, all.lBF.PM, inclprobRB, q, n.iter)
  list2env(summ.Gibbs.list, envir = env)

  if (anyNAvar) {#Pool results for imputed datasets
    imp.array <- imputation.array
    #remove first dummy on each factor, first p0 vars are the fixed ones
    if (mF) imp.array <- imp.array[,-c(indf + p0), , drop = FALSE]
    #Evaluate lm of full model with missings using Rubin's rule
    fit <- list(); mt <- attr(framefull, "terms")
    for (i in 1:n.imp) {
      z <- lm.fit(x = imp.array[,,i], y = y)
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

  result$variables <- depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- q #number of competing vars
  result$k <- q0 #number of fixed vars
  result$HPMbin <- hpm #The binary code for the HPM model
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- depvars

  if (mF) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
    result$modelsrankdefprob <- cbind(all.models.lBF[,-(p+1)], post) # rank deficient models and probs
  }

  #The binary code for all the visited models (after n.thin is applied) and the logBF
  result$modelslogBF <- cf.models.lBF

  result$inclprob <- inclprob #inclusion probability for each variable
  result$inclprobRB <- inclprobRB[n.iter, ] #Rao-Blackwellized inclusion probability
  names(result$inclprobRB) <- depvars

  result$postprobdim <- probdim #vector with the estimated posterior dimension probability
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model
  result$C <- C #estimated normilizing constant
  #Estimation of posterior probabilities based on C
  result$postprobs <- post

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- numeric(q+1)
    priorprobs[1] <- exp(lprior.models(numeric(q))) #prior inclusion prob for dimension 0
    for (i in seq_len(q)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep.int(1, i), rep.int(0, q - i))) + lchoose(q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:q + q0 #dimension prior probability

  if (anyNAvar) {
    #arguments used for imputation
    result$imp.info <- imp.info

    #save the imputed datasets for BMA or sensitivity analysis
    raw.imp.array <- serialize(imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lp.model #function used for model prior
  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models
  result$marginal.factors <- marginal.factors #whether or not factors are marginalized

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}

#' @keywords internal
GM97.Gibbs <- function (X0, X.full, p, namesxnotnull, NAvars, lp.model, lBF.method, BF.approx.method,
                        mF, positions, init.model, n.iter, n.burnin, n.thin) {
  #Gibbs sampling algorithm, originally proposed by George and McCulloch (1997)
  #and further studied by Garcia-Donato and Martinez-Beneito (2013), to explore
  #the model space and approximate the model posterior distribution progress bar for loop
  pb <- txtProgressBar(min = 0, max = 2^p, style = 3, width = 50, char = "=")

  all.models.lBF <- matrix(0, nr = n.iter + n.burnin, nc = p+1) #last column is log(BF_a0)
  all.lBF.PM <- numeric(n.iter + n.burnin) #log(BF_a0*Pr(M))
  #Rao-Blackwellized inclusion probabilities:
  inclprobRB <- matrix(0, nr = n.iter + n.burnin, nc = p)

  current.model <- init.model
  lpm <- lp.model(current.model) #log-model prior
  #lpm of init.model cannot be NA since it has been changed for the non-saturated, if needed

  if (sum(current.model) == 0) { #null
    lBFcurrent <- 0
  } else {
    #check if there are NAs in the model considered to save computation time
    if (any(namesxnotnull[which(current.model == 1)] %in% NAvars)) {
      lBFcurrent <- lBF.method(model = which(current.model == 1))
    } else { #if there are no missings, compute the BF with the method selected
      X.i <- cbind(X0, X.full[, which(current.model == 1)])
      lBFcurrent <- BF.approx.method(k = sum(current.model == 1), X = X.i)
    }
  }
  lBF.PMcurrent <- lBFcurrent + lpm #log(BF_a0*Pr(M))

  #visited models with the corresponding hash, log(BF_a0) and log(BF_a0*Pr(M))
  visited.models <- list()
  visited.models$models <- digest::digest(current.model)
  visited.models$lBF <- lBFcurrent;  visited.models$lBF.PM <- lBF.PMcurrent
  for (i in seq_len(n.iter + n.burnin)){
    setTxtProgressBar(pb, i)
    for (j in seq_len(p)){
      proposal.model <- current.model; proposal.model[j] <- 1 - current.model[j]

      lpm <- lp.model(proposal.model) #log-model prior
      if (is.na(lpm)) next #do not visit saturated or oversaturated models
      hash.proposal <- digest::digest(proposal.model)

      #avoid recomputing BF for models already visited
      already.visited <- which(visited.models$models == hash.proposal)
      if (length(already.visited) > 0) {
        lBFproposal <- visited.models$lBF[already.visited]
        lBF.PMproposal <- visited.models$lBF.PM[already.visited]
      } else {
        #Check if proposal.model is the null model
        if(sum(proposal.model) > 0){
          #check if there are NAs in the model considered to save computation time
          if (any(namesxnotnull[which(proposal.model == 1)] %in% NAvars)) {
            lBFproposal <- lBF.method(model = which(proposal.model == 1))
          } else { #if there are no missings, compute the BF by the method selected
            X.i <- cbind(X0, X.full[, which(proposal.model == 1)])
            lBFproposal <- BF.approx.method(k = sum(proposal.model == 1), X = X.i)
          }
        } else { #null
          lBFproposal <- 0
        }
        lBF.PMproposal <- lBFproposal + lpm #log(BF_a0*Pr(M))
        visited.models$models <- c(visited.models$models, hash.proposal)
        visited.models$lBF <- c(visited.models$lBF, lBFproposal)
        visited.models$lBF.PM <- c(visited.models$lBF.PM, lBF.PMproposal)
      }

      ratio <- exp(lBF.PMproposal - log(exp(lBF.PMproposal) + exp(lBF.PMcurrent)))
      if (runif(1) < ratio) {
        current.model[j] <- proposal.model[j]
        lBFcurrent <- lBFproposal; lBF.PMcurrent <- lBF.PMproposal
      }

      if(i > 1) {
        inclprobRB[i,j] <- inclprobRB[i-1,j] + proposal.model[j]*ratio +
                           (1 - proposal.model[j])*(1 - ratio)
      }
    }

    all.models.lBF[i,] <-  c(current.model, lBFcurrent)
    all.lBF.PM[i] <- lBF.PMcurrent
  }
  cat("\n")
  for(j in seq_len(p)) inclprobRB[,j] <- inclprobRB[,j] / seq(1,(n.iter + n.burnin))

  if (n.burnin > 0) { #remove burnin
    all.models.lBF <- all.models.lBF[-seq_len(n.burnin),]
    all.lBF.PM <- all.lBF.PM[-seq_len(n.burnin)]
  }
  #keep 1 each n.thin iterations
  all.models.lBF <- all.models.lBF[seq(1, n.iter, by = n.thin), ]
  all.lBF.PM <- all.lBF.PM[seq(1, n.iter, by = n.thin)]
  colnames(all.models.lBF) <- c(namesxnotnull, "logBF")

  if (mF) {
    #models matrix at the covariate-factor level
    cf.models.lBF <- all.models.lBF[,seq_len(p)] %*% t(positions)
    cf.models.lBF <- cbind(cf.models.lBF, all.models.lBF[, p+1])
    colnames(cf.models.lBF)[ncol(cf.models.lBF)] <- "logBF"

    inclprobRB <- inclprobRB %*% t(positions)
  } else cf.models.lBF <- all.models.lBF
  #cf.models.lBF is exactly all.models.lBF if there are no factor

  return(list(cf.models.lBF = cf.models.lBF, all.models.lBF = all.models.lBF,
              all.lBF.PM = all.lBF.PM, inclprobRB = inclprobRB))
}

#' @keywords internal
summ.Gibbs <- function (cf.models.lBF, all.lBF.PM, inclprobRB, q, n.iter) {
  #Summ up Gibbs sampling results

  inclprob <- colMeans(cf.models.lBF[,-(q+1)] > 0) #inclusion probabilities except for fixed variables

  nGibbs <- dim(cf.models.lBF)[1]
  dim.tab <- table(c(rowSums(cf.models.lBF[,-(q+1)] > 0), 0:q))
  probdim <- (dim.tab - 1) / nGibbs #posterior probability over the dimension

  #Estimation of the normalizing constant:
  K <- round(nGibbs/2)
  Aset <- sample(x = 1:nGibbs, size = K, replace = F)
  Bset <- (1:nGibbs)[-Aset]
  #Bayes factors multiplied by prior probs of the models in A
  BF.PMAset <- exp(all.lBF.PM)[Aset]
  #Remove duplicates
  BF.PMAset <- unique(BF.PMAset)
  gAset <- sum(BF.PMAset)
  #How many of the models in Bset are in A?
  sumIA <- sum(all.lBF.PM[Bset] %in% all.lBF.PM[Aset])
  #Normalizing constant estimation
  C <- gAset*K / sumIA

  #compute estimated posterior probabilities
  post <- exp(all.lBF.PM - log(C))

  #HPM
  nPmax <- which.max(cf.models.lBF[, q+1])
  hpm <- cf.models.lBF[nPmax, ]

  #MPM
  mpm <- numeric(q)
  mpm[which(inclprobRB[n.iter, ] >= 0.5)] <- 1

  return(list(cf.models.lBF = cf.models.lBF, post = post, C = C,
              inclprob = inclprob, probdim = probdim, hpm = hpm, mpm = mpm))
}

#' @keywords internal
checkGibbsarguments <- function (p, p0, namesnull, namesx, init.model, mF,
                                 positions, positionsfac, l, firstd) {
  #check Gibbs arguments
  #Is there any variable to select from?
  if (p == 0) stop("The number of fixed variables is equal to the number of\n",
                   "regressors in the full model. No model selection can be done.\n")

  if (p <= 20) warning("The number of variables is small enough to visit every model.\n",
                       "Consider using the Bvs version.\n", immediate. = TRUE)

  #check if null model is contained in the full one:
  for (i in 1:p0){
    if (namesnull[i] %notin% namesx) stop("Error in var: ", namesnull[i],
                                          "; null model not nested in full model.\n")
  }

  #Check the initial model:
  if (is.character(init.model) == TRUE) {
    switch (substr(tolower(init.model), 1, 1),
            "n" = {init.model <- numeric(p)}, #null model
            "f" = {init.model <- rep.int(1, p)}, #full model
            "r" = {init.model <- rbinom(n = p, size = 1, prob = .5)},
            stop("Initial model not valid.\n")
    )
  } else {
    init.model <- as.numeric(init.model > 0)
    if (length(init.model) != p) stop("Initial model with incorrect length.\n")
  }

  if (mF) {#change saturated or oversaturated model for c(0,1,...,1)
    gt <- positions %*% init.model > 0 #covariates and/or factors active
    d <- positionsfac %*% init.model #levels of factors
    t <- d > 0 #active factors

    if (sum(t) > 0) {
      M <- t(apply(positionsfac, 1, function(x) x * init.model))
      satf <- (d == l) | ((d == l - 1) & firstd(M))
      if(any(satf)) for (f in which(satf)) {
          init.model[as.logical(positionsfac[f,])] <- c(0, rep.int(1,l[f]-1))
        }
    }
  }

  return(init.model)
}
