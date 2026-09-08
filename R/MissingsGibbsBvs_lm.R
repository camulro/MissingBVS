#' Bayes Factor Averaging for Variable Selection with Missing data in
#' linear regression models using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior model distribution using a
#' Gibbs sampling algorithm to explore the model space. Each posterior model probability
#' is computed following the Bayes Factor Averaging (BFA) framework, using
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
#' \code{\link[MissingBVS]{missingGibbsBVS.lm}} is a heuristic approximation of
#' \code{\link[MissingBVS]{missingBVS.lm}}. See the latter for common details.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param null.model Formula defining which is the simplest (null) model, nested in
#' the full one with possible fixed variables. By default, it is defined to be the one
#' with just the intercept.
#' @param BF.method Method used to compute or approximate data-driven Bayes factors
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
#' @param imp.mice.method Method for \pkg{mice}'s imputation. Can be either a string
#' or a vector of strings of length the number of variables in data, except the response.
#' @param imp.predict.mat Matrix with \code{formula}'s competing variables in rows
#' and some \code{data}'s variables in columns. Each entry equals 1 if the column variable
#' is used as a predictor for the corresponding row variable in the imputation step. Order
#' in columns defines the imputation visit sequence. By default, a shortcut is used to
#' define the most important predictors for each variable based on correlations.
#' @param n.imp Number of imputed datasets for model posterior computation.
#' @param maxit Number of iterations for \pkg{mice}'s imputation. By default, it is 5.
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
#' \item{BF.method}{Method used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
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
#' @seealso Use \code{\link[MissingBVS]{missingBVS.lm}} for an exact computation
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
#' Held L, Sabanés Bové D, Gravestock I (2015).<DOI:10.1214/14-STS510> Approximate
#' Bayesian Model Selection with the Deviance Statistic. Statistical Science. 30.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
#' @examples
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Use a short chain for the example; real analyses need more iterations.
#' #Few imputations for simplicity, real analyses need more.
#' dataS97.mGBVS <- missingGibbsBVS.lm(
#'   formula = gr56092 ~ 1 + lifee060 + gdpsh60l + p60,
#'   data = dataS97, n.iter = 200, n.burnin = 50, n.imp = 2,
#'   Gibbs.seed = 1, imp.seed = 1
#' )
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
#' dataS97.mGBVS$inclprobRB
#' }
#'
missingGibbsBVS.lm <- function (formula,
                                data,
                                null.model = update(as.formula(formula), . ~ 1),
                                BF.method = "BIC",
                                prior.betas = NULL,
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

  #Build matrices and objects needed later on
  matrices <- buildmatrices(formula, null.model, data, marginal.factors)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, matrices$q)

  mF <- matrices$L > 0 & marginal.factors
  positionscov <- if (mF) {
    matrices$positions[matrices$positionsx, , drop = FALSE]
  } else NULL
  positionsfac <- if (mF) matrices$positionsfac else NULL
  l <- if (mF) matrices$l else NULL
  satmodels.repr <- if (mF) matrices$satmodels.repr else NULL

  #Check arguments and compute init.model
  init.model <- checkGibbsarguments(
    matrices$p, matrices$p0, matrices$namesnull, matrices$namesx,
    init.model, mF, positionscov, positionsfac, l
  )

  #Check if factors present and if marginalization of their probabilities. Define model prior
  lp.model <- checkmarg.factorsprior(
    mF, prior.models.dummies, matrices$l, positionscov,
    positionsfac, satmodels.repr, lprior.models
  )

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Check approx method and priors chosen and define the function to be used
  lBF <- checkforprior.betas.lm(
    BF.method, prior.betas, n, matrices$p, matrices$p0, y, SS0
  )

  matrices$X.full <- matrices$X.full[obsnotNA,]

  #check for missings and define competing variables with NAs
  NAvars <- checkformissings(
    y = matrices$framenull[, 1], matrices$framenull[, -1], matrices$X.full
  )

  #Imputation step
  if (anyNAvar <- sum(NAvars) > 0) {
    if (is.null(imp.datasets)) { #if there are no given imputations, build them
      imputation <- buildimputation(
        NAvars, formula, data, imp.predict.mat, n.imp, maxit, n, matrices$q,
        matrices$p0, imp.mice.method, imp.seed, parallelmice, n.core,
        obsnotNA, matrices$ordvars
      )
    } else {
      imputation <- extimputation(
        formula, imp.datasets, n0 = dim(data)[1], matrices$framefull,
        matrices$ordvars, obsnotNA, matrices$p0, NAvars
      )
      n.imp <- imputation$n.imp
    }
  }

  if (anyNAvar && n.imp > 1) {
    #function to compute log(BFa0) for a given model as an average of BF computed
    #by BF.method over the imputed datasets
    lBF.method <- function(model) lBF.av(
      model, imputation.array = imputation$imputation.array,
      lBF = lBF, p0 = matrices$p0, n.imp = n.imp
    )
  } else lBF.method <- function(model) lBF(
    k = length(model),
    X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ]
  )

  #Info:
  cat("Info. . .\n")
  if (mF) {
    cat("Most complex model has a total of", matrices$q + matrices$q0,
        "covariates and/or factors.\n")
  } else cat("Most complex model has a total of", matrices$q + matrices$q0,
             "competing variables.\n")
  if (matrices$q0 == 1) {
    cat("From those 1 is fixed (the intercept) and we should select from the remaining",
        matrices$q, ".\n")
  } else cat("From those", matrices$q0, "are fixed and we should select from the remaining",
             matrices$q, ".\n")
  if (mF) {
    cat("  Numerical covariates:", matrices$depvars[matrices$positionsx], "\n")
    cat("  Factors:", matrices$depvars[!matrices$positionsx], "\n")
  } else cat("  Competing variables:", matrices$depvars, "\n")

  cat("The problem has a total of", 2^matrices$p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  withr::with_seed(Gibbs.seed,
    gibbs <- GM97.Gibbs(
      matrices$X0, matrices$X.full, matrices$p, NAvars, lp.model, lBF.method,
      lBF, mF, matrices$positions, init.model, n.iter, n.burnin, n.thin
    )
  )

  #Summ up Gibbs sampling results
  gibbs.summary <- summ.Gibbs(gibbs, matrices$q, n.iter)

  if (anyNAvar) {#Pool results for imputed datasets
    imp.array <- imputation$imputation.array
    #if marginal factor probs, remove first dummy on each factor, first p0 are the fixed ones
    if (mF) {
      imp.array <- imp.array[, -c(matrices$indf + matrices$p0), , drop = FALSE]
    }
    #Evaluate lm of full model with missings using Rubin's rule
    fit <- list(); mt <- attr(matrices$framefull, "terms")
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

  result$variables <- matrices$depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- matrices$q #number of competing vars
  result$k <- matrices$q0 #number of fixed vars
  result$HPMbin <- gibbs.summary$hpm
  result$MPMbin <- gibbs.summary$mpm
  names(result$MPMbin) <- matrices$depvars

  if (mF) {
    #matrix for the factors index
    result$positions <- matrices$positionsfac
    result$positionsx <- matrices$positionsx
    result$modelsrankdefprob <- cbind(
      gibbs$all.models.lBF[, -(matrices$p + 1)], gibbs.summary$post
    )
  }

  #The binary code for all the visited models (after n.thin is applied) and the logBF
  result$modelslogBF <- gibbs$cf.models.lBF

  result$inclprob <- gibbs.summary$inclprob
  result$inclprobRB <- gibbs$inclprobRB[n.iter, ]
  names(result$inclprobRB) <- matrices$depvars

  result$postprobdim <- gibbs.summary$probdim
  names(result$postprobdim) <- 0:matrices$q + matrices$q0
  result$C <- gibbs.summary$C
  #Estimation of posterior probabilities based on C
  result$postprobs <- gibbs.summary$post

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- numeric(matrices$q + 1)
    priorprobs[1] <- exp(lprior.models(numeric(matrices$q)))
    for (i in seq_len(matrices$q)) {
      priorprobs[i+1] <- exp(lprior.models(
        c(rep.int(1, i), rep.int(0, matrices$q - i))
      ) + lchoose(matrices$q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:matrices$q + matrices$q0

  if (anyNAvar) {
    #arguments used for imputation
    result$imp.info <- imputation$imp.info

    #save the imputed datasets for BMA or sensitivity analysis
    raw.imp.array <- serialize(imputation$imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.method <- BF.method #method used for BF computation
  if (is.null(prior.betas) & BF.method %in% c("gprior", "TBF")) prior.betas <- "gZellner"
  result$prior.betas <- prior.betas
  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models
  result$marginal.factors <- marginal.factors #whether or not factors are marginalized

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}

#' @keywords internal
GM97.Gibbs <- function (X0, X.full, p, NAvars, lp.model, lBF.method, lBF,
                        mF, positions, init.model, n.iter, n.burnin, n.thin) {
  #Gibbs sampling algorithm, originally proposed by George and McCulloch (1997)
  #and further studied by Garcia-Donato and Martinez-Beneito (2013), to explore
  #the model space and approximate the model posterior distribution progress bar for loop
  pb <- txtProgressBar(min = 0, max = n.iter + n.burnin, style = 3, width = 50, char = "=")

  all.models.lBF <- matrix(0, nr = n.iter + n.burnin, nc = p+1) #last column is log(BF_a0)
  all.lBF.PM <- numeric(n.iter + n.burnin) #log(BF_a0*Pr(M))
  #Rao-Blackwellized inclusion probabilities:
  inclprobRB <- matrix(0, nr = n.iter + n.burnin, nc = p)

  current.model <- init.model
  lpm <- lp.model(current.model) #log-model prior
  #lpm of init.model cannot be NA since it has been changed for the non-saturated, if needed

  if (sum(current.model) > 0) {

    #check if there are NAs in the model considered to save computation time
    if (sum(current.model * NAvars) > 0) {
      lBFcurrent <- lBF.method(model = which(current.model == 1))
    } else { #if there are no missings, compute the BF with the method selected
      X.i <- cbind(X0, X.full[, which(current.model == 1)])
      lBFcurrent <- lBF(k = sum(current.model == 1), X = X.i)
    }

  } else lBFcurrent <- 0 #null

  lBF.PMcurrent <- lBFcurrent + lpm #log(BF_a0*Pr(M))

  #visited models with hash, log(BF_a0) and log(BF_a0*Pr(M)) respectively:
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
          if (sum(proposal.model * NAvars) > 0) {
            lBFproposal <- lBF.method(model = which(proposal.model == 1))
          } else { #if there are no missings, compute the BF by the method selected
            X.i <- cbind(X0, X.full[, which(proposal.model == 1)])
            lBFproposal <- lBF(k = sum(proposal.model == 1), X = X.i)
          }

        } else lBFproposal <- 0 #null

        lBF.PMproposal <- lBFproposal + lpm #log(BF_a0*Pr(M))
        #save results:
        visited.models$models <- c(visited.models$models, hash.proposal)
        visited.models$lBF <- c(visited.models$lBF, lBFproposal)
        visited.models$lBF.PM <- c(visited.models$lBF.PM, lBF.PMproposal)
      }

      # ratio <- exp(lBF.PMproposal - log(exp(lBF.PMproposal) + exp(lBF.PMcurrent)))
      ratio <- plogis(lBF.PMproposal - lBF.PMcurrent) #more stable
      if (runif(1) < ratio) { #update current model
        current.model[j] <- proposal.model[j]
        lBFcurrent <- lBFproposal; lBF.PMcurrent <- lBF.PMproposal
      }

      if(i > 1) {
        inclprobRB[i,j] <- inclprobRB[i-1, j] + proposal.model[j] * ratio +
                           (1 - proposal.model[j]) * (1 - ratio)
      }
    }

    all.models.lBF[i,] <-  c(current.model, lBFcurrent)
    all.lBF.PM[i] <- lBF.PMcurrent
  }
  cat("\n")
  for(j in seq_len(p)) inclprobRB[,j] <- inclprobRB[,j] / seq(1,(n.iter + n.burnin))

  if (n.burnin > 0) { #remove burnin
    seqburn <- seq_len(n.burnin)
    all.models.lBF <- all.models.lBF[-seqburn,]
    all.lBF.PM <- all.lBF.PM[-seqburn]
  }

  #keep 1 each n.thin iterations
  seqthin <- seq(1, n.iter, by = n.thin)
  all.models.lBF <- all.models.lBF[seqthin, ]
  all.lBF.PM <- all.lBF.PM[seqthin]
  colnames(all.models.lBF) <- c(names(NAvars), "logBF")

  if (mF) {
    #models matrix at the covariate-factor level with number of active dummies
    cf.models.lBF <- all.models.lBF[,seq_len(p)] %*% t(positions)
    cf.models.lBF <- cbind(cf.models.lBF, all.models.lBF[, p+1])
    colnames(cf.models.lBF)[ncol(cf.models.lBF)] <- "logBF"

    inclprobRB <- inclprobRB %*% t(positions)
  } else cf.models.lBF <- all.models.lBF
  dimnames(cf.models.lBF) <- list(1:nrow(cf.models.lBF), colnames(cf.models.lBF))
  #cf.models.lBF is exactly all.models.lBF if there are no factors

  return(list(cf.models.lBF = cf.models.lBF, all.models.lBF = all.models.lBF,
              all.lBF.PM = all.lBF.PM, inclprobRB = inclprobRB))
}

#' @keywords internal
summ.Gibbs <- function (gibbs, q, n.iter) {
  #Summ up Gibbs sampling results

  cf.models.lBF <- gibbs$cf.models.lBF
  all.lBF.PM <- gibbs$all.lBF.PM
  inclprobRB <- gibbs$inclprobRB

  #inclusion probabilities except for fixed variables:
  inclprob <- colMeans(cf.models.lBF[,-(q+1)] > 0)

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
  nPmax <- which.max(all.lBF.PM)
  hpm <- cf.models.lBF[nPmax, -(q+1)]

  #MPM
  mpm <- numeric(q)
  mpm[which(inclprobRB[n.iter, ] >= 0.5)] <- 1

  return(list(post = post, C = C, inclprob = inclprob, probdim = probdim,
              hpm = hpm, mpm = mpm))
}

#' @keywords internal
checkGibbsarguments <- function (p, p0, namesnull, namesx, init.model, mF,
                                 positionscov, positionsfac, l) {
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

    d <- as.vector(positionsfac %*% init.model) #levels active of factors
    f <- d > 0 #active factors
    cf <- c(positionscov %*% init.model, f) #covariates and/or factors active

    if (any(f)) {
      checkoversat <- which(d == l) #check if oversaturated model
      checksat <- which(d == l - 1) #check if saturated model and not representative

      if (length(c(checkoversat, checksat)) > 0) {
        for (j in c(checkoversat, checksat)) {
          init.model[as.logical(positionsfac[j,])] <- c(0, rep.int(1,l[j]-1))
        }
      }
    }
  }

  return(init.model)
}
