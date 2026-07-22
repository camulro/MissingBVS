#' Bayesian Imputation Averaging for Variable Selection with Missing data in linear regression models
#'
#' Computation and summaries of posterior distribution over the model space in problems
#' of small to moderate size in the presence of (possible) missing data and/or categorical
#' variables in linear models. Each posterior model probability is computed following the
#' Bayesian Imputation Averaging (BIA) framework, using standard priors for model coefficients
#' and the hierarchical approach of García-Donato and Paulo (2022) with factors.
#'
#' The set of competing models is made up by all the possible subsets of regressors
#' specified by \code{formula}: Mi for i in 1,...,2^p, being p the number of potential
#' (non-fixed) regressors in the variable selection problem. It is assumed that the
#' intercept term is present in all models. The simplest one M0, the \code{null.model}
#' nested in the rest, contains the fixed variables, if given, and only the intercept by default.
#' In order to implement BIA, \code{\link[MissingBVS]{missingBVS.lm}} can, either perform
#' \code{n.imp} imputations designed by \code{imp.predict.mat} and \code{imp.mice.method}
#' with the \pck{mice} package, or use user-given imputated datasets by the
#' \code{imp.datasets} argument. Hence, the posterior distribution over the model space
#' is given through Bayes' theorem:
#'
#' Pr(Mi | \code{data}) = Pr(Mi) * AvBi / C,
#'
#' where AvBi is the Average Bayes factor (AvBF) of Mi to M0 under missing data,
#' Pr(Mi) is the prior probability of Mi and C is the normalizing constant.
#' AvBi is an actual Bayes factor (BF) and it is defined as the average of the
#' \code{n.imp} data-driven BFs:
#'
#' AvBi = 1/\code{n.imp} * (Bi(1) + ... + Bi(\code{n.imp})),
#'
#' where Bi(j) corresponds to the BF for model Mi to M0 under the jth imputed dataset.
#' Data-driven BF can be either computed using popular g-prior choices or approximated
#' with the BIC (Schwarz, 1978) or the test-based BF (Held, Gravestock and Sabanés, 2015)
#' with the \code{BF.approx.method} argument.
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and are computed using \pck{BayesVarSel}. The choices currently available are:
#' -"Robust" is the default option and denotes the criteria-based prior of Bayarri,
#' Berger, Forte and Garcia-Donato (2012).
#' -"gZellner" corresponds to the prior in Zellner (1986) with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"ZellnerSiow" is the multivariate Cauchy prior by Zellner and Siow (1980, 1984).
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015).
#' -"IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato,
#' Moreno and Pericchi (2022).
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
#' In the presence of factors, in order to make results do not dependent on codification
#' of factors Pr(Mi) is factorized following García-Donato and Paulo (2022). Each factor is
#' represented by its rank defficient dummy parametrization, making each model with a factor
#' active be the representation of the model space given by the 2^l-l different submodels with
#' at least one dummy active, where l is the number of levels. This approach derives in a
#' hierarchical prior where the part over the covariates and factors is a standard model
#' prior given by \code{prior.models}. The prior over the submodels defined by the dummies,
#' assumes a prior independence between factors and within a factor two options are available:
#' "Constant" assigns the same prior to every submodel and "ScottBerger" to each model size,
#' which is the recommended. Within all these priors, the prior inclusion probabilities
#' of factors and numerical variables are 1/2, which does not happen when the selection is
#' directly over dummies. A non-treatment of factors can be performed through \code{marginal.factors}.
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
#' @param n.keep It can be either the character "all" to return the whole model
#' space or a numeric for the exact number of the most probable models to keep.
#' By default it is set to 10 and automatically adjusted if 10 is greater than
#' the total number of models.
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
#' @param imp.seed Seed for imputation.
#'
#' @return \code{missingBVS.lm} returns an object of class \code{missingBVS}
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
#' \item{MPMbin}{Binary expression of the Median Probability model}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{modelsrankdefprob}{A \code{n.keep} x (p1 * (sum_j l_j - L)+1) matrix which
#' summarizes the \code{n.keep} most probable a posteriori submodels and their
#' probabilities (when relevant)}
#' \item{modelsprob}{A \code{n.keep} x (p+1) matrix which summarizes the \code{n.keep}
#' most probable a posteriori models and their associated probability}
#' \item{inclprob}{Named vector with the inclusion probabilities of p competing variables}
#' \item{postprobdim}{Posterior probabilities over the true model size}
#' \item{C}{The value of the normalizing constant C=sum_i AvBi * Pr(Mi)}
#' \item{call}{The \code{call} to the function}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{imp.info}{List of arguments used for the imputation step and other information
#' (when relevant)}
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
#' \item{method}{String "Full" denoting exhaustive model search}
#'
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGibbsBvs.lm}} for a heuristic
#' approximation based on Gibbs sampling (recommended when p>20).
#'
#' Consider \code{\link[MissingBVS]{plot.MissingBvs}} for graphical summaries of the
#' posterior distribution.
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Garcia-Donato, G. and Paulo, R. (2022)<DOI:10.1080/01621459.2021.1889565>
#' Variable Selection in the Presence of Factors: A Model Selection Perspective.
#' Journal of the American Statistical Association. 117. 1-27.
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
#' Statistics. 6: 461–464.
#'
#' Held, L., Gravestock, I. and Sabanés Bové, D.
#' (2015)<DOI:10.1080/01621459.2014.993077> Objective Bayesian model selection
#' for generalized linear models using test-based Bayes factors. Journal of the
#' American Statistical Association, 110, 1157–1168.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @examples
#'
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Here we keep the 8 competing models:
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(formula = f, data = dataS97, n.keep = 8)
#'
#' #Show the results:
#' dataS97.mBVS
#'
#' #Summ up the results:
#' summary(dataS97.mBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(dataS97.mBVS)
#' }
#'
missingBVS.lm <- function (formula,
                           data,
                           null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                           BF.approx.method = "gprior",
                           prior.betas = "Robust",
                           prior.models = "ScottBerger",
                           prior.models.dummies = "ScottBerger",
                           marginal.factors = TRUE,
                           priorprobs = NULL,
                           n.keep = 10,
                           imp.mice.method = "pmm",
                           imp.predict.mat = NULL,
                           n.imp = 039E1,
                           maxit = 5,
                           parallelmice = NULL,
                           n.core = NULL,
                           imp.datasets = NULL,
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

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(p, p0, namesnull, namesx, n.keep, q)

  #Check model priors chosen and define the functions to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  mF <- L > 0 & marginal.factors
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

  cat("The problem has a total of", 2^q, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #Compute exact posterior distribution and normalizing constant
  posterior.list <- exact.posterior.comput(p, namesxnotnull, NAvars, lBF.method,
                                           lp.model, X0, X.full, BF.approx.method)
  list2env(posterior.list, envir = env)

  #Summ up the posterior distribution
  summ.posterior.list <- summ.posterior(all.models.PM, p, q, mF, positions)
  list2env(summ.posterior.list, envir = env)

  if (anyNAvar) {#Pool results for imputed datasets
    imp.array <- imputation.array
    #if marfinal factor probs, remove first dummy on each factor, first p0 are the fixed ones
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

  ##result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$lmfull <- lmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, lmfull is the lm object for the full model
  result$lmnull <- lmnull # The lm object for the null model (omits NAs)

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
    result$modelsrankdefprob <- all.models.PM # rank deficient models and probs
  }

  #The binary code for the n.keep best models and the correspondent post
  result$modelsprob <- modelsprob[order(modelsprob[,q+1],
                                        decreasing = TRUE)[seq_len(n.keep)],]

  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- depvars

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model
  result$C <- C #normalizing constant

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
  names(result$priorprobs) <- 0:q + q0 #prior dimension probability

  if (anyNAvar) {
    #arguments used for imputation
    result$imp.info <- imp.info

    #save the imputed datasets for BMA or sensitivity analysis
    raw.imp.array <- serialize(imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

    ME(n.imp, imp.seed)
  }

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lp.model #function used for model prior
  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models
  result$marginal.factors <- marginal.factors #whether or not factors are marginalized

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}

#' @keywords internal
exact.posterior.comput <- function (p, namesxnotnull, NAvars, lBF.method, lp.model,
                                    X0, X.full, BF.approx.method) {
  #Compute exact posterior distribution and normalizing constant

  #progress bar for loop
  pb <- txtProgressBar(min = 0, max = 2^p, style = 3, width = 50, char = "=")

  #Posterior computation
  all.models.lPM <- matrix(0, nr = 2^p, nc = p+1) #last column contains log(BF_a0*Pr(M))
  for (i in seq_len(2^p-1)){ # null out of the loop
    setTxtProgressBar(pb, i)

    #transform the number of the model into a binary number
    current.model <- BayesVarSel:::integer.base.b_C(i, p)
    all.models.lPM[i, seq_len(p)] <- current.model

    lpm <- lp.model(current.model) #log-model prior
    if (is.na(lpm)) next #do not visit saturated or oversaturated models

    #check if there are NAs in the model considered to save computation time
    if (any(namesxnotnull[which(current.model == 1)] %in% NAvars)) {
      lBF.PM <- lBF.method(model = which(current.model == 1)) + lpm #log(BF_a0*Pr(M))
    } else { #if there are no missings, compute the BF by the method selected
      X.i <- cbind(X0, X.full[, which(current.model == 1)])
      lBF.PM <- BF.approx.method(k = sum(current.model == 1), X = X.i) + lpm #log(BF_a0*Pr(M))
    }
    all.models.lPM[i, p+1] <- lBF.PM
  }
  setTxtProgressBar(pb, 2^p)
  cat("\n")
  #null model
  all.models.lPM[2^p, seq_len(p)] <- numeric(p)
  all.models.lPM[2^p, p+1] <- lp.model(numeric(p)) #BF = 1 for null model
  all.models.lPM <- na.omit(all.models.lPM) #remove repeated models if dummies

  #renormalize
  C <- sum(exp(all.models.lPM[, p+1]))
  all.models.PM <- all.models.lPM
  all.models.PM[, p+1] <- exp(all.models.lPM[, p+1] - log(C))
  colnames(all.models.PM) <- c(namesxnotnull, "Post")

  return(list(all.models.PM = all.models.PM, C = C))
}

#' @keywords internal
summ.posterior <- function (all.models.PM, p, q, mF, positions) {
  #Summ up the posterior distribution

  if (mF) {
    #compute models matrix at the covariate-factor level
    cf.models.PM <- all.models.PM[,seq_len(p)] %*% t(positions)
    cf.models.PM <- cbind(cf.models.PM, all.models.PM[,p+1])
    colnames(cf.models.PM)[q+1] <- "Post"

    modelsprob <- cbind(t(sapply(seq_len(2^q)-1,
                                 function(j) BayesVarSel:::integer.base.b_C(j, q))), numeric(2^q))
    for (i in seq_len(nrow(cf.models.PM))) {
      modeli <- cf.models.PM[i, seq_len(q)] > 0
      j <- which(apply(modelsprob[, seq_len(q)], 1, function(x) all(x == modeli)))
      modelsprob[j, q + 1] <- modelsprob[j, q + 1] + cf.models.PM[i, q + 1]
    }
    colnames(modelsprob) <- colnames(cf.models.PM)
  } else modelsprob <- all.models.PM

  inclprob <- numeric(q)
  probdim <- numeric(q+1)
  #compute inclusion probabilities (except for fixed variables) and
  #posterior probability of the dimension of the true model
  for (i in seq_len(nrow(modelsprob))) {
    inclprob[which(modelsprob[i, seq_len(q)] == 1)] <-
      inclprob[which(modelsprob[i, seq_len(q)] == 1)] + modelsprob[i, q + 1]

    probdim[sum(modelsprob[i, seq_len(q)] == 1) + 1] <-
      probdim[sum(modelsprob[i, seq_len(q)] == 1) + 1] + modelsprob[i, q + 1]
  }

  #HPM
  nPmax <- which.max(modelsprob[, q+1])
  hpm <- modelsprob[nPmax, ]

  #MPM
  mpm <- numeric(q)
  mpm[which(inclprob >= 0.5)] <- 1

  return(list(modelsprob = modelsprob,
              inclprob = inclprob, probdim = probdim, hpm = hpm, mpm = mpm))
}

#' @keywords internal
buildmatrices <- function (formula, null.model, data, marginal.factors) {
  #Build data matrices and objects to use later

  #Response and fixed vars for imputation
  framenull <- model.frame(null.model, data, na.action = NULL)
  q0 <- dim(framenull)[2] # number of fixed covariates and/or factors,intercept always

  #Missing model matrix of fixed vars
  X0 <- model.matrix(framenull, data)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of fixed covariates or dummies of factors

  #Full design matrix given by formula
  framefull <- model.frame(formula, data, na.action = NULL)
  if (marginal.factors) {
    #Rank deficient fixed model matrix just to remove vars from the full one
    X0rdf <- model.matrix.rankdef(framenull)
    #Rank deficient full model matrix data with missings
    X.full <- model.matrix.rankdef(framefull)
    namesx <- dimnames(X.full)[[2]]
    #Only the non-fixed vars
    namesxnotnull <- setdiff(namesx, dimnames(X0rdf)[[2]])
  } else {
    X.full <- model.matrix(formula, framefull)
    namesx <- dimnames(X.full)[[2]]
    #Only the non-fixed vars
    namesxnotnull <- setdiff(namesx, dimnames(X0)[[2]])
  }
  X.full <- X.full[, namesxnotnull]
  p <- length(namesxnotnull) #Number of covariates and levels of factors to select from

  #the order for the posterior model distribution computation step
  ordvars <- c(namesnull, namesxnotnull) #X0, X.full

  if (marginal.factors) {
    #For factors:
    #covariates and/or factors to select from
    depvars <- setdiff(attr(terms(framefull), "term.labels"),
                       attr(terms(framenull), "term.labels"))

    #positions has number of rows equal to the number of regressors and p columns.
    #A 1 in a row denotes the position in X of a regressor (several positions for
    #the dummies of a factor).
    positions <- t(sapply(depvars, function(var) {
      if(is.factor(data[[var]])) {
        levs <- levels(data[[var]])
        ind <- which(namesxnotnull %in% paste0(var,levs)) #1 if the namelevel matches
      } else ind <- which(namesxnotnull == var) #1 if the name matches

      posi <- numeric(p); posi[ind] <- 1; posi
    }))
    colnames(positions) <- namesxnotnull

    tmp <- colSums(positions %*% t(positions))
    positionsx <- tmp == 1 #vector of length p with TRUE if numeric variable

    L <- sum(!positionsx) #Number of factors to select from
    if (L > 0) {
      #matrix of dim (Lxp) with 1 if dummy variable of the row factor
      positionsfac <- positions[!positionsx, , drop = FALSE]
      l <- tmp[tmp > 1] #Number of levels for each factor

      #vector of length L with the position of the first dummy in each factor to check for repeated models
      indf <- apply(positionsfac, MARGIN = 1, FUN = function(x) head(which(x == 1), n = 1))
      firstd <- ifelse(L > 1, function(M) diag(M[, indf]),  function(M) M[, indf])

      q <- p - sum(l) + L #Number of factors and covariates to select from
      #q = p if there are no factors

      return(list(q0 = q0, p0 = p0, X0 = X0, namesnull = namesnull, framenull = framenull,
                  q = q, p = p, X.full = X.full, namesxnotnull = namesxnotnull, namesx = namesx,
                  framefull = framefull, ordvars = ordvars, depvars = depvars,
                  positions = positions, positionsx = positionsx, positionsfac = positionsfac,
                  L = L, l = l, indf = indf, firstd = firstd))
    }
  }

  return(list(q0 = p0, p0 = p0, X0 = X0, namesnull = namesnull, framenull = framenull,
              q = p, p = p, X.full = X.full, namesxnotnull = namesxnotnull, namesx = namesx,
              framefull = framefull, ordvars = ordvars, depvars = namesxnotnull, L = 0))
}

#' @keywords internal
buildimputation <- function(NAvars, formula, data, imp.predict.mat, n.imp, maxit,
                            n, q, p0, imp.mice.method, imp.seed,
                            parallelmice, n.core, obsnotNA, ordvars, BF.approx.method) {
  #Build imp.predict matrix to imputation, imputed datasets and BF function

  #Impute just competing variables with NAs
  fulldataframe <- model.frame(paste0(formula[[2]], "~."), data, na.action = NULL)
  X.toimp <- fulldataframe[,-1] #full observed design matrix
  #Default prediction matrix by mice:
  quickpredict.mat <- mice::quickpred(X.toimp); Xnames <- colnames(quickpredict.mat)
  if (!is.null(imp.predict.mat)) { #if given by user
    #check predict imputation matrix
    imp.vars <- rownames(imp.predict.mat)
    if (any(NAvars %notin% imp.vars)) { #all formula predictors have to be imputed
      stop("Imputation prediction matrix rows given do not contain all the variables ",
           "given by formula with NAs.", "Make sure to include them all.\n")
    }
    #check row names
    if (any(imp.vars %notin% Xnames)) stop("Row variables in imp.predict.mat not found in data.\n")

    imp.pred.col <- colnames(imp.predict.mat)
    #check column names
    if (any(imp.pred.col %notin% Xnames)) stop("Column variables in imp.predict.mat not found in data.\n")

    #select variables to impute from columns of imp.predict.mat
    quickpredict.mat[imp.vars, imp.pred.col] <- imp.predict.mat
    #set to 0 the ones do not selected by user
    quickpredict.mat[imp.vars, Xnames %notin% imp.pred.col] <- 0
  } else imp.pred.col <- colnames(quickpredict.mat)

  #Do not impute variables with missings that are not in imp.pred.col (also in NAvars)
  not.imp.vars <- setdiff(colnames(X.toimp), imp.pred.col)
  quickpredict.mat[not.imp.vars, ] <- 0

  #visit sequence given by order of rows in imp.predict.mat, if given
  visit.seq <- c(imp.pred.col, not.imp.vars)

  #Check parallel arguments
  if (is.null(parallelmice)) {
    if (n.imp > 120 | (n*q > 50000 & n.imp > 5)) {
      parallelmice <- TRUE #faster
    } else parallelmice <- FALSE
  }

  if (n*q > 10000 | n.imp > 039E1) if (!parallelmice) {
    cat("Do you want to faster imputation running a parallel version of mice? (y/n)\n")
    if (tolower(readline()) == "y") {
      parallelmice <- TRUE
    } else cat("Be aware that imputation could take a while.\n")
  }

  #Imputation of missing data
  cat("Performing imputation of missing data with mice's", imp.mice.method)
  if (parallelmice) cat(" parallel")
  cat(" method.\n", "Please wait . . . \n")

  imput <- mice.imputation(X = X.toimp, formula,
                           n.imp = n.imp,
                           imp.predict.mat = quickpredict.mat,
                           imp.mice.method = imp.mice.method,
                           visit.seq = visit.seq,
                           seed = imp.seed,
                           maxit = maxit,
                           parallel = parallelmice,
                           n.core = n.core)

  #remove observations with missings on the response or fixed vars,
  #select the vars in the order X0, X.full and remove oversaturated for X0 if factors
  imputation.array <- imput$imputation.array[obsnotNA, ordvars, , drop = FALSE]
  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model as an average of BF computed
    #by BF.approx.method over the imputed datasets
    lBF.method <- function (model) lBF.approx(model,
                                              imputation.array = imputation.array,
                                              BF.approx.method = BF.approx.method,
                                              p0 = p0, n.imp = n.imp)
  } else lBF.method <- function (model) BF.approx.method(k = length(model),
                                                         X = imputation.array[,c(1:p0, model+p0),])

  imp.info <- list(loggedEvents = imput$logEvents, parallelmice = parallelmice, imp.mice.method = imp.mice.method,
                   imp.predict.mat = imp.predict.mat, n.imp = n.imp, imp.seed = imp.seed, NAvars = NAvars)

  return(list(imputation.array = imputation.array, lBF.method = lBF.method, imp.info = imp.info))
}

#' @keywords internal
extimputation <- function (formula, imp.datasets, n0, framefull, ordvars, obsnotNA,
                           p0, BF.approx.method, NAvars) {
  X.formula <- as.formula(paste(formula[1], formula[3]))
  isarray <- is.array(imp.datasets)
  islist <- is.list(imp.datasets)
  if (!isarray & !islist) {
    stop("Imputations should be given as an array or list.\n")
  }

  aux <- model.matrix.rankdef(framefull)
  if (isarray) {
    #Check that imputations have the correct format
    if (dim(imp.datasets)[1] != n0) stop("Imputations should be given as an (",
                                         n0, "xn.varsxn.imp) array.\n")

    #Check column names given
    if (is.null(colnames(imp.datasets))) {
      stop("imp.datasets should have dimension names (the variables data names).\n")
    }
    #Check that each var in formula is given by imp.datasets
    varsnotinimp <- colnames(framefull)[-1] %notin% colnames(imp.datasets)
    if (any(varsnotinimp))  stop("Variables: ",
                                 paste0(colnames(framefull)[varsnotinimp], collapse = ", "),
                                 "; not given by imp.datasets.\n")

    n.imp <- dim(imp.datasets)[3] #number of imputed datasets
    #Build rank deficient matrices:
    imputation.array <- array(0, dim = c(n0, ncol(aux), n.imp), #an array with the matrices imputed
                              dimnames = list(seq_len(n0), colnames(aux), seq_len(n.imp)))

    for (s in seq_len(n.imp)) {
      aux.imps <- model.frame(X.formula, data.frame(imp.datasets[,,s]), na.action = NULL)
      imputation.array[,,s] <- model.matrix.rankdef(aux.imps) #build the model matrix
    }
  }

  if (islist) {
    if (dim(imp.datasets[[1]])[1] != n0) stop("Imputations should be given as a list of (",
                                              n0, "xn.vars) matrices.\n")

    #Check column names given
    if (is.null(colnames(imp.datasets[[1]]))) {
      stop("imp.datasets should have dimension names (the variables data names).\n")
    }
    #Check that each var in formula is given by imp.datasets
    varsnotinimp <- colnames(framefull)[-1] %notin% colnames(imp.datasets[[1]])
    if (any(varsnotinimp)) stop("Variables: ",
                                paste0(colnames(framefull)[varsnotinimp], collapse = ", "),
                                "; not given by imp.datasets.\n")

    n.imp <- length(imp.datasets) #number of imputed datasets
    #Build rank deficient matrices:
    imputation.array <- array(0, dim = c(n0, ncol(aux), n.imp), #an array with the matrices imputed
                              dimnames = list(seq_len(n0), colnames(aux), seq_len(n.imp)))

    for (s in seq_len(n.imp)) {
      aux.imps <- model.frame(X.formula, data.frame(imp.datasets[[s]]), na.action = NULL)
      imputation.array[,,s] <- model.matrix.rankdef(aux.imps) #build the model matrix
    }
  }

  #remove observations with missings on the response or fixed vars
  imputation.array <- imputation.array[obsnotNA, ordvars, , drop = FALSE]

  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model as an average of BF computed
    #by BF.approx.method over the imputed datasets
    lBF.method <- function (model) lBF.approx(model,
                                              imputation.array = imputation.array,
                                              BF.approx.method = BF.approx.method,
                                              p0 = p0, n.imp = n.imp)
  } else lBF.method <- function (model) BF.approx.method(k = length(model),
                                                         X = imputation.array[,c(1:p0, model+p0),])

  imp.info <- list(n.imp = n.imp, NAvars = NAvars)

  return(list(imputation.array = imputation.array, lBF.method = lBF.method,
              imp.info = imp.info, n.imp = n.imp))
}


"%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0 #auxiliar function

#' @keywords internal
checkBvsarguments <- function (p, p0, namesnull, namesx, n.keep, q) {
  #check arguments
  #Is there any variable to select from?
  if (p == 0) stop("The number of fixed variables is equal to the number of\n",
                   "regressors in the full model. No model selection can be done.\n")

  #check if the number of regressors is too big.
  if (p > 27) stop("Number of regressors too big. . . Please, use the Gibbs method.\n")
  if (p > 20) warning("Number of regressors too big. . . Consider using the Gibbs method.\n",
                      immediate. = TRUE)

  #check if null model is contained in the full one:
  for (i in 1:p0) {
    if (namesnull[i] %notin% namesx) stop("Error in var: ", namesnull[i],
                                          "; null model not nested in full model.\n")
  }

  #n.keep > 2^q, correct the number of models to keep
  if(is.character(n.keep)) {
    if (n.keep == "all") {
      # n.keep <- 2^(q - L)*prod(2^l - l)
      n.keep <- 2^q
    } else stop("Only n.keep='all' or type the exact number of models to keep instead.\n")
  }
  if (n.keep > 2^q) {
    cat("The number of models to keep (", n.keep,
        ") is larger than the total number of models (", 2^q,
        ") and it has been set to ", 2^q,".\n")
    # n.keep <- 2^(q - L)*prod(2^l - l)
    n.keep <- 2^q
  }
  return(n.keep)
}

#' @keywords internal
checkformissings <- function (y, X0 = NULL, X.full) {
  #checks if there are missings on the response and regressors
  #and returns the name of non-fixed regressors with missings

  ##on the response
  if (sum(is.na(y)) > 0) cat("NA values found on the response variable.",
                             "We are going to omit these observations.\n")

  ##on the fixed vars
  if (sum(is.na(X0)) > 0)  cat("NA values found on the fixed variables.",
                               "We are going to omit these observations.\n")

  ##on the regressors
  if (sum(is.na(X.full)) > 0) {
    O <- 1*(!is.na(X.full))
    #zeros where missing observations on the data without missings on the response
    NAvars <- names(which(colSums(O) < dim(X.full)[1])) #columns with missings
  } else NAvars <- NULL

  return(NAvars)
}

#' @keywords internal
checkforprior.models <- function (prior.models, priorprobs, p) {
  #checks that the model prior given by prior.models is implemented and returns
  #the function to use for model prior computation

  switch (prior.models,
          ScottBerger = {prior.models.f <- function(model) logScottBerger(p = p, model)},
          Constant = {prior.models.f <- function(model) logConstant(p = p)},
          User = {
            if (is.null(priorprobs)) stop("User prior selected but no prior probabilities provided.\n")
            if (length(priorprobs) != (p + 1)) stop("User prior selected but the length of prior",
                                                    "probabilities is not correct (", p+1,").\n")
            if (sum(priorprobs < 0) > 0) stop("Prior probabilities must be positive.\n")

            prior.models.f <- function(model) logUser(p = p, model, priorprobs = priorprobs)
          },
          stop("Only priors 'ScottBerger', 'Constant' and 'User' supported.\n"))

  return(prior.models.f)
}

#' @keywords internal
checkmarg.factorsprior <- function (mF, prior.models.dummies, l, positions, positionsfac,
                                    firstd, lprior.models) {
  #returns the function to use for model prior computation
  if (mF) {
    lprior.models.dummies <- checkforprior.models.dummies(prior.models.dummies, l)

    lp.model <- function (model) {
      gt <- positions %*% model > 0 #covariates and/or factors active
      d <- positionsfac %*% model #levels active of factors
      t <- d > 0 #active factors

      #check if the model is one among the saturated and oversaturated due to the dummies
      if (sum(t) > 0) {
        M <- t(apply(positionsfac, 1, function(x) x * model))
        #to avoid computation of saturated models
        if (any((d == l) | ((d == l - 1) & firstd(M)))) return(NA)
      }

      return(lprior.models(gt) + lprior.models.dummies(d, t))
    }
  } else lp.model <- function (model) lprior.models(model)

  return(lp.model)
}

#' @keywords internal
checkforprior.models.dummies <- function (prior.models.dummies, l) {
  #checks that the model prior given by prior.models.dummies is implemented and
  #returns the function to use for model prior computation

  switch (prior.models.dummies,
          ScottBerger = {prior.models.f <-
            function(delta, tau) logScottBerger.d(delta, tau, l = l)},
          Constant = {prior.models.f <- function(delta, tau) logConstant.d(tau, l = l)},
          stop("Only priors 'ScottBerger' and 'Constant' supported.\n"))

  return(prior.models.f)
}

#' @keywords internal
checkforprior.betas.lm <- function (BF.approx.method, prior.betas, n, p, p0, y, SS0) {
  #checks that the Bayes factor computation method given by BF.approx.method and prior.betas
  #is implemented and returns the function to use for Bayes factor computation on lm
  if (BF.approx.method %notin% c("BIC", "TBF", "gprior")) {
    stop("Only BF approximations 'BIC', 'TBF' and 'gprior' supported.")
  }

  if (BF.approx.method == "gprior") {
    switch (prior.betas, # change the string for the corresponding tag in BayesVarSel code
            gZellner = {prior.betas <- "gBF"}, #fixed g=n
            Robust = {prior.betas <- "RobustBF"}, #random g: criteria-based prior from Bayarri et al (2012)
            Liangetal = {prior.betas <- "LiangBF"}, #random g: hyper-g/n with a=3
            `Zellner-Siow` = {prior.betas <- "ZSBF"}, #random g: cauchy prior
            FLS = {prior.betas <- "flsBF"}, #fixed g Benchmark prior: g=max(n, p*p)
            `intrinsic.MGC` = {prior.betas <- "intrinsicBF"}, #intrinsic prior from Moreno, Giron, Casella (2015)
            IHG = {prior.betas <- "geointrinsicBF"}, #intrinsic hyper-g prior
    stop("prior.betas must be one of 'gZellner', 'Robust', 'Liangetal', 'ZellnerSiow',\n",
         "'FLS', 'intrinsic.MGC' or  'IHG' when using BF.approx.gprior method.\n"))
  }

  switch (BF.approx.method,
          BIC = {BF.approx.method.f <-
            function (k, X) BF.approx.BIC.lm(y = y, X, SS0 = SS0, n = n, k, p0 = p0)},
          TBF = {BF.approx.method.f <-
            function (k, X) BF.approx.TBF.lm(y = y, X, SS0 = SS0, prior.betas = prior.betas,
                                             n = n, k, p0 = p0)},
          gprior = {BF.approx.method.f <-
            ifelse(prior.betas != "flsBF",
                   function (k, X) BF.approx.gprior.lm(y = y, X, SS0 = SS0, prior.betas = prior.betas,
                                                       n = n, k, p0 = p0),
                   function (k, X) BF.approx.FLS.lm(y = y, X, SS0 = SS0, dmax = p + p0,
                                                    n = n, k, p0 = p0))}
  )
  return(BF.approx.method.f)
}

#' @keywords internal
model.matrix.rankdef <- function (model.frame.aux) {
  #internal function to create rank defficient matrices from a given dataframe
  #created from a model.frame call
  if (ncol(model.frame.aux) == 1) { #just the response
    Xnull.def <- cbind(`(Intercept)` = rep.int(1, nrow(model.frame.aux)))
    return(Xnull.def)
  }

  terms <- attr(terms(model.frame.aux), "term.labels")

  Xi.rdef <- sapply(terms, function(var) {
    f <- as.formula(paste0("~ 0 + ", var))
    model.matrix(f, data = model.frame.aux)}, simplify = FALSE)
  Xfull.def <- do.call(cbind, Xi.rdef)
  Xfull.def <- cbind(`(Intercept)` = rep.int(1, nrow(Xfull.def)), Xfull.def)
  return(Xfull.def)
}

#' Print an object of class \code{MissingBvs}
#'
#' Print an object of class \code{MissingBvs}. Top ten models with the highest
#' probabilitis are shown jointly with corresponding log-Bayes factors and
#' posterior probabilities (or an estimation if the object was created
#' by a Gibbs function).
#'
#' @export
#' @param mbvs.object An object of class \code{MissingBvs}.
#' @param ... Additional parameters to be passed.
#'
#' @author Gonzalo Garcia-Donato and Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso  Use \code{\link[MissingBVS]{MissingBvs.lm}},
#' \code{\link[MissingBVS]{MissingGD25}} or \code{\link[MissingBVS]{MissingBvs.glm}}
#' and their Gibbs versions for creating objects of the class \code{MissingBvs}.
#'
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Here we keep the 8 competing models:
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(formula = f, data = dataS97, n.keep = 8)
#'
#' #Show the results:
#' print(dataS97.mBVS)
#' }
#'
print.MissingBvs <- function(mbvs.object,...){

  if (!inherits(mbvs.object, "MissingBvs")){
    warning("An object of class MissingBvs is needed.\n")
  }

  cat("\nCall:\n")
  print(mbvs.object$call)

  if (mbvs.object$method == "Gibbs") {
    p <- mbvs.object$p
    postprob <- mbvs.object$postprobs

    ord <- order(postprob, decreasing = T)
    modelspostprob <- cbind(mbvs.object$modelslogBF[ord,], postprob[ord])
    modelspostprob <- modelspostprob[!duplicated(modelspostprob),]

    n.keep <- min(dim(modelspostprob)[1], 10)
    mod.mat <- as.data.frame(modelspostprob[1:n.keep,])
    colnames(mod.mat) <- c(colnames(mbvs.object$modelslogBF), "Post. prob.")

    cat("\nThe ", n.keep, " most probable models among the visited ones are:\n")
    print(mod.mat)
    cat("---\n")
    cat("Code: Column logBF is the log of Bayes factor and\n")
    cat("column post. prob. is an estimation of posterior probabilities \n")
    cat("based on the normalizing constant.")

  }

  if (mbvs.object$method == "Full") {
    n.keep <- min(dim(mbvs.object$modelsprob)[1], 10)

    cat("\nThe", n.keep, "most probable models and their probabilities are:\n", sep=" ")
    print(mbvs.object$modelsprob[1:n.keep, ])
  }
  cat("\n\n")
}

#' Summary of an object of class \code{MissingBvs}
#'
#' Summary of an object of class \code{MissingBvs}, providing inclusion
#' probabilities and a representation of the Median Probability Model (MPM) and the
#' Highest Posterior probability Model (HPM).
#'
#' @export
#' @param object An object of class \code{MissingBvs}.
#' @param ... Additional parameters to be passed.
#'
#' @author Gonzalo Garcia-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso  Use \code{\link[MissingBVS]{MissingBvs.lm}},
#' \code{\link[MissingBVS]{MissingGD25}} or \code{\link[MissingBVS]{MissingBvs.glm}}
#' and their Gibbs versions for creating objects of the class \code{MissingBvs}.
#'
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Here we keep the 8 competing models:
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(formula = f, data = dataS97, n.keep = 8)
#'
#' #Summ up the results:
#' summary(dataS97.mBVS)
#' }
#'
#' @references Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238>
#' Optimal Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
summary.MissingBvs <- function(object,...){

  if (!inherits(object, "MissingBvs")){
    warning("calling summary.MissingBvs(<fake-MissingBvs-x>) . . . ")
  }

  p <- object$p
  inclprob <- object$inclprob
  HPM <- ifelse(object$HPMbin[1:p], "*", "")
  MPM <- ifelse(object$MPMbin, "*", "")
  summ.missingBvs <- as.data.frame(cbind(round(inclprob, digits = 4), HPM, MPM))
  names(summ.missingBvs) <- c("Incl.prob.", "HPM", "MPM")

  ans <- list()
  ans$summary <- summ.missingBvs
  ans$method <- object$method
  ans$call <- object$call

  cat("\nPosterior Inclusion Probabilities:\n")
  print(ans$summary)
  cat("---\n")
  cat("Code: HPM stands for Highest posterior Probability Model and\n")
  cat("      MPM for Median Probability Model.\n ")
  if (object$method == "Gibbs") {
    cat("Results are estimates based on the visited models.\n")
  }
  class(ans) <- "summary.MissingBvs"
  return(invisible(ans))
}

