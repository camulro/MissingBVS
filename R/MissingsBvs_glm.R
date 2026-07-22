#' Bayesian Imputation Averaging for Variable Selection with Missing data in generalized linear models
#'
#' Computation and summaries of posterior distribution over the model space in problems
#' of small to moderate size in the presence of (possible) missing data and/or categorical
#' variables in generalized linear models. Each posterior model probability is computed
#' following the Bayesian Imputation Averaging (BIA) framework, using standard priors
#' for model coefficients and the hierarchical approach of García-Donato and Paulo (2022)
#' with factors.
#'
#' The set of competing models is made up by all the possible subsets of regressors
#' specified by \code{formula}: Mi for i in 1,...,2^p, being p the number of potential
#' (non-fixed) regressors in the variable selection problem. It is assumed that the
#' intercept term is present in all models. The simplest one M0, the \code{null.model}
#' nested in the rest, contains the fixed variables, if given, and only the intercept by default.
#' In order to implement BIA, \code{\link[MissingBVS]{missingBVS.glm}} can, either perform
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
#' with the \code{BF.approx.method} argument. Approximations can be done through the
#' \pkg{BAS} faster computation if the \code{family} is one of the implemented there:
#' \code{binomial(link = "logit")}, \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and developed for generalized linear models by Li and Clyde (2018). It proceeds
#' using the \pkg{BAS} log-marginal computation (Clyde, 2025) if the \code{family}
#' chosen is ones of the available. Otherwise, method \code{"gprior"} is not provided.
#' The choices currently available are:
#' -"Robust" is the default option and denotes the criteria-based prior of Bayarri,
#' Berger, Forte and Garcia-Donato (2012).
#' -"gZellner" corresponds to the prior in Zellner (1986) with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.WNC" is the intrinsic prior derived by Womack, Novelo and Casella (2014).
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
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link function to be
#' used in the model. If it is one of the implemented families in \pkg{BAS}:
#' \code{binomial(link = "logit")}, \code{poisson(link = "log")} and \code{Gamma(link = "log")};
#' a faster version using \pkg{BAS} log-marginal computation is performed.
#' @param null.model Formula defining which is the simplest (null) model, nested in
#' the full one with possible fixed variables. By default, it is defined to be the one
#' with just the intercept.
#' @param BF.approx.method Method used to compute or approximate data-driven Bayes factors
#' (to be literally specified). Possible choices include "BIC", "TBF" and "gprior"
#' (see details).
#' @param prior.betas Prior distribution for model coefficients if "gprior" method is
#' chosen (to be literally specified). Possible choices are: "Robust", "Liangetal",
#' "gZellner", "FLS" and "intrinsic.WNC" (see details).
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
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#' @param laplace Logical variable to access the Laplace approximation to the
#' marginal likelihood of \pkg{BAS}.
#'
#' @return \code{missingBVS.glm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{Time lasted solving the problem}
#' \item{glmfull}{If missings on the \code{formula} competing variables, combination
#' of the estimates of glm fitted full model over the \code{n.imp} imputed datasets.
#' Otherwise, it is the \code{\link[stats]{glm}} object}
#' \item{glmnull}{The \code{glm} class object that results when \code{null.model}
#' is fitted by \code{\link[stats]{glm}}}
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
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information (when relevant)}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{family}{Family function among \code{\link[stats]{family}} used to specify
#' the error distribution and link function to be used in the model}
#' \item{weights}{Weights vector used in the glm fitting process}
#' \item{offset}{Offset vector used in the glm fitting process}
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
#' @seealso Use \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for a heuristic
#' approximation based on Gibbs sampling (recommended when p>20).
#'
#' Consider \code{\link[MissingBVS]{plot.MissingBvs}} for graphical summaries of the
#' posterior distribution.
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' Garcia-Donato, G. and Paulo, R. (2022)<DOI:10.1080/01621459.2021.1889565>
#' Variable Selection in the Presence of Factors: A Model Selection Perspective.
#' Journal of the American Statistical Association. 117. 1-27.
#'
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6(2): 461–464.
#'
#' Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures
#' of g-Priors in Generalized Linear Models. Journal of the American
#' Statistical Association. 113: 1275–1287.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
#' @examples
#' \dontrun{
#' #Indian Prime Diabetes Data from VIM's package
#'
#' #Here we keep the 32 competing models:
#' f <- Outcome ~ Pregnancies + Glucose + Insulin + BMI + Age
#' diabetes.mBVS <- missingBVS.glm(formula = f, data = VIM::diabetes,
#'   family = binomial(), n.keep = 32, n.imp = 100)
#'
#' #Show the results:
#' diabetes.mBVS
#'
#' #Summ up the results:
#' summary(diabetes.mBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(diabetes.mBVS)
#' }
#'

missingBVS.glm <- function (formula,
                            data,
                            family = binomial(link = "logit"),
                            null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                            BF.approx.method = "BIC",
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
                            imp.seed = runif(1,0,09011975),
                            weights = rep.int(1, dim(data)[1]),
                            offset = rep.int(0, dim(data)[1]),
                            control = glm.control(),
                            laplace = 0L) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #Response in the null model and full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  #select environment to get glm arguments
  environment(formula) <- environment(null.model) <- env <- environment()

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(formula, null.model, data, marginal.factors)
  list2env(buildmatrices.list, envir = env)

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(p, p0, namesnull, namesx, n.keep, q)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  mF <- L > 0 & marginal.factors
  #Check if factors present and if marginalization of their probabilities. Define model prior
  lp.model <- checkmarg.factorsprior(mF, prior.models.dummies, l, positions, positionsfac,
                                     firstd, lprior.models)

  #Evaluate the null model:
  glmnull <- glm(formula = null.model,
                 data,
                 y = TRUE, x = TRUE,
                 family = family,
                 weights = weights,
                 offset = offset,
                 control = control)
  #correct glm arguments
  family <- glmnull$family; weights <- glmnull$prior.weights; offset <- glmnull$offset

  #The response variable
  y <- glmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y) #observations without missings on the response
  devnull <- glmnull$deviance #deviance of the null model
  logmargnull <- as.numeric(-0.5 * devnull) #= log-marginal likelihood of null model - K, K constant

  y <- as.numeric(y); laplace <- as.integer(laplace) #for the C code

  #check whether or not the family chosen is among the options provided by BAS
  inBAS <- checkforfamily(family, BF.approx.method)

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS, n, p, p0, y,
                                              logmargnull, family, devnull,
                                              weights, offset, control, laplace)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define variables with NAs
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
    #remove first dummy on each factor, first p0 vars are the fixed ones
    if (L > 0) imp.array <- imp.array[,-c(indf + p0), , drop = FALSE]
    #Evaluate glm of full model with missings using Rubin's rule
    fit <- list(); mt <- attr(framefull, "terms")
    for (i in 1:n.imp) {
      z <- glm.fit(x = imp.array[,,i], y = y, family = family,
                   weights = weights, offset = offset, control = control)
      z$terms <- mt; class(z) <- "glm"; fit[[i]] <- z
    }
    glmfull <- mice::pool(fit)
    glmfull$call <- NULL #otherwise, Rstudio returns a warning trying to read lmfull$call
  } else glmfull <- glm(formula,
                        data,
                        x = TRUE, y = TRUE,
                        family = family,
                        weights = weights,
                        offset = offset,
                        control = control)

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$glmfull <- glmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, glmfull is the glm object for the full model
  result$glmnull <- glmnull # The glm object for the null model (without NAs)

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
  }

  #glm arguments
  result$family <- family; result$weights <- weights; result$offset <- offset

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
checkforfamily <- function (family, BF.approx.method) {
  #Returns a logical indicating if BAS functions can be used to speed up the process
  #if method is BIC or TBF

  #families implemented in BAS logmarginal computation
  if (family$family %notin% c("binomial", "poisson", "Gamma")) {
    inBAS <- FALSE
  } else {
    if ((family$family == "binomial" & family$link != "logit") |
        (family$family %in% c("poisson", "Gamma") & family$link != "log")) {
      inBAS <- FALSE
    } else inBAS <- TRUE
  }

  if (BF.approx.method == "gprior" & !inBAS) stop("family not implemented in BAS' marginal computation.\n",
                                                  "Try with method 'BIC' or 'TBF'.\n")

  return(inBAS)
}

#' @keywords internal
checkforprior.betas.glm <- function (BF.approx.method, prior.betas, inBAS,
                                     n, p, p0, y, logmargnull,
                                     # null.model,
                                     # data,
                                     family,
                                     devnull, weights, offset, control, laplace) {
  #checks that the Bayes factor computation method given by BF.approx.method and prior.betas
  #is implemented and returns the function to use for Bayes factor computation on glm
  if (BF.approx.method %notin% c("BIC", "TBF", "gprior")) {
    stop("Only BF approximations 'BIC', 'TBF' and 'gprior' supported.")
  }

  #if possible, use faster computation of BAS
  if (inBAS) {
    # c_glm.fit <- function() { #to compute logmarginals from BAS in glm
    #   utils::getFromNamespace("C_glm_deterministic", "BAS")
    # }

    if (BF.approx.method != "BIC") {

      switch (prior.betas, # change the string for the corresponding BAS function
              gZellner = {prior.betas <- BAS::g.prior(g = n)}, #fixed g=n
              Robust = {prior.betas <- BAS::robust(as.numeric(n))}, #random g
              Liangetal = {prior.betas <- BAS::hyper.g.n(alpha = 3, n = n)}, #random g: hyper-g/n with a=3
              # `Zellner-Siow` = {prior.betas <- "ZSBF"}, #random g: cauchy prior, Jeffreys in BAS?
              FLS = {prior.betas <- BAS::g.prior(g = max(n, p^2))}, #fixed Benchmark prior: g=max(n, p*p)
              `intrinsic.WNC` = {prior.betas <- BAS::intrinsic(as.numeric(n))}, #intrinsic prior from Womack, Novelo and Casella (2014)
              # IHG = {prior.betas <- "geointrinsicBF"} #intrinsic hyper-g prior, not available in BAS?
              stop("For now, prior.betas must be one of 'gZellner', 'Robust', 'Liangetal',",
                   "'FLS' or 'intrinsic.WNC' when using TBF or gprior method.\n")
      )
    } else prior.betas <- BAS::bic.prior(n = n)

    switch (BF.approx.method,
            BIC = {BF.approx.method.f <-
              function (k, X) BF.approx.BIC.glm(y = y, X,
                                                family = family,
                                                logmargnull = logmargnull,
                                                k, p0 = p0,
                                                weights = weights,
                                                offset = offset,
                                                control = control,
                                                laplace = laplace)},
            TBF = {BF.approx.method.f <-
              function (k, X) BF.approx.TBF.glm(y = y, X,
                                                family = family,
                                                devnull = devnull,
                                                prior.betas = prior.betas,
                                                k, p0 = p0,
                                                weights = weights,
                                                offset = offset,
                                                control = control,
                                                laplace = laplace)},
            gprior = {BF.approx.method.f <-
              function (k, X) BF.approx.gprior.glm(y = y, X,
                                                   family = family,
                                                   logmargnull = logmargnull,
                                                   prior.betas = prior.betas,
                                                   k, p0 = p0,
                                                   weights = weights,
                                                   offset = offset,
                                                   control = control,
                                                   laplace = laplace)}
    )
  } else { #use slower option that not depends on BAS
    cat("Faster BAS computation cannot be used for chosen arguments.",
        "Be aware that it can take a while.\n")

    switch (BF.approx.method,
            BIC = {BF.approx.method.f <-
              function (k, X) BF.approx.BIC.glm.stats(y = y, X,
                                                      family = family,
                                                      devnull = devnull,
                                                      n, k, p0 = p0,
                                                      weights = weights,
                                                      offset = offset,
                                                      control = control)},
            TBF = {BF.approx.method.f <-
              function (k, X) BF.approx.TBF.glm.stats(y = y, X,
                                                      family = family,
                                                      devnull = devnull,
                                                      n, k, p0 = p0,
                                                      weights = weights,
                                                      offset = offset,
                                                      control = control)}
    )
  }

  return(BF.approx.method.f)
}

