#' Bayesian Imputation Averaging for Variable Selection with Missing data in
#' generalized linear models using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior model distribution using a
#' Gibbs sampling algorithm to explore the model space. Each posterior model probability is computed following the
#' Bayesian Imputation Averaging (BIA) framework, using standard priors for model coefficients
#' and the hierarchical approach of García-Donato and Paulo (2022) with factors.
#'
#' Gibbs sampling search algorithm to avoid exhaustive enumeration of model space
#' when it is unfeasible. It draws from the model posterior distribution and uses
#' frequency of "visits" to construct the estimates. The algorithm was originally
#' proposed by  George and McCulloch (1997). Later, Garcia-Donato and Martinez-Beneito (2013)
#' shown that the sampling strategy in combination with estimates based on frequency of
#' visits provides very reliable results.
#'
#' \code{\link[MissingBVS]{missingGibbsBVS.glm}} is a heuristic approximation of
#' \code{\link[MissingBVS]{missingBVS.glm}}. See the latter for common details.
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
#' @param BF.method Method used to compute or approximate data-driven Bayes factors
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
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#' @param laplace Logical variable to access the Laplace approximation to the
#' marginal likelihood of \pkg{BAS}.
#'
#' @return \code{missingGibbsBVS.glm} returns an object of class \code{missingBVS}
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
#' \item{family}{Family function among \code{\link[stats]{family}} used to specify
#' the error distribution and link function to be used in the model}
#' \item{weights}{Weights vector used in the glm fitting process}
#' \item{offset}{Offset vector used in the glm fitting process}
#' \item{BF.method}{Method used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
#' \item{prior.models}{Two-dimensional vector with \code{prior.models} and
#' \code{prior.models.dummies} chosen. If there are no factors or \code{marginal.factors}
#' is set to FALSE, it saves the only argument used, \code{prior.models}}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' model space probabilities such as García-Donato and Paulo (2022)}
#' \item{method}{String "Gibbs" denoting Gibbs sampling model search}
#'
#' @author  Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{missingBVS.glm}} for an exact computation
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
#' \donttest{
#' # Build a small reproducible binary-response example from airquality.
#' data("airquality")
#' glm_data <- airquality[complete.cases(airquality[, c("Ozone", "Wind",
#'   "Temp", "Solar.R")]), c("Ozone", "Wind", "Temp", "Solar.R")]
#' glm_data$Outcome <- as.integer(glm_data$Ozone > median(glm_data$Ozone))
#' glm_data <- glm_data[, c("Outcome", "Wind", "Temp", "Solar.R")]
#' glm_data$Wind[c(1, 10)] <- NA_real_
#'
#' # Use a short chain for the example; real analyses need more iterations.
#' glm.mGBVS <- missingGibbsBVS.glm(
#'   formula = Outcome ~ Wind + Temp + Solar.R, data = glm_data,
#'   family = binomial(), n.iter = 200, n.burnin = 50, n.imp = 2,
#'   Gibbs.seed = 1, imp.seed = 1
#' )
#'
#' #Show the results:
#' glm.mGBVS
#'
#' #Summ up the results:
#' summary(glm.mGBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(glm.mGBVS)
#' glm.mGBVS$inclprobRB
#' }
#'
missingGibbsBVS.glm <- function (formula,
                                 data,
                                 family = binomial(link = "logit"),
                                 null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
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
                                 imp.seed = runif(1,0,09011975),
                                 weights = rep.int(1, dim(data)[1]),
                                 offset = rep.int(0, dim(data)[1]),
                                 control = glm.control(),
                                 laplace = 0L) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #The response in the null model and in the full model must coincide
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
    matrices$p, matrices$p0, matrices$namesnull, matrices$namesx, init.model,
    mF, positionscov, positionsfac, l
  )

  #Check if factors present and if marginalization of their probabilities. Define model prior
  lp.model <- checkmarg.factorsprior(
    mF, prior.models.dummies, matrices$l, positionscov,
    positionsfac, satmodels.repr, lprior.models
  )

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

  #check whether or not the family chosen is available for BF.method
  checkforfamily(family, BF.method)

  #Check approx method and priors chosen and define the function to be used
  lBF <- checkforprior.betas.glm(
    BF.method, prior.betas, n, matrices$p, matrices$p0, y, glmnull, laplace
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

    switch (as.character(BF.method == "gprior"),
            `TRUE` = {lBF.method <- function(model) lBF.av(
              model, imputation.array = imputation$imputation.array,
              lBF = lBF, p0 = matrices$p0, n.imp = n.imp
            )
            lBFfitnull <- lBF
            },
            `FALSE` = {lBF.method <- function(model) lBF.av.glm.fit(
              model, imputation.array = imputation$imputation.array,
              lBF = lBF, p0 = matrices$p0, n.imp = n.imp, y = y, glmnull = glmnull
            )
            #for posterior computation, if no NAvars active we do not need fitstart
            lBFfitnull <- function(k, X) lBF(k, X, fitstart = NULL)
            }
    )
  } else {
    # If n.imp == 1, we do not need fitstart argument
    if (BF.method != "gprior") lBF <- function(k, X) lBF(k, X, fitstart = NULL)

    lBF.method <- function(model) lBF(
      k = length(model),
      X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ]
    )
    lBFfitnull <- lBF
  }

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
     lBFfitnull, mF, matrices$positions, init.model, n.iter, n.burnin, n.thin
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
    #Evaluate glm of full model with missings using Rubin's rule
    fit <- list()
    mt <- attr(matrices$framefull, "terms")
    for (i in 1:n.imp) {
      z <- glm.fit(x = imp.array[,,i], y = y, family = family,
                   weights = weights, offset = offset, control = control)
      z$terms <- mt; class(z) <- "glm"; fit[[i]] <- z
    }

    glmfull <- mice::pool(fit)
    glmfull$call <- NULL #otherwise, Rstudio returns a warning trying to read glmfull$call
  } else glmfull <- glm(formula,
                        data,
                        x = TRUE, y = TRUE,
                        family = family,
                        weights = weights,
                        offset = offset,
                        control = control)

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the programm to finish
  result$glmfull <- glmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, glmfull is the glm object for the full model
  result$glmnull <- glmnull # The glm object for the null model (without NAs)

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
      priorprobs[i+1] <-
        exp(lprior.models(c(rep.int(1, i), rep.int(0, matrices$q - i))) +
          lchoose(matrices$q, i))
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

  #glm arguments
  result$family <- family; result$weights <- weights; result$offset <- offset

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
