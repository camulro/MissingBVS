#' Bayesian Variable Selection with Missing data for generalized linear
#' regression models
#'
#' Computation and summaries of posterior distribution over the model space
#' in problems of small to moderate size when missingness occurs in generalized
#' linear models.
#'
#' The set of competing models is made up by all the possible subsets of
#' regressors specified by \code{formula}: Mi for i in 1,...,2^p, being p the number
#' of potential (non-fixed) regressors in the variable selection problem. The simplest,
#' nested in all of them, contains only the intercept. \code{MissingBvs} performs
#' \code{n.imp} imputations given by \code{imp.mice.method} with the mice package and
#' computes the posterior distribution over this model space through Bayes' theorem:
#'
#' Pr(Mi | \code{data})=Pr(Mi)*Bi/C,
#'
#' where Bi is an approximation of the Bayes factor of Mi to M0 under missing data,
#' Pr(Mi) is the prior probability of Mi and C is the normalizing constant.
#'
#' Bi is computed as the average of the \code{n.imp} standard Bayes factors, B_i(j)
#' for j in 1,...,\code{n.imp}, for Mi to M0 calculated for the jth imputed data set.
#' If the BF computation method chosen is \code{"gprior"}, BF B_i(j) depends on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and \code{MissingBvs.glm} uses the same choices as
#' \code{[MissingBVS]{missingBVS.lm}} for generalized linear models developed by
#' Li and Clyde (2018). It does so using the \pkg{BAS} logmarginal computation
#' (Clyde, 2025) if the \code{family} chosen is one of the implemented in the package.
#'
#' The BF can also be approximated with the BIC (Schwarz, 1978) or the test-based
#' BF (Held, Gravestock and Sabanés, 2015). It can be done through the \pkg{BAS}
#' faster computation if the \code{family} is one of the implemented in the package.
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
#' model with fixed variables; next p components correspond to the p prior probabilities
#' of the possible model dimensions.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. If it is one of the implemented families in
#' \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}; a faster version
#' using BAS logmarginal computation is performed.
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
#' @param n.keep It can be either the character "all" to return the whole model
#' space or a numeric for the exact number of the most probable models to keep.
#' By default it is set to 10 and automatically adjusted if 10 is greater than
#' the total number of models.
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
#' @param imp.seed Seed for imputation.
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#' @param laplace Logical variable to access the Laplace approximation to the
#' marginal likelihood of \pkg{BAS}. See \code{\link[BAS]{bas.glm}}
#' for more details.
#'
#' @return \code{missingBVS.glm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{glmfull}{If missings on the competing variables, object of class
#' \code{\link[mice]{mipo}} that combines the estimates for the model defined by
#' \code{formula} fitted by \code{\link[stats]{glm}} over
#' the \code{n.imp} (if \code{> 1}) imputed datasets; see \code{\link[mice]{pool}}
#' for details. Otherwise, it is the \code{glm} object for the \code{formula} model}
#' \item{glmnull}{The \code{glm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{glm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model}
#' \item{positions}{\code{matrix} with L rows and p plus the number of dummies
#' resulting from factors columns - L, where L is the number of factors, with 1
#' if the column dummy corresponds to the row factor and 0 otherwise}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate}
#' \item{modelsprob}{A \code{floor(n.iter/n.thin)}x(p+1) \code{matrix} which
#' summaries the keeped models and their associated Bayes factor in logaritmic scale}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables}
#' \item{postprobdim}{Posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the
#' model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{logprior.models.dumm}{Function used to compute the log-prior over the
#' model space of the factor levels}
#' \item{method}{\code{Full}}
#'
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for a heuristic
#' approximation based on Gibbs sampling (recommended when p>20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
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
#' @keywords package
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
                            priorprobs = NULL,
                            n.keep = 10,
                            parallelmice = NULL,
                            n.core = NULL,
                            imp.time.test = TRUE,
                            imp.mice.method = "pmm",
                            imp.predict.mat = NULL,
                            n.imp = 039E1,
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

  #check whether or not the family chosen is among the options provided by BAS
  inBAS <- checkforfamily(family, BF.approx.method)

  #select environment to get glm arguments
  environment(formula) <- environment()
  environment(null.model) <- environment()

  #for the C code
  weights <- as.numeric(weights)
  offset <- as.numeric(offset)
  laplace <- as.integer(laplace)

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(formula, null.model, data)
  list2env(buildmatrices.list, envir = environment())

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(p, p0, namesnull, namesx, n.keep, q)

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  if (L > 0) {
    lprior.models.dummies <- checkforprior.models.dummies(prior.models.dummies, l)
  } else lprior.models.dummies <- function(delta, tau) 0

  #Evaluate the null model:
  glmnull <- glm(formula = null.model,
                 data,
                 y = TRUE, x = TRUE,
                 family = family,
                 weights = weights,
                 offset = offset,
                 control = control)

  #The response variable
  y <- glmnull$y; obsnotNA <- names(y) #without missings
  y <- as.numeric(y)
  n <- length(y) #observations without missings on the response
  devnull <- glmnull$deviance #deviance of the null model

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS,
                                              n, p, p0, y, null.model,
                                              data, family, devnull,
                                              weights, offset, control, laplace)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define variables with NAs
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
  } else cat(paste0("From those ", q0, " are fixed and we should select from the remaining ",
                    q, ".\n"))

  cat("  Numerical covariates:", depvars[positionsx], "\n")
  if (L > 0) cat(" Factors:", depvars[!positionsx], "\n")

  cat("The problem has a total of", 2^q, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #Compute exact posterior distribution and normalizing constant
  posterior.list <- exact.posterior.comput(p, positions, positionsfac, lastd, l,
                                           namesxnotnull, NAvars, lBF.method,
                                           lprior.models, lprior.models.dummies,
                                           X0, X.full, BF.approx.method)
  list2env(posterior.list, envir = environment())

  #Summ up the posterior distribution
  summ.posterior.list <- summ.posterior(L, all.models.PM, p, q, positions)
  list2env(summ.posterior.list, envir = environment())

  if (!is.null(NAvars)) {
    #Pool results for imputed datasets
    if (n.imp > 1) {
      #remove last dummy for each factor, first p0 vars are the fixed ones
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0),]
      #Evaluate glm of full model with missings using Rubin's rule
      fit <- list()
      mt <- attr(framefull, "terms")
      for (i in 1:n.imp) {
        z <- glm.fit(x = imputation.array[,,i], y = y, family = family,
                     weights = weights, offset = offset, control = control)
        z$terms <- mt
        class(z) <- "glm"

        fit[[i]] <- z
      }
      glmfull <- mice::pool(fit)
    } else {
      #compute glm.fit for the unique imputation
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0)]
        glmfull <- glm.fit(x = imputation.array, y = y, family = family,
                           weights = weights, offset = offset, control = control)
    }
  } else glmfull <- glm(formula,
                        data,
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

  if (L > 0) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  result$modelsprob <- modelsprob[order(modelsprob[,q+1],
                                        decreasing = TRUE)[seq_len(n.keep)],]
  #The binary code for the n.keep best models and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- depvars

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, q + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, q))) #prior inclusion prob for dimension 0
    for (i in seq_len(q)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, q - i))) + lchoose(q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:q + q0 #prior dimension probability

  result$C <- C #normalizing constant

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

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}

# "%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0 #auxiliar function

#' @keywords internal
checkforfamily <- function (family, BF.approx.method) {
  #if family is a character redefines it with the corresponding function,
  #if family is a function, calls it and checks the provided argument is among the possible options.
  #Returns a logical indicating if BAS functions can be used to speed up the process
  #if method is BIC or TBF
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()

  #families implemented in BAS logmarginal computation
  if (family$family %notin% c("binomial", "poisson", "Gamma")) {
    inBAS <- FALSE
  } else {
    if ((family$family == "binomial" & family$link != "logit") |
        (family$family %in% c("poisson", "Gamma") & family$link != "log")) {
      inBAS <- FALSE
    } else inBAS <- TRUE
  }

  if (BF.approx.method == "gprior" & !inBAS) {
    stop(paste("family ", family, "not implemented in BAS' marginal computation.\n",
               "Try with method 'BIC' or 'TBF'.\n"))
  }
  return(inBAS)
}

#' @keywords internal
checkforprior.betas.glm <- function (BF.approx.method, prior.betas, inBAS,
                                     n, p, p0, y, null.model,
                                     data, family, devnull, weights, offset, control, laplace) {
  #checks that the Bayes factor computation method given by BF.approx.method and prior.betas
  #is implemented and returns the function to use for Bayes factor computation on glm
  if (BF.approx.method %in% c("BIC", "TBF", "gprior")) {
    if (inBAS) { #use faster computation of BAS
      if (BF.approx.method != "BIC") {
        if (prior.betas %notin% c("gZellner", "Robust", "Liangetal", # "ZellnerSiow",
                                  "FLS", "intrinsic.WNC" # "IHG"
        )) {
          stop(paste0("For now, prior.betas must be one of 'gZellner', 'Robust', 'Liangetal',\n",
                      "'FLS' or 'intrinsic.WNC' when using TBF or gprior method.\n"))
        }

        switch (prior.betas, # change the string for the corresponding BAS function
                gZellner = {prior.betas <- BAS::g.prior(g = n)}, #fixed g=n
                Robust = {prior.betas <- BAS::robust(as.numeric(n))}, #random g
                Liangetal = {prior.betas <- BAS::hyper.g.n(alpha = 3, n = n)}, #random g: hyper-g/n with a=3
                # `Zellner-Siow` = {prior.betas <- "ZSBF"}, #random g: cauchy prior, Jeffreys in BAS?
                FLS = {prior.betas <- BAS::g.prior(g = max(n, p^2))}, #fixed Benchmark prior: g=max(n, p*p)
                `intrinsic.WNC` = {prior.betas <- BAS::intrinsic(as.numeric(n))} #intrinsic prior from Womack, Novelo and Casella (2014)
                # IHG = {prior.betas <- "geointrinsicBF"} #intrinsic hyper-g prior, not available in BAS?
        )

      } else prior.betas <- BAS::bic.prior(n = n)

      logmargnull <- BAS::bas.glm(formula = null.model, # nullmodel
                                  data = data,
                                  n.models = 1,
                                  method = "MCMC",
                                  betaprior = prior.betas,
                                  family = family,
                                  weights = weights,
                                  offset = offset,
                                  control = control,
                                  laplace = laplace)$logmarg

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
      cat("Faster functions from BAS package to compute marginal likelihoods",
          "cannot be used for specified arguments and it can take a while.\n",
          "Do you want to continue? (y/n)\n")
      if (tolower(readline()) != "y") {
        stop("Consider another arguments for your problem.\n")
      }

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

  } else {
    stop("Only approximations 'BIC', 'TBF' and 'gprior' supported.")
  }
}

