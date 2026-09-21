#' Bayes Factor Averaging for Variable Selection with Missing data in generalized linear models
#'
#' Computation and summaries of posterior distribution over the model space in problems
#' of small to moderate size in the presence of (possible) missing data and/or categorical
#' variables in generalized linear models. Each posterior model probability is computed
#' following the Bayes Factor Averaging (BFA) framework, using standard priors
#' for model coefficients and the hierarchical approach of García-Donato and Paulo (2022)
#' with factors.
#'
#' The set of competing models is made up by all the possible subsets of regressors
#' specified by \code{formula}: Mi for i in 1,...,2^p, being p the number of potential
#' (non-fixed) regressors in the variable selection problem. It is assumed that the
#' intercept term is present in all models. The simplest one M0, the \code{null.model}
#' nested in the rest, contains the fixed variables, if given, and only the intercept by default.
#' In order to implement BFA, \code{\link[MissingBVS]{missingBVS.glm}} can, either perform
#' \code{n.imp} imputations designed by \code{imp.predict.mat} and \code{imp.mice.method}
#' with the \pkg{mice} package, or use user-given imputated datasets by the
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
#' with the BIC (Schwarz, 1978), and the default choice, or the test-based BF
#' (Held, Gravestock and Sabanés, 2015) with the \code{BF.method} argument.
#' Approximations can be done through the \pkg{BAS} faster computation if the
#' \code{family} is one of the implemented there: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and developed for generalized linear models by Li and Clyde (2018). It proceeds
#' using the \pkg{BAS} log-marginal computation (Clyde, 2025) if the \code{family}
#' chosen is ones of the available. Otherwise, method \code{"gprior"} is not provided.
#' The choices currently available are:
#' -"Robust" denotes the criteria-based prior of Bayarri, Berger, Forte and
#' Garcia-Donato (2012).
#' -"gZellner" is the default option and corresponds to the prior in Zellner (1986)
#' with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.WNC" is the intrinsic prior derived by Womack, Novelo and Casella (2014).
#'
#' If the BF computation method chosen is \code{"TBF"}, \code{prior.betas} also determines
#' the formula of the BF as before. In this case, the choices available are: "gZellner",
#' "FLS", "Liangetal" and the adapted version of the ZellnerSiow" prior for TBF.
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
#' @param n.keep It can be either the character "all" to return the whole model
#' space or a numeric for the exact number of the most probable models to keep.
#' By default it is set to 10 and automatically adjusted if 10 is greater than
#' the total number of models.
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
#' \item{BF.method}{Method used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
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
#' @seealso Use \code{\link[MissingBVS]{missingGibbsBVS.glm}} for a heuristic
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
#' Held L, Sabanés Bové D, Gravestock I (2015).<DOI:10.1214/14-STS510> Approximate
#' Bayesian Model Selection with the Deviance Statistic. Statistical Science. 30.
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
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + Insulin + BMI + Age
#' #Keep the 32 competing models from five candidate variabeles.
#' #Few imputations for simplicity, real analyses need more.
#'
#' diabetes.mBVS <- missingBVS.glm(
#'   formula = f, data = diabetes, family = binomial(), n.keep = 32,
#'   n.imp = 10, imp.seed = 1
#' )
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
#' diabetes.mBVS$inclprob
#'
#' #Pool of estimates for model given by formula:
#' diabetes.mBVS$glmfull
#'
#' f <- Outcome ~ Pregnancies + Glucose + Insulin + BMI + Age
#'
#' #User given prior probs for model size and Robust g-prior for BFs:
#' diabetes.mBVS.Robust.userprobs <- missingBVS.glm(
#'   formula = f, data = diabetes, family = binomial(),
#'   prior.models = "User", priorprobs = c(0.3, 0.2, 0.2, 0.1, 0.1, 0.1),
#'   #more mass probability over small models
#'   BF.method = "gprior", prior.betas = "Robust",
#'   n.imp = 2, imp.seed = 1
#' )
#'
#' #Other summaries of the posterior distribution:
#' diabetes.mBVS.Robust.userprobs$HPMbin #Highest posterior Probability model
#' diabetes.mBVS.Robust.userprobs$MPMbin #Median Probability model
#' }
#'

missingBVS.glm <- function (formula,
                            data,
                            family = binomial(link = "logit"),
                            null.model = update(as.formula(formula), . ~ 1),
                            BF.method = "BIC",
                            prior.betas = NULL,
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

  formula <- as.formula(formula); environment(formula) <- environment()
  null.model <- as.formula(null.model); environment(null.model) <- environment()

  #Response in the null model and full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  #Build matrices and objects needed later on
  matrices <- buildmatrices(formula, null.model, data, marginal.factors)

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(
    matrices$p, matrices$p0, matrices$namesnull, matrices$namesx, n.keep, matrices$q
  )

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, matrices$q)

  mF <- matrices$L > 0 & marginal.factors
  positionscov <- if (mF) {
    matrices$positions[matrices$positionsx, , drop = FALSE]
  } else NULL
  positionsfac <- if (mF) matrices$positionsfac else NULL

  satmodels.repr <- if (mF) matrices$satmodels.repr else NULL
  l <- if (mF) matrices$l else NULL
  #Check if factors present and if marginalization of their probabilities. Define model prior
  lp.model <- checkmarg.factorsprior(
    mF, prior.models.dummies, matrices$l,
    positionscov, positionsfac, satmodels.repr, lprior.models
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
  y <- as.numeric(y); laplace <- as.integer(laplace) #for the C code

  #check whether or not the family chosen is available for BF.method
  useBAS <- checkforfamily(family, BF.method)

  #Check approx method and priors chosen and define the function to be used
  lBF <- checkforprior.betas.glm(
    BF.method, prior.betas, n, matrices$p, matrices$p0, y, glmnull, useBAS, laplace
  )

  matrices$X.full <- matrices$X.full[obsnotNA,]

  #check for missings and define variables with NAs
  NAvars <- checkformissings(
    y = matrices$framenull[, 1], matrices$framenull[, -1], matrices$X.full
  )

  #Define function to get binary expression for each model
  num2bin.model.fun <- function(x) num2bin.model(x, matrices$p, NAvars)

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

    switch (as.character(useBAS),
            `TRUE` = {lBF.method <- function(model) lBF.av(
                model, imputation.array = imputation$imputation.array,
                lBF = lBF, p0 = matrices$p0, n.imp = n.imp
              )
            },
            `FALSE` = {lBF.method <- function(model) lBF.av.glm.fit(
                model, imputation.array = imputation$imputation.array,
                lBF = lBF, p0 = matrices$p0, n.imp = n.imp, y = y, glmnull = glmnull
              )
            }
    )

  } else {
    #When there are no missings, just compute the BF
    lBF.method <- function(model) lBF(
      k = length(model),
      X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ],
      fitstart = NULL
    )
  }
  #for posterior computation, if no NAvars active we do not need fitstart
  lBFfitnull <- function(k, X) lBF(k, X, fitstart = NULL)

  #Info:
  cat("\nInfo. . .\n")
  if (mF) {
    cat("Most complex model has a total of", matrices$q + matrices$q0,
        "covariates and/or factors.\n")
  } else cat("Most complex model has a total of", matrices$q + matrices$q0,
             "competing variables.\n")
  if (matrices$q0 == 1) {
    cat("From those 1 is fixed (the intercept) and we should select from the remaining",
        matrices$q, "\n")
  } else cat("From those", matrices$q0, "are fixed and we should select from the remaining",
             matrices$q, "\n")
  if (mF) {
    cat("  Numerical covariates:", matrices$depvars[matrices$positionsx], "\n")
    cat("  Factors:", matrices$depvars[!matrices$positionsx], "\n")
  } else cat("  Competing variables:", matrices$depvars, "\n")

  cat("The problem has a total of", 2^matrices$q, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #Compute exact posterior distribution and normalizing constant
  posterior <- exact.posterior.comput(
    matrices, num2bin.model.fun, lBF.method, lp.model, lBFfitnull
  )

  #Summ up the posterior distribution
  posterior.summary <- summ.posterior(posterior, matrices, mF)

  if (anyNAvar) {#Pool results for imputed datasets
    imp.array <- imputation$imputation.array
    #if marginal factor probs, remove first dummy on each factor, first p0 are the fixed ones
    if (mF) {
      imp.array <- imp.array[, -c(matrices$indf + matrices$p0), , drop = FALSE]
    }
    #Evaluate glm of full model with missings using Rubin's rule
    fit <- list(); mt <- attr(matrices$framefull, "terms")
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
  result$time <- Sys.time() - time #The time it took the program to finish
  result$glmfull <- glmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, glmfull is the glm object for the full model
  result$glmnull <- glmnull # The glm object for the null model (without NAs)

  result$variables <- matrices$depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- matrices$q #number of competing vars
  result$k <- matrices$q0 #number of fixed vars
  result$HPMbin <- posterior.summary$hpm
  result$MPMbin <- posterior.summary$mpm
  names(result$MPMbin) <- matrices$depvars

  if (mF) {
    #matrix for the factors index
    result$positions <- matrices$positionsfac
    result$positionsx <- matrices$positionsx
    result$modelsrankdefprob <- posterior$all.models.PM
  }

  #The binary code for the n.keep best models and the correspondent post
  result$modelsprob <- posterior.summary$modelsprob[
    order(posterior.summary$modelsprob[, matrices$q + 1], decreasing = TRUE)[seq_len(n.keep)],
  ]
  dimnames(result$modelsprob) <- list(seq_len(n.keep), c(matrices$depvars, "Post"))

  result$inclprob <- posterior.summary$inclprob
  names(result$inclprob) <- matrices$depvars

  result$postprobdim <- posterior.summary$probdim
  names(result$postprobdim) <- 0:matrices$q + matrices$q0
  result$C <- posterior$C

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

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}

#' Logarithm of the g-prior (g fixed) Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' with g fixed, for complete data in generalized linear models.
#'
#' Computes the approximated expression for Bayes factors under g-priors
#' with g fixed, in logarithmic scale:
#'  lBF_10 = z1/2 + 1/2 log(Jaa0/Jaa1) - k/2 log(1+g) - Q1/(2 * (1+g)).
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param glmnull glm object derived from fitting the null model.
#' @param g Fixed value g for g-prior.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param fitstart Optional starting values for the parameters in the linear
#' predictor. By default, it is \code{NULL}.
#'
#' @return \code{BF.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{missingBVS.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL), n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), y = TRUE)
#' lBF <- MissingBVS:::BF.gprior.glm.fit(y = glmnull$y, X = imp1$imputation.array[,,1],
#'   glmnull, g = length(glmnull$y))
#'
#' @references
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#'
#' @keywords internal
BF.gprior.glm.fit <- function (y, X, #family = binomial(link = "logit"),
                               glmnull, g = n,
                               #default corresponds to gZellner
                               n = length(y), k = ncol(X)-1,
                               weights = rep(1, length(y)),
                               offset = rep(0, length(y)),
                               fitstart = NULL) {

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = glmnull$family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               method = 2) #LLT Cholesky decomposition, faster

  #Get information matrix for model given by X
  W <- diag(fit1$weights)
  J <- crossprod(X, W %*% X) #t(X) %*% W %*% X

  #J = [Jaa Jab; Jba Jbb]
  Jaa <- J[1, 1] #information associated to intercept
  Jab <- J[1, -1, drop = FALSE]
  Jba <- J[-1, 1, drop = FALSE]
  Jbb <- J[-1, -1, drop = FALSE] #information associated with non-fixed terms

  Jbeta <- Jbb - Jba %*% solve(Jaa, Jab) #marginal information matrix for beta

  beta <- fit1$coefficients[-1]
  # Wald statistic
  Q <- as.numeric(crossprod(beta, Jbeta %*% beta))

  #Get information matrix for null model
  devnull <- glmnull$deviance
  X0 <- glmnull$x
  W0 <- diag(glmnull$weights)

  Jaa0 <- crossprod(X0, W0 %*% X0)[1,1] #t(X0) %*% W0 %*% X0

  #change in deviance
  z <- devnull - fit1$deviance

  #Bayes factor derived from approximate marginal likelihood (Li and Clyde, 2018)
  lBFi0 <- z/2 - log(Jaa / Jaa0)/2 - k/2 * log1p(g) - Q/(2 * (1 + g))
  return(lBFi0)
}

#currently only for intercept fixed

#' Logarithm of the hyper g-prior Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given hyper g-prior,
#'  for complete data in generalized linear models.
#'
#' Computes the approximated expression for Bayes factors under hyper g-priors,
#' in logarithmic scale:
#'  lBF_10 = z1/2 - Q1/(v) + 1/2  log(Jaa0/Jaa1) - k/2 * v +
#     log B((a + k)/2, b/2) + log phi1(b/2, r, (a + b + k)/2, (s + Q1)/(2*v), 1 - ka) -
#     log B(a/2, b/2) - log phi1(b/2, r, (a + b)/2, s/(2*v), 1 - ka)
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param glmnull glm object derived from fitting the null model.
#' @param prior.betas.args List of arguments corresponding to a specific hyper
#' g-prior. It must contain a, b, r, s, v and ka.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param fitstart Optional starting values for the parameters in the linear
#' predictor. By default, it is \code{NULL}.
#'
#' @return \code{BF.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{missingBVS.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL), n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), y = TRUE)
#' lBF <- MissingBVS:::BF.hypergprior.glm.fit(y = glmnull$y, X = imp1$imputation.array[,,1],
#'   glmnull)
#'
#' @references
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#'
#' @keywords internal
BF.hypergprior.glm.fit <- function (y, X, # family = binomial(link = "logit"),
                                    glmnull,
                                    prior.betas.args = list(a = 1, b = 2, r = 0,
                                                            s = 0, v = 1, ka = 1),
                                    #default coresponds to liangetal prior
                                    n = length(y), k = ncol(X)-1, #p0 = 1,??
                                    weights = rep(1, length(y)),
                                    offset = rep(0, length(y)),
                                    fitstart = NULL) {

  unlist(prior.betas.args) #contains a, b, r, s, v, ka
  #update arguments that depend on k, if needed
  if (is.function(v)) v <- v(k)
  if (is.function(ka)) ka <- ka(k)

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = glmnull$family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               method = 2) #LLT Cholesky decomposition, faster

  #Get information matrix for model given by X
  # W <- diag(fit1$weights)
  J <- crossprod(X, fit1$weights * X) #t(X) %*% W %*% X

  #J = [Jaa Jab; Jba Jbb]
  Jaa <- J[1, 1] #information associated to intercept
  Jab <- J[1, -1, drop = FALSE]
  # Jba <- J[-1, 1, drop = FALSE]
  Jbb <- J[-1, -1, drop = FALSE] #information associated with non-fixed terms

  # Jbeta <- Jbb - Jba %*% solve(Jaa, Jab) #marginal information matrix for beta
  Jbeta <- Jbb - tcrossprod(Jab) / Jaa

  beta <- fit1$coefficients[-1]
  # Wald statistic
  Q <- as.numeric(crossprod(beta, Jbeta %*% beta))

  #Get information matrix for null model
  devnull <- glmnull$deviance
  # X0 <- glmnull$x
  # W0 <- diag(glmnull$weights)
  #
  # Jaa0 <- crossprod(X0, W0 %*% X0)[1,1] #t(X0) %*% W0 %*% X0
  Jaa0 <- sum(glmnull$weights)

  #change in deviance
  z <- devnull - fit1$deviance

  #Bayes factor derived from approximate marginal likelihood (Li and Clyde, 2018)
  lBFi0 <- z/2 - Q/(2*v) - log(Jaa / Jaa0)/2 - k/2 * log(v) +
    lbeta((a + k)/2, b/2) + log(BAS::phi1(b/2, r, (a + b + k)/2, (s + Q)/(2*v), 1 - ka)) -
    lbeta(a/2, b/2) - log(BAS::phi1(b/2, r, (a + b)/2, s/(2*v), 1 - ka))

  return(lBFi0)
}

#' @keywords internal
checkforfamily <- function (family, BF.method) {
  #Checks if family is among the available ones

  if (BF.method == "gprior") {
    #families implemented in BAS logmarginal computation
    if (family$family %notin% c("binomial", "poisson", "Gamma")) {
      inBAS <- FALSE
    } else {
      if ((family$family == "binomial" & family$link != "logit") |
          (family$family %in% c("poisson", "Gamma") & family$link != "log")) {
        inBAS <- FALSE
      } else inBAS <- TRUE
    }

    # if (!inBAS) stop("family not implemented in BAS' marginal computation.\n",
    #                  "Try with method 'BIC' or 'TBF' instead.\n")
  }

  return(inBAS)
}
