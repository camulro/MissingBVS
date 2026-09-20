#' Bayes Factor Averaging for Variable Selection with Missing data in linear regression models
#'
#' Computation and summaries of posterior distribution over the model space in problems
#' of small to moderate size in the presence of (possible) missing data and/or categorical
#' variables in linear models. Each posterior model probability is computed following the
#' Bayes Factor Averaging (BFA) framework, using standard priors for model coefficients
#' and the hierarchical approach of García-Donato and Paulo (2022) with factors.
#'
#' The set of competing models is made up by all the possible subsets of regressors
#' specified by \code{formula}: Mi for i in 1,...,2^p, being p the number of potential
#' (non-fixed) regressors in the variable selection problem. It is assumed that the
#' intercept term is present in all models. The simplest one M0, the \code{null.model}
#' nested in the rest, contains the fixed variables, if given, and only the intercept by default.
#' In order to implement BFA, \code{\link[MissingBVS]{missingBVS.lm}} can, either perform
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
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and are computed using \pkg{BayesVarSel}. The choices currently available are:
#' -"Robust" denotes the criteria-based prior of Bayarri, Berger, Forte and
#' Garcia-Donato (2012).
#' -"gZellner" is the default option and corresponds to the prior in Zellner (1986)
#' with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"ZellnerSiow" is the multivariate Cauchy prior by Zellner and Siow (1980, 1984).
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015).
#' -"IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato,
#' Moreno and Pericchi (2022).
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
#' @seealso Use \code{\link[MissingBVS]{missingGibbsBVS.lm}} for a heuristic
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
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. A. and Berger, J. O. (2008).
#' Mixtures of g Priors for Bayesian Variable Selection. Journal of the
#' American Statistical Association, 103(481), 410–423.
#'
#' Moreno, E., Giron, J. and Casella, G. (2015) Posterior model consistency
#' in variable selection as the model dimension grows. Statistical Science. 30:
#' 228-241.
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
#' Zellner A, Siow A (1980). Posterior Odds for Selected Regression Hypotheses.
#' In JM Bernardo, MH DeGroot, DV Lindley, AFM Smith (eds.), Bayesian Statistics,
#' pp. 585–603. Valencia University Press.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6: 461–464.
#'
#' Held L, Sabanés Bové D, Gravestock I (2015).<DOI:10.1214/14-STS510> Approximate
#' Bayesian Model Selection with the Deviance Statistic. Statistical Science. 30.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @examples
#'
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' # Keep the eight models generated by three candidate covariates,
#' # BIC approximation for BF by default, few imputations for simplicity, real analyses need more.
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(
#'   formula = f, data = dataS97, n.keep = 8, n.imp = 2, imp.seed = 1
#' )
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
#' dataS97.mBVS$inclprob
#' dataS97.mBVS$postprobdim
#'
#' #Pool of estimates for model given by formula:
#' dataS97.mBVS$lmfull
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#'
#' #User given prior probs for model size and Robust g-prior for BFs:
#' dataS97.mBVS.Robust.userprobs <- missingBVS.lm(
#'   formula = f, data = dataS97, prior.models = "User",
#'   priorprobs = c(0.3, 0.3, 0.25, 0.15), #more mass probability over small models
#'   BF.method = "gprior", prior.betas = "Robust",
#'   n.imp = 2, imp.seed = 1
#' )
#'
#' #Other summaries of the posterior distribution:
#' dataS97.mBVS.Robust.userprobs$HPMbin #Highest posterior Probability model
#' dataS97.mBVS.Robust.userprobs$MPMbin #Median Probability model
#' }
#'
missingBVS.lm <- function (formula,
                           data,
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
                           imp.seed = runif(1,0,09011975)) {
  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #Response in the null model and full model must coincide
  if (formula[[2]] != null.model[[2]]) {
    stop("The response in the full and null model does not coincide.\n")
  }

  #Build matrices and objects needed later on
  matrices <- buildmatrices(formula, null.model, data, marginal.factors)

  #Check arguments and compute n.keep if needed
  n.keep <- checkBvsarguments(
    matrices$p, matrices$p0, matrices$namesnull, matrices$namesx, n.keep, matrices$q
  )

  #Check model priors chosen and define the functions to be used
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
    #by lBF over the imputed datasets
    lBF.method <- function(model) lBF.av(
        model, imputation.array = imputation$imputation.array,
        lBF = lBF, p0 = matrices$p0, n.imp = n.imp
      )
  } else lBF.method <- function(model) lBF(
    k = length(model),
    X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ]
  )

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
    matrices, num2bin.model.fun, lBF.method, lp.model, lBF
  )

  #Summ up the posterior distribution
  posterior.summary <- summ.posterior(posterior, matrices, mF)

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

  ##result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$lmfull <- lmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, lmfull is the lm object for the full model
  result$lmnull <- lmnull # The lm object for the null model (omits NAs)

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

    ME(n.imp, imp.seed)
  }

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

#' @keywords internal
exact.posterior.comput <- function (matrices, num2bin.model.fun, lBF.method,
                                    lp.model, lBF) {
  #Compute exact posterior distribution and normalizing constant
  p <- matrices$p
  cat("\n")

  #progress bar for loop
  pb <- list(i = 0, total = 2^p, start.time = Sys.time())
  pb$tick <- function(i) update.progress(pb$i <- i, pb$total, pb$start.time)

  #Posterior computation
  all.models.lPM <- matrix(0, nr = 2^p, nc = p+1) #last column contains log(BF_a0*Pr(M))
  for (i in seq_len(2^p-1)){ # null out of the loop

    #transform the number of the model into a binary number
    current.model <- num2bin.model.fun(i)
    lpm <- lp.model(current.model["bin",]) #log-model prior
    if (is.na(lpm)) {
      all.models.lPM[i,] <- NA
      next #do not visit saturated or oversaturated models
    }

    all.models.lPM[i, seq_len(p)] <- current.model["bin",]

    #check if there are NAs in the model considered to save computation time
    if (sum(current.model["NA",]) > 0) {
      lBF.PM <- lBF.method(model = which(current.model["bin",] == 1)) + lpm #log(BF_a0*Pr(M))

    } else { #if there are no missings, compute the BF by the method selected
      X.i <- cbind(matrices$X0, matrices$X.full[, which(current.model["bin",] == 1)])
      lBF.PM <- lBF(k = sum(current.model["bin",] == 1), X = X.i) + lpm #log(BF_a0*Pr(M))
    }

    all.models.lPM[i, p+1] <- lBF.PM

    #update bar
    pb$tick(i)
  }

  #null model
  all.models.lPM[2^p, seq_len(p)] <- numeric(p)
  all.models.lPM[2^p, p+1] <- lp.model(numeric(p)) #BF = 1 for null model
  all.models.lPM <- na.omit(all.models.lPM) #remove repeated models if dummies

  #update bar
  pb$tick(i+1)
  cat("\n")

  #renormalize
  logC <- logsumexp.stable(all.models.lPM[, p+1])
  all.models.PM <- all.models.lPM
  all.models.PM[, p+1] <- exp(all.models.lPM[, p+1] - logC)
  colnames(all.models.PM) <- c(colnames(current.model), "Post")

  return(list(all.models.PM = all.models.PM, C = exp(logC)))
}

#' @keywords internal
summ.posterior <- function (posterior, matrices, mF) {
  #Summ up the posterior distribution

  all.models.PM <- posterior$all.models.PM
  p <- matrices$p
  q <- matrices$q
  positions <- if (mF) matrices$positions else NULL

  if (mF) {
    #compute models matrix at the covariate-factor level
    cf.models.PM <- all.models.PM[,seq_len(p)] %*% t(positions)
    cf.models.PM <- cbind(cf.models.PM, all.models.PM[,p+1])
    colnames(cf.models.PM)[q+1] <- "Post"

    # #with hashes:
    # modelhashes <- apply(cf.models.PM[,seq_len(q)] > 0, 1, digest::digest)
    # modelsprob <- buildmodelsmatrix(q, colnames(cf.models.PM)[-q-1])
    # for (i in seq_len(2^q)) {
    #   modeli <- digest::digest(modelsprob[i, seq_len(q)] > 0)
    #
    #   modelsprob[i, q+1] <- sum(cf.models.PM[which(modeli == modelhashes), q + 1])
    # }

    #bin to int codification of models: (not too many in Bvs)
    modelsprob <- buildmodelsmatrix(q)
    pow2 <- 2^(seq_len(q) - 1)
    id_cf <- as.vector((cf.models.PM[,-(q+1)] > 0) %*% pow2)
    id_models <- as.vector(modelsprob[,-(q+1)] %*% pow2)

    modelsprob[, q + 1] <- sapply(id_models,
                                  function (i) sum(cf.models.PM[which(id_cf == i), q + 1]))
  } else modelsprob <- all.models.PM

  #compute inclusion probabilities (except for fixed variables)
  inclprob <- colSums(modelsprob[, -(q+1)] * modelsprob[, q + 1])

  #compute posterior probability of the dimension of the true model
  modeldim <- rowSums(modelsprob[, -(q+1)])
  probdim <- sapply(0:q, function (k) sum(modelsprob[which(modeldim == k), q + 1]))

  #HPM
  nPmax <- which.max(modelsprob[, q+1])
  hpm <- modelsprob[nPmax, -(q+1)]

  #MPM
  mpm <- numeric(q)
  mpm[which(inclprob >= 0.5)] <- 1

  return(list(modelsprob = modelsprob,
              inclprob = inclprob, probdim = probdim, hpm = hpm, mpm = mpm))
}

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
#' @seealso  Use \code{\link[MissingBVS]{missingBVS.lm}},
#' \code{\link[MissingBVS]{missingGD25}} or \code{\link[MissingBVS]{missingBVS.glm}}
#' and their Gibbs versions for creating objects of the class \code{MissingBvs}.
#'
#' @examples
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' # Use the same small model space as in the main example:
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(
#'   formula = f, data = dataS97, n.keep = 8, n.imp = 2, imp.seed = 1
#' )
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
    mod.mat <- as.data.frame(modelspostprob[1:n.keep,], row.names = 1:n.keep)
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
#' @seealso  Use \code{\link[MissingBVS]{missingBVS.lm}},
#' \code{\link[MissingBVS]{missingGD25}} or \code{\link[MissingBVS]{missingBVS.glm}}
#' and their Gibbs versions for creating objects of the class \code{MissingBvs}.
#'
#' @examples
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' # Use the same small model space as in the main example:
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' dataS97.mBVS <- missingBVS.lm(
#'   formula = f, data = dataS97, n.keep = 8, n.imp = 2, imp.seed = 1
#' )
#'
#' #Summ up the results:
#' summary(dataS97.mBVS)
#' summary(dataS97.mBVS)$summary
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
