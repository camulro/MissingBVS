#' Bayes factors and posterior probabilities via Bayesian Imputation Averaging for
#' generalized linear models
#'
#' The model space is build from a list of linear regression models proposed to explain
#' a common response. It returns the Bayes factors and posterior probabilities computed
#' through Bayesian Imputation Averaging (BIA) for generalized linear models in the
#' presence of missing data.
#'
#' Given a list of competing models, the model space is made up by them, assuming that the
#' intercept term is present in every model. The simplest one M0, can be specified (\code{null.model})
#' and must be nested in the rest. In order to implement BIA, \code{\link[MissingBVS]{missingBtest.glm}}
#' can, either perform \code{n.imp} imputations designed by \code{imp.predict.mat} and
#' \code{imp.mice.method} with the \pck{mice} package, or use user-given imputated datasets
#' by the \code{imp.datasets} argument. Hence, the posterior distribution over the model space
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
#' The prior over the model space Pr(Mi) offers three options through \code{prior.models}:
#' -"Constant" assigns the same prior probability to every model, default one.
#' -"ScottBerger" assigns the same prior probability to every different model size.
#' -"User": if chosen, user has to provide a N dimensional vector (where N the number of
#' competing models) with the model prior probabilities of each one in \code{models}
#' through \code{priorprobs}.
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
#' which is the recommended. A non-treatment of factors can be performed through
#' \code{marginal.factors}.
#'
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas, with one
#' nested in all the others. If the list is unnamed, default names are given.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link function to be
#' used in the model. If it is one of the implemented families in \pkg{BAS}:
#' \code{binomial(link = "logit")}, \code{poisson(link = "log")} and \code{Gamma(link = "log")};
#' a faster version using \pkg{BAS} log-marginal computation is performed.
#' @param null.model String for the name of the null model on \code{models}. By default,
#' the names of variables are used to identify the null. If provided, the string
#' must coincide with the one with the largest sum of squared errors and should
#' be the one with the smallest size.
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
#' @param priorprobs A N dimensional vector (being N the number of competing models)
#' defining the prior model probabilities for each one in \code{models} (if
#' \code{prior.models}= "User"; see details).
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
#' marginal likelihood of \pkg{BAS}. See \code{\link[BAS]{bas.glm}}
#' for more details.

#' @return \code{\link[MissingBVS]{MissingBtest.glm}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factors in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{If missings, list of the combined estimates for each model
#' in \code{models} fitted by \code{\link[stats]{glm}} over the \code{n.imp} imputed
#' datasets; or \code{glm} object when there are no missings}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information (when relevant)}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{family}{Family function among \code{\link[stats]{family}} used to specify
#' the error distribution and link function to be used in the model}
#' \item{weights}{Weights vector used in the glm fitting process}
#' \item{offset}{Offset vector used in the glm fitting process}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{prior.models}{Vector with \code{prior.models} and \code{prior.models.dummies}
#' chosen. If there are no factors or \code{marginal.factors} is set to \code{FALSE},
#' it saves the only argument used, \code{prior.models}}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' probabilities such as García-Donato and Paulo (2022)}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#'
#' @author Carolina Mulet and Gonzalo García-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingBvs.glm}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
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
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures
#' of g-Priors in Generalized Linear Models. Journal of the American
#' Statistical Association. 113: 1275–1287.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
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
#' #Indian Prime Diabetes Data from VIM's package
#' data.diabetes <- VIM::diabetes
#'
#' #Default choices are: BIC approximation, Constant prior and 390 imputed
#' #datasets with mice's pmm method.
#' models.list = list(M0 = Outcome ~ 1, M1 = Outcome ~ Pregnancies,
#'   M2 = Outcome ~ Glucose, M3 = Outcome ~ Insulin,
#'   M4 = Outcome ~ Pregnancies + Glucose, M5 = Outcome ~ Pregnancies + Insulin,
#'   M6 = Outcome ~ Pregnancies + Glucose + Insulin)
#'
#' diabetes.mtest <- missingBtest.glm(data = VIM::diabetes,
#'   models = models.list, family = binomial())
#'
#' #Show the results:
#' diabetes.mtest
#' }
#'
missingBtest.glm <- function (data,
                              models,
                              family = binomial(link = "logit"),
                              null.model = NULL,
                              BF.approx.method = "BIC",
                              prior.betas = "Robust",
                              prior.models = "Constant",
                              prior.models.dummies = "ScottBerger",
                              marginal.factors = TRUE,
                              priorprobs = NULL,
                              imp.mice.method = "pmm",
                              imp.predict.mat = NULL,
                              n.imp = 039E1,
                              maxit = 5,
                              parallelmice = NULL,
                              n.core = NULL,
                              imp.datasets = NULL,
                              imp.seed = runif(1,0,09011975),
                              weights = rep.int(1, nrow(data)),
                              offset = rep.int(0, nrow(data)),
                              control = glm.control(),
                              laplace = 0L) {

  #N is the number of models:
  N <- length(models)

  env <- environment()

  #Check Btest given arguments
  Btestarg.list <- checkBtestarguments(models, null.model)
  list2env(Btestarg.list, envir = env)

  #check whether the family chosen is among the options provided by BAS
  inBAS <- checkforfamily(family, BF.approx.method)

  #for the C code
  weights <- as.numeric(weights); offset <- as.numeric(offset)
  laplace <- as.integer(laplace)

  Dev <- numeric(N) #deviances for each model
  Dim <- rep.int(0,N)
  mt <- list() #list of terms for each model

  covar.list <- list() #list that contains the names of the variables in each model
  compvars <- c() #name of original competing vars
  for (i in seq_len(N)) {
    f <- as.formula(models[[i]])
    compvars <- c(compvars, attr(terms(f), "term.labels"))
    environment(f) <- env #select environment to get glm arguments
    temp <- glm(formula = f,
                data = data,
                y = TRUE, x = TRUE,
                family = family,
                weights = weights,
                offset = offset,
                control = control)

    Dev[i] <- temp$deviance

    framei <- model.frame(f, data, na.action = NULL)
    Xi <- model.matrix.rankdef(framei)
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.Dev <- sort(Dev, index.return = TRUE, decreasing = TRUE)
  #Which one acts as null model:
  nullmodel.pos <- ordered.Dev$ix[1]

  #Check null model
  if (!is.null(null.model) & nullmodel.pos != pos.user.null.model) {
      stop("The given null model does not coincide with the one with\n",
           "the largest deviance (and it should).\n")
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(models[[nullmodel.pos]]); environment(null.model) <- env
  competing.models <- seq_len(N)[-nullmodel.pos]

  #Competing vars full formula:
  full.formula <- as.formula(paste0(null.model[[2]], " ~ ",
                                    paste(unique(compvars), collapse = " + ")))

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(full.formula, null.model, data, marginal.factors)
  list2env(buildmatrices.list, envir = env)

  Dim <- Dim - p0 #model dimension (without fixed vars)

  #The response variable
  obsnotNA <- rownames(X0)
  y <- na.omit(model.response(framenull)) #response variable without missings
  n <- length(y)
  devnull <- Dev[nullmodel.pos]

  #is the response a factor?
  if (is.factor(y)) y <- y != levels(y)[1L]

  weights <- weights[as.numeric(obsnotNA)]; offset <- offset[as.numeric(obsnotNA)]

  #Compute model prior
  lprior.models <- priormodels.btest(prior.models, N, Dim, priorprobs)

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS, n, p = max(Dim),
                                              p0, y, null.model, data[obsnotNA,], family, devnull,
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

  mF <- L > 0 & marginal.factors
  #Check if factors present and if marginalization of their probabilities.
  #Define model prior and BF
  lBF.comp <- BFcomp.btest(lprior.models, prior.models.dummies, Dim, mF, positionsfac,
                           namesxnotnull, NAvars, lBF.method, X0, X.full, BF.approx.method, covar.list)

  #Posterior computation of model space defined by models list
  post.btest.list <- posterior.btest(competing.models, namesxnotnull, namesnull, covar.list,
                                     relax.nest, lprior.models, lBF.comp, nullmodel.pos, models)
  list2env(post.btest.list, envir = env)

  #Evaluate glm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesxnotnull %in% covar.list[[j]])
    if (any(namesxnotnull[namesj] %in% NAvars)) {
      fit <- list()
      for (i in 1:n.imp) {
        if (mF) {#remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation.array[,c(1:p0, setdiff(namesj, indf) + p0),i]
        } else Xi <- imputation.array[,c(1:p0, namesj + p0),i]

        z <- glm.fit(x = Xi, y = y, family = family,
                     weights = weights, offset = offset, control = control)
        z$terms <- mt[[j]]; class(z) <- "glm"; fit[[i]] <- z
      }
      modelspool[[j]] <- mice::pool(fit)
      modelspool[[j]]$call <- NULL #otherwise, Rstudio returns a warning trying to read modelspool[[j]]$call
    } else modelspool[[j]] <- glm(models[[j]], data, family = family,
                                  weights = weights, offset = offset, control = control)
  }
  modelspool[[nullmodel.pos]] <- glm(null.model, data, family = family,
                                     weights = weights, offset = offset, control = control)
  names(modelspool) <- names(models)

  result <- list()
  result$lBFi0 <- lBFi0
  result$PostProbi <- PostProbi
  result$models <- models
  result$nullmodel <- names(models)[nullmodel.pos]
  result$modelspool <- modelspool

  if (mF) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  if (!is.null(NAvars)) {
    #arguments used for imputation
    result$imp.info <- imp.info

    #save the imputed datasets for sensitivity analysis
    raw.imp.array <- serialize(imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  #glm arguments
  result$family <- family; result$weights <- weights; result$offset <- offset

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lprior.models #function used for model prior
  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models
  result$marginal.factors <- marginal.factors #whether or not factors are marginalized
  result$priorprobs <- exp(lPriorModels)
  result$call <- match.call()

  class(result) <- "MissingBtest"

  return(result)
}
