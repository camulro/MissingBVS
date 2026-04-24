#' Bayes factors with missing data and posterior probabilities for generalized
#' linear models
#'
#' It computes the Bayes factors and posterior probabilities in the presence
#' of missing data of a list of generalized linear models proposed to explain a
#' common response. See \code{\link[MissingBVS]{MissingBvs.glm}} for more details.
#'
#' The prior over the model space Pr(Mi) offers three possibilities:
#' "Constant" assigns the same prior probability to every model, default choice.
#' "ScottBerger" assigns the same prior probability to every possible model with
#' the same size and, therefore, accounts for multiplicity issues
#' (Scott and Berger 2010).
#' "User" (see below).
#'
#' If \code{prior.models}="User" is chosen, user has to provide a q+1 dimensional
#' parameter vector, where q is the number of diferent sizes among \code{models},
#' with the model dimension prior probabilities through \code{priorprobs}.
#' The first component of \code{priorprobs} must contain the probability of the
#' model with fixed variables; next q components correspond to the q prior probabilities
#' of the possible model dimensions.
#'
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas,
#' with one nested in all the others. If the list is unnamed, default names are
#' given.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. If it is one of the implemented families in
#' \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}; a faster version
#' using BAS logmarginal computation is performed.
#' @param null.model Name of the null model. If \code{NULL}, the names of variables
#' in \code{models} are used to identify the null. If provided, \code{null.model}
#' must coincide with the one with the largest sum of squared errors and should
#' be the one with the smallest dimension.
#' @param BF.approx.method Method used to approximate Bayes factors with missing
#' data (to be literally specified). Possible choices include "BIC", "TBF" and
#' "gprior" (see details).
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices are the ones implemented
#' in \pkg{BayesVarSel}: "Robust", "Liangetal", "gZellner", "ZellnerSiow",
#' "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified).
#' Possible choices are "Constant", "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the model space of the
#' factor levels (to be literally specified). Possible choices are "Constant" and
#' "ScottBerger" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed variables)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details).
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
#' \item{lBFi0}{Bayes factor in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{List of objects of class \code{\link[mice]{mipo}} that
#' combine the estimates for each model on \code{models} fitted by
#' \code{\link[stats]{glm}} over the \code{n.imp} (if \code{> 1}) imputed
#' datasets, see \code{\link[mice]{pool}} for details; and \code{glm} object
#' when there are no missings}
#' \item{positions}{\code{matrix} with L rows and p plus the number of dummies
#' resulting from factors columns - L, where L is the number of factors, with 1
#' if the column dummy corresponds to the row factor and 0 otherwise}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{prior.models}{Function used to compute the prior over the model space}
#' \item{logprior.models.dumm}{Function used to compute the log-prior over the
#' model space of the factor levels}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
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
#' @keywords package
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
                              priorprobs = NULL,
                              imp.mice.method = "pmm",
                              imp.predict.mat = NULL,
                              parallelmice = NULL,
                              n.core = NULL,
                              imp.time.test = TRUE,
                              n.imp = 039E1,
                              imp.seed = runif(1,0,09011975),
                              weights = rep.int(1, nrow(data)),
                              offset = rep.int(0, nrow(data)),
                              control = glm.control(),
                              laplace = 0L) {

  #N is the number of models:
  N <- length(models)

  if (!is.list(models)) stop("Argument models should be a list.\n")

  #If competing models come wihtout a name, give one by default:
  if (is.null(names(models))){
    if (!is.null(null.model)) stop(paste0("Please provide a name for the competing models.\n",
                                          "The null model must be in that list.\n"))
    names(models) <- paste("model", seq_len(N), sep="")
  }

  #Check if the given null model is one of the competing models:
  if (!is.null(null.model)){
    relax.nest = TRUE
    pos.user.null.model <- which(null.model == names(models))
    if (length(pos.user.null.model) == 0) {
      stop("The null model provided is not in the list of competing models.\n")
    }
  } else relax.nest = FALSE

  #check whether the family chosen is among the options provided by BAS
  inBAS <- checkforfamily(family, BF.approx.method)

  #for the C code
  weights <- as.numeric(weights)
  offset <- as.numeric(offset)
  laplace <- as.integer(laplace)

  Dev <- numeric(N) #deviances for each model
  Dim <- rep(0L,N)
  mt <- list() #list of terms for each model

  #list that contains the names of the variables in each model
  covar.list <- list()
  for (i in seq_len(N)) {
    #select environment to get glm arguments
    formulai <- as.formula(models[[i]])
    environment(formulai) <- environment()
    temp <- glm(formula = formulai,
                data = data,
                y = TRUE, x = TRUE,
                family = family,
                weights = weights,
                offset = offset,
                control = control)

    Dev[i] <- temp$deviance

    framei <- model.frame(formulai, data, na.action = NULL)
    Xi <- model.matrix.rankdef(framei)
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.Dev <- sort(Dev, index.return = TRUE, decreasing = TRUE)
  #Which acts as null model:
  nullmodel.pos <- ordered.Dev$ix[1]

  if (!is.null(null.model)){
    if (nullmodel.pos != pos.user.null.model){
      stop(paste0("The given null model does not coincide with the one with\n",
                  "the largest deviance (and it should).\n"))
    }
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(models[[nullmodel.pos]])
  environment(null.model) <- environment()
  competing.models <- seq_len(N)[-nullmodel.pos]

  #Competing vars
  compvars <- unique(unlist(covar.list))[-1] #remove intercept
  full.formula <- as.formula(paste0(null.model[[2]], " ~ ", paste(compvars, collapse = " + ")))

  #Build matrices and objects needed later on
  buildmatrices.list <- buildmatrices(full.formula, null.model, data)
  list2env(buildmatrices.list, envir = environment())

  Dim <- Dim - p0 #model dimension (without fixed vars)

  #Check arguments and define the functions to compute prior model probabilities
  priormodels.list <- priormodels.btest(prior.models, prior.models.dummies,
                                        N, Dim, priorprobs)
  list2env(priormodels.list, envir = environment())

  #The response variable
  obsnotNA <- rownames(X0)
  y <- model.response(framenull) #response variable without missings
  n <- length(y)
  devnull <- Dev[nullmodel.pos]

  #is the response a factor?
  if (is.factor(y)) y <- y != levels(y)[1L]
  y[weights == 0] <- 0

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS,
                                              n, p = max(Dim), p0, y, null.model,
                                              data, family, devnull,
                                              weights, offset, control, laplace)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)

  #Imputation step
  if (!is.null(NAvars)) {
    buildimputation.list <- buildimputation(NAvars, full.formula, data, imp.predict.mat,
                                            n.imp, n, q, p0, imp.time.test, imp.mice.method, imp.seed,
                                            parallelmice, n.core, obsnotNA, ordvars, BF.approx.method)
    list2env(buildimputation.list, envir = environment())
  }

  #Posterior computation of model space defined by models list
  post.btest.list <- posterior.btest(competing.models, namesxnotnull, namesnull, covar.list,
                                     positionsfac, NAvars, lBF.method, lprior.models.dummies,
                                     X0, X.full, prior.models, BF.approx.method,
                                     nullmodel.pos, relax.nest, Dim, models)
  list2env(post.btest.list, envir = environment())

  #Evaluate glm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesxnotnull %in% covar.list[[j]])
    if (any(namesxnotnull[namesj] %in% NAvars)) {
      if (n.imp > 1) {
        fit <- list()
        for (i in 1:n.imp) {
          #remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation.array[,c(1:p0, setdiff(namesj, indf)  + p0),i]
          z <- glm.fit(x = Xi, y = y, family = family,
                       weights = weights, offset = offset, control = control)
          z$terms <- mt[[j]]
          class(z) <- "glm"

          fit[[i]] <- z
        }
        modelspool[[j]] <- mice::pool(fit)
      } else {
        #compute glm.fit for the unique imputation
        Xi <- imputation.array[,c(1:p0, setdiff(namesj, indf)  + p0)]
        modelspool[[j]] <- glm.fit(x = Xi, y = y, family = family,
                                   weights = weights, offset = offset, control = control)
      }
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

  if (L > 0) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  if (!is.null(NAvars)) {
    #arguments used for imputation
    result$imp.args <- imp.args

    #save the imputed datasets for sensitivity analysis
    # raw.imp.array <- serialize(imputation.array, NULL)
    # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$prior.models <- prior.models #function used for model prior
  result$logprior.models.dumm <- lprior.models.dummies
  result$priorprobs <- exp(lPriorModels)

  class(result) <- "MissingBtest"

  return(result)
}
