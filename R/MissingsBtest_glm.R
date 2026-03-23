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
#' model with fixed covariates; next q components correspond to the q prior probabilities
#' of the possible model dimensions.
#'
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas,
#' with one nested in all the others. If the list is unnamed, default names are
#' given.
#' @param null.model Name of the null model. If \code{NULL}, the names of covariates
#' in \code{models} are used to identify the null. If provided, \code{null.model}
#' must coincide with the one with the largest sum of squared errors and should
#' be the one with the smallest dimension.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. If it is one of the implemented families in
#' \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}; a faster version
#' using BAS logmarginal computation is performed.
#' @param BF.approx.method Method used to approximate Bayes factors with missing
#' data (to be literally specified). Possible choices include "BIC", "TBF" and
#' "gprior" (see details).
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices are the ones implemented
#' in \pkg{BayesVarSel}: "Robust", "Liangetal", "gZellner", "ZellnerSiow",
#' "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified).
#' Possible choices are "Constant", "ScottBerger" and "User" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed covariates)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details).
#' @param parallelmice Logital to indicate whether or not to use parallel
#' \code{\link[mice]{mice}} imputation. If \code{NULL}, automatically performs
#' the parallel mice imputation if the number of imputations or the dataset are
#' big enough.
#' @param n.core See \code{\link[mice]{futuremice}} for details.
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 10} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param imp.mice.method Method for mice's imputation.
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

#' @return \code{\link[MissingBVS]{MissingBtest.glm}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factor in logaritmic scale of each model to the null.}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}.}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model.}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{prior.models}{Function used to compute the prior over the model space}
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
#' @examples #To be completed
#'
missingBtest.glm <- function (data,
                              models,
                              null.model = NULL,
                              family = binomial(link = "logit"),
                              BF.approx.method = "BIC",
                              prior.betas = "Robust", #if BF.approx.method = "gprior"
                              prior.models = "Constant",
                              priorprobs = NULL, #needed if prior.models = "User"
                              imp.mice.method = "pmm", #mice's default
                              parallelmice = NULL,
                              n.core = NULL,
                              imp.time.test = TRUE,
                              n.imp = 039E1,
                              imp.seed = runif(1,0,09011975), #seed for the imputation
                              #glm.fit arguments:
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
  inBAS <- checkforfamily(family)

  #for the C code
  weights <- as.numeric(weights)
  offset <- as.numeric(offset)
  laplace <- as.integer(laplace)

  Dev <- rep(0, N) #deviances for each model
  Dim <- rep(0, N)
  lBFi0 <- rep(0, N)
  lPriorModels <- rep(0, N)
  PostProbi <- rep(0, N)

  #list that contains the names of the covariates in each model
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
    Dim[i] <- length(temp$coefficients)
    covar.list[[i]] <- dimnames(temp$x)[[2]]
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

  #Response and fixed vars for imputation
  auxnull <- model.frame(null.model, data, na.action = NULL)
  namesnull.toimp <- dimnames(auxnull)[[2]][-1] #name of fixed variables to imputation

  #Missing model matrix of fixed covariates
  X0 <- model.matrix(null.model, auxnull)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of covariates to select from
  Dim <- Dim - p0 #model dimension (without fixed cov)

  #Full design matrix for imputation
  formula <- as.formula(paste0(null.model[[2]], "~ ."))
  auxfull <- model.frame(formula, data, na.action = NULL)
  namesx.toimp <- dimnames(auxfull)[[2]][-1] #name of variables to imputation
  namesxnotnull.toimp <- namesx.toimp[namesx.toimp %notin% namesnull.toimp]
  X.toimp <- data[,c(namesnull.toimp, namesxnotnull.toimp)] #design matrix with missing data with fixed cov

  #Model matrix data with missings
  X.full <- model.matrix(formula, auxfull)
  namesx <- dimnames(X.full)[[2]]
  namesxnotnull <- namesx[namesx %notin% namesnull]
  X.full <- X.full[, namesxnotnull]
  p <- dim(X.full)[2] #Number of covariates to select from

  #The response variable
  obsnotNA <- rownames(na.omit(auxnull)) #response variable without missings
  y <- as.numeric(auxnull[obsnotNA, 1])
  n <- length(y)
  devnull <- Dev[nullmodel.pos]

  #check for missings and define covariates with NAs
  NAvars <- checkformissings.glm(y = auxnull[,1], X0, X.full, obsnotNA)

  #Check model priors chosen and define the function to be used
  if (prior.models %notin% c("ScottBerger", "Constant", "User")) {
    stop("Only priors 'ScottBerger', 'Constant' and 'User' supported.\n")
  }
  switch (prior.models, #change the string for the corresponding function
          Constant = {prior.models <- function (modeli) 1/length(models)},
          ScottBerger = {prior.models <- function (modeli) 1/length(unique(Dim)) /
                                                           sum(Dim == Dim[modeli])},
          User = {
            if (is.null(priorprobs)) {
              stop("User prior selected but no prior probabilities provided.\n")
            }
            if (length(priorprobs) != N) {
              stop("Vector of prior probabilities with incorrect length.\n")
            }
            if (sum(priorprobs < 0) > 0) {
              stop("Prior probabilities must be positive.\n")
            }
            prior.models <- function(modeli) priorprobs[modeli]}
  )

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS,
                                              n, p = max(Dim), p0, y, null.model,
                                              data, family, devnull,
                                              weights, offset, control, laplace)

  #Imputation of missing data
  if (is.null(parallelmice)) {
    if (n.imp > 120 | n*p > 50000) {
      parallelmice <- TRUE #faster
    } else parallelmice <- FALSE
  }

  if (imp.time.test & (n*p > 10000 | n.imp > 039E1)) {
    #test imputation time
    cat("Time test . . . \n")
    time.test <- mice.imputation(X = X.toimp,
                                 formula,
                                 imp.mice.method = imp.mice.method,
                                 parallel = parallelmice,
                                 n.core = n.core,
                                 time.test = TRUE)

    estim.time <- time.test * n.imp / (60 * 30) #30 imputed datasets used to time
    cat("The whole imputation can take ", estim.time,
        "minutes (approx.) to run.\n Do you want to continue? (y/n)\n")
    if (tolower(readline()) != "y") {
      if (!parallelmice) {
        cat("Do you want to faster imputation running a parallel version of mice? (y/n)\n")
        if (tolower(readline()) == "y") {
          parallelmice <- TRUE
        } else stop("Reduce the number of imputed datasets.\n")
      } else stop("Reduce the number of imputed datasets.\n")
    }
  }

  cat("Performing imputation of missing data with mice's", imp.mice.method)
  if (parallelmice) cat(" parallel")
  cat(" method.\n", "Please wait . . . \n")

  imputation.array <-  mice.imputation(X = X.toimp,
                                       formula,
                                       n.imp = n.imp,
                                       imp.mice.method = imp.mice.method,
                                       seed = imp.seed,
                                       parallel = parallelmice,
                                       n.core = n.core)

  #remove observations with missings on the response
  imputation.array <- imputation.array[obsnotNA,,]
  #function to compute log(BFa0) for a given model as an average of BF computed
  #by BF.approx.method over the imputed datasets
  lBF.method <- function (model) lBF.approx(model,
                                            imputation.array = imputation.array,
                                            BF.approx.method = BF.approx.method,
                                            p0 = p0, n.imp = n.imp)

  for (i in competing.models){
    modeli <- which(namesxnotnull %in% covar.list[[i]])

    #check whether the null is nested in the other ones
    if (!relax.nest & any(namesnull %notin% covar.list[[i]])) {
      stop(paste0("The simplest (null) model may not be nested in all the others.\n",
                  "Please define explicitly the null model if it is the case.\n"))
    }

    #check if there are NAs in the model considered to save computation time
    if (any(covar.list[[i]] %in% NAvars)) {
      lBFi0[i] <- lBF.method(model = modeli)
    } else { #if there are no missings, compute the BF by the method selected
      X.i <- cbind(X0[obsnotNA,], X.full[obsnotNA, modeli])
      lBFi0[i] <- BF.approx.method(k = Dim[i], X = as.matrix(X.i))
    }
    lPriorModels[i] <- log(prior.models(i))
  }

  cat("\n")
  lPriorModels[nullmodel.pos] <- log(prior.models(nullmodel.pos))
  lBFi0[nullmodel.pos] <- 0
  C <- sum(exp(lBFi0 + lPriorModels))
  PostProbi <- exp(lBFi0 + lPriorModels - log(C))

  names(lBFi0) <-
    paste(names(models), ".to.", names(models)[nullmodel.pos], sep = "")
  names(PostProbi) <- names(models)
  names(lPriorModels) <- names(models)

  result <- list()
  result$lBFi0 <- lBFi0
  result$PostProbi <- PostProbi
  result$models <- models
  result$nullmodel <- names(models)[pos.user.null.model]

  #arguments used for imputation
  result$imp.args <- list(parallelmice = parallelmice,
                          imp.mice.method = imp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  #save the imputed datasets for sensitivity analysis
  # raw.imp.array <- serialize(imputation.array, NULL)
  # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$prior.models <- prior.models #function used for model prior
  result$priorprobs <- exp(lPriorModels)

  class(result) <- "MissingBtest"

  return(result)
}
