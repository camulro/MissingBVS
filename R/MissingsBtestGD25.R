#' Bayes factors with missing data and posterior probabilities for lm and normal
#' regressos
#'
#' It computes the Bayes factors and posterior probabilities in the presence
#' of missing data of a list of linear regression models proposed to explain a
#' common response. See \code{\link[MissingBVS]{MissingGD25}} for more details.
#'
#' Null model is fixed to be the one with just the intercept term. If not provided
#' by \code{models}, it is added at the end of the list.
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
#' null model; next q components correspond to the q prior probabilities
#' of the possible model dimensions.
#'
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas,
#' with one nested in all the others. If the list is unnamed, default names are
#' given.
#' @param prior.models Prior distribution over the model space (literally specified).
#' Possible choices are "Constant", "ScottBerger" and "User" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed covariates)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details).
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 30} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param initialimp.mice.method Method for mice's imputation.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param imp.seed Seed for imputation.

#' @return \code{\link[MissingBVS]{MissingGD25Btest}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factor in logaritmic scale of each model to the null.}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}.}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model.}
#' \item{modelspool}{List of objects of class \code{\link[mice]{mipo}} that
#' combine the estimates for each model on \code{models} fitted by
#' \code{\link[stats]{lm}} over the \code{n.imp} (if \code{> 1}) imputed
#' datasets, see \code{\link[mice]{pool}} for details; and \code{lm} object
#' when there are no missings}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{prior.models}{Function used to compute the prior over the model space}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#'
#' @author Carolina Mulet and Gonzalo García-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGD25}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @keywords package
#'
#' @examples
#' \dontrun{
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' #Default choices are: Constant prior and 390 imputed datasets.
#' models.list = list(M0 = Ozone ~ 1, M1 = Ozone ~ Solar.R,
#'   M2 = Ozone ~ Wind, M3 = Ozone ~ Temp, M4 = Ozone ~ Solar.R + Wind,
#'   M5 = Ozone ~ Solar.R + Temp, M6 = Ozone ~ Solar.R + Wind + Temp)
#'
#' airq.mtest <- missingBtestGD25(data = airquality, models = models.list)
#'
#' #Show the results:
#' airq.mtest
#'
#' }
#'
missingBtestGD25 <- function (data,
                              models,
                              prior.models = "Constant",
                              priorprobs = NULL,
                              imp.time.test = TRUE,
                              initialimp.mice.method = "norm",
                              n.imp = 039E1,
                              imp.seed = runif(1,0,09011975)) {

  #N is the number of models:
  N <- length(models)

  if (!is.list(models)) stop("Argument models should be a list.\n")

  namesm <- names(models)
  #If competing models come wihtout a name, give one by default:
  if (is.null(namesm)){
    namesm <- paste("model", 1:N, sep="")
  }

  Dim <- rep(0L, N)
  lBFi0 <- lPriorModels <- PostProbi <- numeric(N)
  mt <- list() #list of terms for each model

  #list that contains the names of the covariates in each model
  covar.list <- list()
  for (i in seq_len(N)) {
    formula <- as.formula(models[[i]])
    #Check for numeric covariates
    aux <- model.frame(formula, data)
    isnum <- sapply(aux, is.numeric)
    # isint <- sapply(aux, is.integer)
    # if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
    if (sum(isnum) < dim(aux)[2]) {
      stop("This method is only for continuous covariates.\n")
    }
    temp <- lm(formula = formula,
               data = data,
               y = TRUE, x = TRUE)

    Dim[i] <- length(temp$coefficients)
    covar.list[[i]] <- dimnames(temp$x)[[2]]
    mt[[i]] <- temp$terms
  }
  nullmodel.pos <- which(covar.list == "(Intercept)")
  #Check if null model is one of the competing models:
  if (length(nullmodel.pos) > 0) {
    competing.models <- (seq_len(N))[-nullmodel.pos]
    null.model <- as.formula(models[[nullmodel.pos]])
    #the null model has to be the one with the intercept
  } else {
    competing.models <- seq_len(N) #null.model is not provided by user
    null.model <- as.formula(paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""))
    #the null model has to be the one with the intercept
    nullmodel.pos <- N + 1
    Dim[nullmodel.pos] <- 1
    models[[nullmodel.pos]] <- null.model
    cat("Null model", paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""),
        "added to the list of competing models.\n")
  }

  #Full design matrix
  formula <- as.formula(paste0(null.model[[2]], "~ ."))
  framefull <- model.frame(formula, data, na.action = NULL)
  X.full <- framefull[,-1] #remove intercept
  namesx <- dimnames(X.full)[[2]]
  p <- length(namesx) #Number of covariates to select from

  #Check arguments and define the functions to compute prior model probabilities
  prior.models <- priormodels.btest(prior.models, NULL, N, Dim, priorprobs)

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Check methods and options
  BF.miss.aux <- function (X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
                                                           y = y, SS0 = SS0,
                                                           n = n, k)

  #check for missings
  NAvars <- checkformissings(y = framefull[,1], X.full = X.full[obsnotNA,])

  #Imputation of missing data
  if (imp.time.test & (p*n > 10000 | n.imp > 039E1)) {
    #test imputation time
    cat("Time test . . . \n")
    time.test <- MC.imputation(X = X.full,
                               time.test = TRUE)

    estim.time <- time.test * n.imp / (60 * 30) #30 imputed datasets used to test
    cat("The whole imputation would take ", estim.time,
        "minutes (approx.) to run.\n Do you want to continue? (y/n)\n")
    if (tolower(readline()) != "y") {
      stop("Reduce the number of imputed datasets.\n")
    }
  }
  cat("Performing imputation of missing data with Garcia-Donato's 2025 method.\n",
      "Please wait . . . \n")
  imputation.list <- MC.imputation(X = X.full,
                                   nMC = n.imp,
                                   seed = imp.seed,
                                   initialimp.mice.method = initialimp.mice.method)

  #remove observations with missings on the response
  imputation.list$rX.imput <- imputation.list$rX.imput[obsnotNA,,]
  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model with García-Donato's 2025 method
    lBF.method <- function (model) lBF.miss(model,
                                            imputation.list = imputation.list,
                                            BF.miss.aux = BF.miss.aux,
                                            n = n, nMC = n.imp)
  } else lBF.method <- function (model) BF.miss.aux(X.center = imputation.list$rX.imput[,model],
                                                    Sigma11 = imputation.list$rSigma[model, model,],
                                                    k = length(model))

  for (i in competing.models){
    modeli <- which(namesx %in% covar.list[[i]])

    lBFi0[i] <- lBF.method(model = modeli)
    lPriorModels[i] <- log(prior.models(i))
  }
  lPriorModels[nullmodel.pos] <- log(prior.models(nullmodel.pos))
  lBFi0[nullmodel.pos] <- 0
  C <- sum(exp(lBFi0 + lPriorModels))
  PostProbi <- exp(lBFi0 + lPriorModels - log(C))

  if (namesm[nullmodel.pos] == "") {
    namesm[nullmodel.pos] <-
      paste("(", as.formula(models[[1]])[[2]], " ~ 1)", sep="")
    names(lBFi0) <-
      paste(namesm, ".to.(",
            paste(as.formula(models[[1]])[[2]], " ~ 1)", sep=""), sep = "")
    names(PostProbi) <- c(namesm[-nullmodel.pos],
                          paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""))
    names(lPriorModels) <- names(PostProbi)
  } else {
    names(lBFi0) <-
      paste(namesm, ".to.", namesm[nullmodel.pos], sep = "")
    names(PostProbi) <- namesm
    names(lPriorModels) <- namesm
  }

  #Evaluate lm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesx %in% covar.list[[j]])
    if (any(namesx[namesj] %in% NAvars)) {
      if (n.imp > 1) {
        fit <- list()
        for (i in 1:n.imp) {
          Xi <- cbind(1, imputation.list$rX.imput[,namesj,i])
          colnames(Xi) <- c("(Intercept)", namesx[namesj])
          z <- lm.fit(x = Xi, y = y)
          z$terms <- mt[[j]]
          class(z) <- "lm"

          fit[[i]] <- z
        }
        modelspool[[j]] <- mice::pool(fit)
      } else {
        #compute lm.fit for the unique imputation
        Xi <- cbind(1, imputation.list$rX.imput[,namesj])
        colnames(Xi) <- c("(Intercept)", namesx[namesj])
        modelspool[[j]] <- lm.fit(x = Xi, y = y)
      }
    } else modelspool[[j]] <- lm(models[[j]], data)
  }
  modelspool[[nullmodel.pos]] <- lm(null.model, data)
  names(modelspool) <- namesm

  result <- list()
  result$lBFi0 <- lBFi0
  result$PostProbi <- PostProbi
  result$models <- models
  result$nullmodel <- namesm[nullmodel.pos]
  result$modelspool <- modelspool

  #arguments used for imputation
  result$imp.args <- list(initialimp.mice.method = initialimp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  #save the imputed datasets for sensitivity analysis
  # raw.imp.array <- serialize(imputation.array, NULL)
  # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$prior.models <- prior.models #function used for model prior
  result$priorprobs <- exp(lPriorModels)

  class(result) <- "MissingBtest"

  return(result)
}
