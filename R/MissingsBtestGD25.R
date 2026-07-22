#' Bayes factors and posterior probabilities for linear regression models with missing data
#'
#' The model space is build from a list of linear regression models proposed to explain
#' a common response. It returns the Bayes factors and posterior probabilities computed
#' in linear models with normally distributed covariates under the g' imputation prior of
#' García-Donato et al. (2025).
#'
#' Given a list of competing models, the model space is made up by them, assuming that the
#' intercept term is present in every model. The simplest one, the null, is fixed to be the
#' one with just the intercept term. If not provided by \code{models}, it is automatically added.
#' \code{\link[MissingBVS]{missingBtestGD25}} performs for the regressors \code{n.imp}
#' imputations from an objective Bayesian multivariate normal model (Sun and Berger, 2006),
#' along with their variance-covariance matrices. Hence, the posterior distribution
#' over the model space is given through Bayes' theorem:
#'
#' Pr(Mi | \code{data}) = Pr(Mi) * Bi / C,
#'
#' where Bi is the Bayes factor (BF) of Mi to M0 derived a the expected value of g'BF of
#' García-Donato et al. (2025) for normal random regressors, Pr(Mi) is the prior
#' probability of Mi and C is the normalizing constant. The g'BFs are computed
#' with \code{\link[MissingBVS]{BF.miss.X}}, which resemble g-prior BFs (Zellner, 1986)
#' for random covariates, for each pair of imputed dataset and  variance-covariance
#' matrix simulated. The integral of g'BF is approximated with a Monte Carlo scheme.
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
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas, with one
#' nested in all the others. If the list is unnamed, default names are given.
#' @param prior.models Model prior distribution over the covariates and/or factors
#' model space (to be literally specified). Possible choices are "Constant",
#' "ScottBerger" and "User" (see details).
#' @param priorprobs A N dimensional vector (being N the number of competing models)
#' defining the prior model probabilities for each one in \code{models} (if
#' \code{prior.models}= "User"; see details).
#' @param n.imp Number of imputed datasets for model posterior computation.
#' @param imp.seed Seed for imputation.

#' @return \code{\link[MissingBVS]{MissingGD25Btest}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factors in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{If missings, list of the combined estimates for each model
#' in \code{models} fitted by \code{\link[stats]{lm}} over the \code{n.imp} imputed
#' datasets; or \code{lm} object when there are no missings}
#' \item{imp.info}{List of arguments used for the imputation step and other information}
#' \item{compress.imp.array}{Compressed array of imputed datasets}
#' \item{logprior.models}{Function used to compute the log-prior over the model space
#' defined by covariates and/or factors}
#' \item{prior.models}{Argument chosen for \code{prior.models}}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{call}{The \code{call} to the function}
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
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' Sun, D. and Berger, J. O. (2006) Objective Bayesian analysis for the multivariate
#' normal model. In Bernardo, J. M., Bayarri, M. J., Berger, J. O., Dawid, A. P.,
#' Heckerman, D., Smith, A. F. M., and West, M. (eds.), Proc. Valencia / ISBA 8th
#' World Meeting on Bayesian statistics. Oxford university Press. MR2433206. 16
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
                              n.imp = 039E1,
                              imp.seed = runif(1,0,09011975)) {

  #N is the number of models:
  N <- length(models)

  #Check Btest given arguments
  Btestarg.list <- checkBtestarguments(models, NULL)
  list2env(Btestarg.list, envir = environment())

  namesm <- names(models)

  Dim <- rep.int(0, N)
  lBFi0 <- lPriorModels <- PostProbi <- numeric(N)
  mt <- list() #list of terms for each model

  covar.list <- list() #list that contains the names of the covariates in each model
  for (i in seq_len(N)) {
    formula <- as.formula(models[[i]])
    #Check for numeric covariates
    aux <- model.frame(formula, data)
    isnum <- sapply(aux, is.numeric)
    # isint <- sapply(aux, is.integer)
    # if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
    if (sum(isnum) < dim(aux)[2]) {
      stop("This method is only for continuous covariates.\nTry missingBtest.lm instead.\n")
    }
    temp <- lm(formula = formula,
               data = data,
               y = TRUE, x = TRUE)

    Dim[i] <- length(temp$coefficients)
    covar.list[[i]] <- dimnames(temp$x)[[2]]
    mt[[i]] <- temp$terms
  }
  nullmodel.pos <- which(covar.list == "(Intercept)")
  response <- as.character(as.formula(models[[1]])[[2]])

  #Check if null model is one of the competing models:
  if (length(nullmodel.pos) > 0) {
    competing.models <- setdiff(seq_len(N), nullmodel.pos)
    #the null model has to be the one with the intercept
    null.model <- as.formula(models[[nullmodel.pos]])
  } else {
    competing.models <- seq_len(N) #null.model is not provided by user
    #the null model has to be the one with the intercept
    null.model <- as.formula(paste(response, "~ 1"))

    nullmodel.pos <- N + 1
    Dim[nullmodel.pos] <- 1
    models[[nullmodel.pos]] <- null.model
    cat("Null model", paste(response, "~ 1"), "added to the list of competing models.\n")
  }

  #Full design matrix
  formula <- as.formula(paste0(null.model[[2]], "~ ."))
  framefull <- model.frame(formula, data, na.action = NULL)
  X.full <- framefull[,-1] #remove intercept
  namesx <- dimnames(X.full)[[2]]
  p <- length(namesx) #Number of covariates to select from

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Compute model prior
  lprior.models <- priormodels.btest(prior.models, N, Dim, priorprobs)

  #Check methods and options
  BF.miss.aux <- function (X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
                                                           y = y, SS0 = SS0,
                                                           n = n, k)

  #check for missings
  NAvars <- checkformissings(y = framefull[,1], X.full = X.full[obsnotNA,])

  #Imputation of missing data
  if (p*n > 10000 | n.imp > 039E1) cat("Imputation step could take a while.\n",
                                       "Consider reducing the number of imputed datasets if that is the case.\n")

  cat("Performing imputation of missing data with Garcia-Donato's 2025 method.\n",
      "Please wait . . . \n")
  imputation.list <- MC.imputation(X = X.full, nMC = n.imp, seed = imp.seed)

  #remove observations with missings on the response
  imputation.list$rX.imput <- imputation.list$rX.imput[obsnotNA,,, drop = FALSE]
  if (n.imp > 1) {
    #function to compute log(BFa0) for a given model with García-Donato's 2025 method
    lBF.method <- function (model) lBF.miss(model,
                                            imputation.list = imputation.list,
                                            BF.miss.aux = BF.miss.aux,
                                            n = n, nMC = n.imp)
  } else lBF.method <- function (model) BF.miss.aux(X.center = imputation.list$rX.imput[,model,],
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

  #Names to save
  null.name <- paste0("(", response, " ~ 1)")
  if (namesm[nullmodel.pos] == "") {
    namesm[nullmodel.pos] <- null.name
    names(lBFi0) <- paste0(namesm, ".to.(", null.name, ")")
    names(PostProbi) <- names(lPriorModels) <- c(namesm[-nullmodel.pos], null.name)
  } else {
    names(lBFi0) <- paste0(namesm, ".to.", namesm[nullmodel.pos])
    names(PostProbi) <- names(lPriorModels) <- namesm
  }

  #Evaluate lm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesx %in% covar.list[[j]])
    if (any(namesx[namesj] %in% NAvars)) {
      fit <- list()
      for (i in 1:n.imp) {
        Xi <- cbind(1, imputation.list$rX.imput[,namesj,i])
        colnames(Xi) <- c("(Intercept)", namesx[namesj])
        z <- lm.fit(x = Xi, y = y)
        z$terms <- mt[[j]]; class(z) <- "lm"; fit[[i]] <- z
      }
      modelspool[[j]] <- mice::pool(fit)
      modelspool[[j]]$call <- NULL #otherwise, Rstudio returns a warning trying to read modelspool[[j]]$call
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
  result$imp.info <- list(n.imp = n.imp, imp.seed = imp.seed)
  #save the imputed datasets for sensitivity analysis
  raw.imp.array <- serialize(imputation.array, NULL)
  result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$logprior.models <- lprior.models #function used for model prior
  result$prior.models <- prior.models #function used for model prior
  result$priorprobs <- exp(lPriorModels)
  result$call <- match.call()

  class(result) <- "MissingBtest"

  return(result)
}
