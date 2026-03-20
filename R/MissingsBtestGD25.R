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
#' of the imputation process with \code{n.imp = 10} if the number of variables or
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
#' @examples #To be completed
#'
missingBtestGD25 <- function (data,
                              models,
                              prior.models = "Constant",
                              priorprobs = NULL, #needed if prior.models = "User"
                              imp.time.test = TRUE,
                              initialimp.mice.method = "pmm", #mice's default
                              n.imp = 039E1,
                              imp.seed = runif(1,0,09011975)) { #seed for the imputation

  #N is the number of models:
  N <- length(models)

  if (!is.list(models)) stop("Argument models should be a list.\n")

  #If competing models come wihtout a name, give one by default:
  if (is.null(names(models))){
    names(models) <- paste("model", 1:N, sep="")
  }

  Dim <- rep(0, N)
  lBFi0 <- rep(0, N)
  lPriorModels <- rep(0, N)
  PostProbi <- rep(0, N)

  #list that contains the names of the covariates in each model
  covar.list <- list()
  for (i in seq_len(N)) {
    formula <- as.formula(models[[i]])
    #Check for numeric covariates
    aux <- model.frame(formula, data)
    isnum <- sapply(aux, is.numeric)
    isint <- sapply(aux, is.integer)
    if (sum(isnum) < dim(aux)[2] | sum(isint) > 0) {
      stop("This method is only for continuous covariates.\n")
    }
    temp <- lm(formula = formula,
               data = data,
               y = TRUE, x = TRUE)

    Dim[i] <- length(temp$coefficients)
    covar.list[[i]] <- dimnames(temp$x)[[2]]
  }
  pos.user.null.model <- which(covar.list == "(Intercept)")
  #Check if null model is one of the competing models:
  if (length(pos.user.null.model) > 0) {
    competing.models <- (seq_len(N))[-pos.user.null.model]
    null.model <- as.formula(paste(as.formula(models[[pos.user.null.model]])[[2]], " ~ 1", sep=""))
    #the null model has to be the one with the intercept
  } else {
    competing.models <- seq_len(N) #null.model is not provided by user
    null.model <- as.formula(paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""))
    #the null model has to be the one with the intercept
    pos.user.null.model <- N + 1
    Dim[pos.user.null.model] <- 1
    models[[pos.user.null.model]] <- null.model
    cat("Null model", paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""),
        "added to the list of competing models.\n")
  }

  #Evaluate the null model:
  lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)

  #The response variable
  y <- lmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y)
  SS0 <- crossprod(lmnull$residuals) #SSE of the null model

  #Full design matrix
  formula <- as.formula(paste0(null.model[[2]], "~ ."))
  auxfull <- model.frame(formula, data, na.action = NULL)
  X.full <- model.matrix(formula, auxfull)[,-1]
  namesx <- dimnames(X.full)[[2]]
  p <- dim(X.full)[2] #Number of covariates to select from

  #check for missings
  checkformissings.lm(y = auxfull[,1], X.full = X.full)

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
            if (length(priorprobs) != length(models)) {
              if (pos.user.null.model == N + 1) {
                stop(paste0("User prior selected but the length of prior probabilities is not correct (", length(models),").\n",
                            "Make sure to provide a prior for the null model, the one with the intercept.\n"))
              } else stop(paste0("User prior selected but the length of prior probabilities is not correct (", length(models),").\n"))
            }
            if (sum(priorprobs < 0) > 0) {
              stop("Prior probabilities must be positive.\n")
            }
            prior.models <- function(modeli) priorprobs[modeli]}
  )

  #Check methods and options
  cat("Be careful, this method is only for normally distributed covariates.\n",
      "Do you want to continue? (y/n)\n")
  if (tolower(readline()) != "y") {
    stop("Try the missingBtest.lm function instead.\n")
  }
  BF.miss.aux <- function (X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
                                                           y = y, SS0 = SS0,
                                                           n = n, k)

  #Imputation of missing data
  if (imp.time.test & (p > 20 | n.imp > 039E1)) {
    #test imputation time
    cat("Time test . . . \n")
    time.test <- MC.imputation(X = X.full,
                               time.test = TRUE)

    estim.time <- time.test * n.imp / (60 * 10) #10 imputed datasets used to test
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
  imputation.list$rX.imput <- imputation.list$rX.imput[names(y),,]
  #function to compute log(BFa0) for a given model with García-Donato's 2025 method
  lBF.method <- function (model) lBF.miss(model,
                                          imputation.list = imputation.list,
                                          BF.miss.aux = BF.miss.aux,
                                          n = n, nMC = n.imp)
  for (i in competing.models){
    modeli <- which(namesx %in% covar.list[[i]])

    lBFi0[i] <- lBF.method(model = modeli)
    lPriorModels[i] <- log(prior.models(i))
  }

  lPriorModels[pos.user.null.model] <- log(prior.models(pos.user.null.model))
  lBFi0[pos.user.null.model] <- 0
  C <- sum(exp(lBFi0 + lPriorModels))
  PostProbi <- exp(lBFi0 + lPriorModels - log(C))

  if (names(models)[pos.user.null.model] == "") {
    names(models)[pos.user.null.model] <-
      paste("(", as.formula(models[[1]])[[2]], " ~ 1)", sep="")
    names(lBFi0) <-
      paste(names(models), ".to.(",
            paste(as.formula(models[[1]])[[2]], " ~ 1)", sep=""), sep = "")
    names(PostProbi) <- c(names(models)[-pos.user.null.model],
                          paste(as.formula(models[[1]])[[2]], " ~ 1", sep=""))
    names(lPriorModels) <- names(PostProbi)
  } else {
    names(lBFi0) <-
      paste(names(models), ".to.", names(models)[pos.user.null.model], sep = "")
    names(PostProbi) <- names(models)
    names(lPriorModels) <- names(models)
  }

  result <- list()
  result$lBFi0 <- lBFi0
  result$PostProbi <- PostProbi
  result$models <- models
  result$nullmodel <- names(models)[pos.user.null.model]

  #arguments used for imputation
  result$imp.args <- list(initialimp.mice.method = initialimp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  result$prior.models <- prior.models #function used for model prior
  result$priorprobs <- exp(lPriorModels)

  class(result) <- "MissingBtest"

  return(result)
}
