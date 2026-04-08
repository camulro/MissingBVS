#' Bayes factors with missing data and posterior probabilities for linear
#' regression models
#'
#' It computes the Bayes factors and posterior probabilities in the presence
#' of missing data of a list of linear regression models proposed to explain a
#' common response. See \code{\link[MissingBVS]{MissingBvs.lm}} for more details.
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
#' @param prior.models Prior distribution over the model space (literally specified).
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
#' @param imp.mice.method Method for mice's imputation.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param imp.seed Seed for imputation.

#' @return \code{\link[MissingBVS]{MissingBtest.lm}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factor in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{List of objects of class \code{\link[mice]{mipo}} that
#' combine the estimates for each model on \code{models} fitted by
#' \code{\link[stats]{lm}} over the \code{n.imp} (if \code{> 1}) imputed
#' datasets, see \code{\link[mice]{pool}} for details; and \code{lm} object
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
#' @seealso Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
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
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6: 461–464.
#'
#' Held, L., Gravestock, I. and Sabanés Bové, D.
#' (2015)<DOI:10.1080/01621459.2014.993077> Objective Bayesian model selection
#' for generalized linear models using test-based Bayes factors. Journal of the
#' American Statistical Association, 110, 1157–1168.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @keywords package
#'
#' @examples #To be completed
#'

missingBtest.lm <- function (data,
                             models,
                             null.model = NULL,
                             BF.approx.method = "gprior",
                             prior.betas = "Robust", #if BF.approx.method = "gprior"
                             prior.models = "Constant",
                             prior.models.dummies = "ScottBerger",
                             priorprobs = NULL, #needed if prior.models = "User"
                             parallelmice = NULL,
                             n.core = NULL,
                             imp.time.test = TRUE,
                             imp.mice.method = "pmm", #mice's default
                             n.imp = 039E1,
                             imp.seed = runif(1,0,09011975)) { #seed for the imputation

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

  SSE <- numeric(N) #SSEs for each model
  Dim <- rep(0L,N)
  lBFi0 <- numeric(N)
  lPriorModels <- numeric(N)
  PostProbi <- numeric(N)
  mt <- list() #list of terms for each model

  #list that contains the names of the variables or factors in each model
  covar.list <- list()
  for (i in seq_len(N)) {
    temp <- lm(formula = as.formula(models[[i]]),
               data = data,
               y = TRUE, x = TRUE)

    SSE[i] <- crossprod(temp$residuals)

    Xi <- model.matrix.rankdef(model.frame(temp))
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.SSE <- sort(SSE, index.return = TRUE, decreasing = TRUE)
  #Which acts as null model:
  nullmodel.pos <- ordered.SSE$ix[1]

  if (!is.null(null.model)){
    if (nullmodel.pos != pos.user.null.model){
      stop(paste0("The given null model does not coincide with the one with the\n",
                  "largest sum of squared error (and it should).\n"))
    }
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(models[[nullmodel.pos]])
  competing.models <- seq_len(N)[-nullmodel.pos]

  #Response and fixed vars for imputation
  framenull <- model.frame(null.model, data, na.action = NULL)
  #Rank deficient fixed model matrix to remove vars from the regressors one
  X0rdf <- model.matrix.rankdef(framenull)

  #Missing model matrix of fixed vars
  X0 <- model.matrix(framenull, data)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of fixed covariates or dummies of factors
  Dim <- Dim - p0 #model dimension (without fixed vars)

  #Full design matrix for imputation
  formula <- as.formula(paste0(null.model[[2]], "~ ."))
  framefull <- model.frame(formula, data, na.action = NULL)
  X.toimp <- data[, attr(terms(framefull), "term.labels")] #design matrix with missing data with fixed vars

  #Rank deficient full model matrix data with missings
  X.full <- model.matrix.rankdef(framefull)
  namesx <- dimnames(X.full)[[2]]
  #Only the non-fixed vars
  namesxnotnull <- namesx[namesx %notin% dimnames(X0rdf)[[2]]]
  X.full <- X.full[, namesxnotnull]
  p <- length(namesxnotnull) #Number of covariates and levels of factors

  #the order for the posterior model distribution computation step
  ordvars <- c(namesnull, namesxnotnull) #X0, X.full

  #covariates and/or factors to select from
  depvars <- setdiff(attr(terms(framefull), "term.labels"),
                     attr(terms(framenull), "term.labels"))

  #Factors: positions has number of rows equal to the number of regressors
  #(factors or numeric covariates) and p columns.
  #A 1 in a row denotes the position in X of a regressor (several positions for
  #the dummies of a factor).
  positions <- t(sapply(depvars, function(var) {
    if(is.factor(data[[var]])) {
      levs <- levels(data[[var]])
      ind <- which(namesxnotnull %in% paste0(var,levs)) #1 if the namelevel matches
    } else ind <- which(namesxnotnull == var) #1 if the name matches

    posi <- rep(0,p); posi[ind] <- 1
    posi
  }))
  colnames(positions) <- namesxnotnull

  tmp <- colSums(positions %*% t(positions))
  positionsx <- tmp == 1 #vector of length p with TRUE if numeric variable

  L <- sum(!positionsx) #Number of factors to select from
  l <- tmp[tmp > 1] #Number of levels for each factor
  q <- p - sum(l) + L #Number of factors and covariates to select from
  #q = p if there are no factors

  if (L > 0) {
    #matrix of dim (Lxp) with 1 if dummy variable of the row factor
    positionsfac <- matrix(positions[!positionsx,], ncol = p, nrow = L)
    rownames(positionsfac) <- depvars[!positionsx]
    colnames(positionsfac) <- namesxnotnull
    #vector of length L with the position of the last dummy for each factor to check for repeated models
    indf <- apply(positionsfac, MARGIN = 1, FUN = function(x) tail(which(x == 1), n = 1))
  } else positionsfac <- indf <- 0

  #The response variable
  obsnotNA <- rownames(X0)
  y <- framenull[obsnotNA, 1] #response variable without missings
  n <- length(y)
  SS0 <- SSE[nullmodel.pos]

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)

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

  #Check model priors for dummies chosen and define the function to be used
  if (prior.models.dummies %notin% c("ScottBerger", "Constant")) {
    stop("Only priors 'ScottBerger' and 'Constant' supported.\n")
  }
  switch (prior.models.dummies, #change the string for the corresponding function
          Constant = {lprior.models.dummies <-
            function (deltai, ltau) {-sum(log(2^(ltau) - 1 - ltau))}},
          ScottBerger = {lprior.models.dummies <-
            function (deltai, ltau) {-sum(mylchoose(ltau, deltai)) - sum(log(ltau - 1))}}
  )

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.lm(BF.approx.method, prior.betas,
                                             n, p = max(Dim), p0, y, SS0)

  if (!is.null(NAvars)) {
    #Imputation of missing data
    if (is.null(parallelmice)) {
      if (n.imp > 120 | (n*q > 50000 & n.imp > 5)) {
        parallelmice <- TRUE #faster
      } else parallelmice <- FALSE
    }

    if (imp.time.test & (n*q > 10000 | n.imp > 039E1)) {
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

    #remove observations with missings on the response or fixed vars,
    #select the vars in the order X0, X.full and remove oversaturated for X0 if factors
    imputation.array <- imputation.array[obsnotNA, ordvars,]
    if (n.imp > 1) {
      #function to compute log(BFa0) for a given model as an average of BF computed
      #by BF.approx.method over the imputed datasets
      lBF.method <- function (model) lBF.approx(model,
                                                imputation.array = imputation.array,
                                                BF.approx.method = BF.approx.method,
                                                p0 = p0, n.imp = n.imp)
    } else lBF.method <- function (model) BF.approx.method(k = length(model),
                                                           X = imputation.array[,c(1:p0, model+p0)])
  }

  for (i in competing.models){
    modeli <- namesxnotnull %in% covar.list[[i]]

    #check whether the null is nested in the other ones
    if (!relax.nest & any(namesnull %notin% covar.list[[i]])) {
      stop(paste0("The simplest (null) model may not be nested in all the others.\n",
                  "Please define explicitly the null model if it is the case.\n"))
    }

    tau <- (positionsfac %*% modeli) > 0; ltau <- (positionsfac %*% modeli)[tau]
    m2 <- sum(tau) #number of factors active
    if (m2 > 0) {
      colsi <- which(colSums(matrix(positionsfac[which(tau),], nrow = m2)) > 0)
      ind <- t(sapply(2:(2^ltau[1])-1, FUN =
                        function(j2) BayesVarSel:::integer.base.b_C(j2, ltau[1])))
      rep <- which((rowSums(ind) == ltau[1]) |
                     ((rowSums(ind) == (ltau[1] - 1)) &  ind[,ltau[1]]))
      mat.ind <- matrix(ind[-rep,], ncol = ltau[1])
      if (m2 > 1) {
        for(j in 2:m2){
          ind <- t(sapply(2:(2^ltau[j])-1, FUN =
                            function(j2) BayesVarSel:::integer.base.b_C(j2, ltau[j])))
          rep <- which((rowSums(ind) == ltau[j]) |
                         ((rowSums(ind) == (ltau[j] - 1)) &  ind[,ltau[j]]))
          ind <- matrix(ind[-rep,], ncol = ltau[j])
          mat.ind <- merge(mat.ind, ind, by = NULL)
        }
      }
      colnames(mat.ind) <- colsi

      lBF <- lpriorM <- numeric(nrow(mat.ind))
      for (j in 1:nrow(mat.ind)) {
        deltaj <- mat.ind[j,]
        deltasumj <- positionsfac[which(tau), colsi] %*% as.integer(deltaj)

        current.model <- as.integer(modeli)
        current.model[colsi] <- as.integer(deltaj)

        #check if there are NAs in the model considered to save computation time
        if (any(namesxnotnull[which(modeli)] %in% NAvars)) {
          lBF[j] <- lBF.method(model = which(current.model == 1)) #log(BF_a0)
          lpriorM[j] <- lprior.models.dummies(deltasumj, ltau) #log(Pr(M_delta))
        } else { #if there are no missings, compute the BF by the method selected
          X.i <- cbind(X0, X.full[, which(current.model == 1)])
          lBF[j] <- BF.approx.method(k = sum(current.model == 1), X = X.i) #log(BF_a0)
          lpriorM[j] <- lprior.models.dummies(deltasumj, ltau) #log(Pr(M_delta))
        }
      }
      lBFi0[i] <- log(sum(exp(lBF + lpriorM)))
      lPriorModels[i] <- log(prior.models(i))
    } else {
      #check if there are NAs in the model considered to save computation time
      if (any(covar.list[[i]] %in% NAvars)) {
        lBFi0[i] <- lBF.method(model = which(modeli))
      } else { #if there are no missings, compute the BF by the method selected
        X.i <- cbind(X0, X.full[,which(modeli)])
        lBFi0[i] <- BF.approx.method(k = Dim[i], X = X.i)
      }
      lPriorModels[i] <- log(prior.models(i))
    }
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

  #Evaluate lm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesxnotnull %in% covar.list[[j]])
    if (any(namesxnotnull[namesj] %in% NAvars)) {
      if (n.imp > 1) {
        fit <- list()
        for (i in 1:n.imp) {
          #remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation.array[,c(1:p0, namesj[namesj %notin% indf] + p0),i]
          z <- lm.fit(x = Xi, y = y)
          z$terms <- mt[[j]]
          class(z) <- "lm"

          fit[[i]] <- z
        }
        modelspool[[j]] <- mice::pool(fit)
      } else {
        #compute lm.fit for the unique imputation
        Xi <- imputation.array[,c(1:p0, namesj[namesj %notin% indf] + p0)]
        modelspool[[j]] <- lm.fit(x = Xi, y = y)
      }
    } else modelspool[[j]] <- lm(models[[j]], data)
  }
  modelspool[[nullmodel.pos]] <- lm(null.model, data)
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
    result$imp.args <- list(parallelmice = parallelmice,
                            imp.mice.method = imp.mice.method,
                            n.imp = n.imp, imp.seed = imp.seed)

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

#' Print an object of class \code{MissingBtest}
#'
#' Print an object of class \code{MissingBtest}
#' @export
#' @param mbtest.object Object of class MissingBtest.
#' @param ... Additional parameters to be passed.
#'
#' @author Gonzalo Garcia-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso See \code{\link[MissingBVS]{MissingBtest.lm}},
#' \code{\link[MissingBVS]{MissingBtest.glm}} and
#' \code{\link[MissingBVS]{MissingBtestGD25}} for creating objects of the class
#' \code{MissingBtest}.
#'
#' @examples #To be completed
#'
print.MissingBtest <- function(mbtest.object,...){
  if (!inherits(mbtest.object, "MissingBtest")){
    warning("An object of class MissingBtest is needed.\n")
  }

  cat("-------\n")
  cat("Competing models:\n")
  print(mbtest.object$models)
  cat("-------\n")
  cat(paste0("log Bayes factors (expressed in relation to ",
            mbtest.object$nullmodel,")\n", sep=""))
  print(round(mbtest.object$lBFi0, 3))
  cat("-------\n")
  cat("Posterior probabilities:\n")
  print(round(mbtest.object$PostProbi,3))
  cat("\n\n")
}
