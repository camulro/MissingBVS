#' Bayes factors and posterior probabilities via Bayesian Imputation Averaging for
#' linear regression models
#'
#' The model space is build from a list of linear regression models proposed to explain
#' a common response. It returns the Bayes factors and posterior probabilities computed
#' through Bayesian Imputation Averaging (BIA) for linear models in the presence of missing
#' data.
#'
#' Given a list of competing models, the model space is made up by them, assuming that the
#' intercept term is present in every model. The simplest one M0, can be specified (\code{null.model})
#' and must be nested in the rest. In order to implement BIA, \code{\link[MissingBVS]{missingBtest.lm}}
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
#' with the \code{BF.approx.method} argument.
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and are computed using \pck{BayesVarSel}. The choices currently available are:
#' -"Robust" is the default option and denotes the criteria-based prior of Bayarri,
#' Berger, Forte and Garcia-Donato (2012).
#' -"gZellner" corresponds to the prior in Zellner (1986) with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"ZellnerSiow" is the multivariate Cauchy prior by Zellner and Siow (1980, 1984).
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015).
#' -"IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato,
#' Moreno and Pericchi (2022).
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
#' @param maxit Number of iterations for \pck{mice}'s imputation. By default, it is 5.
#' @param parallelmice Logical to indicate whether or not to use parallelization on
#' \code{\link[mice]{mice}}'s imputation. By default, automatically performs it if the
#' number of imputations or competing variables given by \code{formula} are big enough.
#' @param n.core Number of cores for parallel imputation.
#' @param imp.datasets Array or list for imputed datasets if given by user. By default
#' it is set to NULL and imputation is performed following other imputation arguments.
#' @param imp.seed Seed for imputation.

#' @return \code{\link[MissingBVS]{MissingBtest.lm}} returns an object of type
#' \code{MissingBtest} with the following elements:
#' \item{lBFi0}{Bayes factors in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{If missings, list of the combined estimates for each model
#' in \code{models} fitted by \code{\link[stats]{lm}} over the \code{n.imp} imputed
#' datasets; or \code{lm} object when there are no missings}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information V}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{BF.approx.method}{Function used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
#' \item{logprior.models}{Function used to compute the log-prior over the model space
#' defined by covariates and/or factors}
#' \item{prior.models}{Two-dimensional vector with \code{prior.models} and
#' \code{prior.models.dummies} chosen. If there are no factors or \code{marginal.factors}
#' is set to FALSE, it saves the only argument used, \code{prior.models}}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' model space probabilities such as García-Donato and Paulo (2022)}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{call}{The \code{call} to the function}
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
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
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
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Default choices are: robust and Constant priors and 390 imputed datasets
#' #with mice's pmm method.
#' models.list = list(M0 = gr56092 ~ 1, M1 = gr56092 ~ lifee060,
#'   M2 = gr56092 ~ gdpsh60l, M3 = gr56092 ~ p60, M4 = gr56092 ~ lifee060 + p60,
#'   M5 = gr56092 ~ lifee060 + gdpsh60l, M6 = gr56092 ~ p60 + gdpsh60l,
#'   M7 = gr56092 ~ lifee060 + gdpsh60l + p60)
#'
#' dataS97.mtest <- missingBtest.lm(data = dataS97, models = models.list)
#'
#' #Show the results:
#' dataS97.mtest
#' }
#'

missingBtest.lm <- function (data,
                             models,
                             null.model = NULL,
                             BF.approx.method = "gprior",
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
                             imp.seed = runif(1,0,09011975)) {

  #N is the number of models:
  N <- length(models)

  env <- environment()

  #Check Btest given arguments
  Btestarg.list <- checkBtestarguments(models, null.model)
  list2env(Btestarg.list, envir = env)

  SSE <- numeric(N) #SSEs for each model
  Dim <- rep.int(0,N)
  mt <- list() #list of terms for each model

  covar.list <- list() #list that contains the names of the variables in each model
  compvars <- c() #name of original competing vars
  for (i in seq_len(N)) {
    f <- as.formula(models[[i]])
    compvars <- c(compvars, attr(terms(f), "term.labels"))
    temp <- lm(formula = f, data = data, y = TRUE, x = TRUE)

    SSE[i] <- crossprod(temp$residuals)

    Xi <- model.matrix.rankdef(model.frame(temp))
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.SSE <- sort(SSE, index.return = TRUE, decreasing = TRUE)
  #Which one acts as null model:
  nullmodel.pos <- ordered.SSE$ix[1]

  #Check null model
  if (!is.null(null.model) & nullmodel.pos != pos.user.null.model){
      stop("The given null model does not coincide with the one with the\n",
           "largest sum of squared error (and it should).\n")
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(models[[nullmodel.pos]])
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
  y <- framenull[obsnotNA, 1] #response variable without missings
  n <- length(y)
  SS0 <- SSE[nullmodel.pos]

  #Compute model prior
  lprior.models <- priormodels.btest(prior.models, N, Dim, priorprobs)

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.lm(BF.approx.method, prior.betas,
                                             n, p = max(Dim), p0, y, SS0)

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)
  #Imputation step
  if (!is.null(NAvars)) {
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

  #Evaluate lm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(namesxnotnull %in% covar.list[[j]])
    if (any(namesxnotnull[namesj] %in% NAvars)) {
      fit <- list()
      for (i in 1:n.imp) {
        if (mF) {#remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation.array[,c(1:p0, setdiff(namesj, indf) + p0),i]
        } else Xi <- imputation.array[,c(1:p0, namesj + p0),i]

        z <- lm.fit(x = Xi, y = y)
        z$terms <- mt[[j]]; class(z) <- "lm"; fit[[i]] <- z
      }
      modelspool[[j]] <- mice::pool(fit)
      modelspool[[j]]$call <- NULL #otherwise, Rstudio returns a warning trying to read modelspool[[j]]$call
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

  if (mF) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  if (!is.null(NAvars)) {
    #arguments used for imputation
    result$imp.info <- imp.info

    # save the imputed datasets for sensitivity analysis
    raw.imp.array <- serialize(imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

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

#' @keywords internal
posterior.btest <- function (competing.models, namesxnotnull, namesnull, covar.list,
                             relax.nest, lprior.models, lBF.comp, nullmodel.pos, models) {
  #posterior computation of model space defined by models list

  lBFi0 <- lPriorModels <- numeric(length(models))
  for (i in competing.models){
    modeli <- namesxnotnull %in% covar.list[[i]]

    #check whether the null is nested in the other ones
    if (!relax.nest & any(namesnull %notin% covar.list[[i]])) {
      stop("The simplest (null) model may not be nested in all the others.\n",
           "Please define explicitly the null model if it is the case.\n")
    }
    lBFi0[i] <- lBF.comp(modeli, i) #log(BF_a0)
    lPriorModels[i] <- lprior.models(i) #log-model prior
  }
  cat("\n")
  lPriorModels[nullmodel.pos] <- lprior.models(nullmodel.pos)
  lBFi0[nullmodel.pos] <- 0
  C <- sum(exp(lBFi0 + lPriorModels))
  PostProbi <- exp(lBFi0 + lPriorModels - log(C))

  names(lBFi0) <-
    paste(names(models), ".to.", names(models)[nullmodel.pos], sep = "")
  names(PostProbi) <- names(models)
  names(lPriorModels) <- names(models)

  return(list(lBFi0 = lBFi0, PostProbi = PostProbi, lPriorModels = lPriorModels))
}

#' @keywords internal
priormodels.btest <- function (prior.models, N, Dim, priorprobs) {
  #Check arguments and define the functions to compute prior model probabilities

  #Check model priors chosen and define the function to be used
  if (prior.models %notin% c("ScottBerger", "Constant", "User")) {
    stop("Only priors 'ScottBerger', 'Constant' and 'User' supported.\n")
  }
  switch (prior.models, #change the string for the corresponding function
          Constant = {lprior.models <- function (modeli) -log(N)},
          ScottBerger = {lprior.models <- function (modeli)
            -log(length(unique(Dim))) - log(sum(Dim == Dim[modeli]))},
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
            lprior.models <- function(modeli) log(priorprobs[modeli])}
  )
  return(lprior.models)
}

#' @keywords internal
BFcomp.btest <- function (lprior.models, prior.models.dummies, Dim, mF,
                          positionsfac,  namesxnotnull, NAvars, lBF.method,
                          X0, X.full, BF.approx.method, covar.list) {
  #Check arguments and define the functions to compute Bayes factors

  if (mF) {
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

    #Define BF to compute the marginal for the dummies
    lBF.comp <- function (modeli, i) {
      pm <- positionsfac %*% modeli
      tau <- pm > 0; ltau <- pm[tau]
      m2 <- sum(tau) #number of factors active
      if (m2 > 0) {
        mats <- lapply(ltau, build_ind)
        mat.ind <- Reduce(function(x, y) merge(x, y, by = NULL), mats)

        cn <- which(colSums(positionsfac[which(tau), , drop = FALSE]) > 0)
        colnames(mat.ind) <- cn

        lBF <- lpriorM <- numeric(nrow(mat.ind))
        for (j in 1:nrow(mat.ind)) {
          dj <- as.integer(mat.ind[j,]); djsum <- positionsfac[which(tau), cn] %*% dj

          current.model <- as.integer(modeli)
          current.model[cn] <- dj

          #check if there are NAs in the model considered to save computation time
          if (any(namesxnotnull[which(modeli)] %in% NAvars)) {
            lBF[j] <- lBF.method(model = which(current.model == 1)) #log(BF_a0)
            lpriorM[j] <- lprior.models.dummies(djsum, ltau) #log(Pr(M_delta))
          } else { #if there are no missings, compute the BF by the method selected
            X.i <- cbind(X0, X.full[, which(current.model == 1)])
            lBF[j] <- BF.approx.method(k = sum(current.model == 1), X = X.i) #log(BF_a0)
            lpriorM[j] <- lprior.models.dummies(djsum, ltau) #log(Pr(M_delta))
          }
        }
        lBFi0 <- log(sum(exp(lBF + lpriorM)))
        # lPriorModels[i] <- lprior.models(i)
      } else {
        #check if there are NAs in the model considered to save computation time
        if (any(covar.list[[i]] %in% NAvars)) {
          lBFi0 <- lBF.method(model = which(modeli))
        } else { #if there are no missings, compute the BF by the method selected
          X.i <- cbind(X0, X.full[,which(modeli)])
          lBFi0 <- BF.approx.method(k = Dim[i], X = X.i)
        }
      }
      return(lBFi0)
    }
  } else { #Define the standard BF
    lBF.comp <- function (modeli, i) {
      #check if there are NAs in the model considered to save computation time
      if (any(covar.list[[i]] %in% NAvars)) {
        lBFi0 <- lBF.method(model = which(modeli))
      } else { #if there are no missings, compute the BF by the method selected
        X.i <- cbind(X0, X.full[,which(modeli)])
        lBFi0 <- BF.approx.method(k = Dim[i], X = X.i)
      }
      return(lBFi0)
    }
  }
  return(lBF.comp)
}

#' @keywords internal
checkBtestarguments <- function (models, null.model) {
  #check arguments
  if (!is.list(models)) stop("Argument models should be a list.\n")

  #If competing models come wihtout a name, give one by default:
  if (is.null(names(models))){
    if (!is.null(null.model)) stop("Please provide a name for the competing models.\n",
                                   "The null model must be in that list.\n")
    names(models) <- paste("model", seq_len(N), sep="")
  }

  #Check if the given null model is one of the competing models:
  if (!is.null(null.model)){
    relax.nest = TRUE
    pos.user.null.model <- which(null.model == names(models))
    if (length(pos.user.null.model) == 0) {
      stop("The null model provided is not in the list of competing models.\n")
    }

    return(list(models = models, relax.nest = relax.nest,
                pos.user.null.model = pos.user.null.model))
  } else relax.nest = FALSE

  return(list(models = models, relax.nest = relax.nest))
}

#' @keywords internal
build_ind <- function(k) {
  ind <- t(sapply(2:2^k - 1,
                  FUN = function(j2) BayesVarSel:::integer.base.b_C(j2,k)))

  rs <- rowSums(ind)
  ind[!(rs == k | (rs == (k - 1) & ind[, k])), , drop = FALSE]
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
#' @examples
#' \dontrun{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' #Default choices are: robust and Constant priors and 390 imputed datasets
#' #with mice's pmm method.
#' models.list = list(M0 = gr56092 ~ 1, M1 = gr56092 ~ lifee060,
#'   M2 = gr56092 ~ gdpsh60l, M3 = gr56092 ~ p60, M4 = gr56092 ~ lifee060 + p60,
#'   M5 = gr56092 ~ lifee060 + gdpsh60l, M6 = gr56092 ~ p60 + gdpsh60l,
#'   M7 = gr56092 ~ lifee060 + gdpsh60l + p60)
#' lifee060 + gdpsh60l + p60
#' dataS97.mtest <- missingBtest.lm(data = dataS97, models = models.list)
#'
#' #Show the results:
#' dataS97.mtest
#' }
#'
print.MissingBtest <- function(mbtest.object,...){
  if (!inherits(mbtest.object, "MissingBtest")){
    warning("An object of class MissingBtest is needed.\n")
  }

  cat("-------\n")
  cat("Competing models:\n")
  print(mbtest.object$models)
  cat("-------\n")
  cat("log Bayes factors (expressed in relation to ",
      mbtest.object$nullmodel,")\n", sep="")
  print(round(mbtest.object$lBFi0, 3))
  cat("-------\n")
  cat("Posterior probabilities:\n")
  print(round(mbtest.object$PostProbi,3))
  cat("\n\n")
}
