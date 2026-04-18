#' Bayesian Variable Selection with Missing data for generalized linear models
#' using Gibbs sampling
#'
#' Approximate computation of summaries of the posterior distribution using a
#' Gibbs sampling algorithm to explore the model space and frequency of
#' "visits" to construct the estimates. It was originally proposed by George
#' and McCulloch (1997) and further studied by Garcia-Donato and
#' Martinez-Beneito (2013). They shown that the sampling strategy in
#' combination with estimates based on frequency of visits provides very
#' reliable results.
#'
#' \code{\link[MissingBVS]{MissingGibbsBvs.glm}} is a heuristic approximation of
#' \code{\link[MissingBVS]{MissingBvs.glm}}. See the later for common details.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data Data frame containing the data.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. If it is one of the implemented families in
#' \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}; a faster version
#' using BAS logmarginal computation is performed.
#' @param null.model Formula defining which is the simplest (null) model, which.
#' should be nested in the full one. By default, it is defined to be the one
#' with just the intercept.
#' @param BF.approx.method Method used to approximate Bayes factors with missing
#' data (to be literally specified). Possible choices include "BIC", "TBF" and
#' "gprior" (see details).
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices are the ones implemented
#' in \pkg{BayesVarSel}: "Robust", "Liangetal", "gZellner", "ZellnerSiow",
#' "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified). Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the model space of the
#' factor levels (to be literally specified). Possible choices are "Constant" and
#' "ScottBerger" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed variables)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details.)
#' @param init.model The model at which the simulation process starts. Options
#' include "Null" (the model only with the variables specified in
#' \code{null.model}), "Full" (the model defined by \code{formula}), "Random" (a
#' randomly selected model) and a vector with p (the number of factors and/or
#' covariates to select from) zeros and ones defining a model.
#' @param n.iter The total number of iterations performed after the burn in
#' process.
#' @param n.burnin Length of burn in, i.e. number of iterations to discard at
#' the beginning.
#' @param n.thin Positive integer defining the thinning rate. Default is 1.
#' Set 'n.thin' > 1 to save memory and computation time if 'n.iter' is large.
#' A large \code{n.thin} can reduce the accuracy because it, along with \code{n.iter},
#' sets the number of simulations used to construct the estimates.
#' @param parallelmice Logital to indicate whether or not to use parallel
#' \code{\link[mice]{mice}} imputation. If \code{NULL}, automatically performs
#' the parallel mice imputation if the number of imputations or the dataset are
#' big enough.
#' @param n.core See \code{\link[mice]{futuremice}} for details.
#' @param imp.time.test Logical to indicate whether to check or not time of performance
#' of the imputation process with \code{n.imp = 30} if the number of variables or
#' the number of imputed datasets are large enough (\code{p>10} or \code{n.imp>390}).
#' @param imp.mice.method Method for mice's imputation.
#' @param imp.predict.mat \code{matrix} with p1 rows and p2 columns, where p1 is
#' the number of independent variables given by \code{formula} with \code{NAs}
#' and p2 is the number of variables used to impute. Each entry equals 1 if the
#' column variable is used as a predictor for the corresponding row variable in the
#' imputation step. By default, the \code{\link[mice]{quickpred}} function
#' is used for the variables with NAs in formula and 0s for the rest of them.
#' @param n.imp Number of imputed data sets used for Bayes factor computation.
#' @param Gibbs.seed Seed for the Gibbs sampler algorithm.
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
#'
#' @return \code{missingGibbsBVS.lm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{glmfull}{If missings on the competing variables, object of class
#' \code{\link[mice]{mipo}} that combines the estimates for the model defined by
#' \code{formula} fitted by \code{\link[stats]{glm}} over
#' the \code{n.imp} (if \code{> 1}) imputed datasets; see \code{\link[mice]{pool}}
#' for details. Otherwise, it is the \code{glm} object for the \code{formula} model}
#' \item{glmnull}{The \code{glm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{glm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model using
#' \code{inclprobRB}}
#' \item{positions}{\code{matrix} with L rows and p plus the number of dummies
#' resulting from factors columns - L, where L is the number of factors, with 1
#' if the column dummy corresponds to the row factor and 0 otherwise}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate}
#' \item{modelsprob}{\code{data.frame} which summaries the \code{n.keep}
#' most probable a posteriori models and their associated Bayes factor in
#' logaritmic scale}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables}
#' \item{inclprobRB}{Rao-Blackwellized inclusion probabilities}
#' \item{postprobdim}{Estimated posterior probabilities over the true model dimension}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the
#' model space)}
#' \item{imp.args}{List of arguments used for the imputation step}
#' \item{BF.approx.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{logprior.models}{Function used to compute the log-prior over the model space}
#' \item{logprior.models.dumm}{Function used to compute the log-prior over the
#' model space of the factor levels}
#' \item{method}{\code{Gibbs}}
#'
#' @author  Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingBvs.glm}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' George E. and McCulloch R. (1997) Approaches for Bayesian variable
#' selection. Statistica Sinica, 7, 339:372.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6(2): 461–464.
#'
#' Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
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
#' @keywords package
#'
#' @examples
#' \dontrun{
#' #Indian Prime Diabetes Data from VIM's package
#'
#' #Default choices are: BIC approximation and ScottBerger prior, 10000 iterations
#' #with 500 of burn in period. Here, 100 imputations with mice's pmm method.
#' diabetes.mGBVS <- missingGibbsBVS.glm(formula = Outcome ~  .,
#'   data = VIM::diabetes, family = binomial(), n.imp = 100)
#'
#' #Show the results:
#' diabetes.mGBVS
#'
#' #Summ up the results:
#' summary(diabetes.mGBVS)
#'
#' #A plot with the posterior inclusion probabilities for each competing variable
#' #and the dimension probability of the true model:
#' plot(diabetes.mGBVS)
#' }
#'
missingGibbsBVS.glm <- function (formula,
                                 data,
                                 family = binomial(link = "logit"),
                                 null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                                 BF.approx.method = "BIC",
                                 prior.betas = "Robust",
                                 prior.models = "ScottBerger",
                                 prior.models.dummies = "ScottBerger",
                                 priorprobs = NULL,
                                 init.model = "Full",
                                 n.iter = 10000,
                                 n.burnin = 500,
                                 n.thin = 1,
                                 parallelmice = NULL,
                                 n.core = NULL,
                                 imp.time.test = TRUE,
                                 imp.mice.method = "pmm",
                                 imp.predict.mat = NULL,
                                 n.imp = 039E1,
                                 Gibbs.seed = runif(1,0,26061970),
                                 imp.seed = runif(1,0,09011975),
                                 weights = rep(1, dim(data)[1]),
                                 offset = rep(0, dim(data)[1]),
                                 control = glm.control(),
                                 laplace = 0L) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #The response in the null model and in the full model must coincide
  if (formula[[2]] != null.model[[2]]){
    stop("The response in the full and null model does not coincide.\n")
  }

  #check whether the family chosen is among the options provided by BAS
  inBAS <- checkforfamily(family, BF.approx.method)

  #select environment to get glm arguments
  environment(formula) <- environment()
  environment(null.model) <- environment()

  #for the C code
  weights <- as.numeric(weights)
  offset <- as.numeric(offset)
  laplace <- as.integer(laplace)

  #Evaluate the null model:
  glmnull <- glm(formula = null.model,
                 data,
                 y = TRUE, x = TRUE,
                 family = family,
                 weights = weights,
                 offset = offset,
                 control = control)

  #Response and fixed vars for imputation
  framenull <- model.frame(null.model, data, na.action = NULL)
  q0 <- dim(framenull)[2] # number of fixed covariates and/or factors,intercept always

  #Rank deficient fixed model matrix to remove vars from the regressors one
  X0rdf <- model.matrix.rankdef(framenull)

  #Missing model matrix of fixed vars
  X0 <- model.matrix(framenull, data)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of fixed covariates or dummies of factors

  #Full design matrix for imputation
  framefull <- model.frame(formula, data, na.action = NULL)
  #Rank deficient full model matrix data with missings
  X.full <- model.matrix.rankdef(framefull)
  namesx <- dimnames(X.full)[[2]]
  #Only the non-fixed vars
  namesxnotnull <- setdiff(namesx, dimnames(X0rdf)[[2]])
  X.full <- X.full[, namesxnotnull]
  p <- length(namesxnotnull) #Number of covariates and levels of factors to select from

  #Is there any variable to select from?
  if (p == 0) {
    stop(paste0("The number of fixed variables is equal to the number of\n",
                "regressors in the full model. No model selection can be done.\n"))
  }

  if (p <= 20) {
    warning(paste0("The number of variables is small enough to visit every model.\n",
                   "Consider using MissingBvs.lm.\n"), immediate. = TRUE)
  }

  #check if null model is contained in the full one:
  for (i in 1:p0){
    if (namesnull[i] %notin% namesx) {
      stop(paste0("Error in var: ", namesnull[i],"; null model not nested in full model.\n"))
    }
  }

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

    posi <- rep(0,p); posi[ind] <- 1; posi
  }))
  colnames(positions) <- namesxnotnull

  tmp <- colSums(positions %*% t(positions))
  positionsx <- tmp == 1 #vector of length p with TRUE if numeric variable

  L <- sum(!positionsx) #Number of factors to select from
  if (L > 0) {
    #matrix of dim (Lxp) with 1 if dummy variable of the row factor
    positionsfac <- matrix(positions[!positionsx,], ncol = p, nrow = L)
    rownames(positionsfac) <- depvars[!positionsx]
    colnames(positionsfac) <- namesxnotnull

    l <- tmp[tmp > 1] #Number of levels for each factor

    #vector of length L with the position of the last dummy for each factor to check for repeated models
    indf <- apply(positionsfac, MARGIN = 1, FUN = function(x) tail(which(x == 1), n = 1))
    checklast <- ifelse(L > 1, function(M) diag(M[, indf]),  function(M) M[, indf])
  } else positionsfac <- l <- 0

  q <- p - sum(l) + L #Number of factors and covariates to select from
  #q = p if there are no factors

  #The response variable
  y <- glmnull$y; obsnotNA <- names(y) #without missings
  y <- as.numeric(y)
  n <- length(y) #observations without missings on the response
  devnull <- glmnull$deviance #deviance of the null model

  X.full <- X.full[obsnotNA,] #remove NA obs from null model

  #check for missings and define competing variables with NAs
  NAvars <- checkformissings(y = framenull[,1], framenull[,-1], X.full)

  if (!is.null(NAvars)) {
    #Impute just competing variables with NAs
    fulldataframe <- model.frame(paste0(formula[[2]], "~."), data, na.action = NULL)
    X.toimp <- fulldataframe[,-1] #full observed design matrix
    quickpredict.mat <- mice::quickpred(X.toimp)
    if (!is.null(imp.predict.mat)) { #if given by user
      #check predict imputation matrix
      imp.vars <- rownames(imp.predict.mat)
      if (any(NAvars %notin% imp.vars)) {
        stop(paste0("Imputation prediction matrix rows given do not contain all the variables ",
                    "given by formula with NAs.", "Make sure to include them all.\n"))
      }
      quickpredict.mat[imp.vars, colnames(imp.predict.mat)] <- imp.predict.mat
      quickpredict.mat[imp.vars, colnames(quickpredict.mat) %notin% colnames(imp.predict.mat)] <- 0
    } else imp.vars <- NAvars

    #Do not impute variables with missings that are not in imp.vars (also in NAvars)
    not.imp.vars <- setdiff(attr(terms(fulldataframe), "term.labels"), imp.vars)
    quickpredict.mat[not.imp.vars, ] <- 0
  }

  #Check the initial model:
  if (is.character(init.model) == TRUE) {
    im <- substr(tolower(init.model), 1, 1)
    if (im %notin% c("n", "f", "r")) stop("Initial model not valid.\n")
    if (im == "n") init.model <- rep(0, p) #null model
    if (im == "f") init.model <- rep(1, p) #full model
    if (im == "r") init.model <- rbinom(n = p,
                                        size = 1,
                                        prob = .5)
  } else {
    init.model <- as.numeric(init.model > 0)
    if (length(init.model) != p) stop("Initial model with incorrect length.\n")
  }

  #Check model priors chosen and define the function to be used
  lprior.models <- checkforprior.models(prior.models, priorprobs, q)

  if (L > 0) {
    lprior.models.dummies <- checkforprior.models.dummies(prior.models.dummies, l)
  } else lprior.models.dummies <- function(delta, tau) 0

  #Check approx method and priors chosen and define the function to be used
  BF.approx.method <- checkforprior.betas.glm(BF.approx.method, prior.betas, inBAS,
                                              n, p, p0, y, null.model,
                                              data, family, devnull,
                                              weights, offset, control, laplace)

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
                                   imp.predict.mat = quickpredict.mat,
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
                                         imp.predict.mat = quickpredict.mat,
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

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", q + q0, "covariates and/or factors.\n")
  if (q0 == 1) {
    cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
               q, ".\n"))
  } else  cat(paste0("From those ", q0, " are fixed and we should select from the remaining ",
                     q, ".\n"))

  cat("  Numerical covariates:", depvars[positionsx], "\n")
  if (L > 0) cat(" Factors:", depvars[!positionsx], "\n")

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these,", n.iter + n.burnin, "are sampled with replacement.\n")
  cat("Then,", floor(n.iter / n.thin), "are kept and used to construct the summaries.\n")

  #George and McCulloch's Gibbs exploration
  gibbs.list <- GM97.Gibbs(y, X0, X.full, p, namesxnotnull, NAvars,
                           lprior.models, lprior.models.dummies, lBF.method, BF.approx.method,
                           positions, positionsfac, l, L, checklast,
                           init.model, n.iter, n.burnin, n.thin, Gibbs.seed)

  cf.models.lPM <- gibbs.list$cf.models.lPM

  inclprob <- colMeans(cf.models.lPM[,-(q+1)]) #inclusion probabilities except for fixed variables

  cf.models.PM <- cf.models.lPM[,seq_len(q)]
  cf.models.PM <- cbind(cf.models.PM, exp(cf.models.lPM[, q + 1]))
  colnames(cf.models.PM)[q+1] <- "BF.PM"

  #Estimation of the normalizing constant:
  K <- round(dim(cf.models.PM)[1]/2)
  Aset <- sample(x = 1:dim(cf.models.PM)[1], size = K, replace = F)
  Bset <- (1:dim(cf.models.PM)[1])[-Aset]
  #Bayes factors multiplied by prior probs of the models in A
  BF.PMAset <- cf.models.PM[Aset, "BF.PM"]
  #Remove duplicates
  BF.PMAset <- unique(BF.PMAset)
  gAset <- sum(BF.PMAset)
  #How many of the models in Bset are in A?
  sumIA <- sum(cf.models.PM[Bset,"BF.PM"] %in% cf.models.PM[Aset,"BF.PM"])
  #Normalizing constant estimation
  C <- gAset*K / sumIA

  #compute estimated posterior probabilities
  cf.models.PM[, q + 1] <- exp(cf.models.lPM[, q + 1] - log(C))
  colnames(cf.models.PM)[q+1] <- "Post"

  probdim <- rep(0, q + 1)
  cf.models.lBF <- cf.models.lPM
  #compute posterior probability of the dimension of the true model and
  #save logBF for each model
  for (i in seq_len(floor(n.iter / n.thin))) {
      probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] <-
        probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] + cf.models.PM[i, q + 1]

    gamma.tau <- cf.models.lPM[i, seq_len(q)] > 0
    deltasum <- cf.models.lPM[i, !positionsx]
    tau <- deltasum > 0
    cf.models.lBF[i, q + 1] <- cf.models.lPM[i, q + 1] -
      lprior.models(gamma.tau) - lprior.models.dummies(deltasum, tau) # lBF
  }
  colnames(cf.models.lBF)[q+1] <- "logBF"

  #HPM
  nPmax <- which.max(cf.models.lBF[, q+1])
  hpm <- cf.models.lBF[nPmax, ]

  #MPM
  mpm <- rep(0,q)
  mpm[which(gibbs.list$inclprobRB[n.iter, ] >= 0.5)] <- 1

  if (!is.null(NAvars)) {
    if (n.imp > 1) {
      #remove last dummy for each factor, first p0 vars are the fixed ones
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0),]
      #Evaluate glm of full model with missings using Rubin's rule
      fit <- list()
      mt <- attr(framefull, "terms")
      for (i in 1:n.imp) {
        z <- glm.fit(x = imputation.array[,,i], y = y, family = family,
                     weights = weights, offset = offset, control = control)
        z$terms <- mt
        class(z) <- "glm"

        fit[[i]] <- z
      }
      glmfull <- mice::pool(fit)
    } else {
      #compute glm.fit for the unique imputation
      if (L > 0) imputation.array <- imputation.array[,-c(indf + p0)]
      glmfull <- glm.fit(x = imputation.array, y = y, family = family,
                         weights = weights, offset = offset, control = control)
    }
  } else glmfull <- glm(formula,
                        data,
                        family = family,
                        weights = weights,
                        offset = offset,
                        control = control)

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the programm to finish
  result$glmfull <- glmfull # If missings, object of class mipo combining the
  # estimates for the n.imp imputed datasets for the fitted full model.
  # Otherwise, glmfull is the glm object for the full model
  result$glmnull <- glmnull # The glm object for the null model (without NAs)

  result$variables <- depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- q #number of competing vars
  result$k <- q0 #number of fixed vars
  result$HPMbin <- hpm #The binary code for the HPM model and its BF.PM
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- depvars

  if (L > 0) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  result$modelsprob <- cf.models.lBF

  #The binary code for all the visited models (after n.thin is applied) and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  result$inclprobRB <- gibbs.list$inclprobRB #Rao-Blackwellized inclusion probability

  result$postprobdim <- probdim/sum(probdim) #vector with the estimated posterior dimension probability
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model

  result$call <- match.call()
  result$C <- C #estimated normilizing constant

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, q + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, q))) #prior inclusion prob for dimension 0
    for (i in seq_len(q)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, q - i))) + lchoose(q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:q + q0 #dimension prior probability

  #Estimation of posterior probabilities based on C
  result$postprobs <- cf.models.PM[,"Post"]

  if (!is.null(NAvars)) {
    #arguments used for imputation
    result$imp.args <- list(parallelmice = parallelmice,
                            imp.mice.method = imp.mice.method,
                            imp.predict.mat = imp.predict.mat,
                            n.imp = n.imp, imp.seed = imp.seed)

    #save the imputed datasets for sensitivity analysis
    # raw.imp.array <- serialize(imputation.array, NULL)
    # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lprior.models #function used for model prior
  if (L > 0) result$logprior.models.dumm <- lprior.models.dummies

  result$method <- "Gibbs"
  class(result) <- "MissingBvs"

  return(result)
}
