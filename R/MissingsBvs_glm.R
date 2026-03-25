#' Bayesian Variable Selection with Missing data for generalized linear
#' regression models
#'
#' Computation and summaries of posterior distribution over the model space
#' in problems of small to moderate size when missingness occurs in generalized
#' linear models.
#'
#' The set of competing models is made up by all the possible subsets of
#' regressors specified by \code{formula}: Mi for i in 1,...,2^p, being p the number
#' of potential (non-fixed) regressors in the variable selection problem. The simplest,
#' nested in all of them, contains only the intercept. \code{MissingBvs} performs
#' \code{n.imp} imputations given by \code{imp.mice.method} with the mice package and
#' computes the posterior distribution over this model space through Bayes' theorem:
#'
#' Pr(Mi | \code{data})=Pr(Mi)*Bi/C,
#'
#' where Bi is an approximation of the Bayes factor of Mi to M0 under missing data,
#' Pr(Mi) is the prior probability of Mi and C is the normalizing constant.
#'
#' Bi is computed as the average of the \code{n.imp} standard Bayes factors, B_i(j)
#' for j in 1,...,\code{n.imp}, for Mi to M0 calculated for the jth imputed data set.
#' If the BF computation method chosen is \code{"gprior"}, BF B_i(j) depends on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}
#' and \code{MissingBvs.glm} uses the same choices as
#' \code{[MissingBVS]{missingBVS.lm}} for generalized linear models developed by
#' Li and Clyde (2018). It does so using the \pkg{BAS} logmarginal computation
#' (Clyde, 2025) if the \code{family} chosen is one of the implemented in the package.
#'
#' The BF can also be approximated with the BIC (Schwarz, 1978) or the test-based
#' BF (Held, Gravestock and Sabanés, 2015). It can be done through the \pkg{BAS}
#' faster computation if the \code{family} is one of the implemented in the package.
#'
#' The prior over the model space Pr(Mi) offers three possibilities:
#' "Constant" assigns the same prior probability to every model.
#' "ScottBerger" is the default choice. It assigns the same prior probability to
#' every possible model dimension and, therefore, accounts for multiplicity issues
#' (Scott and Berger 2010).
#' "User" (see below).
#'
#' If \code{prior.models}="User" is chosen, user has to provide a p+1 dimensional
#' parameter vector with the model dimension prior probabilities through \code{priorprobs}.
#' The first component of \code{priorprobs} must contain the probability of the
#' model with fixed variables; next p components correspond to the p prior probabilities
#' of the possible model dimensions.
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param null.model Formula defining which is the simplest (null) model, which.
#' should be nested in the full one. By default, it is defined to be the one
#' with just the intercept.
#' @param data Data frame containing the data.
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
#' @param prior.models Prior distribution over the model space (to be literally specified). Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the model space of the
#' factor levels (to be literally specified). Possible choices are "Constant" and
#' "ScottBerger" (see details).
#' @param priorprobs A p+1 (being p the number of non-fixed variables)
#' dimensional vector defining the prior model probabilities (used for chosen
#' \code{prior.models}= "User"; see details.)
#' @param n.keep Number of the most probable models kept. By default it is set to
#' 10 and automatically adjusted if 10 is greater than the total number of models.
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
#' @return \code{missingBVS.glm} returns an object of class \code{missingBVS}
#' with the following elements:
#' \item{time}{The internal time consumed in solving the problem}
#' \item{glmfull}{The \code{glm} class object that results when the model
#' defined by \code{formula} is fitted by \code{\link[stats]{glm}}}
#' \item{glmnull}{The \code{glm} class object that results when the null model,
#' the one with just the intercept term, is fitted by \code{\link[stats]{glm}}}
#' \item{variables}{Names of all the potential (non-fixed) explanatory variables}
#' \item{n}{Number of observations}
#' \item{p}{Number of explanatory variables to select from}
#' \item{k}{Number of fixed variables}
#' \item{HPMbin}{Binary expression of the Highest Posterior Probability model}
#' \item{MPMbin}{Binary expression of the Median Probability model}
#' \item{positions}{\code{matrix} with L rows and p plus the number of dummies
#' resulting from factors columns - L, where L is the number of factors, with 1
#' if the column dummy corresponds to the row factor and 0 otherwise}
#' \item{positionsx}{Vector of length p with 1 if the variable is a numerical
#' covariate and 0 otherwise}
#' \item{modelsprob}{\code{data.frame} which summaries the \code{n.keep}
#' most probable a posteriori models and their associated probability}
#' \item{inclprob}{Named vector with the inclusion probabilities of the potential
#' explanatory variables.}
#' \item{postprobdim}{Posterior probabilities over the true model dimension}
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
#' \item{method}{\code{Full}}
#'
#' @author Carolina Mulet and Gonzalo Garcia-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for a heuristic
#' approximation based on Gibbs sampling (recommended when p>20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
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
#' @examples #To be completed
#'

missingBVS.glm <- function (formula,
                            null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                            data,
                            family = binomial(link = "logit"),
                            BF.approx.method = "BIC",
                            prior.betas = "Robust", #if BF.approx.method = "gprior"
                            prior.models = "ScottBerger",
                            prior.models.dummies = "ScottBerger",
                            priorprobs = NULL, #needed if prior.models = User
                            n.keep = 10,
                            parallelmice = NULL,
                            n.core = NULL,
                            imp.time.test = TRUE,
                            imp.mice.method = "pmm", #mice's default
                            n.imp = 039E1, #number of imputed datasets for BF
                            imp.seed = runif(1,0,09011975),
                            weights = rep.int(1, dim(data)[1]),
                            offset = rep.int(0, dim(data)[1]),
                            control = glm.control(),
                            laplace = 0L) {

  time <- Sys.time()

  formula <- as.formula(formula)
  null.model <- as.formula(null.model)

  #Response in the null model and full model must coincide
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
  auxnull <- model.frame(null.model, data, na.action = NULL)
  namesnull.toimp <- dimnames(auxnull)[[2]][-1] #name of fixed variables to imputation
  q0 <- length(namesnull.toimp) + 1 # the intercept

  #Missing model matrix of fixed vars
  X0 <- model.matrix.rankdef(auxnull)
  namesnull <- dimnames(X0)[[2]]
  p0 <- dim(X0)[2] #Number of fixed vars

  #Eval the full model
  glmfull <- glm(formula,
                 data,
                 y = TRUE, x = TRUE,
                 family = family,
                 weights = weights,
                 offset = offset,
                 control = control) #omits NA observations

  #Full design matrix for imputation
  auxfull <- model.frame(formula, data, na.action = NULL)
  namesx.toimp <- dimnames(auxfull)[[2]][-1] #name of variables to imputation
  namesxnotnull.toimp <- namesx.toimp[namesx.toimp %notin% namesnull.toimp]
  X.toimp <- data[,c(namesnull.toimp, namesxnotnull.toimp)] #design matrix with missing data with fixed vars

  #Model matrix data with missings
  X.full <- model.matrix.rankdef(auxfull)
  namesx <- dimnames(X.full)[[2]]
  namesxnotnull <- namesx[namesx %notin% namesnull]
  X.full <- X.full[, namesxnotnull]
  p <- length(namesxnotnull) #Number of covariates and levels of factors to select from

  #Factors: positions has number of rows equal to the number of regressors
  #(factors or numeric covariates) and p columns.
  #A 1 in a row denotes the position in X of a regressor (several positions for
  #the dummies of a factor).
  depvars <- setdiff(attr(terms(auxfull), "term.labels"),
                     attr(terms(auxnull), "term.labels"))

  positions <- t(sapply(depvars, function(var) {
    if(is.factor(data[[var]])) {
      levs <- levels(data[[var]])
      ind <- which(namesxnotnull %in% paste0(var,levs)) #1 in the namelevel matches
    } else ind <- which(namesxnotnull == var) #1 in the name matches

    posi <- rep(0,p); posi[ind] <- 1
    posi
  }))
  colnames(positions) <- namesxnotnull

  #vector of length p with 1 if numeric variable
  positionsx <- as.numeric(colSums(positions %*% t(positions)) == 1)

  L <- sum(!positionsx) #Number of factors to select from
  temp <- rowSums(positions %*% t(positions))
  l <- temp[temp > 1] #Number of levels for each factor

  q <- p - sum(l) + L #Number of factors and covariates to select from
  #q = p if there are no factors

  #matrix of dim (Lxp) with 1 if dummy variable of the row factor
  positionsfac <- positions[!positionsx*1:q,]

  #check if null model is contained in the full one:
  for (i in 1:p0){
    if (namesnull[i] %notin% namesx) {
      stop(paste0("Error in var: ", namesnull[i],"; null model not nested in full model.\n"))
    }
  }

  #Is there any variable to select from?
  if (p == p0) {
    stop(paste0("The number of fixed variables is equal to the number of\n",
                "regressors in the full model. No model selection can be done.\n"))
  }

  #The response variable
  y <- glmnull$y; obsnotNA <- names(y) #without missings
  y <- as.numeric(y)
  n <- length(y) #observations without missings on the response
  devnull <- glmnull$deviance #deviance of the null model

  #check for missings and define variables with NAs
  NAvars <- checkformissings.glm(y = auxnull[,1], X0, X.full, obsnotNA)

  #n.keep > 2^p, the number of models?
  if (n.keep > 2^p) {
    cat(paste0("The number of models to keep (", n.keep,
               ") is larger than the total number of models (",
               2^p, ") and it has been set to ", 2^p,".\n"))
    n.keep <- 2^p
  }

  #check if the number of regressors is too big.
  if (p > 20) {
    warning("Number of regressors too big. . . consider using missingGibbsBvs.lm.\n",
            immediate. = TRUE)
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

  #Imputation of missing data
  if (is.null(parallelmice)) {
    if (n.imp > 120 | n*q > 50000) {
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

  #remove observations with missings on the response or fixed vars
  imputation.array <- imputation.array[obsnotNA,,]
  #function to compute log(BFa0) for a given model as an average of BF computed
  #by BF.approx.method over the imputed datasets
  lBF.method <- function (model) lBF.approx(model,
                                            imputation.array = imputation.array,
                                            BF.approx.method = BF.approx.method,
                                            p0 = p0, n.imp = n.imp)

  #Info:
  cat("Info. . .\n")
  cat("Most complex model has a total of", q + q0, "covariates and/or factors.\n")
  if (q0 == 1) {
    cat(paste0("From those 1 is fixed (the intercept) and we should select from the remaining ",
               q, ".\n"))
  } else cat(paste0("From those ", q0, " are fixed and we should select from the remaining ",
                    q, ".\n"))

  cat("  Numerical covariates:", depvars[positionsx == 1], "\n")
  if (L > 0) cat(" Factors:", depvars[positionsx == 0], "\n")

  cat("The problem has a total of", 2^p, "competing models.\n")
  cat("Of these, the ", n.keep, "most probable (a posteriori) are kept.\n")

  #progress bar for loop
  pb <- txtProgressBar(min = 0,
                       max = 2^p,
                       style = 3,
                       width = 50,
                       char = "=")

  #Posterior computation
  all.models.lPM <- matrix(0, nr = 2^p, nc = p+1) #last column contains log(BF_a0*Pr(M))
  for (i in seq_len(2^p-1)){ # null out of the loop
    setTxtProgressBar(pb, i)

    #transform the number of the model into a binary number
    current.model <- BayesVarSel:::integer.base.b_C(i, p)
    all.models.lPM[i, seq_len(p)] <- current.model

    gamma.tau <- positions %*% current.model > 0 #covariates and/or factors active
    deltasum <- positionsfac %*% current.model #levels of factors
    tau <- deltasum > 0 #factors

    #check if there are NAs in the model considered to save computation time
    if (any(namesxnotnull[which(current.model == 1)] %in% NAvars)) {
      lBF.PM <- lBF.method(model = which(current.model == 1)) +
                lprior.models(gamma.tau) + lprior.models.dummies(deltasum, tau) #log(BF_a0*Pr(M))
    } else { #if there are no missings, compute the BF by the method selected
      X.i <- cbind(X0[obsnotNA,], X.full[obsnotNA, which(current.model == 1)])
      lBF.PM <- BF.approx.method(k = sum(current.model),
                                 X = as.matrix(X.i)) +
        lprior.models(gamma.tau) + lprior.models.dummies(deltasum, tau) #log(BF_a0*Pr(M))
    }
    all.models.lPM[i, p+1] <- lBF.PM
  }
  setTxtProgressBar(pb, 2^p)
  cat("\n")
  #null model
  all.models.lPM[2^p, seq_len(p)] <- rep(0, p)
  all.models.lPM[2^p, p+1] <- lprior.models(rep(0, p)) #BF = 1 for null model

  #renormalize
  C <- sum(exp(all.models.lPM[, p+1]))
  all.models.PM <- all.models.lPM
  all.models.PM[, p+1] <- exp(all.models.lPM[, p+1] - log(C))
  colnames(all.models.PM) <- c(namesxnotnull, "Post")

  #models matrix at the covariate-factor level
  cf.models.PM <- all.models.PM[,seq_len(p)] %*% t(positions)
  cf.models.PM <- cbind(cf.models.PM, all.models.PM[,p+1])
  colnames(cf.models.PM)[q+1] <- "Post"
  #cf.models.PM is exactly all.models.PM if there are no factors

  inclprob <- rep(0, q)
  probdim <- rep(0, q + 1)
  #compute inclusion probabilities (except for fixed variables) and
  #posterior probability of the dimension of the true model
  for (i in seq_len(2^p)) {
    inclprob[which(cf.models.PM[i, seq_len(q)] == 1)] <-
      inclprob[which(cf.models.PM[i, seq_len(q)] == 1)] + cf.models.PM[i, q + 1]
    probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] <-
      probdim[sum(cf.models.PM[i, seq_len(q)] > 0) + 1] + cf.models.PM[i, q + 1]
  }

  #HPM
  nPmax <- which.max(cf.models.PM[, q+1])
  hpm <- cf.models.PM[nPmax, ]

  #MPM
  mpm <- rep(0,q)
  mpm[which(inclprob >= 0.5)] <- 1

  #result
  result <- list()
  result$time <- Sys.time() - time #The time it took the program to finish
  result$glmfull <- glmfull # The lm object for the full model (without NAs)
  result$glmnull <- glmnull # The lm object for the null model

  result$variables <- depvars #The name of the competing variables
  result$n <- n #number of observations
  result$p <- q #number of competing vars
  result$k <- q0 #number of fixed vars
  result$HPMbin <- hpm #The binary code for the HPM model
  names(result$HPMbin) <- C(depvars, "Post")
  result$MPMbin <- mpm #The binary code for the MPM model
  names(result$MPMbin) <- depvars

  if (L > 0) {
    #matrix for the factors index
    result$positions <- positionsfac
    result$positionsx <- positionsx
  }

  result$modelsprob <- cf.models.PM[order(cf.models.PM[,q+1],
                                          decreasing = TRUE)[seq_len(n.keep)],]
  #The binary code for the n.keep best models and the correspondent post
  result$inclprob <- inclprob #inclusion probability for each variable
  names(result$inclprob) <- depvars

  result$postprobdim <- probdim #vector with the dimension probabilities.
  names(result$postprobdim) <- 0:q + q0 #dimension of the true model

  result$call <- match.call()

  if(!identical(lprior.models, logUser)){
    priorprobs <- rep(0, q + 1)
    priorprobs[1] <- exp(lprior.models(rep(0, q))) #prior inclusion prob for dimension 0
    for (i in seq_len(q)) {
      priorprobs[i+1] <- exp(lprior.models(c(rep(1, i), rep(0, q - i))) + lchoose(q, i))
      #prior inclusion probability for each dimension
    }
  }
  result$priorprobs <- priorprobs
  names(result$priorprobs) <- 0:q + q0 #prior dimension probability

  result$C <- C #normalizing constant

  #arguments used for imputation
  result$imp.args <- list(parallelmice = parallelmice,
                          imp.mice.method = imp.mice.method,
                          n.imp = n.imp, imp.seed = imp.seed)

  #save the imputed datasets for sensitivity analysis
  # raw.imp.array <- serialize(imputation.array, NULL)
  # result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")

  result$BF.approx.method <- BF.approx.method #function used for BF computation
  result$prior.betas <- prior.betas
  result$logprior.models <- lprior.models #function used for model prior
  if (L > 0) result$logprior.models.dumm <- lprior.models.dummies

  result$method <- "Full"
  class(result) <- "MissingBvs"

  return(result)
}

"%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0 #auxiliar function

checkforfamily <- function (family, BF.approx.method) {
  #if family is a character redefines it with the corresponding function,
  #if family is a function, calls it and checks the provided argument is among the possible options.
  #Returns a logical indicating if BAS functions can be used to speed up the process
  #if method is BIC or TBF
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()

  #families implemented in BAS logmarginal computation
  if (family$family %notin% c("binomial", "poisson", "Gamma")) {
    inBAS <- FALSE
  } else {
    if ((family$family == "binomial" & family$link != "logit") |
        (family$family %in% c("poisson", "Gamma") & family$link != "log")) {
      inBAS <- FALSE
    } else inBAS <- TRUE
  }

  if (BF.approx.method == "gprior" & !inBAS) {
    stop(paste("family ", family, "not implemented in BAS' marginal computation.\n",
               "Try with method 'BIC' or 'TBF'.\n"))
  }
  return(inBAS)
}

checkformissings.glm <- function (y, X0 = NULL, X.full, obsnotNA = NULL) {
  #checks if there are missings on the response and regressors
  #and returns the name of non-fixed regressors with missings

  ##on the response
  if (sum(is.na(y)) > 0) {
    cat("NA values found on the response variable.",
        "We are going to omit these observations.\n")
  }
  ##on the fixed vars
  if (sum(is.na(X0)) > 0) {
    cat("NA values found on the fixed variables.",
        "We are going to omit these observations.\n")
  }
  ##on the regressors
  if (sum(is.na(X.full)) == 0) {
    stop(paste0("NA values not found for any of the competing variables.\n",
                "We recommend you the BAS package instead.\n"))
  } else {
    O <- 1*(!is.na(X.full[obsnotNA,]))
    #zeros where missing observations on the data without missings on the response
    NAvars <- names(which(colSums(O) < dim(X.full)[1])) #columns with missings
  }
  return(NAvars)
}

checkforprior.betas.glm <- function (BF.approx.method, prior.betas, inBAS,
                                     n, p, p0, y, null.model,
                                     data, family, devnull, weights, offset, control, laplace) {
  #checks that the Bayes factor computation method given by BF.approx.method and prior.betas
  #is implemented and returns the function to use for Bayes factor computation on glm
  if (BF.approx.method %in% c("BIC", "TBF", "gprior")) {
    if (inBAS) { #use faster computation of BAS
      if (BF.approx.method != "BIC") {
        if (prior.betas %notin% c("gZellner", "Robust", "Liangetal", # "ZellnerSiow",
                                  "FLS", "intrinsic.WNC" # "IHG"
        )) {
          stop(paste0("For now, prior.betas must be one of 'gZellner', 'Robust', 'Liangetal',\n",
                      "'FLS' or 'intrinsic.WNC' when using TBF or gprior method.\n"))
        }

        switch (prior.betas, # change the string for the corresponding BAS function
                gZellner = {prior.betas <- BAS::g.prior(g = n)}, #fixed g=n
                Robust = {prior.betas <- BAS::robust(n = n)}, #random g
                Liangetal = {prior.betas <- BAS::hyper.g.n(alpha = 3, n = n)}, #random g: hyper-g/n with a=3
                # `Zellner-Siow` = {prior.betas <- "ZSBF"}, #random g: cauchy prior, Jeffreys in BAS?
                FLS = {prior.betas <- BAS::g.prior(g = max(n, p^2))}, #fixed Benchmark prior: g=max(n, p*p)
                `intrinsic.WNC` = {prior.betas <- BAS::intrinsic(n = n)} #intrinsic prior from Womack, Novelo and Casella (2014)
                # IHG = {prior.betas <- "geointrinsicBF"} #intrinsic hyper-g prior, not available in BAS?
        )

      } else prior.betas <- BAS::bic.prior(n = n)

      logmargnull <- BAS::bas.glm(formula = null.model, # nullmodel
                                  data = data,
                                  n.models = 1,
                                  method = "MCMC",
                                  betaprior = prior.betas,
                                  family = family,
                                  weights = weights,
                                  offset = offset,
                                  control = control,
                                  laplace = laplace)$logmarg

      switch (BF.approx.method,
              BIC = {BF.approx.method.f <-
                function (k, X) BF.approx.BIC.glm(y = y, X,
                                                  family = family,
                                                  logmargnull = logmargnull,
                                                  k, p0 = p0,
                                                  weights = weights,
                                                  offset = offset,
                                                  control = control,
                                                  laplace = laplace)},
              TBF = {BF.approx.method.f <-
                function (k, X) BF.approx.TBF.glm(y = y, X,
                                                  family = family,
                                                  devnull = devnull,
                                                  prior.betas = prior.betas,
                                                  k, p0 = p0,
                                                  weights = weights,
                                                  offset = offset,
                                                  control = control,
                                                  laplace = laplace)},
              gprior = {BF.approx.method.f <-
                function (k, X) BF.approx.gprior.glm(y = y, X,
                                                     family = family,
                                                     logmargnull = logmargnull,
                                                     prior.betas = prior.betas,
                                                     k, p0 = p0,
                                                     weights = weights,
                                                     offset = offset,
                                                     control = control,
                                                     laplace = laplace)}
      )
    } else { #use slower option that not depends on BAS
      cat("Faster functions from BAS package to compute marginal likelihoods",
          "cannot be used for specified arguments and it can take a while.\n",
          "Do you want to continue? (y/n)\n")
      if (tolower(readline()) != "y") {
        stop("Consider another arguments for your problem.\n")
      }

      switch (BF.approx.method,
              BIC = {BF.approx.method.f <-
                function (k, X) BF.approx.BIC.glm.stats(y = y, X,
                                                        family = family,
                                                        devnull = devnull,
                                                        n, k, p0 = p0,
                                                        weights = weights,
                                                        offset = offset,
                                                        control = control)},
              TBF = {BF.approx.method.f <-
                function (k, X) BF.approx.TBF.glm.stats(y = y, X,
                                                        family = family,
                                                        devnull = devnull,
                                                        n, k, p0 = p0,
                                                        weights = weights,
                                                        offset = offset,
                                                        control = control)}
      )
    }

    return(BF.approx.method.f)

  } else {
    stop("Only approximations 'BIC', 'TBF' and 'gprior' supported.")
  }
}

