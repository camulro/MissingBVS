#' \pkg{MissingBVS}'s imputation for normally distributed covariates
#'
#' Performs the multiple imputation developed by García-Donato et al (2025) for
#' continuous and normally distributed covariates to compute Bayesian Variable
#' Selection in the presence of missing data.
#'
#' @export
#' @param X Matrix with missing values to impute.
#' @param nMC Number of samples used to approximate, by MonteCarlo, the integral
#' defining the Bayes factor.
#' @param seed Seed chosen for the gibbs sampler on the imputation step.
#' @param initialimp.mice.method Method used by \code{\link[mice]{mice}} to impute the
#' initial values. See \code{\link[mice]{mice}} for possible choices.
#' @param time.test Logical to indicate whether to check time of performance with
#' \code{nMC = 10} or not.
#'
#' @return \code{MC.imputation} returns an object of class \code{MissingBVS.imputation}
#' with the following elements:
#' \item{rX.imput}{Array of dimension nxpx\code{nMC}, where n is the number of
#' observations and p the number of competing covariates, containing the imputed datasets}
#' \item{rSigma}{Array of dimension pxpx\code{nMC} containing the corresponding
#' covariance matrices}
#' \item{rmu}{Vector of dimension px\code{nMC} containing the corresponding means}
#'
#' @author María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.miss}} for computing the logarithm
#' of the MC approximation of the Bayes factor in linear models.
#' Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20).
#'
#' @examples #To be completed
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
MC.imputation <- function(X, nMC = 039E1,
                          seed = runif(1,0,09011975), initialimp.mice.method = "pmm",
                          time.test = FALSE){
  #(Works for continuous covariates)

  if (time.test) {time <- Sys.time(); nMC <- 10} # to estimate imputation time
  #results:
  rX.imput <- array(0, dim=c(dim(X), nMC))
  rSigma <- array(0, dim=c(dim(X)[2], dim(X)[2], nMC))
  rmu <- matrix(0, nr=dim(X)[2], nc=nMC)

  p <- dim(X)[2]; n <- dim(X)[1]
  if (p < 2) stop("It does not work with p<2")

  O <- 1*(!is.na(X))
  #the following units have at least one NA
  these <- which(rowSums(O) < p)

  #Initial value imputed with mice
  imputed <- mice::mice(X,
                        meth = initialimp.mice.method,
                        m = 1,
                        printFlag = FALSE)
  X.full <- as.matrix(mice::complete(imputed))

  ###Gibbs sampler
  set.seed(seed)
  for(s in seq_len(nMC)) {

    #posterior dist. with Jeffreys independent prior
    Sigma <- LaplacesDemon::rinvwishart(nu=n-1, S=(n-1)*var(X.full))
    mu <- LaplacesDemon::rmvn(1, colMeans(X.full), Sigma/n)

    ###update missing data
    for(i in these) {
      b <- ( O[i,]==0 )
      a <- ( O[i,]==1 )
      solve.Sigma.a.a <- solve(Sigma[a,a])
      thetab.mid.a <- mu[b]+Sigma[b,a]%*%solve.Sigma.a.a%*%(X.full[i,a]-mu[a])
      Sigmab.mid.a <- Sigma[b,b] - as.matrix(Matrix::forceSymmetric(Sigma[b,a]%*%solve.Sigma.a.a%*%Sigma[a,b]))
      X.full[i,b] <- LaplacesDemon::rmvn(1, as.vector(thetab.mid.a), Sigmab.mid.a)
    }

    rX.imput[,,s] <- scale(X.full[, drop = FALSE], center = TRUE, scale = FALSE) # already centered
    rSigma[,,s] <- Sigma
    rmu[,s] <- mu
  }
  if (time.test) return(time <- Sys.time() - time)

  dimnames(rX.imput)[[1]] <- seq_len(n)
  dimnames(rX.imput)[[2]] <- colnames(X)
  dimnames(rSigma)[[1]] <- dimnames(rSigma)[[2]] <- colnames(X)

  imputation.list <- list(rX.imput = rX.imput, rSigma = rSigma, rmu = rmu)
  class(imputation.list) <- "MissingBVS.imputation"
  return(imputation.list)
}
#' \pkg{mice}'s imputation for Missing Bayesian Variable Selection (\pkg{MissingBVS})
#'
#' Performs the Multiple Imputation by Chained Equations of \pkg{mice} and makes
#' an array with the proper form to compute Bayesian Variable Selection in the
#' presence of missing data. It can be applied to any kind of data with the correct
#' imputation method.
#'
#' A parallel version of \pkg{mice} can be performed through the \code{parallelmice}
#' argument, which makes the imputation a lot faster for a big number of imputations
#' or huge datasets.
#'
#' @export
#' @param X Matrix with missing values to impute.
#' @param formula Formula defining the most complex (full) regression model.
#' @param n.imp Number of imputed datasets to compute the average Bayes factor.
#' @param imp.mice.method Method used by \code{\link[mice]{mice}} to impute the
#' values. See \code{\link[mice]{mice}} for possible choices.
#' @param seed Seed chosen for the imputations.
#' @param parallelmice Logical to indicate whether or not to use parallel
#' \code{\link[mice]{mice}} imputation. By default performs parallel mice
#' imputation if the number of imputations is big enough (\code{n.imp > 120}).
#' @param n.core See \code{\link[mice]{futuremice}} for details.
#' @param time.test Logical to indicate whether to check time of performance with
#' \code{n.imp = 10} or not.
#'
#' @return An object of class \code{MissingBVS.imputation}, an array of dimension
#' nxpx\code{n.imp}, where n is the number of observations and p the number of
#' competing covariates, containing the imputed datasets.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} for computing the logarithm
#' of the average Bayes factor among the \code{n.imp} imputed datasets.
#' Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20)
#' in linear models and \code{\link[MissingBVS]{MissingBvs.glm}} for generalized
#' linear models.
#'
#' @examples #To be completed
#'
#' @references van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice:
#' Multivariate Imputation by Chained Equations in R. Journal of Statistical
#' Software. 45(3): 1–67.
#'
#' Volker, T.B. and Vink, G. (2022). futuremice: The future starts today.
#' \link{https://www.gerkovink.com/miceVignettes/futuremice/Vignette_futuremice.html}
#'
mice.imputation <- function(X, formula, n.imp = 039E1,
                            imp.mice.method = "pmm", seed = runif(1,0,09011975),
                            parallel = n.imp > 120, n.core = NULL, time.test = FALSE) {
  formula <- paste(formula)

  if (time.test) {time <- Sys.time(); n.imp <- 30} # to estimate imputation time
  if (!parallel) {
    imput <- mice::mice(X,
                        meth = imp.mice.method,
                        m = n.imp,
                        predictorMatrix = mice::quickpred(X),
                        printFlag = FALSE,
                        seed = seed)
   } else imput <- mice::futuremice(X,
                                    meth = imp.mice.method,
                                    m = n.imp,
                                    predictorMatrix = mice::quickpred(X),
                                    parallelseed = seed,
                                    n.core = n.core)

  imps <- mice::complete(imput, action = "all") #extracts all at once

  #Get final dim of full imputed datasets
  n <- nrow(X)
  X.formula <- as.formula(paste(formula[1], formula[3]))
  aux <- model.matrix.rankdef(model.frame(X.formula, X, na.action = NULL))
  q <- ncol(aux)

  imputation.array <- array(0, dim = c(n, q, n.imp)) #an array with the matrices imputed
  for (s in seq_len(n.imp)) {
    aux.imps <- model.frame(X.formula, imps[[s]])
    imputation.array[, , s] <- model.matrix.rankdef(aux.imps) #build the model matrix
  }
  if (time.test) return(time <- Sys.time() - time)

  dimnames(imputation.array)[[1]] <- seq_len(n)
  dimnames(imputation.array)[[2]] <- colnames(aux)

  class(imputation.array) <- "MissingBVS.imputation"
  return(imputation.array)
}
#' Binary matrix for missing data pattern
#'
#' Generates a matrix with 0 if the original entry was missing and 1 otherwise
#' for observations with missing data and summarizes by variable.
#'
#' @export
#' @param data Data frame containing the data with possible missing values.
#' @param formula Formula defining the most complex (full) regression model.
#' If \code{NULL}, the full data missing entries are chosen.
#' @param show Logical that indicates whether or not to print que matrix pattern.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MC.imputation}},
#' \code{\link[MissingBVS]{mice.imputation}} and
#' \code{\link[MissingBVS]{futuremice.imputation}} for performing multiple
#' imputed dataset for \pkg{MissingBVS}.
#'
#' @examples #To be completed
#'
missing.model <- function (data, formula = NULL, show = TRUE) {
  #data is a matrix or dataframe
  #formula can either be null or a model formula
  data <- data.frame(data)

  #if formula not provided, all variables from data are shown
  if (is.null(formula)) {
    data.model <- data
  } else {
    formula <- as.formula(formula)
    data.model <- model.frame(formula, data, na.action = NULL)
    data.model <- data.model[,which(colnames(data.model) != formula[[2]])] #remove response
  }

  #observations missed
  O <- 1*(!is.na(data.model)) #1=observed, 0=missed
  #the following units have at least one NA
  these <- which(rowSums(O) < ncol(data.model))

  if (length(these) == 0) {
    stop("No missing values for the given model.\n")
  }

  mO <- O[these,]
  rownames(mO) <- these

  Total <- dim(data)[1] - colSums(O) #number of missings per variable
  ord <- order(Total)

  if (show) {
    cat("\nVariables ordered by the number of missings:\n")
    print(mO[,ord])
    cat("---\n")
    cat("Code: 0 = missed, 1 = observed\n")
    cat("  Rows correspond to observations with missings.\n")
    cat("---\nThe total number of missings per variable is:\n")
    print(Total[ord])
    cat("\n")
  }
  mO <- rbind(mO, Total)[,ord]
  return(mO)
}

#' Box-and-whisker plot for numerical observed and imputed data
#'
#' Produce box-and-whisker plots to represent the observed, given by \code{X},
#' vs the imputed, given by \code{imputation}, values for each numerical regressor
#' given by \code{formula}.
#'
#' @export
#' @param X Matrix with missing values to impute.
#' @param imputation Object of class \code{MissingBVS.imputation} containing the
#' imputed datasets.
#' @param formula Formula defining the most complex (full) regression model.
#' @param mfrow Vector of the form \code{c(nr, nc)} for the number of rows and
#' columns respectively. It selects the layout of the figures as an nr-by-nc array.
#' If \code{NULL}, it will be automatically calculated.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MC.imputation}},
#' \code{\link[MissingBVS]{mice.imputation}} and
#' \code{\link[MissingBVS]{futuremice.imputation}} for generating objects of
#' class \code{MissingBVS.imputation}.
#'
#' @examples #To be completed
#'
plot.MissingBVS.imputation <- function (X, imputation, formula, mfrow = NULL) {
  #X is a matrix or dataframe with missing data
  #imputation.array is the array with the whole imputed data
  #formula can either be null or a model formula
  formula <- as.formula(formula)

  if (inherits(imputation, "MissingBVS.imputation")) {
    if (is.list(imputation)) {
      imputation.array <- imputation$rX.imput
    } else imputation.array <- imputation
  } else stop("Only imputation plot available for class MissingBVS.imputation.\n")

  aux <- model.frame(formula, X, na.action = NULL)[,-1]
  #which variables are numerical
  isnum <- sapply(aux, is.numeric)
  X.full <- as.matrix(aux[,isnum])

  missing.matrix <- missing.model(X.full, show = FALSE)
  NAvars <- names(which(missing.matrix["Total",] > 0)) #original vars with missings
  missings <- missing.matrix[-nrow(missing.matrix), NAvars]

  n <- dim(missings)[1]
  p <- dim(missings)[2]
  n.imp <- dim(imputation.array)[3]

  #keep imputed vars
  imputation.array.model <- array(imputation.array[rownames(missings),NAvars,],
                                  dim = c(n, p, n.imp)) #intercept removed
  dimnames(imputation.array.model)[[1]] <- rownames(missings)
  dimnames(imputation.array.model)[[2]] <- NAvars

  #calculate how to grid the composed plot
  if (is.null(mfrow)) {
    nr <- ceiling(sqrt(p))
    nc <- ceiling(p/nr)

    par(mfrow = c(nr, nc), mar = c(2.5,2,0.5,0.5), mgp = c(1.5, 0.5, 0))
  } else if (length(mfrow) != 2) {
    stop("Please provide a vector of the form c(nrow, ncol) to grid the plot.\n")
  } else par(mfrow = mfrow, mar = c(2.5,2,0.5,0.5), mgp = c(1.5, 0.5, 0))

  for (i in NAvars) {
    x.mi <- as.vector(imputation.array.model[which(missings[,i] == 0),i,]) #imputed
    x.oi <- X.full[!is.na(X.full[,i]), i] #observed
    x.i <- c(x.mi, x.oi)

    gr <- c(rep("Imp", length(x.mi)), rep("Obs", length(x.oi)))

    boxplot(x.i ~ gr,
            col = c(rgb(1, 0.3, 0.3, 0.5), rgb(0,0,0,0.5)),
            horizontal = TRUE,
            xlab = paste0(i, ", NA's: ", length(x.mi)/n.imp), ylab = "",)
    stripchart(x.i ~ gr,
               method = "jitter",
               vertical = FALSE,
               pch = 16,
               col = c(rgb(1, 0.3, 0.3, 0.5), rgb(0,0,0,0.5)),
               add = TRUE)
  }
  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0)) #R default values
}
