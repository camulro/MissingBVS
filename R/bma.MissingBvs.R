#' #' Bayesian Model Averaging for regression coefficients in linear models
#' #'
#' #' Samples of model averaged objective posterior distribution of regression
#' #' coefficients:
#' #' \eqn{latex}{sum_M f(\beta | data, M) Pr(M | data),}
#' #' for the retained best (\code{n.keep} in \code{\link[MissingBVS]{MissingBvs.lm}})
#' #' models or the sampled ones with the associated frequencies in the case of
#' #' (\code{\link[MissingBVS]{MissingGibbsBvs.lm}}). \eqn{latex}{f(\beta | data, M)}
#' #' is calculated following Bernardo and Smith (1994) page 442.
#' #'
#' #' @param mbvs.object Object of class \code{MissingBvs}.
#' #' @param n.sim Number of simulations to produce.
#' #' @param method Matrix decomposition method used to determine the matrix root
#' #' of 'sigma' when drawning values from the multivariate t distribution. Default
#' #' method is value decomposition, see \code{\link[mvtnorm]{rmvt}} for more
#' #' details.
#' #'
#' #' @return \code{BMAcoeff.MissingBvs.lm} returns an object of class
#' #' \code{bma.coeffs.MissingBvs}, a \code{n.sim}x\code{n.sim} matrix whose columns
#' #' correspond to the regression coefficients in the full model.
#' #'
#' #' @author Carolina Mulet and Gonzalo García-Donato
#' #' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#' #'
#' #' @export
#' #'
#' #' @seealso See \code{\link[MissingBVS]{MissingBvs.lm}} and
#' #' \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for creating objects of class
#' #' \code{MissingBvs}.
#' #'
#' #' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' #' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' #' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#' #'
#' #' Bernardo, J. M. and Smith, A. F. M. (1994)<DOI:10.1002/9780470316870>
#' #' Bayesian Theory. Chichester: Wiley.
#' #'
#' #' @examples \dontrun{
#' #' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' #' data("dataS97")
#' #'
#' #' #Here we keep the 8 competing models:
#' #' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' #' dataS97.mBVS <- missingBVS.lm(formula = f, data = dataS97, n.keep = 8, n.imp = 5)
#' #'
#' #' #BMA for model regression coefficients
#' #' dataS97.BMAcoef <- BMAcoeff.MissingBvs.lm(dataS97.mBVS)
#' #' }
#' #'
#' #'
#'
#' BMAcoeff.MissingBvs.lm <- function(mbvs.object, n.sim = 1000, method = "svd") {
#'
#'   if (!inherits(mbvs.object, "MissingBvs")){
#'     warning("An object of class MissingBvs is needed.\n")
#'   }
#'
#'   #Check arguments
#'   if (!is.null(mbvs.object$family)) stop("Only BMA for linear models supported.\n")
#'
#'   X0 <- mbvs.object$lmnull$x #fixed vars
#'   p0 <- ncol(X0) #number of fixed vars
#'   #names of fixed vars including the intercept
#'   nullterms <- c("(Intercept)", attr(terms(mbvs.object$lmnull$terms), "term.labels"))
#'   #response
#'   y <- mbvs.object$lmnull$y
#'   n <- length(y)
#'
#'   #to make orthogonal para metrization
#'   P0 <- X0 %*% solve(t(X0) %*% X0) %*% t(X0)
#'   y <- (diag(n) - P0) %*% y
#'
#'   # are factors present?
#'   positionsfac <- mbvs.object$positions; isfac <- sum(positionsfac) > 0
#'
#'   #load stored imputed datasets:
#'   if (is.null(mbvs.object$compress.imp.array) & !is.null(mbvs.object$imp.args)) { #GD25
#'     imp.list <- unserialize(memDecompress(mbvs.object$compress.imp.list, type = "xz"))
#'     imp.array <- imp.list$rX.imput # without intercept
#'     imp.Sigma <- imp.list$rSigma #variance-covariance matrices
#'
#'     NAvars <- mbvs.object$variables #all variables are random
#'
#'     p <- dim(imp.array)[2] #number of competing variablesç
#'     n.imp <- dim(imp.array)[3]
#'
#'     #function to simulate coefficients:
#'     simbeta <- function(covsrMD, howmany) {
#'       rcoeff <- array(0, dim = c(howmany, length(covsrMD), n.imp),
#'                       dimnames = list(1:howmany, covsrMD, 1:n.imp))
#'       for (j in 1:n.imp) {
#'         rcoeff[,,j] <- simbeta.GD25(y, imp.array[,covsrMD,j], imp.Sigma[covsrMD,covsrMD,j],
#'                                     n, howmany, method)
#'       }
#'       rcoeff
#'     }
#'   } else { #mice's imputation or no imputation
#'     if (!is.null(mbvs.object$imp.args)) {
#'
#'       imp.array <- unserialize(memDecompress(mbvs.object$compress.imp.array, type = "xz"))
#'       # if there is only one imputed dataset
#'       if (is.na(dim(imp.array)[3])) imp.array <- array(imp.array, dim = c(dim(imp.array),1),
#'                                                        dimnames = c(dimnames(imp.array),1))
#'
#'       n.imp <- dim(imp.array)[3]
#'       X <- imp.array[,,1] #for the case when no NAs are present on the model
#'     } else X <- mbvs.object$lmfull$x #lm fit for full model when missingness do not occur
#'
#'     p <- ncol(X) - p0 #number of competing variables
#'     NAvars <- mbvs.object$imp.args$NAvars
#'
#'     #function to simulate coefficients:
#'     simbeta <- function(covsrMD, howmany) {
#'       rcoeff <- array(0, dim = c(howmany, length(covsrMD), n.imp),
#'                       dimnames = list(1:howmany, covsrMD, 1:n.imp))
#'       for (j in 1:n.imp) {
#'         rcoeff[,,j] <- simbeta.gprior(y, imp.array[,covsrMD,j], P0, n, howmany, method)
#'       }
#'       rcoeff
#'     }
#'   }
#'
#'   if (mbvs.object$method == "Full") { #for full enumeration
#'     #define model matrix with rank deficient form if factors present
#'     ifelse(isfac, models.matrix <- mbvs.object$modelsrankdefprob,
#'            models.matrix <- mbvs.object$modelsprob)
#'
#'     n.keep <- dim(models.matrix)[1] #n.keep most probable models
#'     vars <- colnames(models.matrix[,-ncol(models.matrix)]) # competing variables
#'     post <- models.matrix[, "Post"] #posterior prob
#'
#'     cat("\nSimulations obtained using the best", n.keep, "models\n")
#'     cat("that accumulate", round(sum(post), 3), "of the total posterior probability.\n")
#'   } else { #Gibbs aproximation
#'     #define model matrix with rank deficient form if factors present
#'     ifelse(isfac, models.matrix <- mbvs.object$modelsrankdeflogBF.PM,
#'            models.matrix <- mbvs.object$modelslogBF)
#'
#'     n.keep <- dim(models.matrix)[1] #n.keep most probable models
#'     vars <- colnames(models.matrix[,-ncol(models.matrix)]) # competing variables
#'     post <- NULL
#'
#'     cat("\nSimulations obtained using the ", n.keep," sampled models.\n")
#'     cat("Their frequencies are taken as the true posterior probabilities.\n")
#'   }
#'
#'   #draw n.sim models with replacement and with weights proportional to post or frequencies
#'   models <- sample(n.keep, n.sim, replace = TRUE, prob = post)
#'   t.models <- table(models)
#'   cs.tmodels <- cumsum(t.models) #number of repeated models
#'
#'   #save BMA coefficients
#'   bma.coeffs <- matrix(0, nrow = n.sim, ncol = p, dimnames = list(1:n.sim, vars))
#'   for (i in 1:length(t.models)) {
#'     #rMD is model drawn (a number between 1 and n.keep)
#'     rMD <- as.numeric(names(t.models)[i])
#'     howmany <- t.models[i]
#'
#'     #active vars in rMD adding the fixed ones (except the intercept term)
#'     activevars <- models.matrix[rMD, 1:p] == 1
#'     covsrMD <- vars[activevars]
#'
#'     if (!any(covsrMD %in% NAvars)) {#if no missing data in model rMD (e.g. null model)
#'       bma.coeffs[(max(cs.tmodels[i - 1], 0) + 1):cs.tmodels[i], covsrMD] <-
#'         simbeta.gprior(y, X[,covsrMD], P0, n, howmany, method)
#'     } else {
#'       #draw over the imputed datasets
#'       rcoeff <- simbeta(covsrMD, howmany)
#'       bma.coeffs[(max(cs.tmodels[i - 1], 0) + 1):cs.tmodels[i], activevars] <-
#'         apply(rcoeff, 1:2, mean) # average of each entry over the diferent imputations
#'     }
#'   }
#'
#'   class(bma.coeffs) <- "bma.coeffs.MissingBvs"
#'   return(bma.coeffs)
#' }
#'
#' #' @keywords internal
#' simbeta.gprior <- function(y, X, P0, n, howmany, method) {
#'   # orthogonal parametrization
#'   V <- (diag(n) - P0) %*% X
#'
#'   #simulated value for the beta:
#'   qrV <- qr(V)
#'   Rinv <- solve(qr.R(qrV))
#'   iVtV <- Rinv %*% t(Rinv)
#'   mean <- iVtV %*% t(V) %*% y
#'
#'   residuals <- y - V %*% mean
#'   Sigma <- sum(residuals^2) * iVtV / (n - qrV$rank)
#'   rcoeff <- mvtnorm::rmvt(n = howmany, sigma = Sigma, df = n - qrV$rank,
#'                           delta = mean, type = "shifted", method = method)
#'   rcoeff
#' }
#'
#' #' @keywords internal
#' simbeta.GD25 <- function(y, X, S, n, howmany, method) {
#'   #X already centered at imputation
#'
#'   #simulated value for the beta:
#'   qrX <- qr(X)
#'   R <- qr.R(qrX)
#'   iXtXplusS <- solve(t(R) %*% R + S) #A
#'   mean <- iXtXplusS %*% t(X) %*% y
#'
#'   residuals <- y - X %*% mean
#'   Sigma <- sum(residuals^2) * iXtXplusS / n - qrX$rank
#'
#'   rcoeff <- mvtnorm::rmvt(n = howmany, sigma = Sigma, df = n - qrX$rank,
#'                           delta = mean, type = "shifted", method = method)
#'   rcoeff
#' }
#'
#' #' Bayesian Model Averaging for regression coefficients in linear models
#' #'
#' #' Samples of model averaged objective posterior distribution of regression
#' #' coefficients:
#' #' \eqn{latex}{sum_M f(\beta | data, M) Pr(M | data),}
#' #' for the retained best (\code{n.keep} in \code{\link[MissingBVS]{MissingBvs.lm}})
#' #' models or the sampled ones with the associated frequencies in the case of
#' #' (\code{\link[MissingBVS]{MissingGibbsBvs.lm}}). \eqn{latex}{f(\beta | data, M)}
#' #' is calculated following Bernardo and Smith (1994) page 442.
#' #'
#' #' @param mbvs.object Object of class \code{MissingBvs}.
#' #' @param n.sim Number of simulations to produce.
#' #' @param method Matrix decomposition method used to determine the matrix root
#' #' of 'sigma' when drawning values from the multivariate t distribution. Default
#' #' method is value decomposition, see \code{\link[mvtnorm]{rmvt}} for more
#' #' details.
#' #'
#' #' @return \code{BMAcoeff.MissingBvs.lm} returns an object of class
#' #' \code{bma.coeffs.MissingBvs}, a \code{n.sim}x\code{n.sim} matrix whose columns
#' #' correspond to the regression coefficients in the full model.
#' #'
#' #' @author Carolina Mulet and Gonzalo García-Donato
#' #' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#' #'
#' #' @export
#' #'
#' #' @seealso See \code{\link[MissingBVS]{MissingBvs.lm}} and
#' #' \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for creating objects of class
#' #' \code{MissingBvs}.
#' #'
#' #' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' #' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' #' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#' #'
#' #' Bernardo, J. M. and Smith, A. F. M. (1994)<DOI:10.1002/9780470316870>
#' #' Bayesian Theory. Chichester: Wiley.
#' #'
#' #' @examples \dontrun{
#' #' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' #' data("dataS97")
#' #'
#' #' #Here we keep the 8 competing models:
#' #' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' #' dataS97.mBVS <- missingBVS.lm(formula = f, data = dataS97, n.keep = 8, n.imp = 5)
#' #'
#' #' #BMA for model regression coefficients
#' #' dataS97.BMAcoef <- BMAcoeff.MissingBvs.lm(dataS97.mBVS)
#' #' }
#' #'
#' BMApredict.MissingBvs.lm <- function(mbvs.object, newdata = NULL, n.sim = 1000) {
#'   #newdata has to be a data.frame object (with possible missings) if not null (full missings)
#'
#'   if (!inherits(mbvs.object, "MissingBvs")){
#'     warning("An object of class MissingBvs is needed.\n")
#'   }
#'
#'   #Check arguments
#'   if (!is.null(mbvs.object$family)) stop("Only BMA for linear models supported.\n")
#'   if (!is.data.frame(newdata)) {
#'     stop("Argument newdata must be a data.frame.\n")#or null?
#'   }
#'
#'   #Add intercept
#'   if (attr(mbvs.object$lmnull$terms, "intercept") == 1) newdata$`(Intercept)` <- 1
#'
#'   X0 <- mbvs.object$lmnull$x #fixed vars
#'
#'   #names of fixed vars including the intercept
#'   nullterms <- c("(Intercept)", attr(terms(mbvs.object$lmnull$terms), "term.labels"))
#'   #response
#'   y <- mbvs.object$lmnull$y
#'   n <- length(y)
#'
#'   # are factors present?
#'   positionsfac <- mbvs.object$positions; isfac <- sum(positionsfac) > 0
#'
#'   if (mbvs.object$method == "Full") { #for full enumeration
#'     #define model matrix with rank deficient form if factors present
#'     ifelse(isfac, models.matrix <- mbvs.object$modelsrankdefprob,
#'            models.matrix <- mbvs.object$modelsprob)
#'
#'     n.keep <- dim(models.matrix)[1] #n.keep most probable models
#'     vars <- colnames(models.matrix[,-ncol(models.matrix)]) # competing variables
#'     post <- models.matrix[, "Post"] #posterior prob
#'
#'     cat("\nSimulations obtained using the best", n.keep, "models\n")
#'     cat("that accumulate", round(sum(post), 3), "of the total posterior probability.\n")
#'   } else { #Gibbs aproximation
#'     #define model matrix with rank deficient form if factors present
#'     ifelse(isfac, models.matrix <- mbvs.object$modelsrankdeflogBF.PM,
#'            models.matrix <- mbvs.object$modelslogBF)
#'
#'     n.keep <- dim(models.matrix)[1] #n.keep most probable models
#'     vars <- colnames(models.matrix[,-ncol(models.matrix)]) # competing variables
#'     post <- NULL
#'
#'     cat("\nSimulations obtained using the ", n.keep," sampled models.\n")
#'     cat("Their frequencies are taken as the true posterior probabilities.\n")
#'   }
#'
#'   #are NAs in newdata?
#'   isNA <- sum(is.na(newdata)) > 0
#'
#'   #load stored imputed datasets:
#'   if (is.null(mbvs.object$compress.imp.array) & !is.null(mbvs.object$imp.args)) { #GD25
#'     imp.list <- unserialize(memDecompress(mbvs.object$compress.imp.list, type = "xz"))
#'     imp.array <- imp.list$rX.imput # without intercept
#'     imp.Sigma <- imp.list$rSigma #variance-covariance matrices
#'     imp.mu <- imp.list$rmu #means
#'
#'     NAvars <- mbvs.object$variables #all variables are random
#'
#'     p <- dim(imp.array)[2] #number of competing variablesç
#'     n.imp <- dim(imp.array)[3]
#'
#'     # #if NAs in newdata
#'     # if (isNA) {
#'     #   for (j in 1:n.imp){
#'     #     new.imp.list <- MC.imputation(imp.array[,,j], nMC = 1)
#'     #   }
#'     #
#'     #
#'     #
#'     # }
#'
#'     #function to simulate predictions:
#'     simpred <- function(covsrMD, howmany) {
#'       covsrMD <- covsrMD[-1]#remove intercept
#'       rpred <- array(0, dim = c(howmany, dim(newdata)[1], n.imp),
#'                      dimnames = list(1:howmany, 1:dim(newdata)[1], 1:n.imp))
#'
#'       if (length(covsrMD) == 0) {#null model
#'         means <-  rep(mean(y), nrow(newdata)) #predictive means
#'         residuals <- (y - mean(y))
#'         sigmas <- as.matrix(sum(residuals * y) * rep(1,dim(newdata)[1]) / (n - 1)) #predictive variance-covariance matrices
#'         for (j in 1:n.imp) {
#'           rpred[,,j] <- mvtnorm::rmvt(n = howmany, delta = means,
#'                                       sigma = diag(as.numeric(sigmas), nrow = length(sigmas)),
#'                                       df = n - 1, type = "shifted")
#'         }
#'         rpred
#'       } else {
#'         newdataMD <- ifelse(dim(newdata)[1] == 1, matrix(newdata[,covsrMD], ncol = length(covsrMD), nrow = 1),
#'                             as.matrix(newdata[,covsrMD]))
#'         for (j in 1:n.imp) {
#'           rpred[,,j] <- simpred.GD25(y, imp.array[,covsrMD,j], imp.Sigma[covsrMD,covsrMD,j],
#'                                      imp.mu[covsrMD,j], n, newdataMD, howmany)
#'         }
#'         rpred
#'       }
#'     }
#'   } else { #mice's imputation or no imputation
#'     if (!is.null(mbvs.object$imp.args)) {
#'
#'       imp.array <- unserialize(memDecompress(mbvs.object$compress.imp.array, type = "xz"))
#'       # if there is only one imputed dataset
#'       if (is.na(dim(imp.array)[3])) imp.array <- array(imp.array, dim = c(dim(imp.array),1),
#'                                                        dimnames = c(dimnames(imp.array),1))
#'
#'       n.imp <- dim(imp.array)[3]
#'       X <- imp.array[,,1] #for the case when no NAs are present on the model
#'     } else X <- mbvs.object$lmfull$x #lm fit for full model when missingness do not occur
#'
#'     p <- ncol(X) - ncol(X0) #number of competing variables
#'     NAvars <- mbvs.object$imp.args$NAvars
#'
#'     #to make orthogonal parametrization
#'     aux <- solve(t(X0) %*% X0) %*% t(X0)
#'     P0 <- X0 %*% aux
#'     Hnew <- as.matrix(newdata[,nullterms]) %*% aux
#'
#'     #if factors, make rank deficient matrix for newdata
#'     newdata <- model.matrix.rankdef(model.frame(as.formula("~."), newdata, na.action = NULL))
#'     #function to simulate predictions:
#'     simpred <- function(covsrMD, howmany) {
#'       rpred <- array(0, dim = c(howmany, dim(newdata)[1], n.imp),
#'                      dimnames = list(1:howmany, 1:dim(newdata)[1], 1:n.imp))
#'       newdataMD <- ifelse(dim(newdata)[1] == 1, matrix(newdata[,covsrMD], ncol = length(covsrMD), nrow = 1),
#'                           as.matrix(newdata[,covsrMD]))
#'       for (j in 1:n.imp) {
#'         rpred[,,j] <- simpred.gprior(y, imp.array[,covsrMD,j], P0, n, newdataMD, Hnew, howmany)
#'       }
#'       rpred
#'     }
#'   }
#'
#'   #draw n.sim models with replacement and with weights proportional to post
#'   models <- sample(n.keep, n.sim, replace = TRUE, prob = post)
#'   t.models <- table(models)
#'   cs.tmodels <- cumsum(t.models) #number of repeated models
#'
#'   #save BMA predictions
#'   ystar <- paste("y",1:dim(newdata)[1], sep = ".")
#'   rpredictions <- matrix(0, nrow = n.sim, ncol = dim(newdata)[1], dimnames = list(1:n.sim, ystar))
#'   for (i in 1:length(t.models)) {
#'     #rMD is model drawn (a number between 1 and n.keep)
#'     rMD <- as.numeric(names(t.models)[i])
#'     howmany <- t.models[i]
#'
#'     #active vars in rMD adding the fixed ones
#'     activevars <- models.matrix[rMD, 1:p] == 1
#'     covsrMD <- c(nullterms, vars[activevars])
#'     if (!any(covsrMD %in% NAvars)) {#if no missing data in model rMD (e.g. null model)
#'       newdataMD <- matrix(newdata[,covsrMD], ncol = length(covsrMD), nrow = dim(newdata)[1])
#'       rpredictions[(max(cs.tmodels[i - 1], 0) + 1):cs.tmodels[i], ] <-
#'         simpred.gprior(y, X[,covsrMD], P0, n, newdata[,covsrMD], Hnew, howmany)
#'     } else {
#'       #draw over the imputed datasets
#'       rpred <- simpred(covsrMD, howmany)
#'
#'       rpredictions[(max(cs.tmodels[i - 1], 0) + 1):cs.tmodels[i], ] <-
#'         apply(rpred, 1:2, mean) # average of each entry over the diferent imputations
#'     }
#'   }
#'
#'   return(rpredictions)
#' }
#'
#' #' @keywords internal
#' simpred.gprior <- function(y, X, P0, n, newdataMD, Hnew, howmany) {
#'   # orthogonal parametrization
#'   V <- (diag(n) - P0) %*% X
#'   qrV <- qr(V)
#'   Rinv <- qr.solve(qr.R(qrV))
#'   iVtV <- Rinv %*% t(Rinv)
#'   mean <- iVtV %*% t(V) %*% y #betahat
#'   residuals <- (diag(n) - P0) %*% y - V %*% mean
#'
#'   newV <- newdataMD - Hnew %*% X
#'   means <- Hnew %*% y + newV %*% mean #predictive means
#'
#'   fnx <- apply(newV, 1, function(x) as.numeric(1 + t(x) %*% iVtV %*% x))
#'   sigmas <- sum(residuals * y) * fnx / (n - qr(P0)$rank - qrV$rank) #predictive variance-covariance matrices
#'
#'   #simulated value for the new y's:
#'   rpred <- mvtnorm::rmvt(n = howmany, delta = means,
#'                          sigma = diag(as.numeric(sigmas), nrow = length(sigmas)),
#'                          df = n - qr(P0)$rank - qrV$rank, type = "shifted")
#'   rpred
#' }
#'
#' #' @keywords internal
#' simpred.GD25 <- function(y, X, S, mu, n, newdataMD, howmany) {
#'   #X already centered at imputation
#'   qrX <- qr(X)
#'   R <- qr.R(qrX)
#'   iXtXplusS <- solve(t(R) %*% R + S) #A
#'   mean <- iXtXplusS %*% t(X) %*% y #betahat
#'   residuals <- (y - mean(y)) - X %*% mean
#'
#'   newX <- as.matrix(sweep(newdataMD, 2, mu, "-"))
#'   means <- mean(y) + newX %*% mean #predictive means
#'
#'   fnx <- apply(newX, 1, function(x) as.numeric(1 + t(x) %*% iXtXplusS %*% x))
#'   sigmas <- as.matrix(sum(residuals^2) * fnx / (n - qrX$rank - 1)) #predictive variance-covariance matrices
#'
#'   rpred <- mvtnorm::rmvt(n = howmany, delta = means,
#'                          sigma = diag(as.numeric(sigmas), nrow = length(sigmas)),
#'                          df = n - qrX$rank - 1, type = "shifted")
#'   rpred
#' }
