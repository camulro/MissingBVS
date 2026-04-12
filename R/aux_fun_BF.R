#' Logarithm of the MC approximation of the Bayes factor in linear models
#'
#' Computes the logarithm of the Bayes factor in the linear regression model
#' when missingness occurs, assuming normally distributed covariates. It does so
#' by approximating the integral via a MonteCarlo scheme with an
#' \code{MC.imputation} object and using the BF computation auxiliary function
#' given by \code{BF.miss.aux}.
#'
#' @param model Vector of indexes in {1,2,...,p} denoting the active variables
#' for a given model with p competing covariates.
#' @param imputation.list Object of class \code{MC.imputation} with the
#' following elements: \code{rX.imput}Array of dimension \code{n}xpx\code{nMC}
#' containing the imputed datasets \code{rSigma}Array of dimension
#' pxpx\code{nMC} containing the corresponding covariance matrices
#' @param BF.miss.aux Auxiliary function to compute the Bayes factor for each
#' entry in \code{imputation.list}.
#' @param n Number of observations.
#' @param nMC Number of samples used to approximate, by MonteCarlo, the integral
#' defining the Bayes factor.
#'
#' @return \code{lBF.miss} returns, in logarithmic scale, the MonteCarlo
#' approximation of the integral defining the Bayes factor for a given \code{model}
#' following García-Donato et al (2025).
#'
#' @author María Eugenia Castellanos and Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MC.imputation}} for computing the
#' \code{MC.imputation} object used in the MonteCarlo approximation and
#' \code{\link[MissingBVS]{BF.miss.X}} for the Bayes factor for each MC step.
#' Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20).
#'
#' @examples
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' Xair = airquality[,c("Ozone", "Wind", "Temp")]
#' imp1 <- MC.imputation(X = Xair, nMC = 1)
#'
#' lmnull <- lm(Solar.R ~ 1, data = airquality, x = T, y = T)
#' BF.fun <- function(X.center, Sigma11, k) BF.miss.X(X.center, Sigma11,
#'   y = lmnull$y, SS0 = crossprod(lmnull$residuals))
#' lBF_f <- lBF.miss(1:3, imp1, BF.miss.aux = BF.fun)
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
lBF.miss <- function(model, imputation.list, BF.miss.aux,
                     n = dim(imputation.list$rX.imput)[1], nMC = dim(imputation.list$rSigma)[3]) {
  k <- length(model)

  imputation.list.model <- list()
  imputation.list.model$rX.imput <- array(imputation.list$rX.imput[,model,],
                                          dim = c(n, k, dim(imputation.list$rSigma)[3]))
  imputation.list.model$rSigma <- array(imputation.list$rSigma[model, model,],
                                        dim = c(k, k, dim(imputation.list$rSigma)[3]))

  lBF.our <- numeric(nMC)
  for(s in 1:nMC) {
    #posterior dist. with Jeffreys independent prior
    Sigma11 <- imputation.list.model$rSigma[,,s]
    # mu <- imputation.list$rmu[,s]
    X.center <- imputation.list.model$rX.imput[,,s] #centered at imputation step

    lBF.our[s] <- BF.miss.aux(X.center, Sigma11, k) #log(BFmodel0) for the sth imputation
  }

  maxlBF.our <- max(lBF.our)
  lBF.our <- lBF.our - maxlBF.our #to avoid infite Bayes factors at the average step

  lBF.miss <- maxlBF.our + log(sum(exp(lBF.our))) - log(nMC) #log(mean(exp(lBF.our)))
  return(lBF.miss)
}

#' Logarithm of the Bayes factor in the MC step
#'
#' Computes the logarithm of the Bayes factor in the linear regression model
#' with complete data, for normally distributed covariates. It employs the imputed
#' data and covariance matrices of an \code{MC.imputation} object.
#'
#' @param X.center Matrix of dimension \code{n}xk containing the imputed data.
#' @param Sigma11 Matrix of dimension kxk containing the corresponding
#' covariance matrix.
#' @param y Response variable in the linear model.
#' @param SS0 Sum of squared error of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#'
#' @return \code{BF.miss.X} returns, in logarithmic scale, the Bayes factor in
#' García-Donato et al (2025) for imputed data given by \code{X.center} and
#' covariance matrix \code{Sigma11}.
#'
#' @author María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{MC.imputation}} for computing the
#' \code{MC.imputation} object containing the \code{X.center} and \code{Sigma11}
#' matrices. Use \code{\link[MissingBVS]{MissingBvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20).
#'
#' @examples
#' #Daily air quality measurements in New York
#' data("airquality")
#'
#' Xair = airquality[,c("Ozone", "Wind", "Temp")]
#' imp1 <- MC.imputation(X = Xair, nMC = 1)
#'
#' lmnull <- lm(Solar.R ~ 1, data = airquality, x = T, y = T)
#' lBF.imp1 <- BF.miss.X(imp1$rX.imput, imp1$rSigma,
#'   y = lmnull$y, SS0 = crossprod(lmnull$residuals))
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
BF.miss.X <- function(X.center, Sigma11, y, SS0, n = length(y), k = ncol(X.center)) {
  tX.center.X.center <- crossprod(X.center)

  lBFi0 <- -.5*(n-1)*log(1-t(y) %*% X.center %*% solve((tX.center.X.center + Sigma11)) %*% t(X.center) %*% y/SS0) -
                .5*determinant(tX.center.X.center %*% solve(Sigma11) + diag(rep(1,k)), log=T)$modulus[1]

  return(lBFi0)
}
#' Logarithm of the Bayes factor with missing data
#'
#' Computes the logarithm of the Bayes factor when missingness occurs, with any
#' type of covariates. It does so by averaging over the \code{n.imp} imputed
#' datasets given by an \code{MissingBVS.imputation} object and using the BF
#' computation auxiliary function given by \code{BF.approx.method}.
#'
#' @param model Vector of indexes in {1,2,...,p} denoting the active variables
#' for a given model with p competing covariates.
#' @param imputation.array Array of dimension \code{n}x(p+p0)x\code{n.imp}
#' containing the imputed datasets, where n is the number of observations.
#' @param BF.approx.method Auxiliary function to compute the Bayes factor for
#' each entry in \code{imputation.list}.
#' @param p0 Number of fixed covariates (including the intercept).
#' @param n.imp Number of imputed datasets by \code{mice.imputation} or
#' \code{futuremice.imputation} functions.
#'
#' @return \code{lBF.approx} returns, in logarithmic scale, the average Bayes
#' factor over the \code{n.imp} imputed datasets.
#'
#' @author María Eugenia Castellanos and Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{mice.imputation}} or
#' \code{\link[MissingBVS]{futuremice.imputation}} for computing the
#' MissingBVS.imputation object used in the average. Use
#' \code{\link[MissingBVS]{MissingBvs.lm}} with linear models or
#' \code{\link[MissingBVS]{MissingBvs.glm}} with generalized linear models for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#' XS97 = dataS97[,c("lifee060", "gdpsh60l", "p60")]
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(X = XS97, formula = f, n.imp = 1)
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, x = T, y = T)
#' BF.fun <- function(X, k) BF.approx.BIC.lm(y = lmnull$y, X,
#'   SS0 = crossprod(lmnull$residuals), k = k)
#' lBF_f <- lBF.approx(1:3, imp1, BF.approx.method = BF.fun)
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
lBF.approx <- function(model, imputation.array, BF.approx.method,
                       p0 = 1, n.imp = dim(imputation.array)[3]) {
  k <- length(model)

  lBF.aux <- numeric(n.imp)
  X1.array <- imputation.array[,c(1:p0, model+p0),] #first p0 columns are fixed
  for(s in 1:n.imp) {
    lBF.aux[s] <- BF.approx.method(k = k, X = X1.array[,,s])
  }

  maxlBF.aux <- max(lBF.aux)
  lBF.aux <- lBF.aux - maxlBF.aux #to avoid infite Bayes factors at the average step

  lBF.approx <- maxlBF.aux + log(sum(exp(lBF.aux))) - log(n.imp) #log(mean(exp(lBF.aux)))

  return(lBF.approx)
}

#' Logarithm of the BIC approximation of the Bayes factor in lm
#'
#' Computes the logarithm of the BIC approximation of Bayes factors for
#' complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#' @return \code{BF.approx.BIC.lm} returns, in logarithmic scale, the BIC
#' approximation of the Bayes factor for a given model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#' XS97 = dataS97[,c("lifee060", "gdpsh60l", "p60")]
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(X = XS97, formula = f, n.imp = 1)
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, x = T, y = T)
#' lBF.imp1 <- BF.approx.BIC.lm(y = lmnull$y, X = imp1[,,1],
#'   SS0 = crossprod(lmnull$residuals))
#'
#' @references Schwarz, G. (1978) Estimating the dimension of a model. The
#' Annals of Statistics. 6(2): 461–464.
#'
BF.approx.BIC.lm <- function(y, X, SS0,
                             n = length(y), k = ncol(X)-p0, p0 = 1L) {
  SSE.model <- crossprod(.lm.fit(y = y, x = X)$residuals)

  # BFi0 <- (SS0/SSE.model)^(n/2) / n^(k/2)
  lBFi0 <- n/2*log(SS0/SSE.model) - k/2*log(n) #exp((BIC0 - BICi)/2)
  return(lBFi0)
}

#' Logarithm of the test-based Bayes factor (TBF) in lm
#'
#' Computes the logarithm of the TBF approximation of Bayes factors for
#' complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param prior.betas Prior distribution for model-specific coefficients in the
#' \pkg{BayesVarSel} codification. Options include "gBF", "RobustBF", "LiangBF",
#' "ZSBF", "flsBF", "intrinsicBF" and "geointrinsicBF". See
#' \code{\link[BayesVarSel]{Bvs}} for more details.
#' @param SS0 Sum of squared error of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#'
#' @return \code{BF.approx.TBF.lm} returns, in logarithmic scale, the TBF
#' approximation of the Bayes factor in linear models for a given model
#' through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#' XS97 = dataS97[,c("lifee060", "gdpsh60l", "p60")]
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(X = XS97, formula = f, n.imp = 1)
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, x = T, y = T)
#' lBF.imp1 <- BF.approx.TBF.lm(y = lmnull$y, X = imp1[,,1],
#'   SS0 = crossprod(lmnull$residuals))
#'
#' @references Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#'
BF.approx.TBF.lm <- function(y, X, SS0, prior.betas = "gBF",
                             n = length(y), k = ncol(X)-p0, p0 = 1L) {
  #fixed g for now
  g <- n

  R2j <- 1 - crossprod(.lm.fit(y = y, x = X)$residuals)/SS0 #for LM
  zj <- -n * log(1 - R2j)
  # BFi0 <- (g + 1)^(-k/2) * exp(g/(g+1) * zj/2)
  lBFi0 <- -k/2 * log(g + 1) + g/(g+1) * zj/2
  return(lBFi0)
}

#' Logarithm of the g-prior Bayes factor in lm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' for complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param prior.betas Prior distribution for model-specific coefficients in the
#' \pkg{BayesVarSel} codification. Options include "gBF", "RobustBF", "LiangBF",
#' "ZSBF", "flsBF", "intrinsicBF" and "geointrinsicBF". See
#' \code{\link[BayesVarSel]{Bvs}} for more details.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#' @return \code{BF.approx.gprior.lm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in linear models for a given model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#' XS97 = dataS97[,c("lifee060", "gdpsh60l", "p60")]
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(X = XS97, formula = f, n.imp = 1)
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, x = T, y = T)
#' lBF.imp1 <- BF.approx.gprior.lm(y = lmnull$y, X = imp1[,,1],
#'   SS0 = crossprod(lmnull$residuals))
#'
#' @references García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
#' Berger, J., Garcıa-Donato, G., Moreno, E., and Pericchi, L. (2022).
#' The intrinsic hyper-g prior for normal linear models. in preparation.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger,J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423
#'
#' Moreno, E., Giron, J. and Casella, G. (2015) Posterior model consistency
#' in variable selection as the model dimension grows. Statistical Science. 30: 228-241.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369> Posterior Odds Ratio
#' for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo,
#' M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia:
#' University Press.
#'
#' Zellner, A. and Siow, A. (1984). Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#'
BF.approx.gprior.lm <- function(y, X, SS0, prior.betas = "RobustBF",
                                n = length(y), k = as.integer(ncol(X)-p0), p0 = 1L) {
  SSE.model <- crossprod(.lm.fit(y = y, x = X)$residuals)

  BFi0 <- .C(prior.betas, n, k + p0, p0, as.double(SSE.model/SS0), 0.0,
             PACKAGE = "BayesVarSel")[5][[1]]
  return(log(BFi0))
}
#' Logarithm of the FLS Bayes factor in lm
#'
#' Computes the logarithm of the Bayes factors derived from the
#' Fernandez, Ley and Steel (2001) Benchmark prior for complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param dmax Maximum model dimension of potential regressors.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#' @return \code{BF.approx.FLS.lm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning the FLS Benchmark g-prior
#' in linear models for a given model through \code{X}. See
#' \code{\link[BayesVarSel]{Bvs}} for more details.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#' XS97 = dataS97[,c("lifee060", "gdpsh60l", "p60")]
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(X = XS97, formula = f, n.imp = 1)
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, x = T, y = T)
#' lBF.imp1 <- BF.approx.FLS.lm(y = lmnull$y, X = imp1[,,1],
#'   SS0 = crossprod(lmnull$residuals), dmax = ncol(XS97))
#'
#' @references García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
BF.approx.FLS.lm <- function(y, X, SS0, dmax,
                             n = length(y), k = as.integer(ncol(X)-p0), p0 = 1L) {
  SSE.model <- crossprod(.lm.fit(y = y, x = X)$residuals)

  BFi0 <- .C("flsBF", dmax - p0, n, k + p0, p0, as.double(SSE.model/SS0), 0.0,
             PACKAGE = "BayesVarSel")[6][[1]]
  return(log(BFi0))
}

c_glm.fit <- utils::getFromNamespace("C_glm_deterministic", "BAS") #to compute logmarginals in glm

#' Logarithm of the BIC approximation of the Bayes factor in glm
#'
#' Computes the logarithm of the BIC approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. For now, only available the implemented
#' families in \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#' @param logmargnull Log-marginal likelihood of the null model considered.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
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
#' @return \code{BF.approx.BIC.glm} returns, in logarithmic scale, the BIC
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.glm}} for
#' an exact computation of the model posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Indian Prime Diabetes Data from VIM's package
#'
#' Xdiab = VIM::diabetes[,c("Pregnancies", "Glucose", "Insulin")]
#' f <- Outcome ~ Pregnancies + Glucose + Insulin
#' imp1 <- mice.imputation(X = Xdiab, formula = f, n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = VIM::diabetes, family = binomial(),
#'   x = T, y = T)
#' lBF.imp1 <- BF.approx.BIC.glm(y = glmnull$y, X = imp1[,,1],
#'   family = binomial(), logmargnull = 0) #returns the logmarginal approx
#'
#' @references Schwarz, G. (1978) Estimating the dimension of a model. The
#' Annals of Statistics. 6: 461–464.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
BF.approx.BIC.glm <- function(y, X, family = binomial(link = "logit"),
                              logmargnull,
                              k = ncol(X)-p0, p0 = 1L,
                              #glm.fit arguments:
                              weights = rep(1, length(y)),
                              offset = rep(0, length(y)),
                              control = glm.control(),
                              laplace = 0L) {
  initprob <- c(rep(1.0, p0), rep(.5, k)) #first p0 columns of X are the fixed covariates
  fit1 <- .Call(c_glm.fit, Y = y, X = X, Roffset = offset, Rweights = weights,
                Rprobinit = initprob, Rmodeldim = 0L, modelprior = BAS::beta.binomial(1, 1),
                betaprior = BAS::bic.prior(length(y)), family = family,
                Rcontrol = control, Rlaplace = laplace)

  # BFi0 <- exp(fit1$logmarg - logmargnull)
  lBFi0 <- fit1$logmarg - logmargnull
  return(lBFi0)
}
#' Logarithm of the BIC approximation of the Bayes factor in glm
#'
#' Computes the logarithm of the BIC approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param devnull Deviance of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#'
#' @return \code{BF.approx.BIC.glm.stats} returns, in logarithmic scale, the BIC
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.glm}} for
#' an exact computation of the model posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Indian Prime Diabetes Data from VIM's package
#'
#' Xdiab = VIM::diabetes[,c("Pregnancies", "Glucose", "Insulin")]
#' f <- Outcome ~ Pregnancies + Glucose + Insulin
#' imp1 <- mice.imputation(X = Xdiab, formula = f, n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = VIM::diabetes, family = binomial(),
#'   x = T, y = T)
#' lBF.imp1 <- BF.approx.BIC.glm.stats(y = glmnull$y, X = imp1[,,1],
#'   family = binomial(), devnull = glmnull$deviance)
#'
#' @references Schwarz, G. (1978) Estimating the dimension of a model. The
#' Annals of Statistics. 6(2): 461–464.
#'
BF.approx.BIC.glm.stats <- function(y, X, family = binomial(link = "logit"),
                                    devnull,
                                    n = length(y), k = ncol(X)-p0, p0 = 1L,
                                    #glm.fit arguments
                                    weights = rep(1, n),
                                    offset = rep(0, n),
                                    control = glm.control()) {
  fit1 <- glm.fit(y = y, x = X, family = family,
                  weights = weights, offset = offset, control = control)

  lBFi0 <- (devnull - fit1$deviance - k * log(n))/2
  return(lBFi0)
}

#' Logarithm of the test-based Bayes factor (TBF) in glm
#'
#' Computes the logarithm of the TBF approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. For now, only available the implemented
#' families in \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#' @param devnull Deviance of the null model considered.
#' @param prior.betas Prior distribution for model-specific coefficients.
#' Options include \code{\link[BAS]{g.prior}}, \code{\link[BAS]{CCH}},
#' \code{\link[BAS]{robust}} and \code{\link[BAS]{intrinsic}} among others.
#' See \code{\link[BAS]{BAS}} for more details.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
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
#' @return \code{BF.approx.TBF.glm} returns, in logarithmic scale, the TBF
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Indian Prime Diabetes Data from VIM's package
#'
#' Xdiab = VIM::diabetes[,c("Pregnancies", "Glucose", "Insulin")]
#' f <- Outcome ~ Pregnancies + Glucose + Insulin
#' imp1 <- mice.imputation(X = Xdiab, formula = f, n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = VIM::diabetes, family = binomial(),
#'   x = T, y = T)
#' lBF.imp1 <- BF.approx.TBF.glm(y = glmnull$y, X = imp1[,,1],
#'   family = binomial(), devnull = glmnull$deviance)
#'
#' @references Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
BF.approx.TBF.glm <- function(y, X, family = binomial(link = "logit"),
                              devnull, prior.betas = BAS::g.prior(g = length(y)),
                              k = ncol(X)-p0, p0 = 1L,
                              #glm.fit arguments
                              weights = rep(1, length(y)),
                              offset = rep(0, length(y)),
                              control = glm.control(),
                              laplace = 0L) {
  initprob <- c(rep(1.0, p0), rep(.5, k)) #first p0 columns of X are the fixed covariates
  fit1 <- .Call(c_glm.fit, Y = y, X = X, Roffset = offset, Rweights = weights,
                Rprobinit = initprob, Rmodeldim = 0L, modelprior = BAS::beta.binomial(1, 1),
                betaprior = prior.betas, family = family, Rcontrol = control, Rlaplace = laplace)
  #fixed g for now
  g <- length(y)
  # BFi0 <- (g + 1)^(-k/2) * exp(g/(g+1) * (devnull - fit1$deviance)/2)
  lBFi0 <- -k/2 * log(g + 1) + g/(g+1) * (devnull - fit1$deviance)/2
  return(lBFi0)
}
#' Logarithm of the test-based Bayes factor (TBF) in glm
#'
#' Computes the logarithm of the TBF approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param devnull Deviance of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#'
#' @return \code{BF.approx.TBF.glm.stats} returns, in logarithmic scale, the TBF
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Indian Prime Diabetes Data from VIM's package
#'
#' Xdiab = VIM::diabetes[,c("Pregnancies", "Glucose", "Insulin")]
#' f <- Outcome ~ Pregnancies + Glucose + Insulin
#' imp1 <- mice.imputation(X = Xdiab, formula = f, n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = VIM::diabetes, family = binomial(),
#'   x = T, y = T)
#' lBF.imp1 <- BF.approx.TBF.glm.stats(y = glmnull$y, X = imp1[,,1],
#'   family = binomial(), devnull = glmnull$deviance)
#'
#' @references Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
BF.approx.TBF.glm.stats <- function(y, X,
                                    family = binomial(link = "logit"),
                                    devnull,
                                    n = length(y), k = ncol(X)-p0, p0 = 1,
                                    #glm.fit arguments
                                    weights = rep(1, n),
                                    offset = rep(0, n),
                                    control = glm.control()) {
  fit1 <- glm.fit(y = y, x = X, family = family,
                  weights = weights, offset = offset, control = control)
  #fixed g for now
  g <- n
  # BFi0 <- (g + 1)^(-k/2) * exp(g/(g+1) * (devnull - fit1$deviance)/2)
  lBFi0 <- -k/2 * log(g + 1) + g/(g+1) * (devnull - fit1$deviance)/2
  return(lBFi0)
}

#' Logarithm of the g-prior Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' for complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. For now, only available the implemented
#' families in \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#' @param logmargnull Log-marginal likelihood of the null model considered.
#' @param prior.betas Prior distribution for model-specific coefficients.
#' Options include \code{\link[BAS]{g.prior}}, \code{\link[BAS]{CCH}},
#' \code{\link[BAS]{robust}} and \code{\link[BAS]{intrinsic}} among others.
#' See \code{\link[BAS]{BAS}} for more details.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
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
#' @return \code{BF.approx.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[MissingBVS]{lBF.approx}} to compute the average Bayes
#' factor for missing data. Use \code{\link[MissingBVS]{MissingBvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examples
#' #Indian Prime Diabetes Data from VIM's package
#'
#' Xdiab = VIM::diabetes[,c("Pregnancies", "Glucose", "Insulin")]
#' f <- Outcome ~ Pregnancies + Glucose + Insulin
#' imp1 <- mice.imputation(X = Xdiab, formula = f, n.imp = 1)
#'
#' glmnull <- glm(Outcome ~ 1, data = VIM::diabetes, family = binomial(),
#'   x = T, y = T)
#' lBF.imp1 <- BF.approx.gprior.glm(y = glmnull$y, X = imp1[,,1],
#'   family = binomial(), logmargnull = 0) #returns the logmarginal
#'
#' @references Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#'
BF.approx.gprior.glm <- function(y, X, family = binomial(link = "logit"),
                                 logmargnull, prior.betas = BAS::robust(n = length(y)),
                                 k = ncol(X)-p0, p0 = 1L,
                                 #glm.fit arguments
                                 weights = rep(1, length(y)),
                                 offset = rep(0, length(y)),
                                 control = glm.control(),
                                 laplace = 0L) {
  initprob <- c(rep(1.0, p0), rep(.5, k)) #first p0 columns of X are the fixed covariates
  fit1 <- .Call(c_glm.fit, Y = y, X = X, Roffset = offset, Rweights = weights,
                Rprobinit = initprob, Rmodeldim = 0L, modelprior = BAS::beta.binomial(1, 1),
                betaprior = prior.betas, family = family, Rcontrol = control, Rlaplace = laplace)

  # BFi0 <- exp(fit1$logmarg - logmargnull)
  lBFi0 <- fit1$logmarg - logmargnull
  return(lBFi0)
}
