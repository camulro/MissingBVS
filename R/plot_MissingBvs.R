#' Summary plot of an object of class \code{MissingBvs}
#'
#' Plots to summarize graphically the results in an object of class
#' \code{MissingBvs}, including prior and posterior probabilities obtained. 
#' Two summaries are available: the model size probability if \code{plotdim}=TRUE 
#' and the inclusion probabilities if  \code{plotpip}=TRUE. By default, 
#' both plots are shown. 
#' 
#' @export
#' @param mbvs.object Object of class \code{MissingBvs}.
#' @param plotdim Logical indicating whether or not to plot model dimension
#' probabilities.
#' @param plotpip Logical indicating whether or not to plot posterior inclusion
#' probabilities.
#' @param ... Additional graphical parameters.
#' 
#' @author Carolina Mulet
#' Maintainer: <carolina.mulet1@@alu.uclm.es>
#' 
#' @seealso See \code{\link[MissingBVS]{MissingBvs.lm}}, 
#' \code{\link[MissingBVS]{MissingBvs.glm}}, 
#' \code{\link[MissingBVS]{MissingGibbsBvs.lm}} and
#' \code{\link[MissingBVS]{MissingGibbsBvs.glm}} for creating objects of class
#' \code{MissingBvs}.
#' 
#' @examples
#' # To be completed

plot.MissingBvs <- function(mbvs.object, plotdim = TRUE, plotpip = TRUE,...) {
    
  if (!inherits(mbvs.object, "MissingBvs")){
    warning("An object of class MissingBvs is needed.\n")
  }
  
  p <- mbvs.object$p #number of covariates
  k <- mbvs.object$k #number of fixed vars
  vars <- mbvs.object$variables
  
  par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(sum(plot), 1))
  if (plot[1]) {
    
    #dimension probabilities
    if (mbvs.object$method == "gibbs") {
      main.1 <- "Estimated Model Size Probabilities"
    } else main.1 <- "Model Size Probabilities"
    
    barplot(
      mbvs.object$priorprobs,
      main = main.1, names.arg = (0:p) + k,
      xlab = "Number of covariates in the true model",
      ylab = "Probability",
      col = rgb(0.8,0.8,0.8,0.8),
      border = NA, ylim = c(0, 1), las = 1,
      cex.names =0.8, cex.axis = 0.8, cex.lab = 0.8
    )
    par(new = TRUE)
    barplot(
      mbvs.object$postprobdim,
      col = rgb(0.3,0.3,0.3,0.8),
      border = NA,
      ylim = c(0, 1),        
      cex.names = 0.8,
      axes = FALSE
    )
    abline(h = seq(0,1,0.2), col = "gray90", lty = "dotted")
    
    legend(
      "topright",
      legend = c("Posterior", "Prior"),
      fill = c(rgb(0.3,0.3,0.3,0.8), rgb(0.8,0.8,0.8,0.8)),
      border = NA,
      bty = "n"
    )
  }

  if (plot[2]) {
    
    #posterior inclusion probabilities
    if (mbvs.object$method == "gibbs") {
      main.2 <- "Estimated Posterior Inclusion Probabilities"
    } else main.2 <- "Posterior Inclusion Probabilities"
    
    barplot(
      mbvs.object$inclprob,
      main = main.2, names.arg = vars,
      col = rgb(0.3,0.3,0.3,0.8),
      border = NA,
      ylim = c(0, 1),
      las = 2,
      cex.names = 0.8, cex.axis = 0.8, cex.lab = 0.8,
      ylab = "Probability",
      ...
    )
    abline(h = 0.5, col = rgb(0.8,0.8,0.8,0.8), lty = 2, lwd = 2)
    abline(h = seq(0,1,0.2), col = "gray90", lty = "dotted")
  }
}