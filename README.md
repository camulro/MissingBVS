### `MissingBVS`: Bayesian Variable Selection with missing data using R

`MissingBVS` is an R package currently on development. It is designed to be highly intuitive, 
similar to the well-known lm and glm functions, in order to appeal a wide range of users from 
different disciplines. The package provides several functions to deal with the dual problem of 
model uncertainty and the presence of missing data by the proper Bayes factor and posterior 
distribution computation.

The process can be done through exact algorithms to perform fast computations in problems of 
small to moderate size and heuristic sampling methods to solve large problems.

## Installing steps

To install the latest development version, you can use `install_github` from the `devtools` package:

```R
## install devtools if necessary for install_github()
install.packages('devtools')
library(devtools)
## get the package from github
install_github('camulro/MissingBVS', dependencies = TRUE)
```
## References

-   Barbieri, M and Berger, J. (2004).  Optimal Predictive Model Selection. 
    *The Annals of Statistics, 32*: 870-897. DOI:
    [10.1214/009053604000000238](https://doi.org/10.1214/009053604000000238)

-   Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G. (2012).
    Criteria for Bayesian Model choice with Application to Variable
    Selection. *Annals of Statistics, 40*: 1550-1577. DOI:
    [10.1214/12-aos1013](http://www.dx.doi.org/10.1214/12-aos1013)
    
-   van Buuren, S. and Groothuis-Oudshoorn, K. (2011). mice:
    Multivariate Imputation by Chained Equations in R. 
    *Journal of Statistical Software, 45*(3): 1–67. DOI:
    [10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)
    
-   Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
    Bayesian Adaptive Sampling. *R package version 2.0.2*
    <https://CRAN.R-project.org/package=BAS>.
    
-   Fernandez, C., Ley, E. and Steel, M.F.J. (2001). Benchmark priors
    for Bayesian model averaging. *Journal of Econometrics, 100*:
    381-427. DOI:
    [10.1016/s0304-4076(00)00076-2](http://www.dx.doi.org/10.1016/s0304-4076(00)00076-2)
    
-   García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
    and Forte, A. (2025). Model Uncertainty and Missing Data: An 
    Objective Bayesian Perspective (with Discussion). *Bayesian Analysis, 20*: 
    1677–1778. DOI:
    [10.1214/25-ba1531](https://doi.org/10.1214/25-BA1531)
    
-   Garcia-Donato, G. and Forte, A. (2018). Bayesian Testing, Variable
    Selection and Model Averaging in Linear Models using R with
    BayesVarSel. *The R Journal, 10*(1): 155–174. Retrieved from
    <https://journal.r-project.org/archive/2018/RJ-2018-021/index.html>
    
-   Garcia-Donato, G. and Martinez-Beneito, M.A.(2013). On sampling
    strategies in Bayesian variable selection problems with large model
    spaces. *Journal of the American Statistical Association, 108*:
    340-352. DOI:
    [10.1080/01621459.2012.742443](http://www.dx.doi.org/10.1080/01621459.2012.742443)
    
-   Held, L., Sabanés Bové, D. and Gravestock, I. (2015). Approximate 
    Bayesian Model Selection with the Deviance Statistic. 
    *Statistical Science, 30*(2): 242-257. DOI:
    [10.1214/14-sts510](https://doi.org/10.1214/14-STS510)
    
-   Li, Y. and Clyde, M. (2018). Mixtures of g-priors in Generalized Linear 
    Models. *Journal of the American Statistical Association, 113*: 1828-1845.
    DOI: [10.1080/01621459.2018.1469992](https://doi.org/10.1080/01621459.2018.1469992)
    
-   Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008).
    Mixtures of g-priors for Bayesian Variable Selection. 
    *Journal of the American Statistical Association, 103*: 410-423. DOI:
    [10.1198/016214507000001337](http://www.dx.doi.org/10.1198/016214507000001337)
    
-   Schwarz, G. (1978). Estimating the dimension of a model. 
    *The Annals of Statistics, 6*(2): 461–464. DOI:
    [10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)
    
-   Scott, J.G. and Berger, J.O. (2010). Bayes and empirical-Bayes multiplicity
    adjustment in the variable-selection problem. 
    *The Annals of Statistics, 38*: 2587–2619. DOI:
    [10.1214/10-AOS792](https://doi.org/10.1214/10-AOS792)
    
-   Zellner, A. and Siow, A. (1980). Posterior Odds Ratio for Selected
    Regression Hypotheses. 
    *Trabajos de Estadística y de Investigación Operativa, 31*: 585. DOI:
    [10.1007/bf02888369](http://www.dx.doi.org/10.1007/bf02888369)
    
-   Zellner, A. and Siow, A. (1984). *Basic Issues in Econometrics*.
    Chicago: University of Chicago.
    
-   Zellner, A. (1986). On Assessing Prior Distributions and Bayesian
    Regression Analysis with g-prior Distributions. In A. Zellner (ed.),
    *Bayesian Inference and Decision techniques: Essays in Honor of Bruno de Finetti*, 
    389-399. Edward Elgar Publishing Limited. DOI:
    [10.2307/2233941](http://www.dx.doi.org/10.2307/2233941)
