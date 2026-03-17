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

