#read example data
data(dataS97)

#example with a few competing covariates
dataS97.MissingBvs <- missingBVS.lm(formula = "gr56092~gdpsh60l+lifee060+p60", 
                                    data = dataS97, n.keep = 2^3, imp.seed = 123)

#print the result
dataS97.MissingBvs

#summary of the result
summary(dataS97.MissingBvs)

#plots for the result
plot(dataS97.MissingBvs)
