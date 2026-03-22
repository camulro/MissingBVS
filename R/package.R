
#' Model Choice And Variable Selection With Bayes Factors With Missing Data
#'
#' Hypothesis testing, model selection and model averaging
#'
#' \tabular{ll}{Package: \tab MissingBVS\cr Type: \tab Package\cr Version:
#' \tab 0.1.0\cr Date: \tab 2026-03-20\cr License: \tab GPL (>= 3)\cr }
#'
#' @name MissingBVS-package
#' @aliases MissingBVS-package MissingBVS
#' @author Carolina Mulet, Gonzalo Garcia-Donato and María Eugenia Castellanos
#'
#' Maintainer: Carolina Mulet \email{Carolina.Mulet1@@alu.uclm.es}
#' @seealso \code{\link[MissingBVS]{missingBtest.lm}},
#' \code{\link[MissingBVS]{missingBvs.lm}},
#' \code{\link[MissingBVSl]{missingGibbsBvs.lm}}
#' @references
#'
#' Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238> Optimal
#' Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1577
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association. 108: 340-352.
#'
#' Held, L., Gravestock, I. and Sabanés Bové, D.
#' (2015)<DOI:10.1080/01621459.2014.993077> Objective Bayesian model selection
#' for generalized linear models using test-based Bayes factors. Journal of the
#' American Statistical Association, 110, 1157–1168.
#'
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures
#' of g-Priors in Generalized Linear Models. Journal of the American
#' Statistical Association. 113: 1275–1287.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6: 461–464.
#'
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369>. Posterior Odds
#' Ratio for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M.
#' Bernardo, M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603.
#' Valencia: University Press.
#'
#' Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' @keywords package
#' @examples
#' # To be completed
#'
"_PACKAGE"


#' Sala-i-Martin 97 data
#'
#' Data set contains 63 variables potentially related with economic Growth (GDP).
#' The dataset is in its original shape with many variables containing missing
#' observations coded as NA. The number of countries gathered is n=140.
#'
#' It was firstly used in Sala-i-Martin (1997) and later analized by Fernández,
#' Ley and Steel (2001) through a BVS perspective deleting missing obeservations.
#' Taken from the Journal of Data Archive - Journal of Applied Econometrics.
#'
#'
#' @name dataS97
#' @docType data
#' @format A data frame with 140 observations and 63 variables
#' \describe{ \item{\code{gr56092}}{Growth rate of GDP per capita between 1960
#' and 1992. Response variable.}
#' \item{\code{Mining}}{Fraction of GDP in mining.}
#' \item{\code{ECORG}}{Degree Capitalism index. Coded with 0,...,5 as follows:
#' 0=statis; 1=mixed statist; 2=mixed capitalist-statist; 3=capitalist-statist;
#' 4=mixed-capitalist; 5=capitalist.}
#' \item{\code{YrsOpen}}{Number of years economy has been open between 1950 and
#' 1990.}
#' \item{\code{EngFrac}}{Fraction of population speaking English.}
#' \item{\code{OthFrac}}{Fraction of population speaking foreign language.}
#' \item{\code{RERD}}{Real exchange rate distortions.}
#' \item{\code{EQINV}}{Equipment Investment.}
#' \item{\code{NONEQINV}}{Non-Equipment Investment.}
#' \item{\code{BMS6087}}{Standard deviation of black market premium.}
#' \item{\code{GDC6089}}{Growth rate of Domestic Credit between 1960 and 1990.}
#' \item{\code{PI6089}}{Average inflation rate between 1960 and 1990.}
#' \item{\code{SCOUT}}{Measure of outward orientation.}
#' \item{\code{STDC6089}}{Standard deviation of Gross Domestic Credit between
#' 1960 and 1989.}
#' \item{\code{STPI6089}}{Standard deviation of the inflation rate between 1960
#' and 1990.}
#' \item{\code{bmp1}}{Log of (1 + Black market premium).}
#' \item{\code{area}}{(Scale effect) Total area of the country.}
#' \item{\code{freeop}}{Free trade oppenness, measure of free trade.}
#' \item{\code{freetar}}{Degree of tariff barriers.}
#' \item{\code{laam}}{Dummy for Latin American countries.}
#' \item{\code{safrica}}{Dummy for Sub-Saharan African countries.}
#' \item{\code{pyr60}}{Average years of primary schooling of total population
#' in 1960.}
#' \item{\code{syr60}}{Average years of secondary schooling of total population
#' in 1960.}
#' \item{\code{hyr60}}{Average years of higher education of total population
#' in 1960.}
#' \item{\code{human60}}{Average years of education of total population in 1960.}
#' \item{\code{p60}}{Enrollment rate in primary education in 1960.}
#' \item{\code{s60}}{Enrollment rate in secondary education in 1960.}
#' \item{\code{h60}}{Enrollment rates in higher education in  1960.}
#' \item{\code{lifee060}}{Life expectancy in 1960.}
#' \item{\code{ggcfc3}}{Average share of expenditures on public investment as
#' fraction of GDP between 1960 and 1965.}
#' \item{\code{geerec1}}{Average share public expenditures on education as
#' fraction of GDP between 1960 and 1965.}
#' \item{\code{gde1}}{Average share public expenditures on defense as fraction
#' of GDP between 1960 and 1965.}
#' \item{\code{assassp2}}{Number of political assassinations.}
#' \item{\code{revcoup}}{Number of revolutions and military coups.}
#' \item{\code{pinstab2}}{Political instability.}
#' \item{\code{wardum}}{Dummy for contries that have been involved in war any
#' time between 1960 and 1990.}
#' \item{\code{prightsb}}{Political rights index.}
#' \item{\code{civlibb}}{Civil liberties index in 1972.}
#' \item{\code{tot1}}{Growth of terms of trade between 1960 and 1990.}
#' \item{\code{ABSLATIT}}{Absolute latitude.}
#' \item{\code{AGE}}{Average age of the population.}
#' \item{\code{BRIT}}{Dummy for former British colony after 1776.}
#' \item{\code{BUDDHA}}{Fraction of population Buddhist in 1960.}
#' \item{\code{CATH}}{Fraction of population Catholic in 1960.}
#' \item{\code{CONFUC}}{Fraction of population Confucian.}
#' \item{\code{DEMOC65}}{Qualitative index of democracy freedom (1965.}
#' \item{\code{FRAC}}{Ethnolinguistic fractionalization (probability two random
#' people in a country do not speak same language).}
#' \item{\code{FRENCH}}{Dummy for former French colonies.}
#' \item{\code{HINDU}}{Fraction of the population Hindu in 1960.}
#' \item{\code{JEW}}{Fraction of the population Jew in 1960.}
#' \item{\code{MUSLIM}}{Fraction of the population Muslim in 1960.}
#' \item{\code{PRIEXP70}}{Fraction of primary exports in total exports in 1970.}
#' \item{\code{PROT}}{Fraction of population Protestant in 1960.}
#' \item{\code{RULELAW}}{Rule of law.}
#' \item{\code{SPAIN}}{Dummy variable for former Spanish colonies.}
#' \item{\code{URB60}}{Fraction of population living in cities.}
#' \item{\code{dpop6090}}{Average growth rate of population between 1960 and 1990.}
#' \item{\code{gdpsh60l}}{Logarithm of GDP per capita in 1960.}
#' \item{\code{humanyl}}{Product of average years of schooling and log of GDP
#' per capita in 1960.}
#' \item{\code{lly1}}{Ratio of liquid liabilities to GDP (a measure of financial
#' development.)}
#' \item{\code{work60l}}{Ratio workers to population.}
#' \item{\code{lforce60}}{Size of labor force (scale effect).}
#' \item{\code{gvxdxe52}}{Public consumption minus education and defense as
#' fraction of GDP.}
#' }
#'
#' @references Sala i Martin, X. (1997) I Just Ran Two Million Regressions.
#' The American Economic Review, 87(2), 178–183.
#'
#' Fernández, C., Ley, E., and Steel., M.F.J. (2001). Model Uncertainty
#' in Cross-Country Growth Regressions. Journal of Applied Econometrics,
#' 16(5), 563–576.
#'
#' @keywords datasets
#'
#' @examples
#' data(dataS97)
#'
"dataS97"


#' OBICE data
#'
#' Dataset corresponding to the OBICE study (Zurriaga et al 2011) where factors
#' associated with childhood obesity are studied. The data were collected in
#' 2007 and 2008 through several questionnaries and n=1188 children were enroled
#' in the study. It contains 155 variables.
#'
#' This is a case and control study with 437 cases (obese) and 751 controls
#' (not obese). Purposedly the dataset is distributed without any post-processing
#' and therefore, many variables contain missing observations coded in different way.
#'
#'
#' @name OBICE
#' @docType data
#' @format A data frame with 1188 entries and 121 variables. The more relevant
#' are described. Contact us if you need specific information of any other.
#' \describe{
#' \item{Acostarse}{does he/she eat before going to bed? (yes/no) }
#' \item{ActFisica}{(physic activity) factor coded as 1 None; 2 less than monthly; 3 less than weekly; 4 less than 2/week; 5 at least 2/week }
#' \item{ActivDepor}{weekly hours devoted to sports activity }
#' \item{Almuerzo}{}
#' \item{Bebida}{(main dring accompanying the main meal) 1 water tap; 2 bottle water; 3 soda; 4 natural juices; 5 bottle juices; 6 Milk (and derivatives); 7 Other }
#' \item{Caso01}{}
#' \item{Cena}{}
#' \item{Chuches}{Sweets and soft drinks weekly consumption (how many times) }
#' \item{CincoComidas}{does he/she have regularly 5 meals per day? (0 is No; 1 is Yes) }
#' \item{clSocEl}{}
#' \item{clSocElla}{}
#' \item{clSocXiquet}{}
#' \item{ComedorEsc}{}
#' \item{Comida}{}
#' \item{Daceite}{}
#' \item{Dcereal}{}
#' \item{Desayuno}{}
#' \item{Descubrimiento}{}
#' \item{Dgalleta}{}
#' \item{Dislipemias}{}
#' \item{DislipeRelacion}{}
#' \item{Dleche}{}
#' \item{Dotros}{}
#' \item{Dpan}{}
#' \item{Dzumoenv}{}
#' \item{Dzumonat}{}
#' \item{Edad}{years old }
#' \item{EntreHoras}{}
#' \item{EstudiosMadre}{}
#' \item{EstudiosMadreSinCon}{}
#' \item{EstudiosPadre}{}
#' \item{EstudiosPadreSinCon}{}
#' \item{Faperitivos}{}
#' \item{Faperitivosmp}{}
#' \item{Farroz}{}
#' \item{Farrozmp}{}
#' \item{Fcarnes}{}
#' \item{Fcarnesmp}{}
#' \item{Fchucherias}{}
#' \item{Fchucheriasmp}{}
#' \item{Fdulces}{}
#' \item{Fdulcesmp}{}
#' \item{Ffiambres}{}
#' \item{Ffiambresmp}{}
#' \item{Ffritos}{}
#' \item{Ffritosmp}{}
#' \item{Ffruta}{}
#' \item{Ffrutamp}{}
#' \item{Fhuevos}{}
#' \item{Fhuevosmp}{}
#' \item{Flacteos}{}
#' \item{Flacteosmp}{}
#' \item{Flegumbres}{}
#' \item{Flegumbresmp}{}
#' \item{Fpan}{}
#' \item{Fpanmp}{}
#' \item{Fpescado}{}
#' \item{Fpescadomp}{}
#' \item{Fprecocina}{}
#' \item{Fprecocinamp}{}
#' \item{Frefrescos}{}
#' \item{Frefrescosmp}{}
#' \item{Fruta}{usual consumption of fruit? (0 is No; 1 is Yes) }
#' \item{FrutaVerdura}{}
#' \item{Fverduras}{}
#' \item{Fverdurasmp}{}
#' \item{HorasPantDia}{}
#' \item{HorasPCDiaPond}{daily hours playing videogames and/or in internet (weekends included) }
#' \item{HorasPCsem1}{}
#' \item{HorasPCsem2}{}
#' \item{HorasTV}{}
#' \item{HorasTVDiaPond}{daily hours watching TV (weekends included) }
#' \item{HorasTVsem1}{}
#' \item{HorasTVsem2}{}
#' \item{HoraSuenyo}{daily hours sleeping }
#' \item{HTA}{}
#' \item{HTARelacion}{}
#' \item{IMC}{}
#' \item{IndEdadComedorEscolar}{}
#' \item{IntolGlucosa}{}
#' \item{IntolRelacion}{}
#' \item{LactMater}{}
#' \item{LactMaterna}{ breast-feeding (1 is Yes; 0 is No) }
#' \item{LactMatMeses}{}
#' \item{LactMatSemanas}{}
#' \item{MadreObesa}{}
#' \item{MadreObesa01}{is the mother obese? (0 is No; 1 is Yes) }
#' \item{Merienda}{Afternoon snack  (1 is Yes; 0 is No) }
#' \item{NumComidas}{}
#' \item{NumContOK}{}
#' \item{NumControles}{}
#' \item{NumHnos}{}
#' \item{NumHnosOb}{}
#' \item{NumPadresEsp02}{}
#' \item{NumPadresObesos}{}
#' \item{OrdenadorDiario}{}
#' \item{OrdenadorFinDe}{}
#' \item{OsteoRelacion}{}
#' \item{OtrosPatol}{}
#' \item{PadreObeso}{is the father obese? (0 is No; 1 is Yes) }
#' \item{PesoActual}{current weight (in kilograms) }
#' \item{PesoNac}{weight born (in grams) }
#' \item{PorcHnosObesos}{}
#' \item{porcHnosObesosOK}{}
#' \item{Postre}{}
#' \item{ProbOsteo}{}
#' \item{ProbPsico}{}
#' \item{ProbResp}{}
#' \item{PsicoRelacion}{}
#' \item{ResoponF01}{}
#' \item{RespRelacion}{}
#' \item{Semlact}{}
#' \item{Sexo}{female (1); male (0) }
#' \item{TallaAct}{current height (in meters) }
#' \item{TallaNac}{height born (in centimeters) }
#' \item{Tipocaso}{}
#' \item{Tipocaso.y}{}
#' \item{TipoObeso}{}
#' \item{TVDiario}{}
#' \item{TVFinSemana}{}
#' \item{Verduras}{usual consumption of vegetables? (0 is No; 1 is Yes) }
#' }
#' @references Zurriaga, O., Perez-Panades, J., Quiles, J. , Gil, M.,
#' Anes, Y., Quiñones, C., Margolles, M., Lopez-Maside, A.,
#' Vega-Alonso, A., Miralles M. and Recent OBICE Research Group (2011)
#' Factors associated with childhood obesity in Spain. The OBICE
#' study: a case–control study based on sentinel networks. Public Health Nutrition
#' 14(6), 1105–1113.
#'
#' @keywords datasets
#'
#' @examples
#' data(OBICE)
#'
"OBICE"
