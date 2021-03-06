\name{thetaEst}
\alias{thetaEst}

\title{Ability estimation (dichotomous and polytomous models)}

\description{
 This command returns the ability estimate for a given response pattern and a given matrix of item parameters, either under the 4PL model or any suitable polytomous IRT model. Available estimators are maximum likelihood, Bayes modal (MAP), expected a posteriori (EAP) and weighted likelihood. 
 }

\usage{
thetaEst(it, x, model = NULL, D = 1, method = "BM", priorDist = "norm", 
 	priorPar = c(0, 1), range = c(-4, 4), parInt = c(-4, 4, 33),
 	constantPatt = NULL, current.th = 0, bRange = c(-2, 2))
 }

\arguments{
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{x}{numeric: a vector of item responses.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"}, \code{"WL"} and \code{"EAP"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}. Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}). Ignored if \code{method=="EAP"}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
\item{constantPatt}{character: the method to estimate ability in case of constant pattern (i.e. only correct or only incorrect responses). Can be either \code{NULL} (default), \code{"BM"}, \code{"EAP"}, \code{"WL"}, \code{"fixed4"}, \code{"fixed7"} or \code{"var"}. \emph{Currently only implemented for dichotomous IRT models}. See \bold{Details}.}
\item{current.th}{numeric: the current ability estimate (default is zero). Required for ability estimation in presence of constant pattern. Ignored if \code{constantPatt} is neither \code{"fixed4"}, \code{"fixed7"} nor \code{"var"}. See \bold{Details}.}
\item{bRange}{numeric: vector of two componentns with the range of difficulty parameters in the parent item bank (default is \code{c(-2,2)}). \emph{Currently only implemented for dichotomous IRT models}. See \bold{Details}.}
}

\value{
 The estimated ability level.
 }
 
\details{
Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

 Four ability estimators are available: the maximum likelihood (ML) estimator (Lord, 1980), the Bayes modal (BM) estimator (Birnbaum, 1969), the expected a posteriori (EAP) estimator (Bock and Mislevy, 1982) and the weighted likelihood (WL) estimator (Warm, 1989). The selected estimator is specified by the \code{method} argument, with values \code{"ML"}, \code{"BM"}, \code{"EAP"} and \code{"WL"} respectively.

 For the BM and EAP estimators, three prior distributions are available: the normal distribution, the uniform distribution and Jeffreys' prior distribution (Jeffreys, 1939, 1946). The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. The \code{priorPar} argument is ignored if \code{method="ML"} or \code{method="WL"}.
 
 The argument \code{priorPar} determines either the prior mean and standard deviation of the normal prior distribution (if 
 \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument
 is ignored if \code{priorDist="Jeffreys"}.

 The \code{parInt} argument sets the range and the number of quadrature points for numerical integration in the EAP process. By default, it takes the vector value (-4, 4, 33), that is, 33 quadrature points on the range [-4; 4] (or, by steps of 0.25). See \code{\link{eapEst}} for further details. 

 The \code{range} argument permits to limit the interval of investigation for the ML, BM and WL ability estimates (in particular, to avoid infinite ability estimates). The default \code{range} is [-4, 4].

Specific ability estimation methods are available in presence of constant patterns (that is with only correct or only incorrect responses) under dichotomous IRT models. These methods are specified by the argument \code{constantPatt}. By default it is set to \code{NULL} and hence ability is estimated with the specified \code{method} (even in presence of constant pattern). Six methods are currently available for constant patterns: \code{"BM"}, \code{"EAP"} and \code{"WL"} that call for Bayes modal, expected a posteriori and weighted likelihood estimation respectively; \code{"fixed4"} and \code{"fixed7"} that perform fixed stepsize adjustment (i.e. increase or decrease of constant magnitude) with step 0.4 and 0.7 respectively; and \code{"var"} for variable stepsize adjustment. Note that in case odf stepsize adjustment, the range of difficultyu parameters must be provided through the \code{bRange} argument, as vector of two components (default value being \code{c(-2,2)}). See Dodd, De Ayala, and Koch (1995) for further details.
 }

\note{
It has been shown that in some cases the weighted likelihood estimator and the Bayes modal estimator with Jeffreys prior return exactly the same ability estimates. This is the case under the 2PL model, and subsequently the 1PL model (Warm, 1989) as well as under all polytomous models currently available (Magis, 2015). Nevertheless, both estimators remain available since (a) Jeffreys prior can also be considered with the EAP estimator, and (b) the 3PL and 4PL models are also available.
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability. \emph{Journal of Mathematical Psychology, 6}, 258-276. doi: 10.1016/0022-2496(69)90005-4

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Applied Psychological Measurement, 6}, 431-444. doi: 10.1177/014662168200600405

Dodd, B. G., De Ayala, R. J., and Koch, W. R. (1995) . Computerized adaptive testing with polytomous items. \emph{Applied Psychological Measurement, 19}, 5-22. doi: 10.1177/014662169501900103

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error.} 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal
 Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

Magis, D. (2015). A note on weighted likelihood and Jeffreys modal estimation of proficiency levels in polytomous item response models. \emph{Psychometrika, 80}, 200-204. doi: 10.1007/S11336-013-9378-5

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. doi: 10.1007/BF02294627
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}
 }

\seealso{
\code{\link{eapEst}}, \code{\link{semTheta}}, \code{\link{genPolyMatrix}}
} 

\examples{

## Dichotomous models ##

 # Generation of an item bank under 3PL with 100 items 
 m.3PL <- genDichoMatrix(100, model = "3PL")
 m.3PL <- as.matrix(m.3PL)


 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.3PL)

 # ML estimation
 thetaEst(m.3PL, x, method = "ML")

 # BM estimation, standard normal prior distribution
 thetaEst(m.3PL, x)

 # BM estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.3PL, x, method = "BM", priorDist = "unif", priorPar = c(-2, 2))

 # BM estimation, Jeffreys' prior distribution  
 thetaEst(m.3PL, x, method = "BM", priorDist = "Jeffreys")

 # EAP estimation, standard normal prior distribution
 thetaEst(m.3PL, x, method = "EAP")

 # EAP estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.3PL, x, method = "EAP", priorDist = "unif", priorPar = c(-2, 2))

 # EAP estimation, Jeffreys' prior distribution  
 thetaEst(m.3PL, x, method = "EAP", priorDist = "Jeffreys")

 # WL estimation
 thetaEst(m.3PL, x, method = "WL")

 # Creation of two constant patterns and estimation with WL, 
 # 'fixed4', 'fixed7' and 'var' stepsize adjustments
 x0 <- rep(0,nrow(m.3PL))
 x1 <- x0 + 1
 thetaEst(m.3PL, x0, constantPatt = "WL") # equivalent to thetaEst(m.3PL, x0, method = "WL")
 thetaEst(m.3PL, x1, constantPatt = "WL") # equivalent to thetaEst(m.3PL, x1, method = "WL")
 thetaEst(m.3PL, x0, constantPatt = "fixed4") 
 thetaEst(m.3PL, x1, constantPatt = "fixed4") 
 thetaEst(m.3PL, x0, constantPatt = "fixed7") 
 thetaEst(m.3PL, x1, constantPatt = "fixed7") 
 thetaEst(m.3PL, x0, constantPatt = "var") 
 thetaEst(m.3PL, x1, constantPatt = "var") 

\dontrun{

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.GRM, model = "GRM")

 # ML estimation
 thetaEst(m.GRM, x, model = "GRM", method = "ML")

 # BM estimation, standard normal prior distribution
 thetaEst(m.GRM, x, model = "GRM")

 # BM estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.GRM, x, model = "GRM", method = "BM", priorDist = "unif", 
          priorPar = c(-2, 2))

 # BM estimation, Jeffreys' prior distribution  
 thetaEst(m.GRM, x, model = "GRM", method = "BM", priorDist = "Jeffreys")

 # EAP estimation, standard normal prior distribution
 thetaEst(m.GRM, x, model = "GRM", method = "EAP")

 # EAP estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.GRM, x, model = "GRM", method = "EAP", priorDist = "unif", 
          priorPar = c(-2, 2))

 # EAP estimation, Jeffreys' prior distribution  
 thetaEst(m.GRM, x, model = "GRM", method = "EAP", priorDist = "Jeffreys")

 # WL estimation
 thetaEst(m.GRM, x, model = "GRM", method = "WL")


 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.PCM, model = "PCM")

 # ML estimation
 thetaEst(m.PCM, x, model = "PCM", method = "ML")

 # BM estimation, standard normal prior distribution
 thetaEst(m.PCM, x, model = "PCM")

 # BM estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.PCM, x, model = "PCM", method = "BM", priorDist = "unif", 
          priorPar = c(-2, 2))

 # BM estimation, Jeffreys' prior distribution  
 thetaEst(m.PCM, x, model = "PCM", method = "BM", priorDist = "Jeffreys")

 # EAP estimation, standard normal prior distribution
 thetaEst(m.PCM, x, model = "PCM", method = "EAP")

 # EAP estimation, uniform prior distribution upon range [-2,2]
 thetaEst(m.PCM, x, model = "PCM", method = "EAP", priorDist = "unif", 
          priorPar = c(-2, 2))

 # EAP estimation, Jeffreys' prior distribution  
 thetaEst(m.PCM, x, model = "PCM", method = "EAP", priorDist = "Jeffreys")

 # WL estimation
 thetaEst(m.PCM, x, model = "PCM", method = "WL")
 }
 }
