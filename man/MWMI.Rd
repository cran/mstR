\name{MWMI}
\alias{MWMI}

\title{Maximum likelihood weighted module information (MLWMI) and maximum posterior weighted module information (MPWMI)}

\description{
 This command returns the value of the likelihood (MLWMI) or the posterior (MPWMI) weighted module information for a given module and an item bank (both under dichotomous and polytomous IRT models).
 }

\usage{
MWMI(itemBank, modules, target.mod, it.given, x, model = NULL, lower = -4, 
    upper = 4, nqp = 33, type = "MLWMI", priorDist = "norm", priorPar = c(0, 1), 
    D = 1)
 }


\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{modules}{a binary matrix that specifies the item membership to th emodules. See \bold{Details}.}
 \item{target.mod}{numeric: the module (referred to as its column number in the \code{modules} matrix) for which the information must be computed.}
 \item{it.given}{numeric: a vector of item indicators for all previously administered items.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items). The length of \code{x} must be equal to the length of \code{it.given}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numerical integration (default is -4).}
 \item{upper}{numeric: the upper bound for numerical integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of information to be computed. Possible values are \code{"MLWMI"} (default) and \code{"MPWMI"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 The required (likelihood or posterior) weighted module information for the selected module.
 }
 
\details{
This function extends the MLWI and the MPWI methods to select the next item in CAT, to the MST framework. This command serves as a subroutine for the \code{\link{nextModule}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set also
 \eqn{I(\theta)} as the information function of the module of interest (specified through \code{target.mod}), made by the sum of all item informations functions from this module, and evaluated at \eqn{\theta}. Set finally \eqn{L(\theta | x_1, ..., x_k)} as the likelihood function evaluated at \eqn{\theta}, given the provisional response pattern. Then, the LWMI for the module is given by
 \deqn{LWMI = \int I(\theta) L(\theta | x_1, ..., x_k) d\theta} and the PWMI by  
 \deqn{PWMI = \int I(\theta) \pi(\theta) L(\theta | x_1, ..., x_k) d\theta} 
 where \eqn{\pi(\theta)} is the prior distribution of the ability level. 

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. But except from this difference, the previous definitions of LWI and PWI remain valid.

 These integrals are approximated by the \code{\link{integrate.mstR}} function. The range of integration is set up 
 by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points. The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

 The argument \code{type} defines the type of information to be computed. The default value, \code{"MLWMI"}, computes the MLWMI value, while the MPWMI value is obtained with \code{type="MPWMI"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. 
 The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that 
 uniform distribution. By default, the standard normal prior distribution is assumed. This argument is ignored whenever \code{method} is not \code{"MPWMI"}.

The provisional response pattern and the related administered items are provided by the vectors \code{x} and \code{it.given} respectively. The target module (for which the maximum information is computed) is given by its column number in the \code{modules} matrix, through the \code{target.mod} argument. 
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 


 }

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{nextModule}}, \code{\link{integrate.mstR}}, \code{\link{genPolyMatrix}}
} 

\examples{

## Dichotomous models ##

 # Generation of an item bank under 2PL, made of 7 successive modules that target
 # different average ability levels and with different lengths
 # (the first generated item parameters hold two modules of 8 items each)
 it <- rbind(genDichoMatrix(16, model = "2PL"),
             genDichoMatrix(6, model = "2PL", bPrior = c("norm", -1, 1)),
             genDichoMatrix(6, model = "2PL", bPrior = c("norm", 1, 1)),
             genDichoMatrix(9, model = "2PL", bPrior = c("norm", -2, 1)),
             genDichoMatrix(9, model = "2PL", bPrior = c("norm", 0, 1)),
             genDichoMatrix(9, model = "2PL", bPrior = c("norm", 2, 1)))
 it <- as.matrix(it)

 # Creation of the 'module' matrix to list item membership in each module
 modules <- matrix(0, 55, 7)
 modules[1:8, 1] <- modules[9:16, 2] <- modules[17:22, 3] <- 1
 modules[23:28, 4] <- modules[29:37, 5] <- modules[38:46, 6] <- 1
 modules[47:55, 7] <- 1

 # Creation of the response pattern for module 1 and true ability level 0
 x <- genPattern(th = 0, it = it[1:8,], seed = 1)

 # MLWMI for module 3
 MWMI(it, modules, target.mod = 3, it.given = 1:8, x = x)

 # MPWMI for module 3
 MWMI(it, modules, target.mod = 3, it.given = 1:8, x = x, type = "MPWMI")

 # MLWMI for for module 3, different integration range
 MWMI(it, modules, target.mod = 3, it.given = 1:8, x = x, lower = -2, upper = 2, nqp = 20)

 # MPWI for module 3, uniform prior distribution on the range [-2,2]
 MWMI(it, modules, target.mod = 3, it.given = 1:8, x = x, type = "MPWMI", 
     priorDist = "unif", priorPar = c(-2, 2))


## Polytomous models ##

 # Same structure as above but parameters are now generated from PCM with at most
 # 4 categories
 it.pol <- genPolyMatrix(55, model = "PCM", nrCat = 4)
 it.pol <- as.matrix(it)

 # Creation of the response pattern for module 1 and true ability level 0
 x <- genPattern(th = 0, it = it.pol[1:8,], seed = 1)

# MLWMI for module 3
 MWMI(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, model = "PCM")

 # MPWMI for module 3
 MWMI(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, type = "MPWMI", 
     model = "PCM")

 # MLWMI for for module 3, different integration range
 MWMI(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, lower = -2, 
     upper = 2, nqp = 20, model = "PCM")

 # MPWI for module 3, uniform prior distribution on the range [-2,2]
 MWMI(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, type = "MPWMI", 
     priorDist = "unif", priorPar = c(-2, 2), model = "PCM")
 }