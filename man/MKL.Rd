\name{MKL}
\alias{MKL}

\title{Module Kullback-Leibler (MKL) and posterior module Kullback-Leibler (MKLP)}

\description{
 This command returns the value of the Kullback-Leibler (MKL) or the posterior Kullback-Leibler (MKLP) weighted likelihood for a given target module and an item bank (both under dichotomous and polytomous IRT models).
 }

\usage{
MKL(itemBank, modules, target.mod, theta = NULL, it.given, x, model = NULL, 
    lower = -4, upper = 4, nqp = 33, type = "MKL", priorDist = "norm", 
    priorPar = c(0, 1), D = 1)
 }


\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{modules}{a binary matrix that specifies the item membership to th emodules. See \bold{Details}.}
 \item{target.mod}{numeric: the module (referred to as its column number in the \code{modules} matrix) for which the information must be computed.}
\item{theta}{either the provisional ability level or \code{NULL} (default). See \bold{Details}.}
 \item{it.given}{numeric: a vector of item indicators for all previously administered items.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items). The length of \code{x} must be equal to the length of \code{it.given}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numerical integration (default is -4).}
 \item{upper}{numeric: the upper bound for numerical integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of Kullback-Leibler information to be computed. Possible values are \code{"MKL"} (default) and \code{"MKLP"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 The required KL (or KLP) weighted module likelihood for the target module.
 }
 
\details{
This function extends the KL and the KLP methods to select the next item in CAT, to the MST framework. This command serves as a subroutine for the \code{\link{nextModule}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

FROM HERE

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set \eqn{\hat{\theta}_k} as the provisional ability estimate (with the first \emph{k} responses). Set \eqn{M} as the number of items in the target module of interest (not yet administered). Set also \eqn{L(\theta | x_1, ..., x_k) } as the likelihood function of the first \eqn{k} items and evaluated at \eqn{\theta}. Set finally \eqn{P_{jt}(\theta)} as the probability of answering response category \emph{t} to item \emph{j} of the target module (\eqn{j = 1, ..., M}) for a given ability level \eqn{\theta}. Then, module Kullack-Leibler (MKL) information is defined as
\deqn{MKL(\theta || \hat{\theta}_k) = \sum_{j=1}^{M} \sum_{t=0}^{g_j} \,P_{jt}(\hat{\theta}_k) \,\log \left( \frac{P_{jt}(\hat{\theta}_k)}{P_{jt}(\theta)}\right).}

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. Set simply \eqn{P_j(\theta)} as the probability of answering item \emph{j} correctly (\eqn{j = 1, ..., M}) for a given ability level \eqn{\theta}. Then, MKL information reduces to
\deqn{MKL(\theta || \hat{\theta}) = \sum_{j=1}^{M} \left\{P_j(\hat{\theta}) \,\log \left( \frac{P_j(\hat{\theta}_k)}{P_j(\theta)}\right) + [1-P_j(\hat{\theta}_k)] \,\log \left( \frac{1-P_j(\hat{\theta}_k)}{1-P_j(\theta)}\right) \right\}.}

The quantity that is returned by this \code{MKL} function is either: the likelihood function weighted by module Kullback-Leibler information (the MKL value):
\deqn{MKL(\hat{\theta}_k) = \int MKL(\theta || \hat{\theta}_k) \, L(\theta | x_1, ..., x_k) \,d\theta}
or the posterior function weighted by module Kullback-Leibler information (the MKLP value):
\deqn{MKLP(\hat{\theta}) = \int MKL(\theta || \hat{\theta}_k) \, \pi(\theta) \,L(\theta | x_1, ..., x_k) \,d\theta}
 where \eqn{\pi(\theta)} is the prior distribution of the ability level.

 These integrals are approximated by the \code{\link{integrate.mstR}} function. The range of integration is set up 
 by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points. The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

 The argument \code{type} defines the type of information to be computed. The default value, \code{"MKL"}, computes the MKL value, while the MKLP value is obtained with \code{type="MKLP"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. 
 The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that 
 uniform distribution. By default, the standard normal prior distribution is assumed. This argument is ignored whenever \code{method} is not \code{"MKLP"}.

The provisional response pattern and the related administered items are provided by the vectors \code{x} and \code{it.given} respectively. The target module (for which the maximum information is computed) is given by its column number in the \code{modules} matrix, through the \code{target.mod} argument. 

The provisioal (ad-interim) ability level can be provided through the \code{theta} argument. If not provided or set to \code{NULL} (default value), it is then internally computed as the ML estimate of ability for the given response pattern \code{x} and the previously administered items \code{it.given}.
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
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

 # MKL for module 3
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x)

 # Same with pre-estimation of ability by ML (same result)
 th <- thetaEst(it[1:8,], x, method = "ML") 
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x, theta = th)

 # Same with pre-estimation of ability by EAP (different result)
 th <- thetaEst(it[1:8,], x, method = "EAP") 
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x, theta = th)

 # MKLP for module 3
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x, type = "MKLP")

 # MKL for for module 3, different integration range
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x, lower = -2, upper = 2, nqp = 20)

 # MKLP for module 3, uniform prior distribution on the range [-2,2]
 MKL(it, modules, target.mod = 3, it.given = 1:8, x = x, type = "MKLP", 
     priorDist = "unif", priorPar = c(-2, 2))


## Polytomous models ##

 # Same structure as above but parameters are now generated from PCM with at most
 # 4 categories
 it.pol <- genPolyMatrix(55, model = "PCM", nrCat = 4)
 it.pol <- as.matrix(it)

 # Creation of the response pattern for module 1 and true ability level 0
 x <- genPattern(th = 0, it = it.pol[1:8,], seed = 1)

# MKL for module 3
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, model = "PCM")

 # Same with pre-estimation of ability by ML (same result)
 th <- thetaEst(it.pol[1:8,], x, method = "ML", model = "PCM") 
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, theta = th, 
     model = "PCM")

 # Same with pre-estimation of ability by EAP (different result)
 th <- thetaEst(it.pol[1:8,], x, method = "EAP", model = "PCM") 
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, theta = th, 
     model = "PCM")

 # MKLP for module 3
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, type = "MKLP", 
     model = "PCM")

 # MKL for for module 3, different integration range
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, lower = -2, 
     upper = 2, nqp = 20, model = "PCM")

 # MKLP for module 3, uniform prior distribution on the range [-2,2]
 MKL(it.pol, modules, target.mod = 3, it.given = 1:8, x = x, type = "MKLP", 
     priorDist = "unif", priorPar = c(-2, 2), model = "PCM")
 }