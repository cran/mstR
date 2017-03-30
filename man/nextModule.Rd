\name{nextModule}
\alias{nextModule}

\title{Selection of the next module in MST}

\description{
 This command selects the next module to be administered in the multistage test, either bases on IRT scoring or on test score and by either providing thresholds or optimally selecting the next module. 
 }

\usage{
nextModule(itemBank, modules, transMatrix, model = NULL, current.module,
   out, x = NULL, cutoff = NULL, theta = 0, criterion = "MFI",
   priorDist = "norm", priorPar = c(0, 1), D = 1, range = c(-4, 4), 
   parInt = c(-4, 4, 33))
 }


\arguments{
 \item{itemBank}{a suitable matrix of item parameters. See \bold{Details}.}
\item{modules}{a binary matrix that specifies the item membership to the modules. See \bold{Details}.}
\item{transMatrix}{a binary squared matrix representing the structure of the MST and the transitions between the moduels and the stages. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{current.module}{integer: the module number (defined as the corresponding column number in the \code{modules} matrix) that indicates the last administered module.}
 \item{out}{numeric: the vector of item indicators (defined as the row numbers in the \code{itemBank} matrix) of previously administered items.}
\item{x}{either a numeric vector of responses to previously administered items or \code{NULL} (default). Ignored if \code{criterion} is either \code{MFI} or \code{random}. See \bold{Details}.}
\item{cutoff}{either a suitable matrix of cut-off values or \code{NULL} (default). See \code{Details}.}
 \item{theta}{numeric: the current ability level for selecting the next module (default is 0). It can also hold the current test score made of the sum of responses to all administered items (when \code{cutoff} is provided). See \bold{Details}.}
 \item{criterion}{character: the method for next item selection. Possible values are \code{"MFI"} (default), \code{"MLWMI"}, \code{"MPWMI"}, \code{"MKL"}, \code{"MKLP"} and \code{"random"}. Ignored if \code{cutoff} is not \code{NULL}. See \bold{Details}.} 
\item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{criterion} is neither \code{"MPWMI"} nor \code{"KLP"}, or if \code{cutoff} is not \code{NULL}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{criterion} is neither \code{"MPWI"} nor \code{"KLP"}, or if \code{cutoff} is not \code{NULL}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).}
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}).}
 \item{parInt}{numeric: a vector of three numeric values, specifying respectively the lower bound, the upper bound and the number of quadrature points for numerical integration (default is \code{c(-4,4,33)}). Ignored if \code{method} is neither \code{"MLWMI"}, \code{"MPWMI"}, \code{"KL"}, nor \code{"KLP"}, or if \code{cutoff} is not \code{NULL}. See \bold{Details}.}
}

\value{
 A list with five arguments:
\item{module}{the selected module (identified by its column number in the \code{modules} argument.}
 \item{items}{the items that belong to the selected module (identified by their number in the item bank).}
 \item{par}{the matrix of item parameters of the selected items (one row per item).}
 \item{info}{either the provisional ability level or score when \code{cutoff} is not \code{NULL}; or \code{NA} when random selection is performed.}
 \item{criterion}{the value of the \code{criterion} argument.}
 }
 

\details{
This function permits to select the next module of the MST. It works with both dichotomous and polytomous item banks.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The \code{modules} argument must be a binary 0/1 matrix with as many rows as the item bank \code{itemBank} and as many columns as the number of modules. Values of 1 indicate to which module(s) the items belong to, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the \emph{i}-th item belongs to the \emph{j}-th module.

The \code{transMatrix} argument must be a binary 0/1 square matrix with as many rows (and columns) as the number of modules. All values of 1 indicate the possible transitions from one module to another, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the MST can move from \emph{i}-th module to \emph{j}-th module.

The two main approaches to select the next module are based on cut-off scores (to be provided) or by optimal module selection.

Optimal module selection is performed by providing an appropriate value to the \code{criterion} argument. Possible methods are: 
\enumerate{
\item{\code{"MFI"} for maximum Fisher information(default);}
\item{\code{"MLWMI"} for maximum likelihood weighted module information;}
\item{\code{"MPWMI"} for maximum posterior weighted module information;}
\item{\code{"MKL"} for module Kullabck-Leibler selection;}
\item{\code{"MKLP"} for module posterior Kullabck-Leibler selection;}
\item{\code{"random"} for random selection.}
}
See \code{\link{MWMI}} and \code{\link{MKL}} for further details.

In case of selection by predefined cut-off scores, the \code{cutoff} argument must be supplied by a matrix with as many rows as modules in the item bank and with two columns. Each row of the \code{cutoff} matrix holds the lower and upper cut-off value and the next module is selected (among those available in the correct stage) such that the current \code{ability} score lies between the two cut-off scores. For instance, if modules 2 and 3 are the two potential modules to be selected, if \code{cutoff[2, 1:2]} is \code{c(0, 5)} and if \code{cutoff[3, 1:2]} is \code{c(5, 10)}, then module 2 will be selected if \code{theta} belongs to \eqn{[0, 5)} and module 3 is selected if \code{theta} belongs to \eqn{[5, 10)}. Other examples are provided in the \bold{Examples} section. All modules from stage 1 must have cut-off scores too for compatibility, but any value (including \code{NA}) is allowed as they will not be used. Note that integer cut-off scores (when \code{theta} is the test score) or numeric values (when \code{theta} is an ability estimate) are allowed in the \code{cutoff} matrix.

By default \code{cutoff} is \code{NULL} and optimal module selection is performed). 
 }


\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error.} 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{genPolyMatrix}}, \code{\link{MWMI}}, \code{\link{MKL}} 
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

 # Creation of the 'modules' matrix to list item membership in each module
 modules <- matrix(0, 55, 7)
 modules[1:8, 1] <- modules[9:16, 2] <- modules[17:22, 3] <- 1
 modules[23:28, 4] <- modules[29:37, 5] <- modules[38:46, 6] <- 1
 modules[47:55, 7] <- 1

 # Creation of the transition matrix to define a 1-2-3 MST
trans <- matrix(0, 7, 7)
trans[1, 3:4] <- trans[2, 3:4] <- trans[3, 5:7] <- trans[4, 5:7] <- 1

 # Module 1 previously administered, provisional ability 0, MFI criterion
 nextModule(it, modules, trans, current.module = 1, out = 1:8)

 # Generation of item responses for module 1
 x <- genPattern(0, it[1:8,])

 # MLWMI criterion 
 nextModule(it, modules, trans, current.module = 1, out = 1:8, x = x, criterion = "MLWMI")

 # MPWMI criterion
 nextModule(it, modules, trans, current.module = 1, out = 1:8, x = x, criterion = "MPWMI")

 # MKL criterion
 nextModule(it, modules, trans, current.module = 1, out = 1:8, x = x, criterion = "MKL")

 # MKLP criterion
 nextModule(it, modules, trans, current.module = 1, out = 1:8, x = x, criterion = "MKLP")

 # Creation of cut-off scores for ability levels: cut score 0 between modules 3 and 4
 # and cut scores -1 and 1 between modules 5, 6 and 7
cut <- matrix(NA, 7, 2)
cut[3,] <- c(-Inf, 0)
cut[4,] <- c(0, Inf)
cut[5,] <- c(-Inf, -1)
cut[6,] <- c(-1, 1)
cut[7,] <- c(1, Inf)

 # Selection by cut-off score, module 1 previously administered, current ability level 0 
 # (=> module 4 is chosen)
 nextModule(it, modules, trans, current.module = 1, out = 1:8, cutoff = cut, theta = 0)

 # Same with current ability level -0.5 (=> module 3 is chosen)
 nextModule(it, modules, trans, current.module = 1, out = 1:8, cutoff = cut, theta = -0.5)

 # Modules 1 and 3 previously administered, current ability level 0 (=> module 6 is chosen)
 nextModule(it, modules, trans, current.module = 3, out = c(1:8, 17:22), cutoff = cut,  
           theta = 0)

 # Same with current ability level 2 (=> module 7 is chosen)
 nextModule(it, modules, trans, current.module = 3, out = c(1:8, 17:22), cutoff = cut, 
           theta = 2)

 # Creation of cut-off scores for test scores: cut score 4 between modules 3 and 4
 # and cut scores 5 and 9 between modules 5, 6 and 7
cut.score <- matrix(NA, 7, 2)
cut.score [3,] <- c(0, 4)
cut.score [4,] <- c(4, 8)
cut.score [5,] <- c(0, 5)
cut.score [6,] <- c(5, 9)
cut.score [7,] <- c(9, 14)

 # Module 1 previously administered, current test score 1 (=> module 3 is chosen)
 nextModule(it, modules, trans, current.module = 1, out = 1:8, cutoff = cut.score, 
           theta = 1)

 # Modules 1 and 3 previously administered, current tes score 6 (=> module 6 is chosen)
 nextModule(it, modules, trans, current.module = 3, out = c(1:8, 17:22), cutoff = cut.score, 
           theta = 6)


## Polytomous models ##

 # Same structure as above but parameters are now generated from PCM with at most
 # 4 categories
 it.pol <- genPolyMatrix(55, model = "PCM", nrCat = 4)
 it.pol <- as.matrix(it)

 # Module 1 previously administered, provisional ability 0, MFI criterion
 nextModule(it.pol, modules, trans, model = "PCM", current.module = 1, out = 1:8)

 # MLWMI criterion 
 nextModule(it.pol, modules, trans, model = "PCM", current.module = 1, out = 1:8, x = x, 
           criterion = "MLWMI")

# MKL criterion 
 nextModule(it.pol, modules, trans, model = "PCM", current.module = 1, out = 1:8, x = x, 
           criterion = "MKL")

# MKLP criterion 
 nextModule(it.pol, modules, trans, model = "PCM", current.module = 1, out = 1:8, x = x, 
           criterion = "MKLP")

# Selection by cut-off score, module 1 previously administered, current ability level 0 
 # (=> module 4 is chosen)
 nextModule(it.pol, modules, trans, model = "PCM", current.module = 1, out = 1:8, 
           cutoff = cut, theta = 0)
  }
