\name{startModule}
\alias{startModule}

\title{Selection of the first module in MST}

\description{
 This command selects the first module of the multistage test, either randomly or on the basis of the module information function.
 }

\usage{
startModule(itemBank, modules, transMatrix, model = NULL, fixModule = NULL, 
  seed = NULL, theta = 0, D = 1) 
 }


\arguments{
 \item{itemBank}{a suitable matrix of item parameters. See \bold{Details}.}
\item{modules}{a binary matrix that specifies the item membership to th emodules. See \bold{Details}.}
\item{transMatrix}{a binary squared matrix representing the structure of the MST and the transitions between the moduels and the stages. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{fixModule}{either a an integer value or \code{NULL} (default). See \bold{Details}.}
 \item{seed}{either a numeric value, \code{NA} or \code{NULL} (default). Ignored if \code{fixModule} is not \code{NULL}. See \bold{Details}.}
 \item{theta}{numeric: the initial ability level for selecting the first module (default is 0). Ignored if either \code{fixModule} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 A list with four arguments:
\item{module}{the selected module (identified by its column number in the \code{modules} argument.}
 \item{items}{the items that belong to the selected module (identified by their number in the item bank).}
 \item{par}{the matrix of item parameters of the selected items (one row per item).}
 \item{thStart}{the starting ability level used for selecting the module or \code{NA} (if not applicable).}
 }
 

\details{
 This function permits to select the first module of the MST. It works with both dichotomous and polytomous item banks.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The \code{modules} argument must be a binary 0/1 matrix with as many rows as the item bank \code{itemBank} and as many columns as the number of modules. Values of 1 indicate to which module(s) the items belong to, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the \emph{i}-th item belongs to the \emph{j}-th module.

The \code{transMatrix} argument must be a binary 0/1 square matrix with as many rows (and columns) as the number of modules. All values of 1 indicate the possible transitions from one module to another, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the MST can move from \emph{i}-th module to \emph{j}-th module.


 The first module of the MST can be selected by one of the following methods.
 \enumerate{
 \item{By specifying the module to be administered. The argument \code{fixModule} then holds the module number as listed in the \code{modules} or \code{transMatrix}. Setting \code{fixModule} to \code{NULL} (default value) disables this method.}
 \item{By selecting it randomly into the \code{modules} matrix. The argument \code{seed} permits to fix the random selection by       specifying the random seed number. Setting \code{seed} to \code{NA} disables the random seed (though items are still picked up randomly in the bank); in other words, successive runs of \code{startModule} with \code{seed=NA} may lead to different module selection. Setting \code{seed} to \code{NULL} (default value) disables this selection method.}
 \item{By selecting the module according to an initial ability value. In this case, the selected module is such that the information fucntion (computed with the items of the module) is maximal for the given initial ability value.}
}

 The third method above will be used if and only if both \code{fixModule} and \code{seed} arguments are fixed to \code{NULL}.
 Otherwise, one of the first two methods will be used.
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
\code{\link{genPolyMatrix}} 
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

 # Creation of the transition matrix to define a 1-2-3 MST
trans <- matrix(0, 7, 7)
trans[1, 3:4] <- trans[2, 3:4] <- trans[3, 5:6] <- trans[4, 6:7] <- 1

 # Selection of module 2 as starting module
 startModule(it, modules, trans, fixModule = 2)
 
\dontrun{

 # Selection of module 3 (not from stage 1 => mistake)
 startModule(it, modules, trans, fixModule = 3)
}

 # Random selection of starting module 
 startModule(it, modules, trans, seed = 1)

 # Selection by maximizing information at ability level 0
 startModule(it, modules, trans, theta = 0)


## Polytomous models ##

 # Same structure as above but parameters are now generated from PCM with at most
 # 4 categories
 it <- genPolyMatrix(55, model = "PCM", nrCat = 4)
 it <- as.matrix(it)

 # Selection of module 2 as starting module
 startModule(it, modules, trans, fixModule = 2, model = "PCM")
 
\dontrun{

 # Selection of module 3 (not from stage 1 => mistake)
 startModule(it, modules, trans, fixModule = 3, model = "PCM")
}

 # Random selection of starting module 
 startModule(it, modules, trans, seed = 1, model = "PCM")

 # Selection by maximizing information at ability level 0
 startModule(it, modules, trans, theta = 0, model = "PCM")
 }
