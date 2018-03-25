\name{randomMST}
\alias{randomMST}
\alias{print.mst}
\alias{plot.mst}


\title{Random generation of multistage tests (dichotomous and polytomous models)}

\description{
 This command generates a response pattern to a multistage test, for a given item bank (with either dichotomous or polytomous models), an MST structure for modules and stages, a true ability level, and several lists of MST parameters.
 }

\usage{
randomMST(trueTheta, itemBank, modules, transMatrix, model = NULL, 
  responses = NULL, genSeed = NULL, start = list(fixModule = NULL, seed = NULL, 
  theta = 0, D = 1), test = list(method = "BM", priorDist = "norm", 
  priorPar = c(0, 1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
  moduleSelect = "MFI", constantPatt = NULL, cutoff = NULL, randomesque = 1,
  random.seed = NULL, score.range = "all"), final = list(method = "BM", 
  priorDist = "norm", priorPar = c(0, 1), range = c(-4, 4), D = 1, 
  parInt = c(-4, 4, 33), alpha = 0.05), allTheta = FALSE, save.output = FALSE, 
  output = c("path", "name", "csv"))
\method{print}{mst}(x, ...)
\method{plot}{mst}(x, show.path = TRUE, border.col = "red", arrow.col = "red",
  module.names = NULL, save.plot = FALSE, save.options = c("path", "name", "pdf"),...)
 }

\arguments{
 \item{trueTheta}{numeric: the value of the true ability level.}
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{modules}{a binary matrix that specifies the item membership to the modules. See \bold{Details}.}
\item{transMatrix}{a binary squared matrix representing the structure of the MST and the transitions between the moduels and the stages. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
\item{responses}{either \code{NULL} (default) or a vector of pre-specified item responses with as many components as the rows of \code{itemBank}. See \bold{Details}.}
\item{genSeed}{either a numeric value to fix the random generation of responses pattern or \code{NULL} (default). Ignored if \code{responses} is not \code{NULL}. See \bold{Details}.}
 \item{start}{a list with the options for starting the multistage test. See \bold{Details}.}
 \item{test}{a list with the options for provisional ability estimation and next module selection. See \bold{Details}.}
 \item{final}{a list with the options for final ability estimation or scoring. See \bold{Details}.}
\item{allTheta}{logical: should all provisional ability estimates and standard errors be computed and returned (even within each module)? Default is \code{FALSE}, meaning that provisional ability estimates and standard errors are computed only at the end of each module administration.}
 \item{save.output}{logical: should the output be saved in an external text file? (default is \code{FALSE}).}
 \item{output}{character: a vector of three components. The first component is either the file path to save the output of \code{"path"} (default), the second component is the name of the output file, and the third component is the file type, either \code{"txt"} or \code{"csv"} (default). See \bold{Details}.}
 \item{x}{either an object of class "mst", typically an output of \code{randomMST} function, or a transition matrix (for \code{plot.mst()} function only).}
\item{show.path}{logical: should the selected path (i.e. set of successive modules) be highlighted in the plot (default is TRUE)?}
\item{border.col}{character: the color for the rectangle border of the path (i.e. selected modules). Default is \code{"red"}. Ignored if \code{show.path} is \code{FALSE}.}
\item{arrow.col}{character: the color for the connecting arrows in the path (i.e. between selected modules). Default is \code{"red"}. Ignored if \code{show.path} is \code{FALSE}.}
\item{module.names}{either \code{NULL} (default) or a vector of character names for the modules. See \bold{Details}.}
 \item{save.plot}{logical: should the plot be saved in an external figure? (default is \code{FALSE}).}
 \item{save.options}{character: a vector of three components. The first component is either the file path or \code{"path"} (default), the second component is the name of the output file or ,\code{"name"} (default), and the third component is the file extension, either \code{"pdf"} (default) or \code{"jpeg"}. Ignored if \code{save.plot} is \code{FALSE}. See \bold{Details}.}
 \item{...}{other generic arguments to be passed to \code{print} and \code{plot} functions.}
}


\value{
 The function \code{randomMST} returns a list of class "mst" with the following arguments:
 \item{trueTheta}{the value of the \code{trueTheta} argument.}
\item{selected.modules}{a vector with the modules (identified by their position in the transition matrix) that were selected for the MST.}
\item{items.per.module}{a vector with the number of items per selected module (in the same order as in \code{selected.modules}).}
\item{transMatrix}{the value of the \code{transMatrix} argument.}
\item{model}{the value of the \code{model} argument.}
 \item{testItems}{a vector with the items that were administered during the test.}
 \item{itemPar}{a matrix with the parameters of the items administered during the test.}
 \item{pattern}{the generated (or selected) response pattern (as vector of 0 and 1 entries for dichotomous items or positive integer values for polytomous items).}
 \item{thetaProv}{a vector with the provisional ability estimates (or test scores if \code{test$method} is  \code{"score"}).}
 \item{seProv}{a vector with the standard errors of the provisional ability estimates (or vector of \code{NA}'s if \code{test$method} is  \code{"score"}).}
 \item{thFinal}{the final ability estimate (or test score if \code{test$method} is  \code{"score"}).}
 \item{seFinal}{the standard error of the final ability estimate (or \code{NA} if \code{test$method} is  \code{"score"}).}
 \item{ciFinal}{the confidence interval of the final ability estimate (or \code{c(NA, NA)} if \code{test$method} is  \code{"score"}).}
\item{genSeed}{the value of the \code{genSeed} argument.}
 \item{startFixModule}{the value of the \code{start$fixModule} argument (or its default value if missing).}
 \item{startSeed}{the value of the \code{start$seed} argument (or its default value if missing).}
 \item{startTheta}{the value of the \code{start$theta} argument (or its default value if missing).}
 \item{startD}{the value of the \code{start$D} argument (or its default value if missing).}
 \item{startThStart}{the starting ability value used for selecting the first module of the test.}
 \item{startSelect}{the value of the \code{start$startSelect} argument (or its default value if missing).}
 \item{provMethod}{the value of the \code{test$method} argument (or its default value if missing).}
 \item{provDist}{the value of the \code{test$priorDist} argument (or its default value if missing).}
 \item{provPar}{the value of the \code{test$priorPar} argument (or its default value if missing).}
 \item{provRange}{the value of the \code{test$range} argument (or its default value if missing).}
 \item{provD}{the value of the \code{test$D} argument (or its default value if missing)or \code{NA} if \code{model} is not \code{NULL}.}
 \item{moduleSelect}{the value of the \code{test$moduleSelect} argument (or its default value if missing).}
\item{constantPattern}{the value of the \code{test$constantPatt} argument (or its default value if missing).}
 \item{cutoff}{the value of the \code{test$cutoff} argument (or its default value if missing).}
 \item{randomesque}{the value of the \code{test$randomesque} argument (or its default value if missing).}
 \item{random.seed}{the value of the \code{test$random.seed} argument (or its default value if missing).}
 \item{score.range}{the value of the \code{test$score.range} argument (or its default value if missing).}
 \item{finalMethod}{the value of the \code{final$method} argument (or its default value if missing).}
 \item{finalDist}{the value of the \code{final$priorDist} argument (or its default value if missing).}
 \item{finalPar}{the value of the \code{final$priorPar} argument (or its default value if missing).}
 \item{finalRange}{the value of the \code{final$range} argument (or its default value if missing).}
 \item{finalD}{the value of the \code{final$D} argument (or its default value if missing), or \code{NA} if \code{model} is not \code{NULL}.}
 \item{finalAlpha}{the value of the \code{final$alpha} argument (or its default value if missing).}
 \item{save.output}{the value of the \code{save.output} argument.}
 \item{output}{the value of the \code{output} argument.}
\item{assigned.responses}{a logical value, being \code{TRUE} if \code{responses} was provided or \code{FALSE} \code{responses} was set to \code{NULL}.}
\item{allTheta}{either a table with all ad-interim ability estimates (even within module, in the CAT spirit) if \code{allTheta} is set to \code{TRUE}, or \code{NULL} if \code{allTheta} is set to \code{FALSE}.}
\item{assigned.responses}{the value of the \code{responses} argument (or its default value if missing).}
 
The function \code{print.mst} returns similar (but differently organized) results. 
 }


\details{
 The \code{randomMST} function generates a multistage test using an item bank specified by arguments \code{itemBank} and \code{model}, an MST structure provided by arguments \code{modules} and \code{transMatrix}, and for a given true ability level specified by argument \code{trueTheta}. 

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). See \code{\link{genDichoMatrix}} for further information.

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The \code{modules} argument must be a binary 0/1 matrix with as many rows as the item bank \code{itemBank} and as many columns as the number of modules. Values of 1 indicate to which module(s) the items belong to, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the \emph{i}-th item belongs to the \emph{j}-th module.

The \code{transMatrix} argument must be a binary 0/1 square matrix with as many rows (and columns) as the number of modules. All values of 1 indicate the possible transitions from one module to another, i.e. a value of 1 on row \emph{i} and column \emph{j} means that the MST can move from \emph{i}-th module to \emph{j}-th module.

By default all item responses will be randomly drawn from parent distribution set by the item bank parameters of the \code{itemBank} matrix (using the \code{\link{genPattern}} function for instance). Moreover, the random generation of the item responses can be fixed (for e.g., replication purposes) by assigning some numeric value to the \code{genSeed} argument. By default this argument is equal to \code{NULL} so the random seed is not fixed (and two successive runs of \code{randomMST} will usually lead to different response patterns).

It is possible, however, to provide a full response pattern of previously recorded responses to each item of the item bank, for instance for post-hoc simulations. This is done by providing to the \code{responses} argument a vector of binary entries (without missing values). By default \code{responses} is set to \code{NULL} and item responses will be drawn from the item bank parameters.

 The test specification is made by means of three lists of options: one list for the selection of the starting module, 
 one list with the options for provisional ability estimation and next module selection, and one list with the options for final ability estimation. These lists are specified respectively by the arguments \code{start}, \code{test} and \code{final}. 

 The \code{start} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{fixModule}: either an integer value, setting the module to be administered as first stage (as its row number in te \code{transMatrix} argument for instance), or \code{NULL} (default) to let the function select the module.}
 \item{\code{seed}: either a numeric value to fix the random seed for module selection, \code{NA} to randomly select the module withour fixing the random seed, or \code{NULL} (default) to make random module selection without fixing the random seed. Ignored if \code{fixModule} is not \code{NULL}.}
 \item{\code{theta}: the initial ability value, used to select the most informative module at this ability level (default is 0). Ignored if either \code{fixModule} or \code{seed} is not \code{NULL}. See \code{startModule} for further details.}
 \item{\code{D}: numeric, the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL} and if \code{startSelect} is not \code{"MFI"}.} 
 }
 These arguments are passed to the function \code{startModule} to select the first module of the multistage test.

The \code{test} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{method}: a character string to specify the method for ability estimation or scoring. Possible values are: \code{"BM"} (default) for Bayesian modal estimation (Birnbaum, 1969), \code{"ML"} for maximum likelihood estimation (Lord, 1980), \code{"EAP"} for expected a posteriori (EAP) estimation (Bock and Mislevy, 1982), and \code{"WL"} for weighted likelihood estimation (Warm, 1989). The \code{method} argument can also take the value \code{"score"}, meaning that module selection is based on the test score from the previously administered modules. The latter works only if \code{cutoff} argument is supplied appropriately, otherwise this leads to an error message.}
 \item{\code{priorDist}: a character string which sets the prior distribution. Possible values are: \code{"norm"} (default) for normal distribution, \code{"unif"} for uniform distribution, and \code{"Jeffreys"} for Jeffreys' noninformative prior distribution (Jeffreys, 1939, 1946). ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}.}
 \item{\code{priorPar}: a vector of two numeric components, which sets the parameters of the prior distribution. If (\code{method="BM"} or \code{method=="EAP"}) and \code{priorDist="norm"}, the components of \code{priorPar} are respectively the mean and the standard deviation of the prior normal density. If (\code{method="BM"} or \code{method="EAP"}) and \code{priorDist="unif"}, the components of \code{priorPar} are respectively the lower and upper bound of the prior uniform density. Ignored in all other cases. By default, \code{priorPar} takes the parameters of the prior standard normal distribution (i.e., \code{priorPar=c(0,1)}). In addition, \code{priorPar} also provides the prior parameters for the comoutation of MLWI and MPWI values for next item selection (see \code{\link{nextModule}} for further details).}
 \item{\code{range}: the maximal range of ability levels, set as a vector of two numeric components. The ability estimate will always lie to this interval (set by default to [-4, 4]). Ignored if \code{method=="EAP"}.}
 \item{\code{D}: the value of the metric constant. Default is \code{D=1} for logistic metric. Setting \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.}
 \item{\code{parInt}: a numeric vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the \code{\link{eapEst}}, \code{\link{eapSem}} and \code{MWI} commands. It specifies the range of quadrature points for numerical integration, and is used for computing the EAP estimate, its standard error, and the MLWI and MPWI values for next item selection. Default vector is (-4, 4, 33), thus setting the range from -4 to 4 by steps of 0.25. Ignored if \code{method} is not \code{"EAP"} and if \code{itemSelect} is neither \code{"MLMWI"} nor \code{"MPMWI"}.}
\item{\code{moduleSelect}: the rule for next module selecion, with possible values:
\enumerate{
\item{\code{"MFI"} (default) for maximum Fisher information criterion;}
\item{\code{"MLMWI"} for maximum likelihood weighted (module) information criterion;}
\item{\code{"MPMWI"} for posterior weighted (module) information criterion;}
\item{\code{"MKL"} for (module) Kullback-Leibler information methods;}
\item{\code{"MKLP"} for posterior (module) Kullback-Leibler information methods;}
\item{\code{"random"} for random selection.}
}
This argument is ignored if \code{cutoff} is supplied appropriately. See \code{nextModule} for further details.}
\item{\code{constantPatt}: the method to estimate ability in case of constant pattern (i.e. only correct or only incorrect responses). Can be either \code{NULL} (default), \code{"BM"}, \code{"EAP"}, \code{"WL"}, \code{"fixed4"} (for fixed stepsize adjustment with step 0.4), \code{"fixed7"} (for fixed stepsize adjustment with step 0.7), or \code{"var"} (for variable stepsize adjustment). This is currently implemented only for dichotomous IRT models and is sgnored if \code{method} is \code{"score"}. See \code{\link{thetaEst}} for further details.}
\item{\code{cutoff}: either \code{NULL} (default) or a suitable matrix of thresholds to select the next module. Thresholds can reflect module selection based on ability estimation (and then \code{method} should hold one of the ability estimation methods) or on provisional test score (and then \code{method} must be set to \code{"score"}. See \code{\link{nextModule}} for further details about suitable definition of the \code{cutoff} matrix (and the examples below).}
\item{\code{randomesque}: a probability value to select the optimal module. Default is one, so the optimal module is always chosen. With a value smaler than one, other elligible modules can be selected.}
\item{\code{random.seed}: either \code{NULL} (default) or a numeric value to fix the random seed of randomesque selection of the module. Ignored if \code{randomesque} is equal to one.}
\item{\code{score.range}: a character value that specifies on which set of modules the provisional test score should be computed. Possible values are \code{"all"} (default) to compute the score with all previously administered modules, or \code{"last"} to compute the score only with the last module. Ignored if \code{method} is not \code{"score"}.}

 }
 These arguments are passed to the functions \code{\link{thetaEst}} and \code{\link{semTheta}} to estimate the ability level and the standard error of this estimate. In addition, some arguments are passed to \code{\link{nextModule}} to select the next module appropriately.

 Finally, the \code{final} list can contain the arguments \code{method}, \code{priorDist}, \code{priorPar}, \code{range}, \code{D} and \code{parInt} of the \code{test} list (with possiblly different values), as well as the additional \code{alpha} argument. The latter specifies the \eqn{\alpha} level of the final confidence interval of ability, which is computed as
 \deqn{[\hat{\theta}-z_{1-\alpha/2} \; se(\hat{\theta}) ; \hat{\theta}+z_{1-\alpha/2} \; se(\hat{\theta})]}
 where \eqn{\hat{\theta}} and \eqn{se(\hat{\theta})} are respectively the ability estimate and its standard error. 

 If some arguments of these lists are missing, they are automatically set to their default value. 

Usually the ability estimates and related standard errors are computed right after the full administration of each module (that is, if current module has \emph{k} items, the \emph{(k-1)} ability levels and standard errors from the first administered \emph{(k-1)} are not computed). This can however be avoided by fixing the argument \code{allTheta} to \code{TRUE} (by default it is \code{FALSE}). In this case, all provisional ability estimates (or test scores) and standard errors (or \code{NA}'s) are computed and returned. 

 The output of \code{randomMST}, as displayed by the \code{print.mst} function, can be stored in a text file provided that \code{save.output} is set to \code{TRUE} (the default value \code{FALSE} does not execute the storage). In this case, the \code{(output} argument mus hold three character values: the path to where the output file must be stored, the name of
 the output file, and the type of output file. If the path is not provided (i.e. left to its default value \code{"path"}), it will be saved in the default working directory. The default name is \code{"name"}, and the default type is \code{"csv"}. Any other value yields a text file. See the \bold{Examples} section for an illustration. 

 The function \code{plot.mst} represents the whole MST structure with as many rectangles as there are available modules, arrows connecting all the modules according to the \code{transMatrix} structure. Each stage is displayed as one horizontal layout with stage 1 on the top and final stage at the bottom of the figure. The selected path (i.e. set of modules) is displayed on the plot when \code{show.path} is \code{TRUE} (which is the default value). Modules from the path and arrows between them are then highlighted in red (by default), and these colors can be modified by settong \code{border.col} and \code{arrow.col} arguments with appropriate color names. By default, modules are labelled as ``module 1", ``module 2" etc., the numbering starting from left module to right module and from stage 1 to last stage. These labels can be modified by providing a vector of character names to argument \code{module.names}. This vector must have as many components as the total number of modules and being ordered identically as described above.

Note that the MST structure can be graphically displayed by only providing (as \code{x} argument) the transition matrix of the MST. In this case, \code{show.path} argument is ignored. This is useful to represent the MST structure set by the transition matrix without running an MST simulation.

 Finally, the plot can be saved in an external file, either as PDF or JPEG format. First, the argument \code{save.plot} must be set to \code{TRUE} (default is \code{FALSE}). Then, the file path for figure storage, the name of the figure and its format are specified through the argument \code{save.options}, all as character strings. See the \code{Examples} section for further information and a practical example.

}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability.
 \emph{Journal of Mathematical Psychology, 6}, 258-276. doi: 10.1016/0022-2496(69)90005-4

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444. doi: 10.1177/014662168200600405

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F. M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. doi: 10.1007/BF02294627
}


\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}

Duanli Yan \cr
Educational Testing Service, Princeton, USA \cr
\email{dyan@ets.org}

 }

\seealso{
\code{\link{thetaEst}}, \code{\link{semTheta}}, \code{\link{eapEst}},
\code{\link{eapSem}}, \code{\link{genPattern}}, \code{\link{genDichoMatrix}} , \code{\link{genPolyMatrix}} ,


\code{\link{nextModule}} 
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

 # Creation of the start list: selection by MFI with ability level 0
 start <- list(theta = 0)

 # Creation of the test list: module selection by MFI, ability estimation by WL,
 # stepsize .4 adjustment for constant pattern
 test <- list(method = "WL", moduleSelect = "MFI", constantPatt = "fixed4")

 # Creation of the final list: ability estimation by ML
 final <- list(method = "ML")

 # Random MST generation for true ability level 1 and all ad-interim ability estimates
 res <- randomMST(trueTheta = 1, itemBank = it, modules = modules, transMatrix = trans,
                   start = start, test = test, final = final, allTheta = TRUE) 

 # Module selection by cut-scores for ability estimates
 # Creation of cut-off scores for ability levels: cut score 0 between modules 3 and 4
 # and cut scores -1 and 1 between modules 5, 6 and 7
 # randomesque selection with probability .8
 cut <- rbind(c(3, 4, 0), c(5, 6, -1), c(6, 7, 1))
 test <- list(method = "WL", constantPatt = "fixed4", cutoff = cut, randomesque = 0.8)
 res <- randomMST(trueTheta = 1, itemBank = it, modules = modules, transMatrix = trans,
                   start = start, test = test, final = final, allTheta = TRUE) 

 # Module selection by cut-scores for test scores
 # Creation of cut-off scores for test scores: cut score 4 between modules 3 and 4
 # and cut scores 5 and 9 between modules 5, 6 and 7
 cut.score <- rbind(c(3, 4, 4), c(5, 6, 5), c(6, 7, 9))
 test <- list(method = "score", cutoff = cut.score)
 final <- list(method = "score")
 res <- randomMST(trueTheta = 1, itemBank = it, modules = modules, transMatrix = trans,
                   start = start, test = test, final = final, allTheta = TRUE) 

 # Modification of cut-scores of stage 3 to use only the last module from stage 2 (6 items):
 # cut scores 2 and 4 between modules 5, 6 and 7
 cut.score2 <- rbind(c(3, 4, 4), c(5, 6, 2), c(6, 7, 4))
 test <- list(method = "score", cutoff = cut.score2, score.range = "last")
 final <- list(method = "score")
 res <- randomMST(trueTheta = 1, itemBank = it, modules = modules, transMatrix = trans,
                   start = start, test = test, final = final, allTheta = TRUE) 

 ## Plot options
 plot(trans)
 plot(res)
 plot(res, show.path = FALSE)
 plot(res, border.col = "blue")
 plot(res, arrow.col = "green")

 ## Saving plot to default folder, in PDF format and with name MST
 plot(res, save.plot = TRUE, save.options = c(tempdir(),"MST", "pdf"))
}
