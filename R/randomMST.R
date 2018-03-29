randomMST<-function (trueTheta=NULL, itemBank, modules, transMatrix, model = NULL, responses = NULL, 
    genSeed = NULL, start = list(fixModule = NULL, seed = NULL, theta = 0, 
        D = 1), test = list(method = "BM", 
        priorDist = "norm", priorPar = c(0, 1), range = c(-4, 
            4), D = 1, parInt = c(-4, 4, 33), moduleSelect = "MFI", 
            constantPatt = NULL, cutoff=NULL,randomesque=1,random.seed=NULL,
score.range="all"), 
    final = list(method = "BM", 
        priorDist = "norm", priorPar = c(0, 1), range = c(-4, 
            4), D = 1, parInt = c(-4, 4, 33), alpha = 0.05), 
    allTheta = FALSE, save.output = FALSE, output = c("path", 
        "name", "csv")) 
{
if (is.null(trueTheta) & is.null(responses)) stop("Either 'trueTheta' or 'responses' argument must be supplied",call.=FALSE)
if (!testListMST(start, type = "start")$test) 
        stop(testListMST(start, type = "start")$message, call. = FALSE)
    if (!testListMST(test, type = "test")$test) 
        stop(testListMST(test, type = "test")$message, call. = FALSE)
    if (!testListMST(final, type = "final")$test) 
        stop(testListMST(final, type = "final")$message, call. = FALSE)
if (!is.null(responses)) 
        assigned.responses <- TRUE
    else assigned.responses <- FALSE
  
    internalMST <- function() {
        startList <- list(fixModule = start$fixModule, seed = start$seed, 
            theta = 0, D = 1)
        startList$theta <- ifelse(is.null(start$theta), 0, start$theta)
        startList$D <- ifelse(is.null(start$D), 1, start$D)
        start <- startList
        testList <- list(method = NULL, priorDist = NULL, priorPar = c(0, 
            1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
            moduleSelect = "MFI", constantPatt = NULL, cutoff=NULL,randomesque=1,
            random.seed=NULL,score.range="all")
        testList$method <- ifelse(is.null(test$method), "BM", 
            test$method)
        testList$priorDist <- ifelse(is.null(test$priorDist), 
            "norm", test$priorDist)
        if (!is.null(test$priorPar)) {
            testList$priorPar[1] <- test$priorPar[1]
            testList$priorPar[2] <- test$priorPar[2]
        }
        if (!is.null(test$range)) {
            testList$range[1] <- test$range[1]
            testList$range[2] <- test$range[2]
        }
        testList$D <- ifelse(is.null(test$D), 1, test$D)
        if (!is.null(test$parInt)) {
            testList$parInt[1] <- test$parInt[1]
            testList$parInt[2] <- test$parInt[2]
            testList$parInt[3] <- test$parInt[3]
        }
        testList$moduleSelect <- ifelse(is.null(test$moduleSelect), 
            "MFI", test$moduleSelect)
              testList$AP <- ifelse(is.null(test$AP), 1, test$AP)
        if (!is.null(test$constantPatt)) 
            testList$constantPatt <- test$constantPatt
if (!is.null(test$cutoff)) 
            testList$cutoff<- test$cutoff
if (!is.null(test$randomesque)) 
            testList$randomesque<- test$randomesque
if (!is.null(test$random.seed)) 
            testList$random.seed<- test$random.seed
if (!is.null(test$score.range)) 
            testList$score.range<- test$score.range
        test <- testList
        finalList <- list(method = NULL, priorDist = NULL, priorPar = c(0, 
            1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
            alpha = 0.05)
        finalList$method <- ifelse(is.null(final$method), "BM", 
            final$method)
        finalList$priorDist <- ifelse(is.null(final$priorDist), 
            "norm", final$priorDist)
        if (is.null(final$priorPar) == FALSE) {
            finalList$priorPar[1] <- final$priorPar[1]
            finalList$priorPar[2] <- final$priorPar[2]
        }
        if (!is.null(final$range)) {
            finalList$range[1] <- final$range[1]
            finalList$range[2] <- final$range[2]
        }
        finalList$D <- ifelse(is.null(final$D), 1, final$D)
        if (!is.null(final$parInt)) {
            finalList$parInt[1] <- final$parInt[1]
            finalList$parInt[2] <- final$parInt[2]
            finalList$parInt[3] <- final$parInt[3]
        }
        finalList$alpha <- ifelse(is.null(final$alpha), 0.05, 
            final$alpha)
        final <- finalList

if (test$method=="score" & is.null(test$cutoff) & test$moduleSelect!="random") stop("'cutoff' argument of 'test' list must be supplied when 'method' is 'score'",call.=FALSE)
        pr0 <- startModule(itemBank = itemBank, modules=modules, transMatrix=transMatrix,
                model = model,  fixModule = start$fixModule, seed = start$seed,
                theta = start$theta, D = start$D)
        ITEMS <- pr0$items
        ITEMS.PER.MOD<-length(ITEMS)
        PAR <- rbind(pr0$par)
MODULE<-pr0$module
            if (!is.null(responses)) 
                PATTERN <- responses[ITEMS]
            else PATTERN <- genPattern(trueTheta, PAR, model = model, 
                D = test$D, seed = genSeed)
        if (test$method!="score") 
            TH <- thetaEst(PAR, PATTERN, model = model, D = test$D, 
                method = test$method, priorDist = test$priorDist, 
                priorPar = test$priorPar, range = test$range, 
                parInt = test$parInt, current.th = start$theta, 
                constantPatt = test$constantPatt, bRange = range(itemBank[, 
                  2]))
        else TH <- sum(PATTERN)
        if (test$method!="score") 
            SETH <- semTheta(TH, PAR, x = PATTERN, model = model, 
                D = test$D, method = test$method, priorDist = test$priorDist, 
                priorPar = test$priorPar, parInt = test$parInt, 
                constantPatt = test$constantPatt)
        else SETH <- NA
        thProv <- TH
        enter.mst <- TRUE
best.mod<-TRUE
if (sum(transMatrix[MODULE,])==0) enter.mst<-FALSE
        if (!enter.mst) {
if (final$method!="score"){
            finalEst <- thetaEst(PAR, PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, range = final$range, 
                parInt = final$parInt)
            seFinal <- semTheta(finalEst, PAR, x = PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, parInt = final$parInt)
            confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                seFinal)
}
else{
finalEst<-sum(PATTERN)
seFinal<-NA
confIntFinal<-c(NA,NA)
}
            RES <- list(trueTheta = trueTheta, selected.modules=MODULE, items.per.module=ITEMS.PER.MOD,
                transMatrix=transMatrix,
                model = model, testItems = ITEMS, itemPar = PAR, 
                pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                genSeed = genSeed, startFixModule = start$fixModule, 
                startSeed = start$seed, startTheta = start$theta, startD = start$D, 
                startThStart = pr0$thStart, startSelect = start$startSelect, 
                provMethod = test$method, provDist = test$priorDist, 
                provPar = test$priorPar, provRange = test$range, 
                provD = test$D, moduleSelect = test$moduleSelect, 
                constantPattern = test$constantPatt, cutoff=test$cutoff,
                randomesque=test$randomesque,random.seed=test$random.seed,
                score.range=test$score.range,best.module=best.mod,
                finalMethod = final$method, finalDist = final$priorDist, 
                finalPar = final$priorPar, finalRange = final$range, 
                finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                output = output,allTheta=NULL,assigned.responses=assigned.responses)
        }
        else {
            repeat {
                pr <- nextModule(itemBank, modules=modules, transMatrix=transMatrix,
                  model = model, current.module= MODULE[length(MODULE)],theta = thProv, 
                  out = ITEMS, x = PATTERN, cutoff=test$cutoff, criterion = test$moduleSelect, 
                  parInt = test$parInt, 
                  priorDist = test$priorDist, priorPar = test$priorPar, 
                  D = test$D, range = test$range,randomesque=test$randomesque,
                  random.seed=test$random.seed)
                ITEMS <- c(ITEMS, pr$items)
best.mod<-c(best.mod,pr$best.module)
ITEMS.PER.MOD<-c(ITEMS.PER.MOD,length(pr$items))
                PAR <- rbind(PAR, pr$par)
                if (!is.null(responses)) 
                  PATTERN <- c(PATTERN, responses[pr$items])
                else PATTERN <- c(PATTERN, genPattern(trueTheta, 
                  pr$par, model = model, D = test$D, seed = genSeed))
if (test$method!="score")
                thProv <- thetaEst(PAR, PATTERN, model = model, 
                  D = test$D, method = test$method, priorDist = test$priorDist, 
                  priorPar = test$priorPar, range = test$range, 
                  parInt = test$parInt, current.th = TH[length(TH)], 
                  constantPatt = test$constantPatt, bRange = range(itemBank[,2]))
else {
if (test$score.range=="all") thProv<-sum(PATTERN)
else{
nr<-length(pr$items)
thProv<-sum(PATTERN[(length(PATTERN)-nr+1):(length(PATTERN))])
}
}
                TH <- c(TH, thProv)
if (test$method!="score")
                seProv <- semTheta(thProv, PAR, x = PATTERN, 
                  model = model, D = test$D, method = test$method, 
                  priorDist = test$priorDist, priorPar = test$priorPar, 
                  parInt = test$parInt, constantPatt = test$constantPatt)
else seProv<-NA
                SETH <- c(SETH, seProv)
MODULE<-c(MODULE,pr$module)
if (sum(transMatrix[MODULE[length(MODULE)],])==0) break
            }
if (final$method!="score"){
            finalEst <- thetaEst(PAR, PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, range = final$range, 
                parInt = final$parInt,current.th = TH[length(TH)], 
                  constantPatt = test$constantPatt, bRange = range(itemBank[,2]))
            seFinal <- semTheta(finalEst, PAR, x = PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, parInt = final$parInt)
            confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                seFinal)
}
else{
finalEst<-sum(PATTERN)
seFinal<-NA
confIntFinal<-c(NA,NA)
}
            RES <- list(trueTheta = trueTheta, selected.modules=MODULE, items.per.module=ITEMS.PER.MOD,
transMatrix=transMatrix,
                model = model, testItems = ITEMS, itemPar = PAR, 
                pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                genSeed = genSeed, startFixModule = start$fixModule, 
                startSeed = start$seed, startTheta = start$theta, 
                startD = start$D, 
                provMethod = test$method, provDist = test$priorDist, 
                provPar = test$priorPar, provRange = test$range, 
                provD = test$D, moduleSelect = test$moduleSelect, 
                constantPattern = test$constantPatt, cutoff=test$cutoff,
                randomesque=test$randomesque,random.seed=test$random.seed,
                score.range=test$score.range,best.module=best.mod,
                finalMethod = final$method, finalDist = final$priorDist, 
                finalPar = final$priorPar, finalRange = final$range, 
                finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                output = output,allTheta=NULL,assigned.responses=assigned.responses)
        }
        if (allTheta) {
           prov.th <- prov.se <- NULL
                for (k in 1:nrow(RES$itemPar)) {
if (test$method!="score"){
                  prov.par <- rbind(RES$itemPar[1:k, ])
                  prov.th[k] <- thetaEst(prov.par, RES$pattern[1:k], 
                    model = model, D = test$D, method = test$method, 
                    priorDist = test$priorDist, priorPar = test$priorPar, 
                    range = test$range, parInt = test$parInt, 
                    constantPatt = test$constantPatt, bRange = range(itemBank[,2]))
                  prov.se[k] <- semTheta(prov.th[k], prov.par, 
                    RES$pattern[1:k], model = model, D = test$D, 
                    method = test$method, priorDist = test$priorDist, 
                    priorPar = test$priorPar, parInt = test$parInt, 
                    constantPatt = test$constantPatt)
}
else{
prov.th[k]<-sum(RES$pattern[1:k])
prov.se[k]<-NA
}
                }
                RES$allTheta <- cbind(prov.th, prov.se)
                colnames(RES$allTheta) <- c("th","se")
        }
        class(RES) <- "mst"
        return(RES)
    }
    resToReturn <- internalMST()
    if (save.output) {
        if (output[1] == "path") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- output[1]
        if (output[3] == "csv") 
            fileName <- paste(wd, output[2], ".csv", sep = "")
        else fileName <- paste(wd, output[2], ".txt", sep = "")
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}


########

print.mst<-function (x, ...) 
{
    if (!x$assigned.responses) {
        cat("Random generation of a MST response pattern", "\n")
        if (is.null(x$genSeed)) 
            cat("  without fixing the random seed", "\n", "\n")
        else cat("  with random seed equal to", x$genSeed, "\n", 
            "\n")
    }
    else cat("Post-hoc simulation of a MST response pattern", 
        "\n", "\n")
    if (is.null(x$model)) {
        if (min(x$itemPar[, 4]) < 1) 
            mod <- "Four-Parameter Logistic model"
        else {
            if (max(x$itemPar[, 3]) > 0) 
                mod <- "Three-Parameter Logistic model"
            else {
                if (length(unique(x$itemPar[, 1])) > 1) 
                  mod <- "Two-Parameter Logistic model"
                else mod <- "One-Parameter Logistic (Rasch) model"
            }
        }
    }
    else {
        if (x$model == "GRM") 
            mod <- "Graded Response Model"
        if (x$model == "MGRM") 
            mod <- "Modified Graded Response Model"
        if (x$model == "PCM") 
            mod <- "Partial Credit Model"
        if (x$model == "GPCM") 
            mod <- "Generalized Partial Credit Model"
        if (x$model == "RSM") 
            mod <- "Rating Scale Model"
        if (x$model == "NRM") 
            mod <- "Nominal Response Model"
    }
    cat(" Item bank calibrated under", mod, "\n", "\n")
    if (!is.null(x$trueTheta)) {
if (x$assigned.responses)
        cat(" True ability level:", round(x$trueTheta, 2), "(not used for post-hoc simulation)", "\n", 
            "\n")
else   cat(" True ability level:", round(x$trueTheta, 2), "\n", 
            "\n")
}
    else cat(" True ability level was not provided", "\n", "\n")

cat(" MST structure:", "\n")
cat("   Number of stages:", length(x$selected.modules),"\n")
# extraction of the number of modules per stages 
nr<-0
tr<-x$transMatrix
nr.st<-NULL
repeat{
nr<-nr+1
ind<-which(colSums(tr)>0)
nr.st[nr]<-ncol(tr)-length(ind)
tr<-tr[ind,ind]
if (sum(tr)==0) break
}
nr.st<-c(nr.st,ncol(tr))
struc<-NULL
for (i in 1:(length(nr.st)-1)) struc<-paste(struc,nr.st[i],"-",sep="")
struc<-paste(struc,nr.st[length(nr.st)],sep="")
cat("   Structure (number of modules per stage):", struc, "\n", "\n")

    cat(" Starting parameters:", "\n")
cat("   Number of available modules at first stage:", sum(colSums(x$transMatrix)==0),"\n")
if (!is.null(x$startFixModule)){
cat("   Selection of the first stage module: chosen by administrator","\n")
cat("   Selected module:",x$startFixModule,"\n")
}
else{
if (!is.null(x$startSeed)) cat("   Selection of the first stage module: by random selection","\n")
else {
cat("   Selection of the first stage module: by maximizing module information","\n")
cat("     for starting ability","\n")
cat("   Starting ability level:", round(x$startTheta,3), "\n")
}
}
cat("\n", "Multistage test parameters:", "\n")
if (is.null(x$cutoff)) {
    itemSel <- switch(x$moduleSelect, MFI = "maximum Fisher information", 
        MLWMI = "Maximum likelihood weighted information (MLWI)", 
        MPWMI = "Maximum posterior weighted information (MPWI)", 
        random = "Random selection", 
        MKL = "Kullback-Leibler (KL) information", 
        MKLP = "Posterior Kullback-Leibler (KLP) information") 
            cat("   Next module selection:", itemSel, "\n")
    if (x$moduleSelect == "MKLP" | x$moduleSelect == "MPWI") {
        met3 <- switch(x$provDist, norm = paste("N(", round(x$provPar[1], 
            2), ",", round(x$provPar[2]^2, 2), ") prior", sep = ""), 
            unif = paste("U(", round(x$provPar[1], 2), ",", round(x$provPar[2], 
                2), ") prior", sep = ""), Jeffreys = "Jeffreys' prior")
        cat("     Prior ability distribution for", ifelse(x$moduleSelect=="MKLP", "KLP", "MPWI"), 
            "method:", met3, "\n")
    }
}
else {
cat("   Next module selection: by pre-specified cut scores","\n")
cat("     (random selection among allowed modules)","\n")
}
    met2 <- switch(x$provMethod, BM = "Bayes modal (MAP) estimator", 
        WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
        EAP = "Expected a posteriori (EAP) estimator",score="Test score (sum-score) computation")
    if (x$provMethod == "BM" | x$provMethod == "EAP") {
        met3 <- switch(x$provDist, norm = paste("N(", round(x$provPar[1], 
            2), ",", round(x$provPar[2]^2, 2), ") prior", sep = ""), 
            unif = paste("U(", round(x$provPar[1], 2), ",", round(x$provPar[2], 
                2), ") prior", sep = ""), Jeffreys = "Jeffreys' prior")
    }
    if (x$provMethod == "ML") 
        ra1 <- paste("[", round(x$provRange[1], 2), ",", round(x$provRange[2], 
            2), "]", sep = "")
    cat("   Provisional ability estimator:", met2, "\n")
    if (x$provMethod == "BM" | x$provMethod == "EAP") 
        cat("     Provisional prior ability distribution:", met3, 
            "\n")
    if (x$provMethod == "ML") 
        cat("   Provisional range of ability values:", ra1, "\n")
    if (x$provMethod == "score") {
if (x$score.range=="all") type1<-"all previous modules)"
else type1<-"last module only)"
        cat("     (provisional score computed on", type1, "\n")
}
    if (!is.null(x$model) | is.null(x$constantPattern) | x$provMethod=="score") 
        adj <- "none"
    else adj <- switch(x$constantPattern, fixed4 = "fixed .4 stepsize", 
        fixed7 = "fixed .7 stepsize", var = "variable stepsize")
    cat("   Ability estimation adjustment for constant pattern:", 
        adj, "\n")
    if (x$randomesque==1)  cat("   Randomesque selection of optimal module: no","\n")
else {
cat("   Randomesque selection of optimal module: yes","\n")
cat("     Probability to select optimal module:",x$randomesque,"\n")
}
    cat("\n", "Multistage test details:", "\n")
binf<-1
for (co in 1:length(x$selected.modules)){
cat("\n","  Stage ",co,":","\n",sep="")
cat("    Module administered:",x$selected.modules[co],"\n")
cat("    Number of items in module ",x$selected.modules[co],": ",x$items.per.module[co]," items","\n",sep="")
boo<-ifelse(x$best.module[co],"yes","no")
cat("    Optimal module:",boo,"\n")
cat("    Items and responses:","\n")
bsup<-binf+x$items.per.module[co]-1
its<-x$testItems[binf:bsup]
mat<-rbind(as.character(1:x$items.per.module[co]),as.character(its),as.character(x$pattern[binf:bsup]))
rownames(mat)<-c("Nr", "Item", "Resp.")
colnames(mat)<-rep("",ncol(mat))
if (!is.null(x$allTheta)){
mat<-rbind(mat,as.character(round(x$allTheta[binf:bsup,1],3)),as.character(round(x$allTheta[binf:bsup,2],3)))
rownames(mat)[4:5]<-c("Est.","SE")
for (tt in 1:(ncol(mat)-1)) mat[4,tt]<-paste("(",mat[4,tt],")",sep="")
for (tt in 1:(ncol(mat)-1)) mat[5,tt]<-paste("(",mat[5,tt],")",sep="")
}
binf<-bsup+1
print(format(mat, justify = "right"), quote = FALSE)
cat("\n","    Provisional ability estimate (SE) after stage ",co,": ",round(x$thetaProv[co],3)," (",round(x$seProv[co],3),")","\n",sep="")
}
cat("\n")
    cat(" Final results:", "\n")
    met <- switch(x$finalMethod, BM = "Bayes modal (MAP) estimator", 
        WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
        EAP = "Expected a posteriori (EAP) estimator", score="Total test score")
    if (x$finalMethod == "BM" | x$finalMethod == "EAP") {
        met2 <- switch(x$finalDist, norm = paste("N(", round(x$finalPar[1], 
            2), ",", round(x$finalPar[2]^2, 2), ") prior", sep = ""), 
            unif = paste("U(", round(x$finalPar[1], 2), ",", 
                round(x$finalPar[2], 2), ") prior", sep = ""), 
            Jeffreys = "Jeffreys' prior")
    }
    if (x$finalMethod == "ML") 
        ra1 <- paste("[", round(x$finalRange[1], 2), ",", round(x$finalRange[2], 
            2), "]", sep = "")
    cat("   Total length of multistage test:", length(x$testItems), "items", 
        "\n")
    cat("   Final ability estimator:", met, "\n")
    if (x$finalMethod == "BM" | x$finalMethod == "EAP") 
        cat("   Final prior distribution:", met2, "\n")
    if (x$finalMethod == "ML") 
        cat("   Final range of ability values:", ra1, "\n")
    cat("   Final ability estimate (SE):", round(x$thFinal, 3), 
        paste("(", round(x$seFinal, 3), ")", sep = ""), "\n")
    if (x$finalMethod!="score") cat(paste("   ", (1 - x$finalAlpha) * 100, "% confidence interval: [", 
        round(x$ciFinal[1], 3), ",", round(x$ciFinal[2], 3), 
        "]", sep = ""), "\n")
    if (!x$save.output) 
        cat("\n", "Output was not captured!", "\n")
    else {
        if (x$output[1] == "path") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- x$output[1]
        if (x$output[3] == "csv") 
            fileName <- paste(wd, x$output[2], ".csv", sep = "")
        else fileName <- paste(wd, x$output[2], ".txt", sep = "")
        cat("\n", "Output was captured and saved into file", 
            "\n", " '", fileName, "'", "\n", "on ", as.character(Sys.Date()), 
            "\n", "\n", sep = "")
    }
}


#####

plot.mst<-function (x, show.path = TRUE, border.col = "red", arrow.col = "red", 
    module.names = NULL, save.plot = FALSE, save.options = c("path", 
        "name", "pdf"), ...) 
{
    internalMST <- function() {
        nr <- 0
if (class(x)!="mst" & class(x)!="matrix") stop("'x' must be a list of class 'mst' or a transition matrix",call.=FALSE)
  if (class(x)=="matrix") {
x<-list(transMatrix=x)
SHOW.path<-FALSE
}
else SHOW.path<-show.path
tr <- x$transMatrix
        nr.st <- NULL
        repeat {
            nr <- nr + 1
            ind <- which(colSums(tr) > 0)
            nr.st[nr] <- ncol(tr) - length(ind)
            tr <- tr[ind, ind]
            if (sum(tr) == 0) 
                break
        }
        nr.st <- c(nr.st, ncol(tr))
        height <- 2 * length(nr.st) + 3 * (length(nr.st) - 1)
        width <- 2 * (max(nr.st) * 2 - 1)
        xl <- c(-width/2, width/2)
        yl <- c(-height/2, height/2)
        plot(0, 0, xlim = xl, ylim = yl, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "", bty = "n", col = "white")
        xcenter <- ycenter <- NULL
        allleft <- allright <- allup <- alldown <- NULL
        for (NR in 1:length(nr.st)) {
            up <- yl[2] - 5 * (NR - 1)
            down <- up - 2
            left <- seq(from = -nr.st[NR] * 2 + 1, length = nr.st[NR], 
                by = 4)
            right <- left + 2
            rect(left, down, right, up)
            xcenter <- c(xcenter, (left + right)/2)
            ycenter <- c(ycenter, rep((up + down)/2, nr.st[NR]))
            allleft <- c(allleft, left)
            allright <- c(allright, right)
            allup <- c(allup, rep(up, nr.st[NR]))
            alldown <- c(alldown, rep(down, nr.st[NR]))
        }
        for (i in 1:nrow(x$transMatrix)) {
            for (j in 1:ncol(x$transMatrix)) {
                if (x$transMatrix[i, j] == 1) 
                  arrows(xcenter[i], ycenter[i] - 1.1, xcenter[j], 
                    ycenter[j] + 1.1, length = 0.1, angle = 20)
            }
        }
        for (i in 1:length(xcenter)) {
            if (is.null(module.names)) 
                text(xcenter[i], ycenter[i], paste("Module", 
                  i))
            else text(xcenter[i], ycenter[i], module.names[i])
        }
        if (SHOW.path) {
            for (i in 1:length(x$selected.modules)) {
                ind <- x$selected.modules[i]
                rect(allleft[ind], alldown[ind], allright[ind], 
                  allup[ind], lwd = 2, border = border.col)
            }
            for (i in 1:(length(x$selected.modules) - 1)) {
                ind <- x$selected.modules[i]
                ind2 <- x$selected.modules[i + 1]
                arrows(xcenter[ind], ycenter[ind] - 1.1, xcenter[ind2], 
                  ycenter[ind2] + 1.1, length = 0.1, angle = 20, 
                  lwd = 2, col = arrow.col)
            }
        }
    }
    internalMST()
    if (save.plot) {
        plotype <- NULL
        if (save.options[3] == "pdf") 
            plotype <- 1
        if (save.options[3] == "jpeg") 
            plotype <- 2
        if (is.null(plotype)) 
            cat("Invalid plot type (should be either 'pdf' or 'jpeg').", 
                "\n", "The plot was not captured!", "\n")
        else {
            if (save.options[1] == "path") 
                wd <- paste(getwd(), "/", sep = "")
            else wd <- save.options[1]
            nameFile <- paste(wd, save.options[2], switch(plotype, 
                `1` = ".pdf", `2` = ".jpg"), sep = "")
            if (plotype == 1) {
                {
                  pdf(file = nameFile)
                  internalMST()
                }
                dev.off()
            }
            if (plotype == 2) {
                {
                  jpeg(filename = nameFile)
                  internalMST()
                }
                dev.off()
            }
            cat("The plot was captured and saved into", "\n", 
                " '", nameFile, "'", "\n", "\n", sep = "")
        }
    }
    else cat("The plot was not captured!", "\n", sep = "")
}