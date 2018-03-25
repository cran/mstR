nextModule<-function (itemBank, modules, transMatrix, model = NULL, current.module, 
    out, x = NULL, cutoff = NULL, theta = 0, criterion = "MFI", 
    priorDist = "norm", priorPar = c(0, 1), D = 1, range = c(-4, 
        4), parInt = c(-4, 4, 33),randomesque=1,random.seed = NULL) 
{
    crit <- switch(criterion, MFI = "MFI", MLWMI = "MLWMI", MPWMI = "MPWMI", 
        MKL = "MKL", MKLP = "MKLP", random = "random")
    if (is.null(cutoff) & is.null(crit)) 
        stop("invalid 'criterion' name", call. = FALSE)
    if (is.null(cutoff) & !is.null(model)) {
        mod <- switch(model, GRM = 1, MGRM = 2, PCM = 3, GPCM = 4, 
            RSM = 5, NRM = 6)
        if (is.null(mod)) 
            stop("invalid 'model' type!", call. = FALSE)
    }
    pot.mods <- which(transMatrix[current.module, ] == 1)
    sel.stage <- NULL
    for (i in 1:length(pot.mods)) {
        items <- which(modules[, pot.mods[i]] == 1)
        if (length(out) + length(items) == length(unique(c(out, 
            items)))) 
            sel.stage <- c(sel.stage, pot.mods[i])
    }
    if (is.null(sel.stage)) 
        stop("No available module without overlap with administered items", 
            call. = FALSE)
    if (!is.null(cutoff)) {
        thr <- NULL
        for (i in 1:(length(sel.stage) - 1)) thr <- c(thr, cutoff[cutoff[, 
            1] == sel.stage[i] & cutoff[, 2] == sel.stage[i + 
            1], 3])
        thr <- c(-Inf, thr, Inf)
        if (sum(thr == theta) == 1) {
            ind <- which(thr == theta)
        }
        else {
            n <- length(thr)
            ind <- which(thr[1:(n - 1)] < theta & thr[2:n] > 
                theta)
        }
if (length(sel.stage)>1){
probs<-rep((1-randomesque)/(length(sel.stage)-1),length(sel.stage))
probs[ind]<-randomesque
if (!is.null(random.seed)) set.seed(random.seed)
ind<-which(c(rmultinom(1,1,probs))==1)
}
        final.module <- sel.stage[ind]
        select <- which(modules[, final.module] == 1)
        res <- list(module = final.module, items = select, par = itemBank[select, 
            ], info = theta, criterion = "cutoff")
    }
    else {
        if (criterion == "MFI") {
            infos <- NULL
            for (i in 1:length(sel.stage)) {
                items <- which(modules[, sel.stage[i]] == 1)
                infos[i] <- sum(Ii(theta, itemBank[items, ], 
                  model = model, D = D)$Ii)
            }
            maxinfo <- which(infos == max(infos))
            if (length(maxinfo) > 1) 
                maxinfo <- sample(maxinfo, 1)
if (length(sel.stage)>1){
probs<-rep((1-randomesque)/(length(sel.stage)-1),length(sel.stage))
probs[maxinfo]<-randomesque
if (!is.null(random.seed)) set.seed(random.seed)
maxinfo<-which(c(rmultinom(1,1,probs))==1)
}
            final.module <- sel.stage[maxinfo]
            select <- which(modules[, final.module] == 1)
            res <- list(module = final.module, items = select, 
                par = itemBank[select, ], info = max(infos), 
                criterion = "MFI")
        }
        if (criterion == "MLWMI" | criterion == "MPWMI") {
            infos <- NULL
            for (i in 1:length(sel.stage)) {
                infos[i] <- MWMI(itemBank, modules, target.mod = sel.stage[i], 
                  it.given = out, x = x, lower = parInt[1], upper = parInt[2], 
                  nqp = parInt[3], type = criterion, priorDist = priorDist, 
                  priorPar = priorPar, D = D)
            }
            maxinfo <- which(infos == max(infos))
            if (length(maxinfo) > 1) 
                maxinfo <- sample(maxinfo, 1)
if (length(sel.stage)>1){
probs<-rep((1-randomesque)/(length(sel.stage)-1),length(sel.stage))
probs[maxinfo]<-randomesque
if (!is.null(random.seed)) set.seed(random.seed)
maxinfo<-which(c(rmultinom(1,1,probs))==1)
}
            final.module <- sel.stage[maxinfo]
            select <- which(modules[, final.module] == 1)
            res <- list(module = final.module, items = select, 
                par = itemBank[select, ], info = max(infos), 
                criterion = criterion)
        }
        if (criterion == "MKL" | criterion == "MKLP") {
            infos <- NULL
            for (i in 1:length(sel.stage)) {
                infos[i] <- MKL(itemBank, modules, target.mod = sel.stage[i], 
                  it.given = out, x = x, theta = theta, lower = parInt[1], 
                  upper = parInt[2], nqp = parInt[3], type = criterion, 
                  priorDist = priorDist, priorPar = priorPar, 
                  D = D)
            }
            maxinfo <- which(infos == max(infos))
            if (length(maxinfo) > 1) 
                maxinfo <- sample(maxinfo, 1)
if (length(sel.stage)>1){
probs<-rep((1-randomesque)/(length(sel.stage)-1),length(sel.stage))
probs[maxinfo]<-randomesque
if (!is.null(random.seed)) set.seed(random.seed)
maxinfo<-which(c(rmultinom(1,1,probs))==1)
}
            final.module <- sel.stage[maxinfo]
            select <- which(modules[, final.module] == 1)
            res <- list(module = final.module, items = select, 
                par = itemBank[select, ], info = max(infos), 
                criterion = criterion)
        }
        if (criterion == "random") {
            final.module <- sample(sel.stage, 1)
            select <- which(modules[, final.module] == 1)
            res <- list(module = final.module, items = select, 
                par = itemBank[select, ], info = NA, criterion = "random")
        }
    }
set.seed(NULL)
    return(res)
}
