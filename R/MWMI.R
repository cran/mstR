MWMI<-function (itemBank, modules, target.mod, it.given, x, model = NULL, lower = -4, 
    upper = 4, nqp = 33, type = "MLWMI", priorDist = "norm", priorPar = c(0, 
        1), D = 1) 
{
    if (type != "MLWMI" & type != "MPWMI") 
        stop("'type' must be either 'MLWMI' or 'MPWMI'", call. = FALSE)
    if (is.null(model)) {
        L <- function(th, x, par) prod(Pi(th, par, D = D)$Pi^x * 
            (1 - Pi(th, par, D = D)$Pi)^(1 - x))
        X <- seq(from = lower, to = upper, length = nqp)
        lik <- sapply(X, L, x, itemBank[it.given,])
        items<-which(modules[,target.mod]==1)
        Iprov <- function(t) sum(Ii(t, itemBank, D = D)$Ii[items])
        info <- sapply(X, Iprov)
        crit.value <- lik * info
        if (type == "MPWMI") {
            pd <- NULL
            for (k in 1:length(X)) pd[k] <- switch(priorDist, 
                norm = dnorm(X[k], priorPar[1], priorPar[2]), 
                unif = dunif(X[k], priorPar[1], priorPar[2]))
            crit.value <- crit.value * pd
        }
    }
    else {
        LL <- function(th, it.given, x, model, D = 1) {
            if (dim(it.given)[1] == 0) 
                res <- 1
            else {
                prob <- Pi(th, it.given, model = model, D = D)$Pi
                res <- 1
                for (i in 1:length(x)) res <- res * prob[i, x[i] + 
                  1]
            }
            return(res)
        }
        X <- seq(from = lower, to = upper, length = nqp)
        lik <- sapply(X, LL, itemBank[it.given,], x, model, D = D)
items<-which(modules[,target.mod]==1)
        Iprov <- function(t) sum(Ii(t, itemBank[items, ], model = model, 
            D = D)$Ii)
        info <- sapply(X, Iprov)
        crit.value <- lik * info
        if (type == "MPWMI") {
            pd <- NULL
            for (k in 1:length(X)) pd[k] <- switch(priorDist, 
                norm = dnorm(X[k], priorPar[1], priorPar[2]), 
                unif = dunif(X[k], priorPar[1], priorPar[2]))
            crit.value <- crit.value * pd
        }
    }
    RES <- integrate.mstR(X, crit.value)
    return(RES)
}
