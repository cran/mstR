MKL <- function (itemBank, modules, target.mod, theta=NULL, it.given, x, model = NULL, lower = -4, upper = 4, nqp = 33, 
                type = "MKL", priorDist="norm", priorPar = c(0, 1), D=1) 
{
  if (type != "MKL" & type != "MKLP") 
    stop("'type' must be either 'MKL' or 'MKLP'", call. = FALSE)
if (is.null(theta)) theta<-thetaEst(itemBank[it.given,],x,D=D,model=model, method="ML")
 KLF<-NULL
 par <- rbind(itemBank[modules[,target.mod]==1,])
X<-seq(from=lower,to=upper,length=nqp)
par.given<-itemBank[it.given,]
if (is.null(model)){
L <- function(th, r, param) prod(Pi(th, param,D=D)$Pi^r*(1-Pi(th,param,D=D)$Pi)^(1-r))
lik<-sapply(X,L,x,par.given)
for (t in 1:nqp) KLF[t] <- sum(Pi(theta,par,D=D)$Pi*log(Pi(theta,par,D=D)$Pi/Pi(X[t],par,D=D)$Pi)+(1-Pi(theta,par,D=D)$Pi)*log((1-Pi(theta,par,D=D)$Pi)/(1-Pi(X[t],par,D=D)$Pi)))
crit.value <- lik*KLF 
if (type=="MKLP") {
pd<-switch(priorDist,norm=dnorm(X,priorPar[1],priorPar[2]),unif=dunif(X,priorPar[1],priorPar[2]))
crit.value <- crit.value*pd 
}
}
else{
LL <- function(th, param, r, model,D=1) {
prob <- Pi(th, param, model = model,D=D)$Pi
res <- 1
for (i in 1:length(r)) res <- res * prob[i, r[i] + 1]
return(res)
}
lik<-sapply(X,LL,par.given,x,model=model,D=D)
pi<-Pi(theta,par,model=model,D=D)$Pi
for (i in 1:length(X)){
pri<-Pi(X[i],par,model=model,D=D)$Pi
KLF[i]<-sum(pi*log(pi/pri),na.rm=TRUE)
}
crit.value <- lik*KLF 
if (type=="MKLP") {
pd<-switch(priorDist,norm=dnorm(X,priorPar[1],priorPar[2]),unif=dunif(X,priorPar[1],priorPar[2]))
crit.value <- crit.value*pd 
}
}
RES <- integrate.mstR(X, crit.value)
return(RES)
}
