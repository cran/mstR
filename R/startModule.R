
startModule<-function (itemBank, modules, transMatrix, model = NULL, fixModule = NULL, seed = NULL,
 theta = 0, D = 1) 
{
if (!is.null(fixModule)) {
if (sum(transMatrix[,fixModule])>0) stop("Selected module is not from stage 1!",call.=FALSE)
items<-which(modules[,fixModule]==1)
par <- itemBank[items, ]
thStart <- NA
res <- list(module=fixModule,items = items, par = par, thStart = thStart)
}
else{
if (!is.null(seed)){
if (!is.na(seed)) set.seed(seed)
mod<-sample(which(colSums(transMatrix)==0),1)
items<-which(modules[,mod]==1)
par <- itemBank[items, ]
thStart <- NA
module<-mod
set.seed(NULL)
}
else{
mods<-which(colSums(transMatrix)==0)
info<-NULL
for (i in 1:length(mods)){
items<-which(modules[,mods[i]]==1)
info[i]<-sum(Ii(theta,itemBank[items,],model=model,D=D)$Ii)
}
keep<-min(which(info==max(info)))
items<-which(modules[,mods[keep]]==1)
par <- itemBank[items, ]
thStart <- theta
module<-mods[keep]
}
res <- list(module=module,items = items, par = par, thStart = thStart)
}
return(res)
}
