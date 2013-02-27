pkgname <- "c060"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('c060')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("plot.peperr.curves")
### * plot.peperr.curves

flush(stderr()); flush(stdout())

### Name: plot.peperr.curves
### Title: Plot method for prediction error curves of a peperr object
### Aliases: plot.peperr.curves
### Keywords: models regression survival

### ** Examples

## Not run: 
##D 
##D # example from glmnet package
##D set.seed(10101)
##D library(glmnet)
##D library(survival)
##D library(peperr)
##D 
##D N=1000;p=30
##D nzc=p/3
##D x=matrix(rnorm(N*p),N,p)
##D beta=rnorm(nzc)
##D fx=x[,seq(nzc)]##D 
##D hx=exp(fx)
##D ty=rexp(N,hx)
##D tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
##D y=Surv(ty,1-tcens)
##D 
##D peperr.object <- peperr(response=y, x=x, 
##D                         fit.fun=fit.glmnet, args.fit=list(family="cox"), 
##D                         complexity=complexity.glmnet,  
##D                         args.complexity=list(family="cox",nfolds=10),
##D                         indices=resample.indices(n=N, method="sub632", sample.n=10))
##D 
##D plot.peperr.curves(peperr.object)
##D # peperr
##D plot(peperr.object)
## End(Not run)



cleanEx()
nameEx("plotstabpath")
### * plotstabpath

flush(stderr()); flush(stdout())

### Name: plotstabpath
### Title: function to plot a stability path
### Aliases: plotstabpath
### Keywords: stability selection

### ** Examples

## Not run: 
##D #gaussian
##D set.seed(1234)
##D x=matrix(rnorm(100*1000,0,1),100,1000)
##D y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
##D res <- stability.path(y,x,weakness=1,mc.cores=2)
##D plotstabpath(res,fwer=.5)
## End(Not run)



cleanEx()
nameEx("stability.path")
### * stability.path

flush(stderr()); flush(stdout())

### Name: stability.path
### Title: Stability path for glmnet models
### Aliases: stability.path
### Keywords: stability selection

### ** Examples
## Not run: 
##D #gaussian
##D set.seed(1234)
##D x <- matrix(rnorm(100*1000,0,1),100,1000)
##D y <- x[1:100,1:1000]%*% c(rep(2,5),rep(-2,5),rep(.1,990))
##D res <- stability.path(y,x,weakness=1,mc.cores=2)
##D plotstabpath(res)
##D 
##D #binomial
##D y=sample(1:2,100,replace=TRUE)
##D res <- stability.path(y,x,weakness=1,mc.cores=2,family="binomial")
##D plotstabpath(res)
##D     
##D #multinomial
##D y=sample(1:4,100,replace=TRUE)
##D res <- stability.path(y,x,weakness=1,mc.cores=2,family="multinomial")
##D plotstabpath(res)
##D     
##D #poisson
##D N=100; p=1000
##D nzc=5
##D x=matrix(rnorm(N*p),N,p)
##D beta=rnorm(nzc)
##D f = x[,seq(nzc)]%*%beta
##D mu=exp(f)
##D y=rpois(N,mu)
##D res <- stability.path(y,x,weakness=1,mc.cores=2,family="poisson")
##D plotstabpath(res)
##D 
##D #Cox
##D library(survival)
##D set.seed(10101)
##D N=100;p=1000
##D nzc=p/3
##D x=matrix(rnorm(N*p),N,p)
##D beta=rnorm(nzc)
##D fx=x[,seq(nzc)]%*%beta/3
##D hx=exp(fx)
##D ty=rexp(N,hx)
##D tcens=rbinom(n=N,prob=.3,size=1)
##D y=cbind(time=ty,status=1-tcens)
##D res <- stability.path(y,x,weakness=1,mc.cores=2,family="cox")
##D plotstabpath(res)
## End(Not run)



cleanEx()
nameEx("stability.selection")
### * stability.selection

flush(stderr()); flush(stdout())

### Name: stability.selection
### Title: function to estimate a stable set of variables
### Aliases: stability.selection
### Keywords: stability selection

### ** Examples

## Not run: 
##D #gaussian
##D set.seed(1234)
##D x=matrix(rnorm(100*1000,0,1),100,1000)
##D y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
##D res <- stability.path(y,x,weakness=1,mc.cores=2)
##D stability.selection(res,fwer=.5)
## End(Not run)


### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
