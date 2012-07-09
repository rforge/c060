### R code from vignette source 'c060_vignette.Rnw'

###################################################
### code chunk number 1: setup
###################################################
rm(list=ls())

# load libraries
library(cacheSweave)
library(Biobase)
library(genefilter)
library(ggplot2)
library(GEOquery)
library(limma)
library(xtable)
library(survival)

library(glmnet)
require(penalizedSVM)
library(parallel) #for stabilityselection.R
library(peperr) #for peperr_glmnet.R
library(tgp) #for EPSGO.R
library(mlegp) #for EPSGO.R
library(pamr) #for tune.glmnet.interval.R
library(lattice) # for plotting 

library(c060)

dir.create("figures",showWarnings = FALSE)

# set cache folder; needs different driver, to run as:
# Sweave("c060_vignette.Rnw",driver=cacheSweaveDriver)
setCacheDir("cache")


###################################################
### code chunk number 2: GEOdataGSE12417
###################################################
if (file.exists("Rdata/expressionSet.Rdata")) {
    load("Rdata/expressionSet.Rdata") 
} else {
    library(GEOquery)
    
    dir.create("GEOdata",showWarnings = FALSE)
    dir.create("Rdata",showWarnings = FALSE)
    
    # download processed data files from url or if available get data from locally saved copies
    geo_exprs_data  <- getGEO("GSE12417", destdir="GEOdata",AnnotGPL=T)[[1]]
    annotation(geo_exprs_data) <- "hgu133plus2"
    
    #clinical characteristics including survival data are in "characteristics_ch1" in the phenoData 
    clin.df <- as.data.frame(strsplit2(as.character(pData(geo_exprs_data)$characteristics_ch1), split=";"))
    clin.df[,1] <- strsplit2(clin.df[,1]," ")[,7]
    clin.df[,2] <- as.numeric(sub("=", "", strsplit2(clin.df[,2]," ")[,3]))
    clin.df[,3] <- as.numeric(strsplit2(clin.df[,3]," ")[,4])
    clin.df[,4] <- as.numeric(strsplit2(clin.df[,4]," ")[,4])
    colnames(clin.df) <- c("FAB", "age", "os", "os_status")
    rownames(clin.df) <- rownames(pData(geo_exprs_data))
    pData(geo_exprs_data) <- clin.df

    save(geo_exprs_data, file="Rdata/expressionSet.Rdata")  
    unlink("GEOdata", recursive=TRUE) #delete GEOdata
} 


###################################################
### code chunk number 3: preprocessData
###################################################
# unspecific filtering: select most varying probe sets
nfeatures <- 10000
rowvars    <- rowVars(exprs(geo_exprs_data))
topprobes   <- which(rowvars>=sort(rowvars,decreasing=T)[nfeatures])
eset       <- geo_exprs_data[topprobes,]     


###################################################
### code chunk number 4: glmnet
###################################################
set.seed(1234)
cvres <- cv.glmnet(y=Surv(pData(eset)$os, pData(eset)$os_status), x=t(exprs(eset)), family="cox", nfolds=10)
res <- cvres$glmnet.fit
plot(cvres)

cof <- coef(res, s=cvres$lambda.min)
names.cof <- rownames(cof)
cofn <- cof[which(cof!=0)]
names(cofn) <- names.cof[which(cof!=0)]


###################################################
### code chunk number 5: glmnet0
###################################################
fit <- glmnet(y=Surv(pData(eset)$os, pData(eset)$os_status), 
              x=t(exprs(eset)), family="cox")


###################################################
### code chunk number 6: glmnet2
###################################################
print(cofn)
cofn.lasso <- cofn
#cofn.ix <- predict(res, type="nonzero", s=cvres$lambda.min)


###################################################
### code chunk number 7: glmnet3
###################################################
bet <- res$beta[match(names(cofn), rownames(res$beta)),]
par(mar=c(4,4,2.5,1), mgp=c(2.5,1,0), mfrow=c(2,2))

plot(res, xvar="lambda", col="gray")
glmnet:::plotCoef(bet, lambda = res$lambda, df = res$df, dev = res$dev.ratio, xvar = "lambda", add=TRUE, col="red")
abline(v=log(cvres$lambda.min), lty=3)
abline(v=log(cvres$lambda.1se), lty=3)

glmnet:::plotCoef(bet, lambda = res$lambda, df = res$df, dev = res$dev.ratio, xvar = "lambda", add=FALSE, col="red")
abline(v=log(cvres$lambda.min), lty=3)
abline(v=log(cvres$lambda.1se), lty=3)

norm <- apply(abs(res$beta), 2, sum)
plot(res, xvar="norm", col="gray")
glmnet:::plotCoef(bet, xvar = "norm", add=TRUE, col="red",
                  norm = norm, lambda = res$lambda, df = res$df, dev = res$dev.ratio)
abline(v=norm[match(cvres$lambda.min, cvres$lambda)], lty=3)
abline(v=norm[match(cvres$lambda.1se, cvres$lambda)], lty=3)

plot(res, xvar="dev", col="gray")
glmnet:::plotCoef(bet, lambda = res$lambda, df = res$df, dev = res$dev.ratio, xvar = "dev", add=TRUE, col="red")
abline(v=res$dev.ratio[match(cvres$lambda.min, cvres$lambda)], lty=3)
abline(v=res$dev.ratio[match(cvres$lambda.1se, cvres$lambda)], lty=3)


###################################################
### code chunk number 8: peperrMandatoryParallel
###################################################
obj  <- peperr(response=Surv(pData(eset)$os, pData(eset)$os_status),
         x=data.frame(eset$age,t(exprs(eset))),
         fit.fun=fit.glmnet, 
         args.fit=list(standardize=F, family="cox", 
         penalty.factor=rep(0:1, times=c(1,dim(eset)[1]))),
         complexity=complexity.glmnet,
         args.complexity=list(standardize=F, nfolds=10, 
         family="cox", 
         penalty.factor=rep(0:1, times=c(1,dim(eset)[1]))),
         trace=F, RNG="fixed", seed=0815, 
         cpus=3, parallel=T, clustertype="SOCK", 
         load.list=list(packages=c("c060")),
         indices=resample.indices(n=dim(eset)[2], sample.n=100, 
         method="sub632"))


###################################################
### code chunk number 9: peperrPlot
###################################################
plot.peperr.curves(obj, at.risk=T)


###################################################
### code chunk number 10: stabilitySelection
###################################################
set.seed(1234)
y <- cbind(time=pData(eset)$os, status=pData(eset)$os_status)
cores <- 2
fwer <- 0.5
pi.thr <- 0.6
#calculate stability path
res <- stability.path(y=y, x=t(exprs(eset)), weakness=1, mc.cores=cores, family="cox")
par(mar=c(4,4,2,1), mgp=c(2.5,1,0))
#plot stability path and perform stability selection controlling a fwer of 0.5
sel <- plotstabpath(res, fwer=fwer, pi_thr=pi.thr, xvar="lambda", col.all="gray")


###################################################
### code chunk number 11: stabilitySelection2
###################################################
print(sel$stable)


###################################################
### code chunk number 12: interval_search_cox_setup
###################################################
  x <- t(exprs(eset))
  y <- cbind(time=pData(geo_exprs_data)$os,status=pData(geo_exprs_data)$os_status)
  
  Q.func <- "tune.glmnet.interval"
  
  # bounds for alpha [0,1]
  bounds <- t(data.frame(alpha=c(0, 1)))
  colnames(bounds) <- c("lower", "upper")  
  parms.coding <- "none"
  seed <- 1234
  type.measure <- 'deviance'
    
  family <- "cox"
  nfolds <- 10
  show <- "none"
  fminlower <- -100  
  type.min <- "lambda.1se"
  verbose <- TRUE


###################################################
### code chunk number 13: interval_search_cox_pre
###################################################
print(bounds)


###################################################
### code chunk number 14: interval_search_cox
###################################################
  # fix folds for each model
  set.seed(seed)
  foldid <- my.balanced.folds(class.column.factor=y[,2], cross.outer=nfolds)
  #table(foldid,y[,2])

  print("start interval search")
  fit <- EPSGO(Q.func, 
             bounds=bounds, 
             parms.coding=parms.coding, 
             show=show, N=NULL,   
             seed=seed, 
             fminlower=fminlower,
             # Q.func arguments
    				 x=x, y=y, family=family, 
             #nfolds=nfolds,  
             foldid=foldid,
             type.min=type.min,
             type.measure=type.measure,
             verbose=verbose)
  names(fit) 
 

  print("chose the model with min num of FS ")      
  print("# FS: ")
  sel.models <- sapply(fit$model.list, "[", "model") [fit$Ytrain == fit$fmin ]
  sel.alpha <- fit$xmin
  sel.error <- fit$fmin  	
  		
  out <- list(model=sel.models, alpha=sel.alpha, error=sel.error) 
  save(fit, out, file="Rdata/fit_interval_search.RData") 



###################################################
### code chunk number 15: interval_search_cox_extract
###################################################
# 
alphas <- fit$Xtrain[,1]
lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
deviances <- fit$Ytrain
# round problems!!!! take first from the fit 
# number of selected features in the models; dfs
tmp.models<-sapply(sapply(sapply(fit$model, "[", "model"), "[", "cvreg"), "[", "glmnet.fit")

n.features<-mapply( function(List, lam) List$df[which(List$lambda %in% lam)], tmp.models, lambdas)
             
opt.alpha <- out$alpha
opt.lambda <- out$model[[1]]$lambda

summary_int_search <- as.data.frame(cbind(alpha=alphas,lambda=lambdas,deviance=deviances,n.features=n.features))
rownames(summary_int_search) <- c(1:nrow(summary_int_search))


###################################################
### code chunk number 16: interval_search_cox_fit_opt
###################################################
#select the model with optimal parameters from the object fit
# 
f <- fit$model.list[unlist(sapply(fit$model.list,"[",1)) == fit$fmin]
if (length(f)>1) print("more than one optimal model!")
#take the first model
f1 <- f[[1]]$model$cvreg
res<-f1$glmnet.fit
cof <- coef(res, s=opt.lambda)
names.cof <- rownames(cof)
cofn <- cof[which(cof!=0)]
names(cofn) <- names.cof[which(cof!=0)] 
bet <- res$beta[match(names(cofn), rownames(res$beta)),]


###################################################
### code chunk number 17: interval_search_cox_summary
###################################################
summary(cofn)
head(sort(cofn))
tail(sort(cofn))


###################################################
### code chunk number 18: interval_search_cox_check_overlap_stability
###################################################
names(cofn.lasso) %in% names.cof
'206932_at' %in% names.cof


###################################################
### code chunk number 19: interval_search_cox_out_plot1
###################################################
 
# breaks <- do.breaks(range(summary_int_search$deviance), 20)
# n_init <- 20 # number of initial alpha values at iteration zero
# summary_int_search$cols <- level.colors(summary_int_search$deviance,at=breaks, col.regions = cm.colors(n=20,alpha=0.5))
# print(xyplot(alpha ~ lambda, data = summary_int_search, groups = cols, cex = 2, col = "black",   jitter.y=T, amount=0.01,
#              scales=list(x=list(log=T, equispaced.log = FALSE)), # x axis on log-scale
#              panel = function(x, y, groups, ..., subscripts) { 
#              fill <- groups[subscripts] 
#              panel.grid(h = -1, v = -1)
#              panel.xyplot(x, y, pch = rep(c(22,21),c(n_init,nrow(summary_int_search)-n_init)), fill = fill, ...) ;
#              ltext(x=x, y=y, labels=n.features, pos=ifelse(y<0.1,3,4), offset=3, cex=0.75,col=4)},
#              legend = list(top = list(fun = draw.colorkey, args = list(key = list(space = "top", col = cm.colors(n=20,alpha=0.5), at = breaks), draw = FALSE))),
#              main="cross-validated model deviance",sub="number of selected features are given in blue \n rectangles show initial alpha values"))

summary_int_search$log_lambda <- log(summary_int_search$lambda)
breaks                        <- do.breaks(range(summary_int_search$deviance), 20)
n_init                        <- 21 # number of initial alpha values at iteration zero
summary_int_search$cols       <- level.colors(summary_int_search$deviance,at=breaks, col.regions = gray.colors)

print(xyplot(alpha ~ log_lambda, data = summary_int_search, groups = cols, cex = 2, col = "black",  
             jitter.y=T, amount=0.01, xlab=expression(paste("log ",lambda)),ylab=expression(alpha),
#             scales=list(x=list(log=T, equispaced.log = FALSE)), # x axis on log-scale
             panel = function(x, y, groups, ..., subscripts) { 
             fill <- groups[subscripts] 
             panel.grid(h = -1, v = -1)
             panel.xyplot(x, y, pch = rep(c(22,21),c(n_init,nrow(summary_int_search)-n_init)), fill = fill, ...) ;
             ltext(x=x, y=y, labels=n.features, pos=ifelse(y<0.1,3,4), offset=3, cex=1,col=1)},
             legend = list(top = list(fun = draw.colorkey, args = list(key = list(space = "top", col = gray.colors, at = breaks), draw = FALSE))),
             main="Cross-validated partial log likelihood deviance",
             #sub="number of selected features are printed next to symbol \n rectangles show initial alpha values"
             ))


###################################################
### code chunk number 20: interval_search_cox_out_plot2
###################################################
niter<-nrow(summary_int_search) - 21
iter<-c(rep(0,21), c(1:niter))
plot(summary_int_search$alpha, iter, xlab=expression(alpha), ylab="Iteration")
grid(NA, niter+1, lwd=2)
abline(v=opt.alpha, col="red")


