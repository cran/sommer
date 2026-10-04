## -----------------------------------------------------------------------------
library(sommer)
data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

Ai <- solve(A)
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE

ansSingle <- mmes(Yield~1,
              random= ~ vsm(ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)



## -----------------------------------------------------------------------------

ansMain <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)



## -----------------------------------------------------------------------------

ansDG <- mmes(Yield~Env,
              random= ~ vsm(dsm(Env),ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)


## -----------------------------------------------------------------------------
E <- diag(length(unique(DT$Env)));rownames(E) <- colnames(E) <- unique(DT$Env)
Ei <- solve(E)
Ai <- solve(A)
EAi <- kronecker(Ei,Ai, make.dimnames = TRUE)
Ei <- as(as(as( Ei,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
EAi <- as(as(as( EAi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE
attr(EAi, "inverse")=TRUE
ansCS <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=Ai) + vsm(ism(Env:Name), Gu=EAi),
              rcov= ~ units, 
              data=DT, verbose = FALSE)
summary(ansCS)


## -----------------------------------------------------------------------------

ansUS <- mmes(Yield~Env,
              random= ~ vsm(usm(Env),ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)



## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))

ansRR <- mmes(Yield~Env,
              random= ~ vsm(dsm(leg(EnvN,1)),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))

ansRR <- mmes(Yield~Env,
              random= ~ vsm(usm(leg(EnvN,1)),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------

ansAR1 <- mmes(Yield~Env,
              random= ~ vsm(csm(Env),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansAR1)


## -----------------------------------------------------------------------------

data(DT_h2, package="enhancer")
DT <- DT_h2

## build the environmental index
ei <- aggregate(y~Env, data=DT,FUN=mean)
colnames(ei)[2] <- "envIndex"
ei$envIndex <- ei$envIndex - mean(ei$envIndex,na.rm=TRUE) # center the envIndex to have clean VCs
ei <- ei[with(ei, order(envIndex)), ]

## add the environmental index to the original dataset
DT2 <- merge(DT,ei, by="Env")

# numeric by factor variables like envIndex:Name can't be used in the random part like this
# they need to come with the vsm() structure
DT2 <- DT2[with(DT2, order(Name)), ]
mix2 <- mmes(y~ envIndex, 
             random=~ Name + vsm(ism(envIndex),ism(Name)), data=DT2,
             rcov=~vsm(dsm(Name),ism(units)),
             nIters = 50, verbose = FALSE
)
# summary(mix2)$varcomp

b=mix2$uList$`vsm(ism(envIndex), ism(Name` # adaptability (b) or genotype slopes
mu=mix2$uList$`vsm(ism(Name`# general adaptation (mu) or main effect
e=sqrt(summary(mix2)$varcomp[-c(1:2),"estimate"]) # error variance for each individual

## general adaptation (main effect) vs adaptability (response to better environments)
plot(mu[,1]~b[,1], ylab="general adaptation", xlab="adaptability")
text(y=mu[,1],x=b[,1], labels = rownames(mu), cex=0.5, pos = 1)

## prediction across environments
Dt <- mix2$Dtable
Dt[1,"average"]=TRUE
Dt[2,"include"]=TRUE
Dt[3,"include"]=TRUE

mix2 <- postPEV(mix2, mode=2)
pp <- predict(mix2,Dtable = Dt, D="Name")
preds <- pp$pvals
# preds[with(preds, order(-predicted.value)), ]
## performance vs stability (deviation from regression line)
plot(preds[,2]~e, ylab="performance", xlab="stability")
text(y=preds[,2],x=e, labels = rownames(mu), cex=0.5, pos = 1)


## -----------------------------------------------------------------------------

data(DT_h2, package="enhancer")
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

# factor analytic model with 2 factors
ansFA2 <- mmes(y~Env, 
               random=~vsm( fam(Env, 2) , ism(Name)) ,
               rcov=~units,
               nIters = 100, verbose = FALSE,
               data=DT)

Dt <- ansFA2$Dtable; Dt
Dt[1:2,"average"]=TRUE
Dt[3,c("average","include")]=TRUE

ppfa <- predict(ansFA2, D="Env:Name", Dtable = Dt)
head(ppfa$pvals)


## -----------------------------------------------------------------------------

# reduced rank model with 2 factors
ansRR2 <- mmes(y~Env, henderson=TRUE,
              random=~vsm( rrm(Env, 2) , ism(Name)) + # rr
                vsm(dsm(Env), ism(Name)), # diag
              rcov=~units,
              nIters = 100, verbose = FALSE,
              data=DT)


Dt <- ansRR2$Dtable; Dt
Dt[1:2,"average"]=TRUE
Dt[3:4,c("average","include")]=TRUE

pprr <- predict(ansRR2, D="Env:Name", Dtable = Dt)
head(pprr$pvals)

# compare 
plot(ppfa$pvals[,"predicted.value"],pprr$pvals[,"predicted.value"])


## -----------------------------------------------------------------------------

# fit diagonal model first to produce H matrix
ansDG <- mmes(y~Env, henderson=TRUE,
              random=~ vsm(dsm(Env), ism(Name)),
              rcov=~units, nIters = 100,
              data=DT, verbose = FALSE)

H0 <- ansDG$uList$`vsm(dsm(Env), ism(Name))` # GxE table

# # reduced rank model
# ansFA <- mmes(y~Env, henderson=TRUE,
#               random=~vsm( usm(rrmat(Env, H = H0, nPC = 2)) , ism(Name)) + # rr
#                 vsm(dsm(Env), ism(Name)), # diag
#               rcov=~units,
#               # we recommend giving more iterations to these models
#               nIters = 100, verbose = FALSE,
#               # we recommend giving more EM iterations at the beggining
#               data=DT)
# 
# vcFA <- ansFA$theta[[1]]
# vcDG <- ansFA$theta[[2]]
# 
# loadings=with(DT, rrmat(Env, nPC = 2, H = H0, returnGamma = TRUE) )$Gamma
# scores <- ansFA$uList[[1]]
# 
# vcUS <- loadings %*% vcFA %*% t(loadings)
# G <- vcUS + vcDG
# # colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# # hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)
# 
# uFA <- scores %*% t(loadings)
# uDG <- ansFA$uList[[2]]
# u <- uFA + uDG
# 
# plot(ppfa$pvals[,"predicted.value"],as.vector(t(u)))
# plot(as.vector(t(u)),pprr$pvals[,"predicted.value"])


## -----------------------------------------------------------------------------

##########
## stage 1
## use mmes for dense field trials
##########
data(DT_h2, package="enhancer")
DT <- DT_h2
head(DT)
envs <- unique(DT$Env)
BLUEL <- list()
XtXL <- list()
for(i in 1:length(envs)){
  ans1 <- mmes(y~Name-1,
               random=~Block,
               verbose=FALSE,
               computeCi = 2,
               data=droplevels(DT[which(DT$Env == envs[i]),])
  )
  ans1$Beta$Env <- envs[i]
  
  BLUEL[[i]] <- data.frame( Effect=factor(rownames(ans1$b)), 
                            Estimate=ans1$b[,1], 
                            Env=factor(envs[i]))
  # to be comparable to 1/(se^2) = 1/PEV = 1/Ci = 1/[(X'X)inv]
  XtXL[[i]] <- solve(ans1$Ci[1:nrow(ans1$b),1:nrow(ans1$b)]) 
}

DT2 <- do.call(rbind, BLUEL)
OM <- Reduce(adiag1,lapply(XtXL,as.matrix))

##########
## stage 2
## use mmes for sparse equation
##########
m <- matrix(1/var(DT2$Estimate, na.rm = TRUE))

ans2 <- mmes(Estimate~Env, henderson=TRUE,
             random=~ Effect + Env:Effect, 
             rcov = ~ vsm(
               ism(units),
               sigma2 = 1,
               fixedSigma2 = TRUE
             ),
             W=OM, 
             verbose=FALSE,
             data=DT2
)
summary(ans2)$varcomp


