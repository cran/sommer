## -----------------------------------------------------------------------------
library(sommer)
data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

ansSingle <- mmes(Yield~1,
              random= ~ vsm(ism(Name), Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)

# if setting henderson=TRUE provide the inverse
# Ai <- solve(A)
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, "inverse")=TRUE


## -----------------------------------------------------------------------------

ansMain <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)

# if setting henderson=TRUE provide the inverse
# Ai <- solve(A)
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, "inverse")=TRUE


## -----------------------------------------------------------------------------

ansDG <- mmes(Yield~Env,
              random= ~ vsm(dsm(Env),ism(Name), Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)

# if setting henderson=TRUE provide the inverse
# Ai <- solve(A)
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, "inverse")=TRUE


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
              random= ~ vsm(usm(Env),ism(Name), Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)
# if setting henderson=TRUE provide the inverse
Ai <- solve(A)
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE



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

E <- AR1(DT$Env) # can be AR1() or CS(), etc.
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=A) + vsm(ism(Env:Name), Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)


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
mix2 <- mmes(y~ envIndex, henderson=TRUE,
             random=~ Name + vsm(dsm(envIndex),ism(Name)), data=DT2,
             rcov=~vsm(dsm(Name),ism(units)),
             tolParConvNorm = .0001,
             nIters = 50, verbose = FALSE
)
# summary(mix2)$varcomp

b=mix2$uList$`vsm(dsm(envIndex), ism(Name))` # adaptability (b) or genotype slopes
mu=mix2$uList$`sommer::vsm( sommer::ism( Name ) ` # general adaptation (mu) or main effect
e=sqrt(summary(mix2)$varcomp[-c(1:2),1]) # error variance for each individual

## general adaptation (main effect) vs adaptability (response to better environments)
plot(mu[,1]~b[,1], ylab="general adaptation", xlab="adaptability")
text(y=mu[,1],x=b[,1], labels = rownames(mu), cex=0.5, pos = 1)

## prediction across environments
Dt <- mix2$Dtable
Dt[1,"average"]=TRUE
Dt[2,"include"]=TRUE
Dt[3,"include"]=TRUE
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

# fit diagonal model first to produce H matrix
ansDG <- mmes(y~Env, henderson=TRUE,
              random=~ vsm(dsm(Env), ism(Name)),
              rcov=~units, nIters = 100,
              data=DT, verbose = FALSE)

H0 <- ansDG$uList$`vsm(dsm(Env), ism(Name))` # GxE table

# reduced rank model
ansFA <- mmes(y~Env, henderson=TRUE,
              random=~vsm( usm(rrm(Env, H = H0, nPC = 3)) , ism(Name)) + # rr
                vsm(dsm(Env), ism(Name)), # diag
              rcov=~units,
              # we recommend giving more iterations to these models
              nIters = 100, verbose = FALSE,
              # we recommend giving more EM iterations at the beggining
              data=DT)

vcFA <- ansFA$theta[[1]]
vcDG <- ansFA$theta[[2]]

loadings=with(DT, rrm(Env, nPC = 3, H = H0, returnGamma = TRUE) )$Gamma
scores <- ansFA$uList[[1]]

vcUS <- loadings %*% vcFA %*% t(loadings)
G <- vcUS + vcDG
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)

uFA <- scores %*% t(loadings)
uDG <- ansFA$uList[[2]]
u <- uFA + uDG


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
               data=droplevels(DT[which(DT$Env == envs[i]),]
               )
  )
  ans1$Beta$Env <- envs[i]
  
  BLUEL[[i]] <- data.frame( Effect=factor(rownames(ans1$b)), 
                            Estimate=ans1$b[,1], 
                            Env=factor(envs[i]))
  # to be comparable to 1/(se^2) = 1/PEV = 1/Ci = 1/[(X'X)inv]
  XtXL[[i]] <- solve(ans1$Ci[1:nrow(ans1$b),1:nrow(ans1$b)]) 
}

DT2 <- do.call(rbind, BLUEL)
OM <- do.call(adiag1,XtXL)

##########
## stage 2
## use mmes for sparse equation
##########
m <- matrix(1/var(DT2$Estimate, na.rm = TRUE))
ans2 <- mmes(Estimate~Env, henderson=TRUE,
             random=~ Effect + Env:Effect, 
             rcov=~vsm(ism(units,thetaC = matrix(3), theta = m)),
             W=OM, 
             verbose=FALSE,
             data=DT2
)
summary(ans2)$varcomp


