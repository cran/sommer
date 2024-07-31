## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example

ansSingle <- mmer(Yield~1,
              random= ~ vsr(Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)

# or
Ai <- as(solve(A), Class="dgCMatrix")
ansSingle <- mmec(Yield~1,
              random= ~ vsc(isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)


## -----------------------------------------------------------------------------

ansMain <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)

# or 

Ai <- as(solve(A), Class="dgCMatrix")
ansMain <- mmec(Yield~Env,
              random= ~ vsc(isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)


## -----------------------------------------------------------------------------

ansDG <- mmer(Yield~Env,
              random= ~ vsr(dsr(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)

# or
Ai <- as(solve(A), Class="dgCMatrix")
ansDG <- mmec(Yield~Env,
              random= ~ vsc(dsc(Env),isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)


## -----------------------------------------------------------------------------
E <- diag(length(unique(DT$Env)))
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A) + vsr(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)

## or
E <- diag(length(unique(DT$Env)));rownames(E) <- colnames(E) <- unique(DT$Env)
Ei <- solve(E)
Ai <- solve(A)
EAi <- kronecker(Ei,Ai, make.dimnames = TRUE)
Ei <- as(Ei, Class="dgCMatrix")
Ai <- as(Ai, Class="dgCMatrix")
EAi <- as(EAi, Class="dgCMatrix")
ansCS <- mmec(Yield~Env,
              random= ~ vsc(isc(Name), Gu=Ai) + vsc(isc(Env:Name), Gu=EAi),
              rcov= ~ units, 
              data=DT, verbose = FALSE)
summary(ansCS)


## -----------------------------------------------------------------------------

ansUS <- mmer(Yield~Env,
              random= ~ vsr(usr(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)
# adjust variance BLUPs by adding covariances
# ansUS$U[1:6] <- unsBLUP(ansUS$U[1:6])

# or
Ai <- solve(A)
Ai <- as(Ai, Class="dgCMatrix")
ansUS <- mmec(Yield~Env,
              random= ~ vsc(usc(Env),isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)


## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))
ansRR <- mmer(Yield~Env,
              random= ~ vsr(leg(EnvN,1),Name),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)

# or

ansRR <- mmec(Yield~Env,
              random= ~ vsc(dsc(leg(EnvN,1)),isc(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))
ansRR <- mmer(Yield~Env,
              random= ~ vsr(usr(leg(EnvN,1)),Name),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)

# or

ansRR <- mmec(Yield~Env,
              random= ~ vsc(usc(leg(EnvN,1)),isc(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------

E <- AR1(DT$Env) # can be AR1() or CS(), etc.
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A) + vsr(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)


## -----------------------------------------------------------------------------

data(DT_h2)
DT <- DT_h2

## build the environmental index
ei <- aggregate(y~Env, data=DT,FUN=mean)
colnames(ei)[2] <- "envIndex"
ei$envIndex <- ei$envIndex - mean(ei$envIndex,na.rm=TRUE) # center the envIndex to have clean VCs
ei <- ei[with(ei, order(envIndex)), ]

## add the environmental index to the original dataset
DT2 <- merge(DT,ei, by="Env")

# numeric by factor variables like envIndex:Name can't be used in the random part like this
# they need to come with the vsc() structure
DT2 <- DT2[with(DT2, order(Name)), ]
mix2 <- mmec(y~ envIndex,
             random=~ Name + vsc(dsc(envIndex),isc(Name)), data=DT2,
             rcov=~vsc(dsc(Name),isc(units)),
             tolParConvNorm = .0001,
             nIters = 50, verbose = FALSE
)
# summary(mix2)$varcomp

b=mix2$uList$`vsc(dsc(envIndex), isc(Name))` # adaptability (b) or genotype slopes
mu=mix2$uList$`vsc( isc( Name ) )` # general adaptation (mu) or main effect
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

data(DT_h2)
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

# fit diagonal model first to produce H matrix
ansDG <- mmec(y~Env,
              random=~ vsc(dsc(Env), isc(Name)),
              rcov=~units, nIters = 100,
              # we recommend giving more EM iterations at the beggining
              emWeight = c(rep(1,10),logspace(10,1,.05), rep(.05,80)),
              data=DT, verbose = FALSE)

H0 <- ansDG$uList$`vsc(dsc(Env), isc(Name))` # GxE table

# reduced rank model
ansFA <- mmec(y~Env,
              random=~vsc( usc(rrc(Env, H = H0, nPC = 3)) , isc(Name)) + # rr
                vsc(dsc(Env), isc(Name)), # diag
              rcov=~units,
              # we recommend giving more iterations to these models
              nIters = 100, verbose = FALSE,
              # we recommend giving more EM iterations at the beggining
              emWeight = c(rep(1,10),logspace(10,1,.05), rep(.05,80)),
              data=DT)

vcFA <- ansFA$theta[[1]]
vcDG <- ansFA$theta[[2]]

loadings=with(DT, rrc(Env, nPC = 3, H = H0, returnGamma = TRUE) )$Gamma
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
## use mmer for dense field trials
##########
data(DT_h2)
DT <- DT_h2
head(DT)
envs <- unique(DT$Env)
BLUEL <- list()
XtXL <- list()
for(i in 1:length(envs)){
  ans1 <- mmer(y~Name-1,
                random=~Block,
                verbose=FALSE,
                data=droplevels(DT[which(DT$Env == envs[i]),]
               )
  )
  ans1$Beta$Env <- envs[i]
  
  BLUEL[[i]] <- ans1$Beta
  # to be comparable to 1/(se^2) = 1/PEV = 1/Ci = 1/[(X'X)inv]
  XtXL[[i]] <- solve(ans1$VarBeta) 
}

DT2 <- do.call(rbind, BLUEL)
OM <- do.call(adiag1,XtXL)

##########
## stage 2
## use mmec for sparse equation
##########
m <- matrix(1/var(DT2$Estimate, na.rm = TRUE))
ans2 <- mmec(Estimate~Env,
             random=~Effect + Env:Effect, 
             rcov=~vsc(isc(units,thetaC = matrix(3), theta = m)),
             W=OM, 
             verbose=FALSE,
             data=DT2
             )
summary(ans2)$varcomp


