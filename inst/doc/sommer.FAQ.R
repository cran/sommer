## ------------------------------------------------------------------------
# iteration    LogLik     wall    cpu(sec)   restrained
#     1      -224.676   18:11:23      3           0
# Sistem is singular. Stopping the job
# matrix multiplication: incompatible matrix dimensions: 0x0 and ...x...

## ------------------------------------------------------------------------
library(sommer)
## rrBLUP for makers
data(DT_cpdata)
mix.rrblup <- mmer(fixed=color~1,
                   random=~vs(GT,Gtc=unsm(1)) + vs(Rowf,Gtc=diag(1)),
                   rcov=~vs(units,Gtc=unsm(1)), getPEV = FALSE,
                   data=DT, verbose = FALSE)
summary(mix.rrblup)
## GBLUP for individuals
A <- A.mat(GT)
mix.gblup <- mmer(fixed=color~1,
                  random=~vs(id,Gu=A, Gtc=unsm(1)) + vs(Rowf,Gtc=diag(1)),
                  rcov=~vs(units,Gtc=unsm(1)),
                  data=DT, verbose = FALSE)
summary(mix.gblup)
## Equivalence
plot(GT%*%mix.rrblup$U$`u:GT`$color, mix.gblup$U$`u:id`$color)

## ------------------------------------------------------------------------
library(sommer)

data(DT_cpdata)
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
#### look at the data and fit the model
set.seed(12)
DT2 <- droplevels(DT[sample(1:nrow(DT),100),]) # we simulate a dataset with only 100 animals
nrow(DT2); length(levels(DT2$id))

# we fit a model with the reduced datatset where only 100 blups will be returned since only
# 100 levels exist in the "id" column
mix1 <- mmer(Yield~1,
              random=~vs(id,Gu=A)
                      + Rowf + Colf,
              rcov=~units,
              data=DT2, verbose = FALSE)
summary(mix1)
length(mix1$U$`u:id`$Yield) # only 100 levels

# we add additional levels to the "id" column and also provide them in the relationship matrix
levels(DT2$id) <- c(levels(DT2$id), setdiff(levels(DT$id), levels(DT2$id)))
mix2 <- mmer(Yield~1,
             random=~vs(id,Gu=A)
             + Rowf + Colf,
             rcov=~units,
             data=DT2, verbose = FALSE)
summary(mix2)
length(mix2$U$`u:id`$Yield) # now 363 levels

## ------------------------------------------------------------------------
library(sommer)
data(DT_cpdata)
mix1 <- mmer(Yield~1,
              random=~ Rowf + Colf,
              rcov=~units,
              data=DT, verbose = FALSE)
summary(mix1)$varcomp

## ------------------------------------------------------------------------
library(sommer)
data(DT_cpdata)
mixAR1row <- mmer(Yield~1,
             random=~ vs(Rowf, Gu=AR1(Rowf, rho=0.3)) + Colf,
             rcov=~units,
             data=DT, verbose = FALSE)
summary(mixAR1row)$varcomp

## ------------------------------------------------------------------------
library(sommer)
data(DT_cpdata)
mixAR1col <- mmer(Yield~1,
             random=~ Rowf + vs(Colf, Gu=AR1(Colf, rho=0.3)),
             rcov=~units,
             data=DT, verbose = FALSE)
summary(mixAR1col)$varcomp

## ------------------------------------------------------------------------
library(sommer)
data(DT_cpdata)
mixAR1rowcol <- mmer(Yield~1,
                  random=~ vs(Rowf:Colf, 
                              Gu=kronecker(AR1(Rowf, rho=0.3),AR1(Colf, rho=0.3),make.dimnames = TRUE) 
                              ),
                  rcov=~units,
                  data=DT, verbose = FALSE)
summary(mixAR1rowcol)$varcomp

## ------------------------------------------------------------------------
library(sommer)
data(DT_example)
M <- matrix(rep(0,41*1000),1000,41)
for (i in 1:41) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
tM <- t(M)

## ------------------------------------------------------------------------
## GWAS for main term in CS model
ansx <- GWAS(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT,
             M=tM,
             gTerm = "Name", verbose = FALSE)

ms <- as.data.frame(t(ansx$scores))
plot(ms$`Yield score`, ylim=c(0,8))

## ------------------------------------------------------------------------
## GWAS for the interaction term in CS model
E <- matrix(1,nrow = length(unique(DT$Env)));E
EtM <- kronecker(E,tM)
ansx <- GWAS(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT,
             M=EtM,
             gTerm = "Env:Name", verbose = FALSE)

ms <- as.data.frame(t(ansx$scores))
plot(ms$`Yield score`, ylim=c(0,8))

## ------------------------------------------------------------------------
## GWAS for the interaction term in DIAG model
E <- matrix(1,nrow = length(unique(DT$Env)));E
EtM <- kronecker(E,tM)
ansx <- GWAS(Yield~Env,
             random= ~Name + vs(ds(Env),Name),
             rcov= ~ vs(ds(Env),units),
             data=DT,
             M=tM,
             gTerm = "CA.2011:Name", verbose = FALSE )

ms <- as.data.frame(t(ansx$scores))
plot(ms$`Yield score`, ylim=c(0,8))

## ------------------------------------------------------------------------
## GWAS for main term in US model
ansx <- GWAS(Yield~Env,
             random= ~vs(us(Env),Name),
             rcov= ~ vs(us(Env),units),
             data=DT,
             M=tM,
             gTerm = "CA.2011:Name", verbose = FALSE)

ms <- as.data.frame(t(ansx$scores))
plot(ms$`Yield score`, ylim=c(0,8))

## ------------------------------------------------------------------------
## GWAS for main term in multitrait DIAG model
ansx <- GWAS(cbind(Weight,Yield)~Env,
             random= ~vs(ds(Env),Name, Gtc=unsm(2)),
             rcov= ~ vs(ds(Env),units, Gtc=diag(2)),
             data=DT,
             M=tM,
             gTerm = "CA.2011:Name", verbose = FALSE)

ms <- as.data.frame(t(ansx$scores))
plot(ms$`Yield score`, ylim=c(0,8))

