## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
## solving for r records
ans1r <- mmer(Yield~Env,
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ans1r)$varcomp
## solving for c coefficients
ans1c <- mmec(Yield~Env,
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ans1c)$varcomp


## -----------------------------------------------------------------------------

data(DT_example)
DT <- DT_example
 
ans2r <- mmer(Yield~Env,
              random= ~Name + vsr(dsr(Env),Name),
              rcov= ~ vsr(dsr(Env),units),
              data=DT, verbose = FALSE)
summary(ans2r)$varcomp

DT=DT[with(DT, order(Env)), ]
ans2c <- mmec(Yield~Env,
              random= ~Name + vsc(dsc(Env),isc(Name)),
              rcov= ~ vsc(dsc(Env),isc(units)),
              data=DT, verbose = FALSE)
summary(ans2c)$varcomp


## -----------------------------------------------------------------------------

data(DT_example)
DT <- DT_example

ans3r <- mmer(Yield~Env,
             random=~ vsr(usr(Env),Name),
             rcov=~vsr(dsr(Env),units), 
             data=DT, verbose = FALSE)
summary(ans3r)$varcomp

DT=DT[with(DT, order(Env)), ]
ans3c <- mmec(Yield~Env,
             random=~ vsc(usc(Env),isc(Name)),
             rcov=~vsc(dsc(Env),isc(units)), 
             data=DT, verbose = FALSE)
summary(ans3c)$varcomp


## -----------------------------------------------------------------------------

data(DT_example)
DT <- DT_example
DT$EnvName <- paste(DT$Env,DT$Name)

DT$Yield <- as.vector(scale(DT$Yield))
DT$Weight <- as.vector(scale(DT$Weight))

ans4r <- mmer(cbind(Yield, Weight) ~ Env,
             random= ~ vsr(Name, Gtc=unsm(2)),
             rcov= ~ vsr(units, Gtc=diag(2)),
             data=DT, verbose = FALSE)
summary(ans4r)$varcomp

DT2 <- reshape(DT, idvar = c("Name","Env","Block"), varying = list(6:7),
        v.names = "y", direction = "long", timevar = "trait", times = colnames(DT)[6:7])
DT2$trait <- as.factor(DT2$trait)

# ans4c <- mmec(y ~ Env:trait,
#              random= ~ vsc(usc(trait),isc(Name)),
#              rcov= ~ vsc(dsc(trait),isc(units)), returnParam = T, 
#              data=DT2, verbose = T)
# summary(ans4c)$varcomp


## -----------------------------------------------------------------------------
unsm(4)

## -----------------------------------------------------------------------------
fixm(4)

## -----------------------------------------------------------------------------
fcm(c(1,0,1,0))

## -----------------------------------------------------------------------------
library(orthopolynom)
data(DT_legendre)
DT <- DT_legendre

mRR2r<-mmer(Y~ 1 + Xf
           , random=~ vsr(usr(leg(X,1)),SUBJECT)
           , rcov=~vsr(units)
           , data=DT, verbose = FALSE)
summary(mRR2r)$varcomp

mRR2c<-mmec(Y~ 1 + Xf
           , random=~ vsc(usc(leg(X,1)),isc(SUBJECT))
           , rcov=~vsc(isc(units))
           , data=DT, verbose = FALSE)
summary(mRR2c)$varcomp


## -----------------------------------------------------------------------------

data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
#### look at the data and fit the model
head(DT,3)
head(MP,3)
GT[1:3,1:4]
mix1 <- GWAS(color~1,
             random=~vsr(id,Gu=A)
             + Rowf + Colf,
             rcov=~units,
             data=DT, nIters=3,
             M=GT, gTerm = "u:id",
             verbose = FALSE)

ms <- as.data.frame(mix1$scores)
ms$Locus <- rownames(ms)
MP2 <- merge(MP,ms,by="Locus",all.x = TRUE);
manhattan(MP2, pch=20,cex=.5, PVCN = "color")


## -----------------------------------------------------------------------------

data("DT_halfdiallel")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modhr <- mmer(sugar~1, 
             random=~vsr(overlay(femalef,malef)) 
             + genof, data=DT,verbose = FALSE)

modhc <- mmec(sugar~1, 
             random=~vsc(isc(overlay(femalef,malef, sparse = TRUE))) 
             + genof,data=DT,verbose = FALSE)


## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### mimic two fields
A <- A.mat(GT)
mix <- mmer(Yield~1,
            random=~vsr(id, Gu=A) +
              vsr(Rowf) +
              vsr(Colf) +
              spl2Da(Row,Col),
            rcov=~vsr(units), nIters=3,
            data=DT,verbose = FALSE)
summary(mix)

## -----------------------------------------------------------------------------

data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata

#### look at the data and fit the model
mix1 <- mmer(Yield~1,
              random=~vsr(list(GT)),
              rcov=~units, nIters=3,
              data=DT,verbose = FALSE)


## -----------------------------------------------------------------------------
data(DT_wheat)
DT <- DT_wheat
GT <- GT_wheat
colnames(DT) <- paste0("X",1:ncol(DT))
DT <- as.data.frame(DT);DT$id <- as.factor(rownames(DT))
# select environment 1
rownames(GT) <- rownames(DT)
K <- A.mat(GT) # additive relationship matrix
colnames(K) <- rownames(K) <- rownames(DT)
# GBLUP pedigree-based approach
set.seed(12345)
y.trn <- DT
vv <- sample(rownames(DT),round(nrow(DT)/5))
y.trn[vv,"X1"] <- NA

## GBLUP
ans <- mmer(X1~1,
            random=~vsr(id,Gu=K), 
            rcov=~units, nIters=3,
            data=y.trn,verbose = FALSE) # kinship based
ans$U$`u:id`$X1 <- as.data.frame(ans$U$`u:id`$X1)
rownames(ans$U$`u:id`$X1) <- gsub("id","",rownames(ans$U$`u:id`$X1))
cor(ans$U$`u:id`$X1[vv,],DT[vv,"X1"], use="complete")

## rrBLUP
ans2 <- mmer(X1~1,
             random=~vsr(list(GT), buildGu = FALSE), 
             rcov=~units, getPEV = FALSE, nIters=3,
             data=y.trn,verbose = FALSE) # kinship based

u <- GT %*% as.matrix(ans2$U$`u:GT`$X1) # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP
# same can be achieved with the mmec function (see ?DT_wheat)

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### mimic two fields
A <- A.mat(GT)

mix1 <- mmer(Yield~1,
            random=~vsr(id, Gu=A) +
              vsr(Rowf) +
              vsr(Colf),
            rcov=~vsr(units), nIters=3,
            data=DT, verbose = FALSE)

## -----------------------------------------------------------------------------
mix2 <- mmer(Yield~1,
            random=~vsr(id, Gu=A) +
              vsr(Rowf) +
              vsr(Colf) +
              spl2Da(Row,Col),
            rcov=~vsr(units), nIters=3,
            data=DT,verbose = FALSE)

## -----------------------------------------------------------------------------
lrt <- anova(mix1, mix2)

## -----------------------------------------------------------------------------

data(DT_example)
DT <- DT_example

DT$EnvName <- paste(DT$Env,DT$Name)
modelBase <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(Name, Gtc=diag(2)), # here is diag()
              rcov= ~ vsr(units, Gtc=unsm(2)), nIters=3,
              data=DT,verbose = FALSE)

modelCov <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(usr(Env),Name, Gtc=unsm(2)), # here is unsm()
              rcov= ~ vsr(dsr(Env),units, Gtc=unsm(2)), nIters=3,
              data=DT,verbose = FALSE)

lrt <- anova(modelBase, modelCov)

## -----------------------------------------------------------------------------
library(sommer)
data(DT_yatesoats)
DT <- DT_yatesoats
m3 <- mmer(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT, verbose=FALSE)
summary(m3)$varcomp

## -----------------------------------------------------------------------------
Dt <- m3$Dtable; Dt
# first fixed effect just average
Dt[1,"average"] = TRUE
# second fixed effect include
Dt[2,"include"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N")
pp$pvals

