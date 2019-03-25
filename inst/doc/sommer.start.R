## ------------------------------------------------------------------------
library(sommer)
data(DT_example)
head(DT)

ans1 <- mmer(Yield~Env,
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=DT)
summary(ans1)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
ans2 <- mmer(Yield~Env,
              random= ~Name + vs(ds(Env),Name),
              rcov= ~ vs(ds(Env),units),
              data=DT)
summary(ans2)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
ans3 <- mmer(Yield~Env,
             random=~ vs(us(Env),Name),
             rcov=~vs(us(Env),units), 
             data=DT)
summary(ans3)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
DT$EnvName <- paste(DT$Env,DT$Name)
ans4 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vs(Name, Gtc=unsm(2)) + vs(EnvName, Gtc=unsm(2)),
              rcov= ~ vs(units, Gtc=unsm(2)),
              data=DT)
summary(ans4)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
DT$EnvName <- paste(DT$Env,DT$Name)
ans5 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vs(Name, Gtc=unsm(2)) + vs(ds(Env),Name, Gtc=unsm(2)),
              rcov= ~ vs(ds(Env),units, Gtc=unsm(2)),
              data=DT)
summary(ans5)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
DT$EnvName <- paste(DT$Env,DT$Name)
ans6 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vs(us(Env),Name, Gtc=unsm(2)),
              rcov= ~ vs(ds(Env),units, Gtc=unsm(2)),
              data=DT)
summary(ans6)


## ------------------------------------------------------------------------
unsm(4)

## ------------------------------------------------------------------------
uncm(4)

## ------------------------------------------------------------------------
fixm(4)

## ------------------------------------------------------------------------
fcm(c(1,0,1,0))

## ------------------------------------------------------------------------
library(orthopolynom)
data(DT_legendre)
head(DT)
mRR2<-mmer(Y~ 1 + Xf
           , random=~ vs(us(leg(X,1)),SUBJECT)
           , rcov=~vs(units)
           , data=DT)
summary(mRR2)$varcomp


## ------------------------------------------------------------------------

data(DT_cpdata)
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
#### look at the data and fit the model
head(DT,3)
head(MP,3)
GT[1:3,1:4]
mix1 <- GWAS(color~1,
             random=~vs(id,Gu=A)
             + Rowf + Colf,
             rcov=~units,
             data=DT,
             M=GT, gTerm = "u:id")

ms <- as.data.frame(t(mix1$scores))
ms$Locus <- rownames(ms)
MP2 <- merge(MP,ms,by="Locus",all.x = TRUE);
manhattan(MP2, pch=20,cex=.5, PVCN = "color score")


## ------------------------------------------------------------------------

data("DT_halfdiallel")
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modh <- mmer(sugar~1, 
             random=~vs(overlay(femalef,malef)) 
             + genof,
             data=DT)


## ------------------------------------------------------------------------
data("DT_cpdata")
### mimic two fields
A <- A.mat(GT)
mix <- mmer(Yield~1,
            random=~vs(id, Gu=A) +
              vs(Rowf) +
              vs(Colf) +
              vs(spl2D(Row,Col)),
            rcov=~vs(units),
            data=DT)
summary(mix)

## ------------------------------------------------------------------------

data(DT_cpdata)
GT[1:4,1:4]
#### look at the data and fit the model
mix1 <- mmer(Yield~1,
              random=~vs(list(GT)),
              rcov=~units,
              data=DT)


## ------------------------------------------------------------------------
data("DT_wheat"); 
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
head(y.trn)
## GBLUP
ans <- mmer(X1~1,
            random=~vs(id,Gu=K), 
            rcov=~units, 
            data=y.trn) # kinship based
ans$U$`u:id`$X1 <- as.data.frame(ans$U$`u:id`$X1)
rownames(ans$U$`u:id`$X1) <- gsub("id","",rownames(ans$U$`u:id`$X1))
cor(ans$U$`u:id`$X1[vv,],DT[vv,"X1"], use="complete")

## rrBLUP
ans2 <- mmer(X1~1,
             random=~vs(list(GT)), 
             rcov=~units,
             data=y.trn) # kinship based

u <- GT %*% as.matrix(ans2$U$`u:GT`$X1) # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP

