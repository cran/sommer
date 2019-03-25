## ------------------------------------------------------------------------
library(sommer)
data(DT_example)
head(DT)

ans1 <- mmer(Yield~1,
             random= ~ Name + Env + Env:Name + Env:Block,
             rcov= ~ units,
             data=DT)
summary(ans1)$varcomp
(n.env <- length(levels(DT$Env)))
pin(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

## ------------------------------------------------------------------------
data("DT_cpdata")
DT$idd <-DT$id; DT$ide <-DT$id
### look at the data
A <- A.mat(GT) # additive relationship matrix
D <- D.mat(GT) # dominance relationship matrix
E <- E.mat(GT) # epistatic relationship matrix
ans.ADE <- mmer(color~1, 
                 random=~vs(id,Gu=A) + vs(idd,Gu=D), 
                 rcov=~units,
                 data=DT)
(summary(ans.ADE)$varcomp)
pin(ans.ADE, h2 ~ (V1) / ( V1+V3) )
pin(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) )

## ---- fig.show='hold'----------------------------------------------------
data("DT_cornhybrids")
### fit the model
modFD <- mmer(Yield~1, 
               random=~ vs(at(Location,c("3","4")),GCA2), 
               rcov= ~ vs(ds(Location),units),
               data=DT)
summary(modFD)

## ------------------------------------------------------------------------
data("DT_cornhybrids")
GT[1:4,1:4]
### fit the model
modFD <- mmer(Yield~1, 
              random=~ vs(at(Location,c("3","4")),GCA2,Gu=GT), 
              rcov= ~ vs(ds(Location),units),
              data=DT)
summary(modFD)

## ------------------------------------------------------------------------
data("DT_cpdata")
### look at the data
A <- A.mat(GT) # additive relationship matrix
ans <- mmer(color~1, 
                random=~vs(id,Gu=A), 
                rcov=~units,
                data=DT)
(summary(ans.ADE)$varcomp)
pin(ans, h2 ~ (V1) / ( V1+V2) )

## ------------------------------------------------------------------------
data("DT_cornhybrids")

modFD <- mmer(Yield~Location, 
               random=~GCA1+GCA2+SCA, 
               rcov=~units,
               data=DT)
(suma <- summary(modFD)$varcomp)
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

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
summary(modh)$varcomp

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

## ------------------------------------------------------------------------
data("DT_technow")
# RUN THE PREDICTION MODEL
y.trn <- DT
vv1 <- which(!is.na(DT$GY))
vv2 <- sample(vv1, 100)
y.trn[vv2,"GY"] <- NA
anss2 <- mmer(GY~1, 
               random=~vs(dent,Gu=Ad) + vs(flint,Gu=Af), 
               rcov=~units,
               data=y.trn) 
summary(anss2)$varcomp

zu1 <- model.matrix(~dent-1,y.trn) %*% anss2$U$`u:dent`$GY
zu2 <- model.matrix(~flint-1,y.trn) %*% anss2$U$`u:flint`$GY
u <- zu1+zu2+anss2$Beta[1,"Estimate"]
cor(u[vv2,], DT$GY[vv2])

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
data("DT_cpdata")
A <- A.mat(GT)
ans.m <- mmer(cbind(Yield,color)~1,
               random=~ vs(id, Gu=A, Gtc=unsm(2))
               + vs(Rowf,Gtc=diag(2))
               + vs(Colf,Gtc=diag(2)),
               rcov=~ vs(units, Gtc=unsm(2)),
               data=DT)

## ------------------------------------------------------------------------
cov2cor(ans.m$sigma$`u:id`)

