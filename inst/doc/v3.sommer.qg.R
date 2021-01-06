## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example

ans1 <- mmer(Yield~1,
             random= ~ Name + Env + Env:Name + Env:Block,
             rcov= ~ units,
             data=DT, verbose = FALSE)
summary(ans1)$varcomp
(n.env <- length(levels(DT$Env)))
vpredict(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
DT$idd <-DT$id; DT$ide <-DT$id
### look at the data
A <- A.mat(GT) # additive relationship matrix
D <- D.mat(GT) # dominance relationship matrix
E <- E.mat(GT) # epistatic relationship matrix
ans.ADE <- mmer(color~1, 
                 random=~vs(id,Gu=A) + vs(idd,Gu=D), 
                 rcov=~units,
                 data=DT,verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans.ADE, h2 ~ (V1) / ( V1+V3) )
vpredict(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) )

## ---- fig.show='hold'---------------------------------------------------------
data(DT_cornhybrids)
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids
### fit the model
modFD <- mmer(Yield~1, 
               random=~ vs(at(Location,c("3","4")),GCA2), 
               rcov= ~ vs(ds(Location),units),
               data=DT, verbose = FALSE)
summary(modFD)

## -----------------------------------------------------------------------------
data(DT_cornhybrids)
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids
GT[1:4,1:4]
### fit the model
modFD <- mmer(Yield~1, 
              random=~ vs(at(Location,c("3","4")),GCA2,Gu=GT), 
              rcov= ~ vs(ds(Location),units),
              data=DT, verbose = FALSE)
summary(modFD)

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### look at the data
A <- A.mat(GT) # additive relationship matrix
ans <- mmer(color~1, 
                random=~vs(id,Gu=A), 
                rcov=~units,
                data=DT, verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans, h2 ~ (V1) / ( V1+V2) )


## -----------------------------------------------------------------------------
data(DT_btdata)
DT <- DT_btdata
mix3 <- mmer(cbind(tarsus, back) ~ sex,
               random = ~ vs(dam, Gtc=unsm(2)) + vs(fosternest,Gtc=diag(2)),
               rcov=~vs(units,Gtc=unsm(2)),
               data = DT, verbose = FALSE)
summary(mix3)
#### calculate the genetic correlation
vpredict(mix3, gen.cor ~ V2 / sqrt(V1*V3))

## -----------------------------------------------------------------------------
data(DT_cornhybrids)
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids

modFD <- mmer(Yield~Location, 
               random=~GCA1+GCA2+SCA, 
               rcov=~units,
               data=DT, verbose = FALSE)
(suma <- summary(modFD)$varcomp)
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

## -----------------------------------------------------------------------------
data("DT_halfdiallel")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modh <- mmer(sugar~1, 
             random=~vs(overlay(femalef,malef)) 
             + genof,
             data=DT, verbose = FALSE)
summary(modh)$varcomp

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
head(y.trn)
## GBLUP
ans <- mmer(X1~1,
            random=~vs(id,Gu=K), 
            rcov=~units, 
            data=y.trn, verbose = FALSE) # kinship based
ans$U$`u:id`$X1 <- as.data.frame(ans$U$`u:id`$X1)
rownames(ans$U$`u:id`$X1) <- gsub("id","",rownames(ans$U$`u:id`$X1))
cor(ans$U$`u:id`$X1[vv,],DT[vv,"X1"], use="complete")

## rrBLUP
ans2 <- mmer(X1~1,
             random=~vs(list(GT)), 
             rcov=~units,
             data=y.trn, verbose = FALSE) # kinship based

u <- GT %*% as.matrix(ans2$U$`u:GT`$X1) # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP

## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
Md <- Md_technow
Mf <- Mf_technow
Ad <- Ad_technow
Af <- Af_technow
# RUN THE PREDICTION MODEL
y.trn <- DT
vv1 <- which(!is.na(DT$GY))
vv2 <- sample(vv1, 100)
y.trn[vv2,"GY"] <- NA
anss2 <- mmer(GY~1, 
               random=~vs(dent,Gu=Ad) + vs(flint,Gu=Af), 
               rcov=~units,
               data=y.trn, verbose = FALSE) 
summary(anss2)$varcomp

zu1 <- model.matrix(~dent-1,y.trn) %*% anss2$U$`u:dent`$GY
zu2 <- model.matrix(~flint-1,y.trn) %*% anss2$U$`u:flint`$GY
u <- zu1+zu2+anss2$Beta[1,"Estimate"]
cor(u[vv2,], DT$GY[vv2])

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### mimic two fields
A <- A.mat(GT)
mix <- mmer(Yield~1,
            random=~vs(id, Gu=A) +
              vs(Rowf) +
              vs(Colf) +
              vs(spl2D(Row,Col)),
            rcov=~vs(units),
            data=DT, verbose = FALSE)
summary(mix)
# make a plot to observe the spatial effects found by the spl2D()
W <- with(DT,spl2D(Row,Col)) # 2D spline incidence matrix
DT$spatial <- W%*%mix$U$`u:Row`$Yield # 2D spline BLUPs
lattice::levelplot(spatial~Row*Col, data=DT) # plot the spatial effect by row and column

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
A <- A.mat(GT)
ans.m <- mmer(cbind(Yield,color)~1,
               random=~ vs(id, Gu=A, Gtc=unsm(2))
               + vs(Rowf,Gtc=diag(2))
               + vs(Colf,Gtc=diag(2)),
               rcov=~ vs(units, Gtc=unsm(2)),
               data=DT, verbose = FALSE)

## -----------------------------------------------------------------------------
cov2cor(ans.m$sigma$`u:id`)

## -----------------------------------------------------------------------------
library(sommer)
data("DT_cpdata")
DT <- DT_cpdata
M <- GT_cpdata

################
# MARKER MODEL
################
mix.marker <- mmer(color~1,
                   random=~Rowf+vs(M),
                   rcov=~units,data=DT,
                   verbose = FALSE)


me.marker <- mix.marker$U$`u:M`$color

################
# PARTITIONED GBLUP MODEL
################

MMT<-M%*%t(M) ## additive relationship matrix 
MMTinv<-solve(MMT) ## inverse
MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-

mix.part <- mmer(color~1,
                 random=~Rowf+vs(id, Gu=MMT),
                 rcov=~units,data=DT,
                 verbose = FALSE)

#convert BLUPs to marker effects me=M'(M'M)- u
me.part<-MTMMTinv%*%matrix(mix.part$U$`u:id`$color,ncol=1)

# compare marker effects between both models
plot(me.marker,me.part)


## -----------------------------------------------------------------------------

data("DT_wheat")
rownames(GT_wheat) <- rownames(DT_wheat)
G <- A.mat(GT_wheat)
Y <- data.frame(DT_wheat)

# make the decomposition
UD<-eigen(G) # get the decomposition: G = UDU'
U<-UD$vectors
D<-diag(UD$values)# This will be our new 'relationship-matrix'
rownames(D) <- colnames(D) <- rownames(G)
X<-model.matrix(~1, data=Y) # here: only one fixed effect (intercept)
UX<-t(U)%*%X # premultiply X and y by U' 
UY <- t(U) %*% as.matrix(Y) # multivariate

# dataset for decomposed model
DTd<-data.frame(id = rownames(G) ,UY, UX =UX[,1])
DTd$id<-as.character(DTd$id)

modeld <- mmer(cbind(X1,X2) ~ UX - 1, 
              random = ~vs(id,Gu=D), 
              rcov = ~vs(units),
              data=DTd, verbose = FALSE)

# dataset for normal model
DTn<-data.frame(id = rownames(G) , DT_wheat)
DTn$id<-as.character(DTn$id)

modeln <- mmer(cbind(X1,X2) ~ 1, 
              random = ~vs(id,Gu=G), 
              rcov = ~vs(units),
              data=DTn, verbose = FALSE)

## compare regular and transformed blups
plot(x=(solve(t(U)))%*%modeld$U$`u:id`$X2[colnames(D)], 
     y=modeln$U$`u:id`$X2[colnames(D)], xlab="UDU blup",
     ylab="blup")


## -----------------------------------------------------------------------------

data(DT_expdesigns)
DT <- DT_expdesigns$car1
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
levelplot(yield~male*female|set, data=DT, main="NC design I")
##############################
## Expected Mean Square method
##############################
mix1 <- lm(yield~ setf + setf:repf + femalef:malef:setf + malef:setf, data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=2
nfem=2
Vfm <- (ms2-mse)/nrep
Vm <- (ms1-ms2)/(nrep*nfem)

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg
##############################
## REML method
##############################
mix2 <- mmer(yield~ setf + setf:repf,
            random=~femalef:malef:setf + malef:setf, 
            data=DT, verbose = FALSE)
vc <- summary(mix2)$varcomp; vc
Vfm <- vc[1,"VarComp"]
Vm <- vc[2,"VarComp"]

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg


## -----------------------------------------------------------------------------
DT <- DT_expdesigns$car2
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
levelplot(yield~male*female|set, data=DT, main="NC desing II")
head(DT)

N=with(DT,table(female, male, set))
nmale=length(which(N[1,,1] > 0))
nfemale=length(which(N[,1,1] > 0))
nrep=table(N[,,1])
nrep=as.numeric(names(nrep[which(names(nrep) !=0)]))

##############################
## Expected Mean Square method
##############################

mix1 <- lm(yield~ setf + setf:repf + 
             femalef:malef:setf + malef:setf + femalef:setf, data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef","Mean Sq"]
ms3 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=length(unique(DT$rep))
nfem=length(unique(DT$female))
nmal=length(unique(DT$male))
Vfm <- (ms3-mse)/nrep; 
Vf <- (ms2-ms3)/(nrep*nmale); 
Vm <- (ms1-ms3)/(nrep*nfemale); 

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg

##############################
## REML method
##############################

mix2 <- mmer(yield~ setf + setf:repf ,
            random=~femalef:malef:setf + malef:setf + femalef:setf, 
            data=DT, verbose = FALSE)
vc <- summary(mix2)$varcomp; vc
Vfm <- vc[1,"VarComp"]
Vm <- vc[2,"VarComp"]
Vf <- vc[3,"VarComp"]

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg


## -----------------------------------------------------------------------------

data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
#### look at the data and fit the model
mix1 <- mmer(Yield~1,
              random=~vs(id,Gu=A),
              rcov=~units,
              data=DT, verbose = FALSE)

####=========================================####
#### adding dominance and forcing the other VC's
####=========================================####
DT$idd <- DT$id;
D <- D.mat(GT) # dominance relationship matrix
mm <- matrix(3,1,1) ## matrix to fix the var comp

mix2 <- mmer(Yield~1,
              random=~vs(id, Gu=A, Gti=mix1$sigma_scaled$`u:id`, Gtc=mm)
                      + vs(idd, Gu=D, Gtc=unsm(1)),
              rcov=~vs(units,Gti=mix1$sigma_scaled$units, Gtc=mm),
              data=DT, verbose = FALSE)

# analyze variance components
summary(mix1)$varcomp
summary(mix2)$varcomp


