## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example

ans1 <- mmec(Yield~1,
             random= ~ Name + Env + Env:Name + Env:Block,
             rcov= ~ units, nIters=3,
             data=DT, verbose = FALSE)
summary(ans1)$varcomp
(n.env <- length(levels(DT$Env)))
# vpredict(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

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
                 random=~vsr(id,Gu=A) + vsr(idd,Gu=D), 
                 rcov=~units, nIters=3,
                 data=DT,verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans.ADE, h2 ~ (V1) / ( V1+V3) ) # narrow sense
vpredict(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense

## ----fig.show='hold'----------------------------------------------------------
# data(DT_cornhybrids)
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- GT_cornhybrids
# ### fit the model
# modFD <- mmec(Yield~1, 
#               random=~ vsc(atc(Location,c("3","4")),isc(GCA2)), 
#               rcov= ~ vsc(dsc(Location),isc(units)), nIters=3,
#               returnParam = F,
#               data=DT, verbose = FALSE)
# summary(modFD)

## -----------------------------------------------------------------------------
# data(DT_cornhybrids)
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- as(GT_cornhybrids, Class = "dgCMatrix")
# GT[1:4,1:4]
# DT=DT[with(DT, order(Location)), ]
# ### fit the model
# modFD <- mmec(Yield~1, 
#               random=~ vsc(atc(Location,c("3","4")),isc(GCA2),Gu=GT), 
#               rcov= ~ vsc(dsc(Location),isc(units)), nIters=3,
#               data=DT, verbose = FALSE)
# summary(modFD)

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### look at the data
A <- A.mat(GT) # additive relationship matrix
ans <- mmer(color~1, 
                random=~vsr(id,Gu=A), 
                rcov=~units, nIters=3,
                data=DT, verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans, h2 ~ (V1) / ( V1+V2) )


## -----------------------------------------------------------------------------
## just silenced to avoid too much time building the vignettes
# data(DT_btdata)
# DT <- DT_btdata
# mix3 <- mmer(cbind(tarsus, back) ~ sex,
#                random = ~ vsr(dam, Gtc=unsm(2)) + vsr(fosternest,Gtc=diag(2)),
#                rcov=~vsr(units,Gtc=unsm(2)), nIters=3,
#                data = DT, verbose = FALSE)
# summary(mix3)
# #### calculate the genetic correlation
# vpredict(mix3, gen.cor ~ V2 / sqrt(V1*V3))

## -----------------------------------------------------------------------------
data(DT_cornhybrids)
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids

modFD <- mmec(Yield~Location, 
              random=~GCA1+GCA2+SCA, 
              rcov=~units, nIters=3,
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
modh <- mmec(sugar~1, 
             random=~vsc(isc(overlay(femalef,malef)) )
             + genof, nIters=3,
             data=DT, verbose = FALSE)
summary(modh)$varcomp

## -----------------------------------------------------------------------------
data(DT_wheat)
DT <- DT_wheat
GT <- GT_wheat[,1:200]
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
            random=~vsr(id,Gu=K), 
            rcov=~units,nIters=3,
            data=y.trn, verbose = FALSE) # kinship based
ans$U$`u:id`$X1 <- as.data.frame(ans$U$`u:id`$X1)
rownames(ans$U$`u:id`$X1) <- gsub("id","",rownames(ans$U$`u:id`$X1))
cor(ans$U$`u:id`$X1[vv,],DT[vv,"X1"], use="complete")

## rrBLUP
ans2 <- mmer(X1~1,
             random=~vsr(list(GT), buildGu = FALSE), 
             rcov=~units, getPEV = FALSE, nIters=3,
             data=y.trn, verbose = FALSE) # kinship based

u <- GT %*% as.matrix(ans2$U$`u:GT`$X1) # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP

## -----------------------------------------------------------------------------
# data(DT_ige)
# DT <- DT_ige
# Af <- A_ige
# An <- A_ige

## Direct genetic effects model
# modDGE <- mmec(trait ~ block,
#                random = ~ focal,
#                rcov = ~ units, nIters=3,
#                data = DT, verbose=FALSE)
# summary(modDGE)$varcomp


## -----------------------------------------------------------------------------
# data(DT_ige)
# DT <- DT_ige
# A <- A_ige
# 
# ## Indirect genetic effects model
# modIGE <- mmec(trait ~ block, dateWarning = FALSE,
#                random = ~ focal + neighbour, verbose = FALSE,
#                rcov = ~ units, nIters=100,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

# ### Indirect genetic effects model
# modIGE <- mmec(trait ~ block, dateWarning = FALSE,
#                random = ~ covc( vsc(isc(focal)), vsc(isc(neighbour)) ),
#                rcov = ~ units, nIters=100, verbose = FALSE,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

### Indirect genetic effects model
# Ai <- as( solve(A_ige + diag(1e-5, nrow(A_ige),nrow(A_ige) )), Class="dgCMatrix")
# # Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
# modIGE <- mmec(trait ~ block, dateWarning = FALSE,
#                random = ~ covc( vsc(isc(focal), Gu=Ai), vsc(isc(neighbour), Gu=Ai) ),
#                rcov = ~ units, nIters=100, verbose = FALSE,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
Md <- (Md_technow*2) - 1
Mf <- (Mf_technow*2) - 1
Ad <- A.mat(Md)
Af <- A.mat(Mf)
Adi <- as(solve(Ad + diag(1e-4,ncol(Ad),ncol(Ad))), Class="dgCMatrix")
Afi <- as(solve(Af + diag(1e-4,ncol(Af),ncol(Af))), Class="dgCMatrix")
# RUN THE PREDICTION MODEL
y.trn <- DT
vv1 <- which(!is.na(DT$GY))
vv2 <- sample(vv1, 100)
y.trn[vv2,"GY"] <- NA
anss2 <- mmec(GY~1, 
              random=~vsc(isc(dent),Gu=Adi) + vsc(isc(flint),Gu=Afi), 
              rcov=~units, nIters=15,
              data=y.trn, verbose = FALSE) 
summary(anss2)$varcomp

# zu1 <- model.matrix(~dent-1,y.trn) %*% anss2$uList$`vsc(isc(dent), Gu = Adi)`
# zu2 <- model.matrix(~flint-1,y.trn) %*% anss2$uList$`vsc(isc(flint), Gu = Afi)`
# u <- zu1+zu2+as.vector(anss2$b)
# cor(u[vv2,], DT$GY[vv2])

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
              spl2Da(Row,Col), nIters=3,
            rcov=~vsr(units),
            data=DT, verbose = FALSE)
summary(mix)
# make a plot to observe the spatial effects found by the spl2D()
W <- with(DT,spl2Da(Row,Col)) # 2D spline incidence matrix
DT$spatial <- W$Z$`A:all`%*%mix$U$`A:all`$Yield # 2D spline BLUPs
# lattice::levelplot(spatial~Row*Col, data=DT) # plot the spatial effect by row and column

## -----------------------------------------------------------------------------
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# A <- A.mat(GT)
# ans.m <- mmer(cbind(Yield,color)~1,
#                random=~ vsr(id, Gu=A, Gtc=unsm(2))
#                + vsr(Rowf,Gtc=diag(2))
#                + vsr(Colf,Gtc=diag(2)),
#                rcov=~ vsr(units, Gtc=unsm(2)), nIters=3,
#                data=DT, verbose = FALSE)

## -----------------------------------------------------------------------------
# cov2cor(ans.m$sigma$`u:id`)

## -----------------------------------------------------------------------------
library(sommer)
data("DT_cpdata")
DT <- DT_cpdata
M <- GT_cpdata

################
# MARKER MODEL
################
mix.marker <- mmer(color~1,
                   random=~Rowf+vsr(M),
                   rcov=~units,data=DT, 
                   verbose = FALSE)


me.marker <- mix.marker$U$`u:M`$color

################
# PARTITIONED GBLUP MODEL
################

MMT <-tcrossprod(M) ## MM' = additive relationship matrix 
MMTinv<-solve(MMT) ## inverse
MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-

mix.part <- mmer(color~1,
                 random=~Rowf+vsr(id, Gu=MMT),
                 rcov=~units,data=DT,
                 verbose = FALSE)

#convert BLUPs to marker effects me=M'(M'M)- u
me.part<-MTMMTinv%*%matrix(mix.part$U$`u:id`$color,ncol=1)

# compare marker effects between both models
plot(me.marker,me.part)


## -----------------------------------------------------------------------------

# data("DT_wheat")
# rownames(GT_wheat) <- rownames(DT_wheat)
# G <- A.mat(GT_wheat)
# Y <- data.frame(DT_wheat)
# 
# # make the decomposition
# UD<-eigen(G) # get the decomposition: G = UDU'
# U<-UD$vectors
# D<-diag(UD$values)# This will be our new 'relationship-matrix'
# rownames(D) <- colnames(D) <- rownames(G)
# X<-model.matrix(~1, data=Y) # here: only one fixed effect (intercept)
# UX<-t(U)%*%X # premultiply X and y by U' 
# UY <- t(U) %*% as.matrix(Y) # multivariate
# 
# # dataset for decomposed model
# DTd<-data.frame(id = rownames(G) ,UY, UX =UX[,1])
# DTd$id<-as.character(DTd$id)
# 
# modeld <- mmer(cbind(X1,X2) ~ UX - 1, 
#               random = ~vsr(id,Gu=D), 
#               rcov = ~vsr(units),
#               data=DTd, verbose = FALSE)
# 
# # dataset for normal model
# DTn<-data.frame(id = rownames(G) , DT_wheat)
# DTn$id<-as.character(DTn$id)
# 
# modeln <- mmer(cbind(X1,X2) ~ 1, 
#               random = ~vsr(id,Gu=G), 
#               rcov = ~vsr(units),
#               data=DTn, verbose = FALSE)
# 
# ## compare regular and transformed blups
# plot(x=(solve(t(U)))%*%modeld$U$`u:id`$X2[colnames(D)], 
#      y=modeln$U$`u:id`$X2[colnames(D)], xlab="UDU blup",
#      ylab="blup")


## -----------------------------------------------------------------------------

data(DT_expdesigns)
DT <- DT_expdesigns$car1
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
#levelplot(yield~male*female|set, data=DT, main="NC design I")
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
            random=~femalef:malef:setf + malef:setf, nIters=3,
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
#levelplot(yield~male*female|set, data=DT, main="NC desing II")
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
            nIters=3,
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
GT <- GT_cpdata[,1:200]
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
n <- nrow(DT) # to be used for degrees of freedom
k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)

## -----------------------------------------------------------------------------
###########################
#### Regular GWAS/EMMAX approach
###########################
mix2 <- GWAS(color~1,
             random=~vsr(id, Gu=A) + Rowf + Colf,
             rcov=~units, M=GT, gTerm = "u:id",
             verbose = FALSE, nIters=3,
             data=DT)

## -----------------------------------------------------------------------------
###########################
#### GWAS by RRBLUP approach
###########################
Z <- GT[as.character(DT$id),]
mixRRBLUP <- mmer(color~1,
              random=~vsr(Z) + Rowf + Colf,
              rcov=~units, nIters=3,
              verbose = FALSE,
              data=DT)

a <- mixRRBLUP$U$`u:Z`$color # marker effects
se.a <- sqrt(diag(kronecker(diag(ncol(Z)),mixRRBLUP$sigma$`u:Z`) - mixRRBLUP$PevU$`u:Z`$color)) # SE of marker effects
t.stat <- a/se.a # t-statistic
pvalRRBLUP <- dt(t.stat,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
###########################
#### GWAS by GBLUP approach
###########################
M<- GT
MMT <-tcrossprod(M) ## MM' = additive relationship matrix
MMTinv<-solve(MMT + diag(1e-6, ncol(MMT), ncol(MMT))) ## inverse of MM'
MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
mixGBLUP <- mmer(color~1,
             random=~vsr(id, Gu=MMT) + Rowf + Colf,
             rcov=~units, nIters=3,
             verbose = FALSE,
             data=DT)
a.from.g <-MTMMTinv%*%matrix(mixGBLUP$U$`u:id`$color,ncol=1)
var.g <- kronecker(MMT,mixGBLUP$sigma$`u:id`) - mixGBLUP$PevU$`u:id`$color
var.a.from.g <- t(M)%*%MMTinv%*% (var.g) %*% t(MMTinv)%*%M
se.a.from.g <- sqrt(diag(var.a.from.g))
t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
###########################
#### Compare results
###########################
# plot(mix2$scores[,1], main="GWAS")
plot(-log(pvalRRBLUP), main="GWAS by RRBLUP/SNP-BLUP") 
plot(-log(pvalGBLUP), main="GWAS by GBLUP")


