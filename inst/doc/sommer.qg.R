## -----------------------------------------------------------------------------
library(sommer)
data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

ans1 <- mmes(Yield~1,
             random= ~ Name + Env + Env:Name + Env:Block,
             rcov= ~ units, nIters=10,
             data=DT, verbose = FALSE)
summary(ans1)$varcomp
(n.env <- length(levels(DT$Env)))
vpredict(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

## -----------------------------------------------------------------------------
data(DT_cpdata, package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
DT$idd <-DT$id; DT$ide <-DT$id
### look at the data
A <- A.mat(GT) # additive relationship matrix
D <- D.mat(GT) # dominance relationship matrix
E <- E.mat(GT) # epistatic relationship matrix
ans.ADE <- mmes(color~1, 
                 random=~vsm(ism(id),Gu=A) + vsm(ism(idd),Gu=D), 
                 rcov=~units, nIters=10,
                 data=DT,verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans.ADE, h2 ~ (V1) / ( V1+V3) ) # narrow sense
vpredict(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense

## ----fig.show='hold'----------------------------------------------------------
# data(DT_cornhybrids, package="enhancer")
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- GT_cornhybrids
# ### fit the model
# modFD <- mmes(Yield~1, 
#               random=~ vsm(atr(Location,c("3","4")),ism(GCA2)), 
#               rcov= ~ vsm(dsm(Location),ism(units)), nIters=10,
#               returnParam = F,
#               data=DT, verbose = FALSE)
# summary(modFD)

## -----------------------------------------------------------------------------
# data(DT_cornhybrids, package="enhancer")
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- as(as(as( GT_cornhybrids,  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
# GT[1:4,1:4]
# DT=DT[with(DT, order(Location)), ]
# ### fit the model
# modFD <- mmes(Yield~1, 
#               random=~ vsm(atr(Location,c("3","4")),ism(GCA2),Gu=GT), 
#               rcov= ~ vsm(dsm(Location),ism(units)), nIters=10,
#               data=DT, verbose = FALSE)
# summary(modFD)

## -----------------------------------------------------------------------------
data(DT_cpdata, package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### look at the data
A <- A.mat(GT) # additive relationship matrix
ans <- mmes(color~1, 
                random=~vsm(ism(id),Gu=A), 
                rcov=~units, nIters=10,
                data=DT, verbose = FALSE)
(summary(ans.ADE)$varcomp)
vpredict(ans, h2 ~ (V1) / ( V1+V2) )


## -----------------------------------------------------------------------------
# data(DT_cornhybrids, package="enhancer")
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- GT_cornhybrids
# 
# modFD <- mmes(Yield~Location, 
#               random=~GCA1+GCA2+SCA, 
#               rcov=~units, nIters=10,
#               data=DT, verbose = FALSE)
# (suma <- summary(modFD)$varcomp)
# Vgca <- sum(suma[1:2,1])
# Vsca <- suma[3,1]
# Ve <- suma[4,1]
# Va = 4*Vgca
# Vd = 4*Vsca
# Vg <- Va + Vd
# (H2 <- Vg / (Vg + (Ve)) )
# (h2 <- Va / (Vg + (Ve)) )

## -----------------------------------------------------------------------------
data("DT_halfdiallel", package="enhancer")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modh <- mmes(sugar~1, 
             random=~vsm(ism(overlay(femalef,malef)) )
             + genof, nIters=10,
             data=DT, verbose = FALSE)
summary(modh)$varcomp

## -----------------------------------------------------------------------------
data(DT_wheat, package="enhancer")
DT <- DT_wheat
GT <- apply(GT_wheat,2,as.numeric)
rownames(GT) <- rownames(GT_wheat)

colnames(DT) <- paste0("X",1:ncol(DT))
DT <- as.data.frame(DT);DT$id <- as.factor(rownames(DT))
# select environment 1
K <- A.mat(GT) # additive relationship matrix
colnames(K) <- rownames(K) <- rownames(DT)
# GBLUP pedigree-based approach
set.seed(12345)
y.trn <- DT
vv <- sample(rownames(DT),round(nrow(DT)/5))
y.trn[vv,"X1"] <- NA
head(y.trn)
## GBLUP
ans <- mmes(X1~1,
            random=~vsm(ism(id),Gu=K), 
            rcov=~units,nIters=10,
            data=y.trn, verbose = FALSE) # kinship based
cor(ans$u[vv,] ,DT[vv,"X1"], use="complete")

## rrBLUP
ans2 <- mmes(X1~1,
             random=~vsm(ism(GT), buildGu = FALSE), 
             rcov=~units, getPEV = TRUE, nIters=10, # you do more iterations
             data=y.trn, verbose = FALSE) # kinship based

u <- GT %*% as.matrix(ans2$uList$`vsm(ism(GT), buildGu = FALSE`) # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP

## -----------------------------------------------------------------------------
# data(DT_ige, package="enhancer")
# DT <- DT_ige
# Af <- A_ige
# An <- A_ige

## Direct genetic effects model
# modDGE <- mmes(trait ~ block,
#                random = ~ focal,
#                rcov = ~ units, nIters=30,
#                data = DT, verbose=FALSE)
# summary(modDGE)$varcomp


## -----------------------------------------------------------------------------
# data(DT_ige, package="enhancer")
# DT <- DT_ige
# A <- A_ige
# 
# ## Indirect genetic effects model
# modIGE <- mmes(trait ~ block, dateWarning = FALSE,
#                random = ~ focal + neighbour, verbose = FALSE,
#                rcov = ~ units, nIters=100,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

# ### Indirect genetic effects model
# modIGE <- mmes(trait ~ block, dateWarning = FALSE,
#                random = ~ covm( vsm(ism(focal)), vsm(ism(neighbour)) ),
#                rcov = ~ units, nIters=100, verbose = FALSE,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

### Indirect genetic effects model
# Ai <- solve(A_ige + diag(1e-5, nrow(A_ige),nrow(A_ige) ))
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# # Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
# modIGE <- mmes(trait ~ block, dateWarning = FALSE,
#                random = ~ covm( vsm(ism(focal), Gu=Ai), vsm(ism(neighbour), Gu=Ai) ),
#                rcov = ~ units, nIters=100, verbose = FALSE,
#               data = DT)
# summary(modIGE)$varcomp


## -----------------------------------------------------------------------------
data(DT_technow, package="enhancer")
DT <- DT_technow

Md <- apply(Md_technow,2,as.numeric)
rownames(Md) <- rownames(Md_technow)
Mf <- apply(Mf_technow,2,as.numeric)
rownames(Mf) <- rownames(Mf_technow)

Md <- (Md*2) - 1
Mf <- (Mf*2) - 1
Ad <- A.mat(Md)
Af <- A.mat(Mf)
Adi <- solve(Ad + diag(1e-4,ncol(Ad),ncol(Ad)))
Adi <- as(as(as( Adi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Adi, 'inverse')=TRUE
Afi <- solve(Af + diag(1e-4,ncol(Af),ncol(Af)))
Afi <- as(as(as( Afi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Afi, 'inverse')=TRUE
# RUN THE PREDICTION MODEL
y.trn <- DT
vv1 <- which(!is.na(DT$GY))
vv2 <- sample(vv1, 100)
y.trn[vv2,"GY"] <- NA
anss2 <- mmes(GY~1,  henderson=TRUE,
              random=~vsm(ism(dent),Gu=Adi) + vsm(ism(flint),Gu=Afi), 
              rcov=~units, nIters=15,
              data=y.trn, verbose = FALSE) 
summary(anss2)$varcomp

# zu1 <- model.matrix(~dent-1,y.trn) %*% anss2$uList$`vsm(ism(dent), Gu = Adi)`
# zu2 <- model.matrix(~flint-1,y.trn) %*% anss2$uList$`vsm(ism(flint), Gu = Afi)`
# u <- zu1+zu2+as.vector(anss2$b)
# cor(u[vv2,], DT$GY[vv2])

## -----------------------------------------------------------------------------
# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# traits <- c("color","Yield")
# DT[,traits] <- apply(DT[,traits],2,scale)
# DTL <- reshape(DT[,c("id", traits)],
#                idvar = c("id"),
#                varying = traits,
#                v.names = "value", direction = "long",
#                timevar = "trait", times = traits )
# DTL <- DTL[with(DTL, order(trait)), ]
# head(DTL)
# 
# A <- A.mat(GT) # additive relationship matrix
# # if using mmes=TRUE you need to provide the inverse
# Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, 'inverse')=TRUE
# #### be patient take some time
# ansm <- mmes( value ~ trait, # henderson=TRUE,
#                random=~ vsm(usm(trait), ism(id), Gu=A), # Ai if henderson
#                rcov=~ vsm(dsm(trait), ism(units)),
#                data=DTL)
# cov2cor(ansm$theta[[1]])

## -----------------------------------------------------------------------------
library(sommer)
data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
M <- GT_cpdata

################
# MARKER MODEL
################
mix.marker <- mmes(color~1,
                   random=~Rowf+vsm(ism(M)),
                   rcov=~units,data=DT, 
                   verbose = FALSE)


me.marker <- mix.marker$uList$`vsm(ism(M`

################
# PARTITIONED GBLUP MODEL
################

MMT <-tcrossprod(M) ## MM' = additive relationship matrix 
MMTinv<-solve(MMT) ## inverse
MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-

mix.part <- mmes(color~1,
                 random=~Rowf+vsm(ism(id), Gu=MMT),
                 rcov=~units,data=DT,
                 verbose = FALSE)

#convert BLUPs to marker effects me=M'(M'M)- u
me.part<-MTMMTinv%*%matrix(mix.part$uList$`vsm(ism(id), Gu = MMT`,ncol=1)

# compare marker effects between both models
plot(me.marker,me.part)


## -----------------------------------------------------------------------------

# data("DT_wheat", package="enhancer")
# rownames(GT_wheat) <- rownames(DT_wheat)
# GT <- apply(GT_wheat,2,as.numeric)
# rownames(GT) <- rownames(GT_wheat)
# G <- A.mat(GT)
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
# modeld <- mmes(cbind(X1,X2) ~ UX - 1, 
#               random = ~vsm(id,Gu=D), 
#               rcov = ~vsm(units),
#               data=DTd, verbose = FALSE)
# 
# # dataset for normal model
# DTn<-data.frame(id = rownames(G) , DT_wheat)
# DTn$id<-as.character(DTn$id)
# 
# modeln <- mmes(cbind(X1,X2) ~ 1, 
#               random = ~vsm(id,Gu=G), 
#               rcov = ~vsm(units),
#               data=DTn, verbose = FALSE)
# 
# ## compare regular and transformed blups
# plot(x=(solve(t(U)))%*%modeld$U$`u:id`$X2[colnames(D)], 
#      y=modeln$U$`u:id`$X2[colnames(D)], xlab="UDU blup",
#      ylab="blup")


## -----------------------------------------------------------------------------

data(DT_expdesigns, package="enhancer")
DT <- DT_expdesigns$car1
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
# lattice::levelplot(yield~male * malef:femalef|set, data=DT, main="NC design I")
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
mix2 <- mmes(yield~ setf + setf:repf,
            random=~femalef:malef:setf + malef:setf, nIters=10,
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

mix2 <- mmes(yield~ setf + setf:repf ,
            random=~femalef:malef:setf + malef:setf + femalef:setf, 
            nIters=10,
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
data(DT_cpdata, package="enhancer")
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
# mix2 <- GWAS(color~1,
#              random=~vsm(ism(id), Gu=A) + Rowf + Colf,
#              rcov=~units, M=GT, gTerm = "u:id",
#              verbose = FALSE, nIters=10,
#              data=DT)

## -----------------------------------------------------------------------------
###########################
#### GWAS by RRBLUP approach
###########################
# Z <- GT[as.character(DT$id),]
# mixRRBLUP <- mmes(color~1,
#               random=~vsm(ism(Z)) + Rowf + Colf,
#               rcov=~units, nIters=10,
#               verbose = FALSE,
#               data=DT)
# 
# a <- mixRRBLUP$uList$`vsm(ism(Z`# marker effects
# start=mixRRBLUP$partitions[[1]][1]
# end=mixRRBLUP$partitions[[1]][2]
# se.a <- sqrt( diag(kronecker(diag(ncol(Z)),mixRRBLUP$theta[[1]]) - mixRRBLUP$Ci[start:end,start:end] ) ) # SE of marker effects
# t.stat <- a/se.a # t-statistic
# pvalRRBLUP <- dt(t.stat[,1],df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
# ###########################
# #### GWAS by GBLUP approach
# ###########################
# M<- GT
# MMT <-tcrossprod(M) ## MM' = additive relationship matrix
# MMTinv<-solve(MMT + diag(1e-6, ncol(MMT), ncol(MMT))) ## inverse of MM'
# MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
# mixGBLUP <- mmes(color~1,
#              random=~vsm(ism(id), Gu=MMT) + Rowf + Colf,
#              rcov=~units, nIters=10,
#              verbose = FALSE,
#              data=DT)
# a.from.g <-MTMMTinv%*%matrix(mixGBLUP$uList$`vsm(ism(id), Gu = MMT`,ncol=1)
# start=mixGBLUP$partitions[[1]][1]
# end=mixGBLUP$partitions[[1]][2]
# var.g <- kronecker(MMT,mixGBLUP$theta[[1]]) - mixGBLUP$Ci[start:end,start:end]
# var.a.from.g <- t(M)%*%MMTinv%*% (var.g) %*% t(MMTinv)%*%M
# se.a.from.g <- sqrt(diag(var.a.from.g))
# t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
# pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
###########################
#### Compare results
###########################
# plot(mix2$scores[,1], main="GWAS")
# plot(-log(pvalRRBLUP), main="GWAS by RRBLUP/SNP-BLUP") 
# plot(-log(pvalGBLUP), main="GWAS by GBLUP")


