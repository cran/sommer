## ------------------------------------------------------------------------
library(sommer)
data(h2)
head(h2)

ans1 <- mmer2(y~1, 
              random = ~Name + Env + Name:Env + Block,
              rcov = ~units,
              data=h2, silent = TRUE)

suma <- summary(ans1)
n.env <- length(levels(h2$Env))
pin(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )


## ------------------------------------------------------------------------
library(sommer)
data(h2)
head(h2)
Z1 <- model.matrix(~Name-1, h2)
Z2 <- model.matrix(~Env-1, h2)
Z3 <- model.matrix(~Env:Name-1, h2)
Z4 <- model.matrix(~Block-1, h2)
ETA <- list(name=list(Z=Z1),env=list(Z=Z2),name.env=list(Z=Z3),block=list(Z=Z4))
y <- h2$y
ans1 <- mmer(Y=y, Z=ETA, silent = TRUE)
vc <- ans1$var.comp

## ------------------------------------------------------------------------
data(CPdata)
CPpheno$idd <-CPpheno$id; CPpheno$ide <-CPpheno$id
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
## fit a model including additive and dominance effects
A <- A.mat(CPgeno) # additive relationship matrix
D <- D.mat(CPgeno) # dominance relationship matrix
E <- E.mat(CPgeno) # epistatic relationship matrix

ans.ADE <- mmer2(color~1, 
                 random=~g(id) + g(idd) + g(ide), 
                 rcov=~units,
                 G=list(id=A,idd=D,ide=E), 
                 silent = TRUE, data=CPpheno)
suma <- summary(ans.ADE)$var.comp.table
(H2 <- sum(suma[1:3,1])/sum(suma[,1]))
(h2 <- sum(suma[1,1])/sum(suma[,1]))


## ------------------------------------------------------------------------
data(CPdata)
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
## fit a model including additive and dominance effects
Z1 <- model.matrix(~id-1, CPpheno); colnames(Z1) <- gsub("id","",colnames(Z1))
A <- A.mat(CPgeno) # additive relationship matrix
D <- D.mat(CPgeno) # dominance relationship matrix
E <- E.mat(CPgeno) # epistatic relationship matrix
y <- CPpheno$color

ETA <- list(id=list(Z=Z1,K=A),idd=list(Z=Z1,K=D),ide=list(Z=Z1,K=E))
ans.ADE <- mmer(Y=y, Z=ETA, silent = TRUE)
ans.ADE$var.comp

## ---- fig.show='hold'----------------------------------------------------
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
head(hybrid2)
### fit the model
modFD <- mmer2(Yield~1, 
               random=~ at(Location,c("3","4")):GCA2, 
               rcov= ~ at(Location):units,
               data=hybrid2, silent = TRUE)
summary(modFD)

## ------------------------------------------------------------------------
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
## get the covariance structure for GCA2
A <- cornHybrid$K
## fit the model
modFD <- mmer2(Yield~1, 
               random=~ g(GCA2) + at(Location):g(GCA2), 
               rcov= ~ at(Location):units,
               data=hybrid2, G=list(GCA2=A),
               silent = TRUE, draw=FALSE)
summary(modFD)

## ------------------------------------------------------------------------
data(CPdata)
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
mix1 <- mmer2(color~1,
              random=~g(id), 
              rcov=~units,
              G=list(id=A), data=CPpheno, silent=TRUE)
summary(mix1)
#### run the pin function
pin(mix1, h2 ~ V1 / ( V1 + V2 ) )

## ------------------------------------------------------------------------
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
head(hybrid2)

modFD <- mmer2(Yield~Location, 
               random=~GCA1+GCA2+SCA, 
               rcov=~units,
               data=hybrid2,silent = TRUE, draw=FALSE)
(suma <- summary(modFD))
Vgca <- sum(suma$var.comp.table[1:2,1])
Vsca <- suma$var.comp.table[3,1]
Ve <- suma$var.comp.table[4,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

## ------------------------------------------------------------------------
data(HDdata)
head(HDdata)
HDdata$geno <- as.factor(HDdata$geno)
HDdata$male <- as.factor(HDdata$male)
HDdata$female <- as.factor(HDdata$female)
# Fit the model
modHD <- mmer2(sugar~1, 
               random=~male + and(female) + geno, 
               rcov=~units,
               data=HDdata, silent = TRUE)
summary(modHD)
suma <- summary(modHD)$var.comp.table
Vgca <- suma[1,1]
Vsca <- suma[2,1]
Ve <- suma[3,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve/2)) ) # 2 technical reps
(h2 <- Va / (Vg + (Ve/2)) )

## ------------------------------------------------------------------------
  data(HDdata)
  head(HDdata)
  #### GCA matrix for half diallel using male and female columns
  #### use the 'overlay' function to create the half diallel matrix
  Z1 <- overlay(HDdata[,c(3:4)])
  #### Obtain the SCA matrix
  Z2 <- model.matrix(~as.factor(geno)-1, data=HDdata)
  #### Define the response variable and run
  y <- HDdata$sugar
  ETA <- list(list(Z=Z1), list(Z=Z2)) # Zu component
  modHD <- mmer(Y=y, Z=ETA, draw=FALSE, silent=TRUE)
  summary(modHD)

## ------------------------------------------------------------------------
data(wheatLines); 
X <- wheatLines$wheatGeno; X[1:5,1:4]; dim(X)
Y <- data.frame(wheatLines$wheatPheno); Y$id <- rownames(Y); head(Y);
rownames(X) <- rownames(Y)
# select environment 1
K <- A.mat(X) # additive relationship matrix
# GBLUP pedigree-based approach
set.seed(12345)
y.trn <- Y
vv <- sample(rownames(Y),round(dim(Y)[1]/5))
y.trn[vv,"X1"] <- NA
ans <- mmer2(X1~1,
             random=~g(id), 
             rcov=~units,
             G=list(id=K), method="NR", 
             data=y.trn, silent = TRUE) # kinship based
cor(ans$u.hat$`g(id)`[vv,],Y[vv,"X1"])

## ------------------------------------------------------------------------
data(Technow_data)

A.flint <- Technow_data$AF # Additive relationship matrix Flint
A.dent <- Technow_data$AD # Additive relationship matrix Dent

pheno <- Technow_data$pheno # phenotypes for 1254 single cross hybrids
head(pheno);dim(pheno) 
# CREATE A DATA FRAME WITH ALL POSSIBLE HYBRIDS
DD <- kronecker(A.dent,A.flint,make.dimnames=TRUE)
hybs <- data.frame(sca=rownames(DD),yield=NA,matter=NA,gcad=NA, gcaf=NA)
hybs$yield[match(pheno$hy, hybs$sca)] <- pheno$GY
hybs$matter[match(pheno$hy, hybs$sca)] <- pheno$GM
hybs$gcad <- as.factor(gsub(":.*","",hybs$sca))
hybs$gcaf <- as.factor(gsub(".*:","",hybs$sca))
head(hybs)
# RUN THE PREDICTION MODEL
y.trn <- hybs
vv1 <- which(!is.na(hybs$yield))
vv2 <- sample(vv1, 100)
y.trn[vv2,"yield"] <- NA
anss2 <- mmer2(yield~1, 
               random=~g(gcad) + g(gcaf), 
               rcov=~units,
               G=list(gcad=A.dent, gcaf=A.flint), 
               method="NR", silent=TRUE, data=y.trn) 
summary(anss2)
cor(anss2$fitted.y[vv2], hybs$yield[vv2])

## ------------------------------------------------------------------------
data(CPdata)
### look at the data
head(CPpheno);CPgeno[1:5,1:4]
## fit a model including additive effects
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ans.A <- mmer2(cbind(color,Yield)~1, 
               random=~us(trait):g(id),
               rcov=~us(trait):units,
               G=list(id=A),
               data=CPpheno, silent = TRUE)
summary(ans.A)

## ------------------------------------------------------------------------
## genetic variance covariance
gvc <- ans.A$var.comp$`g(id)`
## extract variances (diagonals) and get standard deviations
sd.gvc <- as.matrix(sqrt(diag(gvc))) 
## get possible products sd(Vgi) * sd(Vgi')
prod.sd <- sd.gvc %*% t(sd.gvc)
## genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')]
(gen.cor <- gvc/prod.sd)
## heritabilities
(h2 <- diag(gvc) / diag(cov(CPpheno[,names(diag(gvc))], use = "complete.obs")))

