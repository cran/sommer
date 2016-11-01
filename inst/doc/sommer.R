## ------------------------------------------------------------------------
library(sommer)
data(h2)
head(h2)

ans1 <- mmer2(y~1, random=~Name + Env + Name:Env + Block,data=h2, silent = TRUE)
vc <- ans1$var.comp
V_E <- vc[2,1];V_GE <- vc[3,1];V_G <- vc[1,1];Ve <- vc[5,1]

n.env <- length(levels(h2$Env))
h2 <- V_G/(V_G + V_GE/n.env + Ve/(2*n.env)) #the 2 is a reference for block
h2

## ------------------------------------------------------------------------
data(CPdata)
CPpheno <- CPdata$pheno
CPgeno <- CPdata$geno
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
## fit a model including additive and dominance effects
y <- CPpheno$color
Za <- diag(length(y)); Zd <- diag(length(y)); Ze <- diag(length(y))
A <- A.mat(CPgeno) # additive relationship matrix
D <- D.mat(CPgeno) # dominance relationship matrix
E <- E.mat(CPgeno) # epistatic relationship matrix

ETA.ADE <- list(add=list(Z=Za,K=A),dom=list(Z=Zd,K=D),epi=list(Z=Ze,K=E))
ans.ADE <- mmer(Y=y, Z=ETA.ADE,silent = TRUE)
(H2 <- sum(ans.ADE$var.comp[1:3,1])/sum(ans.ADE$var.comp[,1]))
(h2 <- sum(ans.ADE$var.comp[1,1])/sum(ans.ADE$var.comp[,1]))

## ---- fig.show='hold'----------------------------------------------------
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
head(hybrid2)

modFD <- mmer2(Yield~Location, random=~GCA1+GCA2+SCA, data=hybrid2,silent = TRUE, draw=FALSE)
summary(modFD)
Vgca <- sum(modFD$var.comp[1:2,1])
Vsca <- modFD$var.comp[3,1]
Ve <- modFD$var.comp[4,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

## ---- fig.show='hold'----------------------------------------------------
data(HDdata)
head(HDdata)
# GCA matrix for half diallel using male and female columns
Z1 <- hdm(HDdata[,c(3:4)])
# SCA matrix
Z2 <- model.matrix(~as.factor(geno)-1, data=HDdata)
# Fit the model
y <- HDdata$sugar
ETA <- list(GCA=list(Z=Z1), SCA=list(Z=Z2)) # Zu component
modHD <- mmer(Y=y, Z=ETA,silent = TRUE)
summary(modHD)
Vgca <- modHD$var.comp[1,1]
Vsca <- modHD$var.comp[2,1]
Ve <- modHD$var.comp[3,1]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve/2)) ) # 2 technical reps
(h2 <- Va / (Vg + (Ve/2)) )
plot(density(randef(modHD)$GCA[,1]), col="cadetblue", 
     lwd=10, main="GCA BLUPs")

## ------------------------------------------------------------------------
data(CPdata)
CPpheno <- CPdata$pheno
CPgeno <- CPdata$geno
my.map <- CPdata$map
### look at the data
head(CPpheno); CPgeno[1:5,1:4]

y <- CPpheno$color # response
Za <- diag(length(y)) # inidence matrix for random effect
A <- A.mat(CPgeno) # additive relationship matrix
ETA.A <- list(add=list(Z=Za,K=A)) # create random component
ans.A <- mmer(Y=y, Z=ETA.A, W=CPgeno, silent=TRUE) # fit the model
### if you have a genetic map you can use it
ans.B <- mmer(Y=y, Z=ETA.A, W=CPgeno, silent=TRUE, map=my.map) # fit the model

## ------------------------------------------------------------------------
data(PolyData)
genotypes <- PolyData$PGeno
phenotypes <- PolyData$PPheno
## convert markers to numeric format
numo <- atcg1234(data=genotypes, ploidy=4, silent = TRUE); numo[1:5,1:4]; dim(numo)
# get only plants with both genotypes and phenotypes
common <- intersect(phenotypes$Name,rownames(numo))
marks <- numo[common,]; marks[1:5,1:4]
phenotypes2 <- phenotypes[match(common,phenotypes$Name),];
phenotypes2[1:5,1:4]
# Additive relationship matrix, specify ploidy
yy <- phenotypes2$tuber_shape
K1 <- A.mat(marks, ploidy=4)
Z1 <- diag(length(yy))
ETA <- list( list(Z=Z1, K=K1)) # random effects for genotypes
# run the model
models <- c("additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
ans2 <- mmer(Y=yy, Z=ETA, W=marks, method="EMMA", 
              ploidy=4, models=models[1], silent = TRUE)
summary(ans2)

## ------------------------------------------------------------------------
data(wheatLines)
X <- wheatLines$wheatGeno; X[1:5,1:4]; dim(X)
Y <- wheatLines$wheatPheno
rownames(X) <- rownames(Y)
# select environment 1
y <- Y[,1] # response grain yield
Z1 <- diag(length(y)) # incidence matrix
K <- A.mat(X) # additive relationship matrix
# GBLUP pedigree-based approach
set.seed(12345)
y.trn <- y
vv <- sample(1:length(y),round(length(y)/5))
y.trn[vv] <- NA
ETA <- list(g=list(Z=Z1, K=K))
ans <- mmer(Y=y.trn, Z=ETA, method="EMMA", silent = TRUE) # kinship based
cor(ans$u.hat$g[vv],y[vv])
## maximum prediction value that can be achieved
sqrt(ans$var.comp[1,1]/sum(ans$var.comp[,1]))

## ------------------------------------------------------------------------
data(Technow_data)

A.flint <- Technow_data$AF # Additive relationship matrix Flint
A.dent <- Technow_data$AD # Additive relationship matrix Dent
M.flint <- Technow_data$MF # Marker matrix Flint
M.dent <- Technow_data$MD # Marker matrix Dent

pheno <- Technow_data$pheno # phenotypes for 1254 single cross hybrids
pheno$hy <- paste(pheno$dent, pheno$flint, sep=":");head(pheno);dim(pheno) 
# CREATE A DATA FRAME WITH ALL POSSIBLE HYBRIDS
DD <- kronecker(A.dent,A.flint,make.dimnames=TRUE)

hybs <- data.frame(sca=rownames(DD),yield=NA,matter=NA,gcad=NA, gcaf=NA)
hybs$yield[match(pheno$hy, hybs$sca)] <- pheno$GY
hybs$matter[match(pheno$hy, hybs$sca)] <- pheno$GM
hybs$gcad <- as.factor(gsub(":.*","",hybs$sca))
hybs$gcaf <- as.factor(gsub(".*:","",hybs$sca))
head(hybs)
# CREATE INCIDENCE MATRICES
Z1 <- model.matrix(~gcad-1, data=hybs)
Z2 <- model.matrix(~gcaf-1, data=hybs)
# SORT INCIDENCE MATRICES ACCORDING TO RELATIONSHIP MATRICES, REAL ORDERS
real1 <- match(  colnames(A.dent), gsub("gcad","",colnames(Z1)))
real2 <- match(  colnames(A.flint), gsub("gcaf","",colnames(Z2)))
Z1 <- Z1[,real1]
Z2 <- Z2[,real2]
# RUN THE PREDICTION MODEL
y.trn <- hybs$yield
vv1 <- which(!is.na(hybs$yield))
vv2 <- sample(vv1, 100)
y.trn[vv2] <- NA
ETA2 <- list(GCA1=list(Z=Z1, K=A.dent), GCA2=list(Z=Z2, K=A.flint)) 
anss2 <- mmer(Y=y.trn, Z=ETA2, method="EM", silent=TRUE) 
summary(anss2)
cor(anss2$fitted.y[vv2], hybs$yield[vv2])

## ------------------------------------------------------------------------
data(CPdata)
CPpheno <- CPdata$pheno[,-c(1:4)]
CPgeno <- CPdata$geno
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
## fit a model including additive and dominance effects
Y <- CPpheno
Za <- diag(dim(Y)[1])
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(add=list(Z=Za,K=A))
ans.A <- mmer(Y=Y, Z=ETA.A, MVM=TRUE, method="EMMA")
summary(ans.A)

## ------------------------------------------------------------------------
## genetic variance covariance
gvc <- ans.A$var.comp$Vu
## extract variances (diagonals) and get standard deviations
sd.gvc <- as.matrix(sqrt(diag(gvc))) 
## get possible products sd(Vgi) * sd(Vgi')
prod.sd <- sd.gvc %*% t(sd.gvc)
## genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')]
(gen.cor <- gvc/prod.sd)
## heritabilities
(h2 <- diag(gvc) / diag(cov(Y, use = "complete.obs")))

## ------------------------------------------------------------------------
data(CPdata)
CPpheno <- CPdata$pheno[,-c(1:4)]
CPgeno <- CPdata$geno
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
## fit a model including additive and dominance effects
Y <- CPpheno
Za <- diag(dim(Y)[1])
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(add=list(Z=Za,K=A))
ans.A <- mmer(Y=Y, Z=ETA.A, W=CPgeno, MVM=TRUE, EIGEND=TRUE,
              map=CPdata$map, silent=TRUE,gwas.plots = FALSE)
summary(ans.A)

