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
              random= ~ vs(Name) + vs(EnvName),
              rcov= ~ vs(units),
              data=DT)
summary(ans4)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
DT$EnvName <- paste(DT$Env,DT$Name)
ans5 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vs(Name) + vs(ds(Env),Name),
              rcov= ~ vs(ds(Env),units),
              data=DT)
summary(ans5)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
DT$EnvName <- paste(DT$Env,DT$Name)
ans6 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vs(us(Env),Name),
              rcov= ~ vs(ds(Env),units),
              data=DT)
summary(ans6)


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

data(DT_example)
head(DT)
ans2 <- mmer(Yield~Env,
              random= ~ vs(ds(Env),Name, Gu=A),
              rcov= ~ vs(ds(Env),units),
              data=DT)
summary(ans2)


## ------------------------------------------------------------------------

data(DT_example)
head(DT)
ans2 <- mmer(cbind(Yield,Weight)~Env,
              random= ~ vs(ds(Env),Name, Gu=A, Gtc=unsm(2)),
              rcov= ~ vs(ds(Env),units, Gtc=diag(2)),
              data=DT)
summary(ans2)


## ------------------------------------------------------------------------

data(DT_cpdata)
GT[1:4,1:4]
#### look at the data and fit the model
mix1 <- mmer(Yield~1,
              random=~vs(list(GT)),
              rcov=~units,
              data=DT)


## ------------------------------------------------------------------------

data("DT_halfdiallel")
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modh <- mmer(sugar~1, 
             random=~vs(overlay(DT$femalef,DT$malef)) 
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
unsm(4)

## ------------------------------------------------------------------------
uncm(4)

## ------------------------------------------------------------------------
fixm(4)

## ------------------------------------------------------------------------
fcm(c(1,0,1,0))

## ------------------------------------------------------------------------
data(DT_example)
ansf <- mmer(cbind(Yield,Weight)~vs(Env,Gtc=fcm(c(0,1))),
             random= ~ vs(ds(Env),Name),
             rcov= ~ vs(ds(Env),units),
             data=DT)
summary(ansf)

## ------------------------------------------------------------------------
data(DT_example)
ans.uns <- mmer(cbind(Yield,Weight)~Env,
             random= ~ vs(Name,Gtc=unsm(2)),
             rcov= ~ vs(units,Gtc=unsm(2)),
             data=DT)
summary(ans.uns)

ans.diag <- mmer(cbind(Yield,Weight)~Env,
             random= ~ vs(Name,Gtc=diag(2)),
             rcov= ~ vs(units,Gtc=diag(2)),
             data=DT)
summary(ans.diag)

## ------------------------------------------------------------------------
# Generate some fake data: 
# 100 males and 100 females
# Two traits are measured on each male, and two traits on each female
# 20 individuals per sex are measured for each of 5 different genotypes 
set.seed(3434)
df <- data.frame(
  sex = rep(c("female", "male"), each = 100),
  female_trait_1 = c(rnorm(100), rep(NA, 100)),
  female_trait_2 = c(rnorm(100), rep(NA, 100)),
  male_trait_1 = c(rep(NA, 100), rnorm(100)),
  male_trait_2 = c(rep(NA, 100), rnorm(100)),
  genotype = rep(rep(1:5, each = 20), 2),
  individual = 1:200
)
df$genotype <- as.factor(df$genotype)
df$individual <- as.factor(df$individual)

mm <- adiag1(unsm(2),unsm(2));mm
# mix <- mmer(cbind(female_trait_1, 
#                   female_trait_2,
#                   male_trait_1,
#                   male_trait_2) ~ 1,
#             random=~vs(genotype,Gtc=unsm(4)) + vs(individual,Gtc=mm),
#             rcov=~vs(units), na.method.Y = "include",
#             data=df)
# summary(mix)

