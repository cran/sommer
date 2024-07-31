## -----------------------------------------------------------------------------
# install.packages("lme4")
# library(lme4)
library(sommer)
data(DT_sleepstudy)
DT <- DT_sleepstudy
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (1 | Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Subject  (Intercept) 1378.2   37.12   
#  Residual              960.5   30.99   
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject, 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

# fm2 <- mmec(Reaction ~ Days,
#             random= ~ Subject, 
#             data=DT, tolParInv = 1e-6, verbose = FALSE)
# summary(fm2)$varcomp



## -----------------------------------------------------------------------------
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (Days || Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  Subject   (Intercept) 627.57   25.051  
#  Subject.1 Days         35.86    5.988  
#  Residual              653.58   25.565  
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject + vsr(Days, Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

# fm2 <- mmec(Reaction ~ Days,
#             random= ~ Subject + vsc(dsc(Days), isc(Subject)),
#             data=DT, tolParInv = 1e-6, verbose = FALSE)
# summary(fm2)$varcomp


## -----------------------------------------------------------------------------
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=DT)
# summary(fm1) # or # vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name        Variance Std.Dev. Corr
#  Subject  (Intercept) 612.10   24.741       
#           Days         35.07    5.922   0.07
#  Residual             654.94   25.592       
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
## using mmec() the model would be
fm2 <- mmec(Reaction ~ Days,
            random= ~ covc( vsc(isc(Subject)) , vsc(isc(Days), isc(Subject)) ), 
            nIters = 200, data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
cov2cor(fm2$theta[[1]])


## -----------------------------------------------------------------------------
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (0 + Days | Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name Variance Std.Dev.
#  Subject  Days  52.71    7.26   
#  Residual      842.03   29.02   
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(Days, Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

# fm2 <- mmec(Reaction ~ Days,
#             random= ~ vsc(dsc(Days), isc(Subject)), 
#             data=DT, tolParInv = 1e-6, verbose = FALSE)
# summary(fm2)$varcomp


## -----------------------------------------------------------------------------
library(orthopolynom)
## diagonal model
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(dsr(Daysf), Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## unstructured model
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(usr(Daysf), Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## random regression (legendre polynomials)
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(leg(Days,1), Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## unstructured random regression (legendre)
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(usr(leg(Days,1)), Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

# same can be done with the mmec function

