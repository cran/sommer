## -----------------------------------------------------------------------------
# install.packages("lme4")
# library(lme4)
library(sommer)
data(DT_sleepstudy)
DT <- DT_sleepstudy

## lme4
# fm1 <- lmer(Reaction ~ Days + (1 | Subject), data=DT)
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject, 
            data=DT, tolparinv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))


## -----------------------------------------------------------------------------

## lme4
# fm1 <- lmer(Reaction ~ Days + (Days || Subject), data=DT)
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject + vs(Days, Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)

# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
summary(fm2)$varcomp


## -----------------------------------------------------------------------------

## lme4
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=DT)
## sommer
## no equivalence in sommer to find the correlation between the 2 vc
## this is the most similar which is equivalent to (intercept || slope)
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject + vs(Days, Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)

# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
summary(fm2)$varcomp


## -----------------------------------------------------------------------------

## lme4
# fm1 <- lmer(Reaction ~ Days + (0 + Days | Subject), data=DT)
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ vs(Days, Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)

# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
summary(fm2)$varcomp


## -----------------------------------------------------------------------------
library(orthopolynom)
## diagonal model
fm2 <- mmer(Reaction ~ Days,
            random= ~ vs(ds(Daysf), Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## unstructured model
fm2 <- mmer(Reaction ~ Days,
            random= ~ vs(us(Daysf), Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## random regression (legendre polynomials)
fm2 <- mmer(Reaction ~ Days,
            random= ~ vs(leg(Days,1), Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp
## unstructured random regression (legendre)
fm2 <- mmer(Reaction ~ Days,
            random= ~ vs(us(leg(Days,1)), Subject), 
            data=DT, tolparinv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp



