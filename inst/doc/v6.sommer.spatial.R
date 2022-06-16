## -----------------------------------------------------------------------------
library(sommer)
data(DT_yatesoats)
DT <- DT_yatesoats
DT$row <- as.numeric(as.character(DT$row))
DT$col <- as.numeric(as.character(DT$col))
DT$R <- as.factor(DT$row)
DT$C <- as.factor(DT$col)

# SPATS MODEL
# m1.SpATS <- SpATS(response = "Y",
#                   spatial = ~ PSANOVA(col, row, nseg = c(14,21), degree = 3, pord = 2),
#                   genotype = "V", fixed = ~ 1,
#                   random = ~ R + C, data = DT,
#                   control = list(tolerance = 1e-04))
# 
# summary(m1.SpATS, which = "variances")
# 
# Spatial analysis of trials with splines 
# 
# Response:                   Y         
# Genotypes (as fixed):       V         
# Spatial:                    ~PSANOVA(col, row, nseg = c(14, 21), degree = 3, pord = 2)
# Fixed:                      ~1        
# Random:                     ~R + C    
# 
# 
# Number of observations:        72
# Number of missing data:        0
# Effective dimension:           17.09
# Deviance:                      483.405
# 
# Variance components:
#                   Variance            SD     log10(lambda)
# R                 1.277e+02     1.130e+01           0.49450
# C                 2.673e-05     5.170e-03           7.17366
# f(col)            4.018e-15     6.339e-08          16.99668
# f(row)            2.291e-10     1.514e-05          12.24059
# f(col):row        1.025e-04     1.012e-02           6.59013
# col:f(row)        8.789e+01     9.375e+00           0.65674
# f(col):f(row)     8.036e-04     2.835e-02           5.69565
# 
# Residual          3.987e+02     1.997e+01 

# SOMMER MODEL
m1.sommer <- mmer(Y~1+V+spl2Db(col,row, nsegments = c(14,21), degree = c(3,3), penaltyord = c(2,2), what = "base"), 
                  random = ~R+C+spl2Db(col,row, nsegments = c(14,21), degree = c(3,3), penaltyord = c(2,2), what="bits"),
                  data=DT, tolParConv = 1e-6, verbose = FALSE)
summary(m1.sommer)$varcomp
# get the fitted values for the spatial kernel and plot
# ff <- fitted.mmer(m1.sommer)
# DT$fit <- as.matrix(Reduce("+",ff$Zu[-c(1:2)])) 
# lattice::levelplot(fit~row*col,data=DT)


## ---- fig.show='hold'---------------------------------------------------------

# SOMMER MODEL
m2.sommer <- mmer(Y~1+V, 
                  random = ~R+C+spl2Da(col,row, nsegments = c(14,21), degree = c(3,3), penaltyord = c(2,2)),
                  data=DT, tolParConv = 1e-6, verbose = FALSE)
summary(m1.sommer)$varcomp
# get the fitted values for the spatial kernel and plot
# ff <- fitted.mmer(m2.sommer)
# DT$fit <- as.matrix(Reduce("+",ff$Zu[-c(1:2)])) 
# lattice::levelplot(fit~row*col,data=DT)


## ---- fig.show='hold'---------------------------------------------------------

DT2 <- rbind(DT,DT)
DT2$Y <- DT2$Y + rnorm(length(DT2$Y))
DT2$trial <- c(rep("A",nrow(DT)),rep("B",nrow(DT)))
head(DT2)
# SOMMER MODEL
m3.sommer <- mmer(Y~1+V, 
                  random = ~vsr(dsr(trial),R)+vsr(dsr(trial),C)+
                    spl2Da(col,row, nsegments = c(14,21), degree = c(3,3), penaltyord = c(2,2), at.var = trial),
                  rcov = ~vsr(dsr(trial),units),
                  data=DT2, tolParConv = 1e-6, verbose = FALSE)
summary(m3.sommer)$varcomp
# get the fitted values for the spatial kernel and plot
# ff <- fitted.mmer(m3.sommer)
# DT2$fit <- as.matrix(Reduce("+",ff$Zu[-c(1:4)])) 
# lattice::levelplot(fit~row*col|trial,data=DT2)


