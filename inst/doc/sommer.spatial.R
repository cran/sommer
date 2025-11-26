## -----------------------------------------------------------------------------
library(sommer)
data(DT_yatesoats, package="enhancer")
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
M <- spl2Dmats(x.coord.name = "col", y.coord.name = "row", data=DT, 
               nseg =c(14,21), degree = c(3,3), penaltyord = c(2,2) 
               )
mix <- mmes(Y~V, henderson = TRUE,
            random=~ R + C + vsm(ism(M$fC)) + vsm(ism(M$fR)) + 
              vsm(ism(M$fC.R)) + vsm(ism(M$C.fR)) +
              vsm(ism(M$fC.fR)),
            rcov=~units, verbose=FALSE,
            data=M$data)
summary(mix)$varcomp


## ----fig.show='hold'----------------------------------------------------------

# SOMMER MODEL
mix <- mmes(Y~V,
            random=~ R + C +
              spl2Dc(row,col),
            rcov=~units, verbose=FALSE,
            data=DT)
summary(mix)$varcomp


## ----fig.show='hold'----------------------------------------------------------

DT2 <- rbind(DT,DT)
DT2$Y <- DT2$Y + rnorm(length(DT2$Y))
DT2$trial <- c(rep("A",nrow(DT)),rep("B",nrow(DT)))
head(DT2)
# SOMMER MODEL
mix <- mmes(Y~V,
            random=~ R + C +
              spl2Dc(row,col, at.var = trial),
            rcov=~units, verbose=FALSE,
            data=DT2)
summary(mix)$varcomp


