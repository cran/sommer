## ------------------------------------------------------------------------
# iteration    LogLik     wall    cpu(sec)   restrained
#     1      -224.676   18:11:23      3           0
# Sistem is singular. Stopping the job
# matrix multiplication: incompatible matrix dimensions: 0x0 and ...x...

## ------------------------------------------------------------------------
library(sommer)
## rrBLUP for makers
data(DT_cpdata)
mix.rrblup <- mmer(fixed=cbind(color,Yield)~1,
                   random=~vs(GT,Gtc=unsm(2)) + vs(Rowf,Gtc=diag(2)),
                   rcov=~vs(units,Gtc=unsm(2)),
                   data=DT)
summary(mix.rrblup)
## GBLUP for individuals
A <- A.mat(GT)
mix.gblup <- mmer(fixed=cbind(color,Yield)~1,
                  random=~vs(id,Gu=A, Gtc=unsm(2)) + vs(Rowf,Gtc=diag(2)),
                  rcov=~vs(units,Gtc=unsm(2)),
                  data=DT)
summary(mix.gblup)

