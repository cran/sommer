## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example

ansSingle <- mmer(Yield~1,
              random= ~ vs(Name, Gu=A),
              rcov= ~ units,
              data=DT)
summary(ansSingle)


## -----------------------------------------------------------------------------

ansMain <- mmer(Yield~Env,
              random= ~ vs(Name, Gu=A),
              rcov= ~ units,
              data=DT)
summary(ansMain)


## -----------------------------------------------------------------------------

ansDG <- mmer(Yield~Env,
              random= ~ vs(ds(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT)
summary(ansDG)


## -----------------------------------------------------------------------------
E <- diag(length(unique(DT$Env)))
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vs(Name, Gu=A) + vs(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT)
summary(ansCS)


## -----------------------------------------------------------------------------

ansMain <- mmer(Yield~Env,
              random= ~ vs(Name, Gu=A) + vs(ds(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT)
summary(ansMain)


## -----------------------------------------------------------------------------

ansUS <- mmer(Yield~Env,
              random= ~ vs(us(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT)
summary(ansUS)
# adjust variance BLUPs by adding covariances
# ansUS$U[1:6] <- unsBLUP(ansUS$U[1:6])

## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))
ansRR <- mmer(Yield~Env,
              random= ~ vs(leg(EnvN,1),Name),
              rcov= ~ units,
              data=DT)
summary(ansRR)


## -----------------------------------------------------------------------------

E <- CS(DT$Env)
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vs(Name, Gu=A) + vs(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT)
summary(ansCS)


