## -----------------------------------------------------------------------------
library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example

ansSingle <- mmer(Yield~1,
              random= ~ vsr(Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)

# or
Ai <- as(solve(A), Class="sparseMatrix")
ansSingle <- mmec(Yield~1,
              random= ~ vsc(isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansSingle)


## -----------------------------------------------------------------------------

ansMain <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)

# or 

Ai <- as(solve(A), Class="sparseMatrix")
ansMain <- mmec(Yield~Env,
              random= ~ vsc(isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansMain)


## -----------------------------------------------------------------------------

ansDG <- mmer(Yield~Env,
              random= ~ vsr(dsr(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)

# or
Ai <- as(solve(A), Class="sparseMatrix")
ansDG <- mmec(Yield~Env,
              random= ~ vsc(dsc(Env),isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansDG)


## -----------------------------------------------------------------------------
E <- diag(length(unique(DT$Env)))
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A) + vsr(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)

## or
E <- diag(length(unique(DT$Env)));rownames(E) <- colnames(E) <- unique(DT$Env)
Ei <- solve(E)
Ai <- solve(A)
EAi <- kronecker(Ei,Ai, make.dimnames = TRUE)
Ei <- as(Ei, Class="sparseMatrix")
Ai <- as(Ai, Class="sparseMatrix")
EAi <- as(EAi, Class="sparseMatrix")
ansCS <- mmec(Yield~Env,
              random= ~ vsc(isc(Name), Gu=Ai) + vsc(isc(Env:Name), Gu=EAi),
              rcov= ~ units, 
              data=DT, verbose = FALSE)
summary(ansCS)


## -----------------------------------------------------------------------------

ansUS <- mmer(Yield~Env,
              random= ~ vsr(usr(Env),Name, Gu=A),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)
# adjust variance BLUPs by adding covariances
# ansUS$U[1:6] <- unsBLUP(ansUS$U[1:6])

# or
Ai <- solve(A)
Ai <- as(Ai, Class="sparseMatrix")
ansUS <- mmec(Yield~Env,
              random= ~ vsc(usc(Env),isc(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansUS)


## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))
ansRR <- mmer(Yield~Env,
              random= ~ vsr(leg(EnvN,1),Name),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)

# or

ansRR <- mmec(Yield~Env,
              random= ~ vsc(dsc(leg(EnvN,1)),isc(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))
ansRR <- mmer(Yield~Env,
              random= ~ vsr(usr(leg(EnvN,1)),Name),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)

# or

ansRR <- mmec(Yield~Env,
              random= ~ vsc(usc(leg(EnvN,1)),isc(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansRR)


## -----------------------------------------------------------------------------

E <- AR1(DT$Env) # can be AR1() or CS(), etc.
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vsr(Name, Gu=A) + vsr(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)


