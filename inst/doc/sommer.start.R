## ------------------------------------------------------------------------
library(sommer)
data(example)
head(example)

ans1 <- mmer2(Yield~Env, 
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=example, silent = TRUE)
summary(ans1)


## ------------------------------------------------------------------------

data(example)
head(example)
ans1 <- mmer2(Yield~Env, 
              random= ~Name + at(Env):Name,
              rcov= ~ at(Env):units,
              data=example, silent = TRUE)
summary(ans1)


## ------------------------------------------------------------------------

data(example)
head(example)
ans1 <- mmer2(cbind(Yield, Weight) ~ Env, 
              random= ~ us(trait):Name + us(trait):Env:Name,
              rcov= ~ us(trait):units,
              data=example, silent = TRUE)
summary(ans1)


## ------------------------------------------------------------------------

data(example)
head(example)
ans1 <- mmer2(cbind(Yield, Weight) ~ Env, 
              random= ~ us(trait):Name + us(trait):at(Env):Name,
              rcov= ~ us(trait):at(Env):units,
              data=example, silent = TRUE)
summary(ans1)


## ------------------------------------------------------------------------

data(example)
head(example)
K[1:4,1:4]
ans1 <- mmer2(Yield ~ Env, 
              random= ~ g(Name) + at(Env):g(Name),
              rcov= ~ at(Env):units,
              G=list(Name=K),
              data=example, silent = TRUE)
summary(ans1)


## ------------------------------------------------------------------------

data(example)
head(example)
K[1:4,1:4]
ans1 <- mmer2(cbind(Yield, Weight) ~ Env, 
              random= ~ us(trait):g(Name) + us(trait):at(Env):g(Name),
              rcov= ~ us(trait):at(Env):units,
              G=list(Name=K),
              data=example, silent = TRUE)
summary(ans1)


