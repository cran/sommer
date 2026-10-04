## -----------------------------------------------------------------------------
library(sommer)
data(DT_h2, package="enhancer")
DT <- DT_h2
DT <- DT[with(DT, order(Env)), ]
length(unique(DT$Env))
length(unique(DT$Name))
head(DT)

## -----------------------------------------------------------------------------
fitFA <- mmes(y ~ Env,
              random = ~ vsm(fam(Env, 2), ism(Name)),
              rcov = ~ units,
              nIters = 150, verbose = FALSE,
              data = DT)
summary(fitFA)$varcomp

## -----------------------------------------------------------------------------
fitRR <- mmes(y ~ Env,
              random = ~ vsm(rrm(Env, 2), ism(Name)),
              rcov = ~ units,
              nIters = 150, verbose = FALSE,
              data = DT)
summary(fitRR)$varcomp

## -----------------------------------------------------------------------------
c(AIC_FA = fitFA$AIC, AIC_RR = fitRR$AIC)
anova.mmes(fitFA, fitRR)

## -----------------------------------------------------------------------------
faInfo <- loadings_mmes(fitFA)
round(faInfo$loadings, 3)
round(faInfo$specific, 3)
faInfo$sigma2

## -----------------------------------------------------------------------------
rrInfo <- loadings_mmes(fitRR)
round(rrInfo$specific, 3)

## -----------------------------------------------------------------------------
faScores <- scores_mmes(fitFA)
head(faScores)

## -----------------------------------------------------------------------------
varPerFactor <- colSums(faInfo$loadings^2)
totalVar <- sum(faInfo$loadings^2) + sum(faInfo$specific)
propExplained <- varPerFactor / totalVar
round(100 * propExplained, 1)

## ----fig.show='hold'----------------------------------------------------------
barplot(t(faInfo$loadings), beside = TRUE,
        col = c("steelblue4", "tomato"),
        las = 2, cex.names = 0.7,
        ylab = "Loading",
        main = "Factor-analytic loadings by environment")
legend("topright", legend = colnames(faInfo$loadings),
       fill = c("steelblue4", "tomato"), bty = "n")

## ----fig.show='hold'----------------------------------------------------------
plot(faScores[,2] ~ faScores[,1],
     xlab = "Factor 1 score", ylab = "Factor 2 score",
     main = "Genotype scores")
text(faScores[,2] ~ faScores[,1], labels = rownames(faScores),
     cex = 0.6, pos = 1)
abline(h = 0, v = 0, lty = 3)

## ----fig.show='hold'----------------------------------------------------------
Sigma <- faInfo$sigma2 * (
  faInfo$loadings %*% t(faInfo$loadings) + diag(faInfo$specific)
)
corMat <- cov2cor(Sigma)
heatmap(corMat, symm = TRUE,
        main = "Fitted genetic correlation among environments")

