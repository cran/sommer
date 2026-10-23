## ----eval=FALSE---------------------------------------------------------------
# fit <- mmes(
#   response ~ locale + cell,
#   data = observations,
#   control = mmes_control(representation = list(fixedRank = "drop"))
# )
# fit$fixedRank[c("rank", "retained", "omitted")]
# summary(fit)$betas

## ----eval=FALSE---------------------------------------------------------------
# cells <- predict(fit, D = "locale:cell", observed = TRUE)
# cells$observed

## ----eval=FALSE---------------------------------------------------------------
# rules <- fit$Dtable
# rules$include[rules$type == "random"] <- TRUE
# conditional <- predict(fit, D = "locale:cell", observed = TRUE, Dtable = rules)

