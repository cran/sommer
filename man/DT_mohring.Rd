\name{DT_mohring}
\alias{DT_mohring}
\docType{data}
\title{
Full diallel data for corn hybrids
}
\description{
This dataset contains phenotpic data for 36 winter bean hybrids, coming from a full diallel design and evaluated for 9 traits. The column male and female origin columns are included as well.

}
\usage{data("DT_mohring")}
\format{
  The format is:
 chr "DT_mohring"
}
\source{
This data was generated by a winter bean study and originally included in the agridat package.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}
\examples{

# ####=========================================####
# #### For CRAN time limitations most lines in the
# #### examples are silenced with one '#' mark,
# #### remove them and run the examples
# ####=========================================####
# data(DT_mohring)
# DT <- DT_mohring
# head(DT)
# DT2 <- add.diallel.vars(DT,par1="Par1", par2="Par2")
# head(DT2)
# # is.cross denotes a hybrid (1)
# # is.self denotes an inbred (1)
# # cross.type denotes one way (-1, e.g. AxB) and reciprocal (1, e.g., BxA) and no cross (0)
# # cross.id denotes the name of the cross (same name for direct & reciprocal)
#
# # GRIFFING MODEL 2 with reciprocal effects ###########################
#
# mod1hb <- mmes(Ftime ~ 1, data=DT2,
#                random = ~ Block
#                # GCA male & female overlayed
#                + vsm(ism(overlay(Par1, Par2)))
#                # SCA effects (includes cross and selfs)
#                + cross.id
#                # SCAR reciprocal effects (remaining variance in crosses; 
#                # if zero there's no reciprocal effects)
#                + vsm(dsm(cross.type), ism(cross.id)) )
# summary(mod1hb)$varcomp
#
# ##                                    VarComp VarCompSE   Zratio
# ## Block.Ftime-Ftime                  0.00000   9.32181 0.000000
# ## overlay(Par1, Par2).Ftime-Ftime 1276.73089 750.17269 1.701916
# ## cross.id.Ftime-Ftime            1110.99090 330.16921 3.364914
# ## cross.id:cross.type.Ftime-Ftime   66.02295  49.26876 1.340057
# ## units.Ftime-Ftime                418.47949  74.56442 5.612321
# ##
# # GRIFFING MODEL 2, no reciprocal effects ##############################
#
# mod1h <- mmes(Ftime ~ Block + is.cross, data=DT2, nIters = 50,
#               random = ~
#                 # GCA effects for all (hybrids and inbreds)
#                 vsm(ism(overlay(Par1, Par2)))
#               # GCA effect (calculated only in hybrids; remaining variance)
#               + vsm(ism(is.cross),ism(overlay(Par1, Par2)))
#               # SCA effect (calculated in hybrids only)
#               + vsm(ism(is.cross), ism(cross.id))
#               )
# summary(mod1h)$varcomp
#
# ##                                           VarComp  VarCompSE   Zratio
# ## overlay(Par1, Par2).Ftime-Ftime          2304.1781 1261.63193 1.826347
# ## overlay(Par1, Par2):is.cross.Ftime-Ftime  613.6040  402.74347 1.523560
# ## cross.id:is.cross.Ftime-Ftime             340.7030  148.56225 2.293335
# ## units.Ftime-Ftime                         501.6275   74.36075 6.745864
# ##
# # GRIFFING MODEL 3, no reciprocal effects ###############################
#
# mod1h <- mmes(Ftime ~ Block + is.cross, data=DT2, nIters = 100,
#               random = ~
#                 # GCAC (only for hybrids)
#                 vsm(ism(is.cross),ism(overlay(Par1, Par2)))
#               # male GCA (only for inbreds)
#               + vsm(ism(is.self),ism(Par1))
#               # SCA (for hybrids only)
#               + vsm(ism(is.cross), ism(cross.id))
#               )
# summary(mod1h)$varcomp
# ##                                           VarComp  VarCompSE   Zratio
# ## overlay(Par1, Par2):is.cross.Ftime-Ftime  927.7895  537.91218 1.724797
# ## Par1:is.self.Ftime-Ftime                 9960.9247 5456.58188 1.825488
# ## cross.id:is.cross.Ftime-Ftime             341.4567  148.53667 2.298804
# ## units.Ftime-Ftime                         498.5974   73.92066 6.745035
# ##
# # GRIFFING MODEL 2, with reciprocal effects #############################
# # In Mohring: mixed model 3 reduced
#
# mod1h <- mmes(Ftime ~ Block + is.cross, data=DT2, nIters = 100,
#               random = ~
#                 # GCAC (for hybrids only)
#                 vsm(ism(is.cross),ism(overlay(Par1, Par2)))
#               # male GCA (for selfs only)
#               + vsm(ism(is.self),ism(Par1))
#               # SCA (for hybrids only)
#               + vsm(ism(is.cross), ism(cross.id))
#               # SCAR reciprocal effects (remaning SCA variance)
#               + vsm(ism(cross.type), ism(cross.id))
#               )
# summary(mod1h)$varcomp
#
# ##                                             VarComp  VarCompSE   Zratio
# ## overlay(Par1, Par2):is.cross.Ftime-Ftime   927.78742  537.89981 1.724833
# ## Par1:is.self.Ftime-Ftime                 10001.78854 5456.47578 1.833013
# ## cross.id:is.cross.Ftime-Ftime              361.89712  148.54264 2.436318
# ## cross.id:cross.type.Ftime-Ftime             66.43695   49.24492 1.349113
# ## units.Ftime-Ftime                          416.82960   74.27202 5.612203
# ##
# # GRIFFING MODEL 3, with RGCA + RSCA ####################################
# # In Mohring: mixed model 3
#
# mod1h <- mmes(Ftime ~ Block + is.cross, data=DT2,nIters = 100,
#               random = ~
#                 # GCAC (for hybrids only)
#                 vsm(ism(is.cross),ism(overlay(Par1, Par2)))
#               # RGCA: exclude selfs (to identify reciprocal GCA effects)
#               + vsm(ism(cross.type),ism(overlay(Par1, Par2)))
#               # male GCA (for selfs only)
#               + vsm(ism(is.self),ism(Par1))
#               # SCA (for hybrids only)
#               + vsm(ism(is.cross), ism(cross.id))
#               # SCAR: exclude selfs (if zero there's no reciprocal SCA effects)
#               + vsm(ism(cross.type), ism(cross.id))
#               )
# summary(mod1h)$varcomp
#
# ##                                            VarComp  VarCompSE    Zratio
# ## overlay(Par1, Par2):is.cross.Ftime-Ftime   927.7843  537.88164 1.7248857
# ## Par1:is.self.Ftime-Ftime                 10001.7570 5456.30125 1.8330654
# ## cross.id:is.cross.Ftime-Ftime              361.8958  148.53670 2.4364068
# ## overlay(Par1, Par2):cross.type.Ftime-Ftime  17.9799   19.92428 0.9024114
# ## cross.id:cross.type.Ftime-Ftime             30.9519   46.43908 0.6665054
# ## units.Ftime-Ftime                         416.09922  447.2101 0.93043333



}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmes}} }
