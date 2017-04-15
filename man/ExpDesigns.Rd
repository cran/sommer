\name{ExpDesigns}
\alias{ExpDesigns}
\docType{data}
\title{
  Data for different experimental designs
}
\description{
  The following data is a list containing data frames for different type of experimental designs relevant in plant breeding:
  
  1) Augmented designs (2 examples)
  
  2) Incomplete block designs (1 example)
  
  3) Split plot design (2 examples)
  
  4) Latin square designs (1 example)
  
  5) North Carolina designs I,II and III
  
  How to fit each is shown at the Examples section. This may help you get introduced to experimental designs relevant to plant breeding. Good luck.
}

\format{
 Different based on the design.
}

\source{
  Datasets and more detail about them can be found in the agricolae package. Here we just show the datasets and how to analyze them using the \code{\link{sommer}} package.
}

\references{
 Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}

\examples{

#### =================================== ####
#### ===== Augmented Block Design 1 ==== ####
#### =================================== ####
data(ExpDesigns)
data1 <- ExpDesigns$au1
head(data1)
## response variable: "yield"
## check indicator: "entryc" ('nc' for all unreplicated, but personal.name for checks)
## blocking factor: "block"
## treatments, personal names for replicated and non-replicated: "trt"
## check no check indicator: "new"
mix1 <- mmer2(yield~entryc, random=~block+trt, data=data1)
summary(mix1)

# ## compare raw unreplicated measure with adjusted blup
# library(plyr)
# avers <- ddply(data1, .(trt), summarize,mean = round(mean(yield), 2))
# plot(avers$mean[-c(1:4)], mix1$u.hat$trt[-c(1:4),1],ylab="BLUP",xlab="Raw", 
#      pch=20, cex=2, col="blue")
# ## if you have row and column information you can see the design
# library(agridat)
# data1$row <- c(sample(1:7),sample(1:6),sample(1:7))
# data1$col <- c(1:3)[data1$block]
# d1 <- desplot(yield ~ row*col, data1, main="ABD1",
#               text=trt,strip.cex=3, cex=1)
# print(d1)

#### =================================== ####
#### ===== Augmented Block Design 2 ==== ####
#### =================================== ####
data(ExpDesigns)
data2 <- ExpDesigns$au2
head(data2)
## response variable: "TSW"
## check indicator: "entryc"
## blocking factor: "Block"
## treatments, replicated and non-replicated: "Entry"
## check no check indicator: "new"
## this is also known as Federer's unreplicated design
mix2<- mmer2(TSW ~ entryc, random=~Block+Entry, data=data2)
summary(mix2)

#### =================================== ####
#### ===== Incomplete block design  ==== ####
#### =================================== ####
data(ExpDesigns)
data.ibd <- ExpDesigns$ibd$book
head(data.ibd)
ExpDesigns$ibd$sketch
## response variable: "yield"
## 2 replications (r)
## 30 genotypes (trt)
## 10 incomplete blocks (s) with 3 trts each (k)
## design was an alpha design
## agricolae::design.alpha(trt=paste("gen",1:30,sep=""),k=3,r=2,seed=5)$sketch
mix.ibd <- mmer2(yield~Genotype,random=~replication+replication:block,
                 data=data.ibd)
summary(mix.ibd)
# rownames(a)[1] <-"geno1"
#a[-1] <- a[-1]+a[1]
#plot(density(mix.ibd$beta.hat))

## map of the field
# library(agridat)
# data.ibd$block <- as.numeric(as.character(data.ibd$block))
# data.ibd$cols <- as.numeric(as.character(data.ibd$cols))
# d1 <- desplot(yield ~ block*cols, data.ibd, main="IBD",
#               text=Genotype,strip.cex=3, cex=1)
# print(d1)

#### =================================== ####
#### ======= Split Plot Design  ======== ####
#### =================================== ####
data(ExpDesigns)
data.spd <- ExpDesigns$spd
head(data.spd)
## response variable: "yield"
## 3 blocks or reps (r)
## 2 whole plot treatment (A)
## 3 small plot treatments (B)
##
##        i.e BLOCK 1
##[]======================[]
##[] A1(B1) A1(B2) A1(B3) []
##[] A2(B1) A2(B2) A2(B3) []
##[]======================[]
##
## more replication in whole plot treatments (A)
## less replication in sub plot treatments (B)
mix.split <- mmer2(yield ~block + A + B ,random=~ A:B, data=data.spd)
summary(mix.split)

#### =================================== ####
#### ==== Split-Split Plot Design  ===== ####
#### =================================== ####
data(ExpDesigns)
data.sspd <- ExpDesigns$sspd
head(data.sspd)
## response variable: "yield"
## 5 levels of nitrogen (N) main plot
## 3 levels of management (M) sub-plot
## 3 varieties (B) sub-sub-plot
##
##        i.e BLOCK 1
##[]==================================[]
##[] N1(M1(V1)) N1(M2(V1)) N1(M3(V1)) []
##[] N2(M1(V1)) N2(M2(V1)) N2(M3(V1)) []
##[] N3(M1(V1)) N3(M2(V1)) N3(M3(V1)) []
##[] N4(M1(V1)) N4(M2(V1)) N4(M3(V1)) []
##[] N5(M1(V1)) N5(M2(V1)) N5(M3(V1)) []
##[]==================================[]
##
head(data.sspd)
mix.sspd <- mmer2(yield ~1,random=~ block + nitrogen + management + 
                    variety + nitrogen:management + variety:nitrogen +
                    variety:management + variety:nitrogen:management,
                  data=data.sspd)
summary(mix.sspd)

#### =================================== ####
#### ======= Latin Square Design  ====== ####
#### =================================== ####
data(ExpDesigns)
data.lsd <- ExpDesigns$lsd
head(data.lsd)
## response variable: "yield"
## 4 columns (c)
## 4 rows (r)
## 4 varieties (V)
##
##   c1 c2 c3 c4
##[]=============[]
##[] V1 V4 V2 V3 [] row 1
##[] V2 V3 V4 V1 [] row 2
##[] V3 V2 V4 V1 [] row 3
##[] V4 V1 V3 V2 [] row 4
##[]=============[]
##   c1 c2 c3 c4
##
mix.lsd <- mmer2(yield ~ variety ,random=~ row + col, data=data.lsd)
summary(mix.lsd)

# library(agridat)
# desplot(yield ~ row*col, data.lsd, main="LSD",
#               strip.cex=3, cex=1, text=variety)

#### =================================== ####
#### =====  North Carolina Design I ==== ####
#### =================================== ####
data(ExpDesigns)
data.car1 <- ExpDesigns$car1
head(data.car1)
## response variable: "yield"
## male indicator: "male"
## female indicator: "female"
## replication: "rep"
## set of males: "set"
mix.car1 <- mmer2(yield~set,random=~ set:rep + set:male 
                  +set:male:female + set:male:female:rep, data=data.car1)
summary(mix.car1)
(Var.A <- 4*mix.car1$var.comp[2,1])
(Var.D <- 4*mix.car1$var.comp[3,1] - 4*mix.car1$var.comp[2,1])

#### =================================== ####
#### ===== North Carolina Design II ==== ####
#### =================================== ####
data(ExpDesigns)
data.car2 <- ExpDesigns$car2
head(data.car2)
## response variable: "yield"
## male indicator: "male"
## female indicator: "female"
## replication: "rep"
## set of males: "set"
mix.car2 <- mmer2(yield ~ 1, random=~ set + set:rep + set:male
                  + set:female + set:male:female, data=data.car2)
summary(mix.car2)
(Var.Am <- 4*mix.car2$var.comp[3,1])
(Var.Af <- 4*mix.car2$var.comp[4,1])
(Var.D <- 4*mix.car2$var.comp[5,1])

#### =================================== ####
#### ==== North Carolina Design III ==== ####
#### =================================== ####
data(ExpDesigns)
data.car3 <- ExpDesigns$car3
head(data.car3)
## response variable: "yield"
## male indicator: "male"
## female indicator: "female"
## replication: "rep"
## set of males: "set"
mix.car3 <- mmer2(yield ~ set + set:rep, random=~ set:male
                  + set:female + set:male:female, data=data.car3)
summary(mix.car3)
(Var.A <- 4*mix.car3$var.comp[1,1]) # var males
(Var.D <- 2*mix.car3$var.comp[3,1]) # var females in males


}

