phase.F1 <- function(genos, silent=FALSE, start=7){
  ########### convert to rqtl function
  ConvertToRqtl<-function(sup1,p=6, silent=FALSE){
    markerT<-c('l','m','n','p','h','k','e','f','g','a','b','c','d')
    markerT2<-c(1,1,2,2,3,3,4,4,4,5,5,5,5)
    matrixK<-matrix(NA,nrow=nrow(sup1),ncol=ncol(sup1)-p)
    ######## progress bar
    if(!silent){
      count <- 0
      tot <- nrow(sup1)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
    
    for (i in 1:nrow(sup1)){
      
      ### fill bar
      if(!silent){
        count <- count + 1
      }
      ####
      
      for (ii in 1:(ncol(sup1)-p)){
        geno1<-sup1[i,ii+p]
        if (geno1!='--'){
          markerType<-markerT2[which(substr(sup1$Segregation[i],start=2,stop=2)==markerT)]
          ## lmxll
          if (markerType==1){
            markerPhase<-numeric()
            markerPhase<-as.numeric(substr(sup1$Phase[i],start=2,stop=2))
            if (markerPhase==0){
              if (geno1=='lm'){
                matrixK[i,ii]<-6
              }else{
                matrixK[i,ii]<-5
              }
            }else{
              if (geno1=='lm'){
                matrixK[i,ii]<-5
              }else{
                matrixK[i,ii]<-6
              }
            }
          }
          ## nnxnp
          if (markerType==2){
            markerPhase<-numeric()
            markerPhase<-as.numeric(substr(sup1$Phase[i],start=3,stop=3))
            if (markerPhase==0){
              if (geno1=='np'){
                matrixK[i,ii]<-8
              }else{
                matrixK[i,ii]<-7
              }
            }else{
              if (geno1=='np'){
                matrixK[i,ii]<-7
              }else{
                matrixK[i,ii]<-8
              }
            }
          }
          ##hkxhk
          if (markerType==3){
            markerPhase<-numeric()
            for (j in 1:2){
              markerPhase[j]<-as.numeric(substr(sup1$Phase[i],start=j+1,stop=j+1))
            }
            
            if (markerPhase[1]==0 & markerPhase[2]==0){
              if (geno1=='hh'){
                matrixK[i,ii]<-1
              }
              if (geno1=='hk'){
                matrixK[i,ii]<-10
              }
              if (geno1=='kk'){
                matrixK[i,ii]<-4
              }
            }
            if (markerPhase[1]==0 & markerPhase[2]==1){
              if (geno1=='hh'){
                matrixK[i,ii]<-3
              }
              if (geno1=='hk'){
                matrixK[i,ii]<-9
              }
              if (geno1=='kk'){
                matrixK[i,ii]<-2
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==0){
              if (geno1=='hh'){
                matrixK[i,ii]<-2
              }
              if (geno1=='hk'){
                matrixK[i,ii]<-9
              }
              if (geno1=='kk'){
                matrixK[i,ii]<-3
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==1){
              if (geno1=='hh'){
                matrixK[i,ii]<-4
              }
              if (geno1=='hk'){
                matrixK[i,ii]<-10
              }
              if (geno1=='kk'){
                matrixK[i,ii]<-1
              }
            }
            
          }  # end of phase hkhk
          ##efxeg
          if (markerType==4){
            markerPhase<-numeric()
            for (j in 1:2){
              markerPhase[j]<-as.numeric(substr(sup1$Phase[i],start=j+1,stop=j+1))
            }
            if (markerPhase[1]==0 & markerPhase[2]==0){
              if (geno1=='ee'){
                matrixK[i,ii]<-1
              }
              if (geno1=='ef'){
                matrixK[i,ii]<-2
              }
              if (geno1=='eg'){
                matrixK[i,ii]<-3
              }
              if (geno1=='fg'){
                matrixK[i,ii]<-4
              }
            }
            if (markerPhase[1]==0 & markerPhase[2]==1){
              if (geno1=='ee'){
                matrixK[i,ii]<-3
              }
              if (geno1=='ef'){
                matrixK[i,ii]<-4
              }
              if (geno1=='eg'){
                matrixK[i,ii]<-1
              }
              if (geno1=='fg'){
                matrixK[i,ii]<-2
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==0){
              if (geno1=='ee'){
                matrixK[i,ii]<-2
              }
              if (geno1=='ef'){
                matrixK[i,ii]<-1
              }
              if (geno1=='eg'){
                matrixK[i,ii]<-4
              }
              if (geno1=='fg'){
                matrixK[i,ii]<-3
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==1){
              if (geno1=='ee'){
                matrixK[i,ii]<-4
              }
              if (geno1=='ef'){
                matrixK[i,ii]<-3
              }
              if (geno1=='eg'){
                matrixK[i,ii]<-2
              }
              if (geno1=='fg'){
                matrixK[i,ii]<-1
              }
            }
            
          }  # end of phase efeg
          ##abxcd
          if (markerType==5){
            markerPhase<-numeric()
            for (j in 1:2){
              markerPhase[j]<-as.numeric(substr(sup1$Phase[i],start=j+1,stop=j+1))
            }
            if (markerPhase[1]==0 & markerPhase[2]==0){
              if (geno1=='ac'){
                matrixK[i,ii]<-1
              }
              if (geno1=='bc'){
                matrixK[i,ii]<-2
              }
              if (geno1=='ad'){
                matrixK[i,ii]<-3
              }
              if (geno1=='bd'){
                matrixK[i,ii]<-4
              }
            }
            if (markerPhase[1]==0 & markerPhase[2]==1){
              if (geno1=='ac'){
                matrixK[i,ii]<-3
              }
              if (geno1=='ad'){
                matrixK[i,ii]<-1
              }
              if (geno1=='bc'){
                matrixK[i,ii]<-4
              }
              if (geno1=='bd'){
                matrixK[i,ii]<-2
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==0){
              if (geno1=='ac'){
                matrixK[i,ii]<-2
              }
              if (geno1=='ad'){
                matrixK[i,ii]<-4
              }
              if (geno1=='bc'){
                matrixK[i,ii]<-1
              }
              if (geno1=='bd'){
                matrixK[i,ii]<-3
              }
            }
            if (markerPhase[1]==1 & markerPhase[2]==1){
              if (geno1=='ac'){
                matrixK[i,ii]<-4
              }
              if (geno1=='ad'){
                matrixK[i,ii]<-2
              }
              if (geno1=='bc'){
                matrixK[i,ii]<-3
              }
              if (geno1=='bd'){
                matrixK[i,ii]<-1
              }
            }
          }  # end of phase ebcd
        }
      }
      
      ###########
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      ############
    }
    
    matrixK2<-as.data.frame(cbind(sup1[,1:p],matrixK))
    colnames(matrixK2)<-colnames(sup1)
    return(matrixK2)
  }
  ####################################
  
  
  
  ###################################
  ## impute hks function
  imputeHKs<-function(SNPdata, start=5, cM=3, silent=FALSE){
    chr <- max(SNPdata$LG)
    ######################## pick function
    PickHKtype<-function(x,y,z){
      #genotype number i nr/QTL
      geno<-numeric()
      #Matrix of possible genotype combinations#
      pickMatrix<-matrix(NA,8,8)
      pickMatrix[1,]<-c(1,NA,NA,NA,1,NA,1,NA)
      pickMatrix[2,]<-c(NA,2,NA,NA,NA,2,2,NA)
      pickMatrix[3,]<-c(NA,NA,3,NA,3,NA,NA,3)
      pickMatrix[4,]<-c(NA,NA,NA,4,NA,4,NA,4)
      pickMatrix[5,]<-c(1,NA,3,NA,11,NA,1,3)
      pickMatrix[6,]<-c(NA,2,NA,4,NA,12,2,4)
      pickMatrix[7,]<-c(1,2,NA,NA,1,3,13,NA)
      pickMatrix[8,]<-c(NA,NA,3,4,3,4,NA,14)
      rownames(pickMatrix)<-c("ac","bc","ad","bd","a-","b-","-c","-d")
      colnames(pickMatrix)<-rownames(pickMatrix)
      #genotype 9 is ac(1) or bd(4)
      #genotype 10 is ad(3) or bc(2)
      if(is.na(x)| is.na(y)|x==9|x==10|y==9|y==10){geno<-z
      }else{
        if(is.na(pickMatrix[x,y])){
          geno<-z
        } 
        else if(pickMatrix[x,y]==11){
          if(z==9){
            geno<-1
          } else {geno<-3}
        }
        else if(pickMatrix[x,y]==12){
          if(z==9){
            geno<-4
          } else {geno<-2}
        }
        else if(pickMatrix[x,y]==13){
          if(z==9){
            geno<-1
          } else {geno<-2}
        }
        else if(pickMatrix[x,y]==14){
          if(z==9){
            geno<-4
          } else {geno<-3}
        }
        else{
          geno<-pickMatrix[x,y]
        }
      }
      return(geno)
    }
    ########################
    imputedHKs<-SNPdata
    imputedHKs[nrow(SNPdata),start:ncol(SNPdata)]<-NA
    for(g in 1:chr){
      playMatrix<-SNPdata[which(SNPdata$LG==g),]
      tracker<-playMatrix
      tracker[1:nrow(playMatrix),start:ncol(playMatrix)]<-NA
      #### progress bar initiation
      
      if(!silent){
        count <- 0
        tot <- nrow(playMatrix)
        pb <- txtProgressBar(style = 3)
        setTxtProgressBar(pb, 0)
      }
      ########
      for(i in 1:nrow(playMatrix)) {
        ### fill bar
        if(!silent){
          count <- count + 1
        }
        ####
        if (playMatrix$Segregation[i] == "<hkxhk>") {
          for (j in start:ncol(playMatrix)) {
            if (is.na(playMatrix[i,j])){
              
            } else if (playMatrix[i,j] == 9 | playMatrix[i,j] == 10 ) {
              counttop <- i
              countbottom <- i
              while (counttop >= 1) {
                counttop<-counttop-1
                if(is.na(playMatrix[counttop+1,j])) next;
                counttop<-counttop+1;
                if (playMatrix[counttop,j] != 9 & playMatrix[counttop,j] != 10)
                  break;
                counttop <- counttop - 1
                
              }
              while (countbottom < nrow(playMatrix)) {
                countbottom<-countbottom+1
                if(is.na(playMatrix[countbottom-1,j])) next;
                countbottom<-countbottom-1;
                if (playMatrix[countbottom,j] != 9 & playMatrix[countbottom,j] != 10)
                  break;
                countbottom <- countbottom + 1
                #print(c("row:",i))
                #print(c("column:",j))
                #print(c("linkage group:",g))
                #print(c("top genotype:",SNPdata[counttop,j]))
                #print(c("bottom genotype:",SNPdata[countbottom,j]))
              }
              if (counttop == 0) {
                counttop <- 1
              }
              if(abs(playMatrix$Position[countbottom]-playMatrix$Position[counttop])<=cM){
                playMatrix[i,j]<-PickHKtype(playMatrix[counttop,j],playMatrix[countbottom,j],playMatrix[i,j])
                tracker[i,j]<-playMatrix[i,j]
              } else {tracker[i,j] <- playMatrix[i,j]}
            }
          }
        }
        ###########
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        ############
      }
      SNPdata[which(SNPdata$LG==g),]<-playMatrix
      imputedHKs[which(imputedHKs$LG==g),]<-tracker
      
      cat(paste("\nChromosome",g,"imputed\n"))
    }
    totalSNPs<-nrow(imputedHKs)*(ncol(imputedHKs[,start:ncol(imputedHKs)]))
    total9sand10s<-totalSNPs-sum(is.na(imputedHKs))
    overallTotal<-total9sand10s-length(which(imputedHKs[,start:ncol(imputedHKs)]==9))-length(which(imputedHKs[,start:ncol(imputedHKs)]==10))
    send<-c(totalSNPs,total9sand10s,overallTotal)
    names(send)<-c("TotalGenotype Data","Total hks","Total hks imputed")
    result<-list(imputed=SNPdata,checkTracker=imputedHKs,explanation=send)
    return(result)
  }
  ########################################
  
  ###Converts from r/qtl to JoinMap Parental Form
  ConvertToJoinMapParentals<-function(rQTLData, silent=FALSE){
    Maternal<-rQTLData
    Paternal<-rQTLData
    
    ######## progress bar
    if(!silent){
      count <- 0
      tot <- nrow(rQTLData)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
    
    for (i in 1:nrow(rQTLData)){
      ### fill bar
      if(!silent){
        count <- count + 1
      }
      ####
      for (j in 7:ncol(rQTLData)){
        x<-rQTLData[i,j]
        if (is.na(x)==F){
          if (x==1){
            Maternal[i,j]<-'ll'
            Paternal[i,j]<-'nn'
          }
          if (x==2){
            Maternal[i,j]<-'lm'
            Paternal[i,j]<-'nn'
          }
          if (x==3){
            Maternal[i,j]<-'ll'
            Paternal[i,j]<-'np'
          }
          if (x==4){
            Maternal[i,j]<-'lm'
            Paternal[i,j]<-'np'
          }
          if (x==5){
            Maternal[i,j]<-'ll'
          }
          if (x==6){
            Maternal[i,j]<-'lm'
          }
          if (x==7){
            Paternal[i,j]<-'nn'
          }
          if (x==8){
            Paternal[i,j]<-'np'
          }
          if (x==9){
            Maternal[i,j]<-'--'
            Paternal[i,j]<-'--'
          }
          if (x==10){
            Maternal[i,j]<-'--'
            Paternal[i,j]<-'--'
          }
          
        }else{
          Maternal[i,j]<-'--'
          Paternal[i,j]<-'--'
        }
        
      }
      if (rQTLData$Segregation[i]=='<efxeg>'|rQTLData$Segregation[i]=='<abxcd>'|rQTLData$Segregation[i]=='<hkxhk>'){
        Maternal$Segregation[i]<-'<lmxll>'
        Paternal$Segregation[i]<-'<nnxnp>'
      }
      Maternal$Phase[i]<-"{0-}"
      Paternal$Phase[i]<-"{-0}"
      
      ###########
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      ############
    }
    Paternal<-Paternal[-which(Paternal$Segregation!="<nnxnp>"),]
    Maternal<-Maternal[-which(Maternal$Segregation!="<lmxll>"),]
    return(list(maternal=Maternal,paternal=Paternal))
  }
  
  
  ########################################
  ###Converts from r/qtl to JoinMap Cp form
  ConvertToJoinMapCP<-function(matrixKK){
    if(is.data.frame(matrixKK)){
      matrixKK$Segregation<-as.character(matrixKK$Segregation)
      
      ######## progress bar
      if(!silent){
        count <- 0
        tot <- nrow(matrixKK)
        pb <- txtProgressBar(style = 3)
        setTxtProgressBar(pb, 0)
      }
      #####################
      
      for (i in 1:nrow(matrixKK)){
        
        ### fill bar
        if(!silent){
          count <- count + 1
        }
        ####
        
        for (j in 1:ncol(matrixKK)){
          x<-matrixKK[i,j]
          if (is.na(x)==F){
            if (x==1){
              matrixKK[i,j]<-'ac'
            }
            if (x==2){
              matrixKK[i,j]<-'bc'
            }
            if (x==3){
              matrixKK[i,j]<-'ad'
            }
            if (x==4){
              matrixKK[i,j]<-'bd'
            }
            if (x==5){
              matrixKK[i,j]<-'ll'
            }
            if (x==6){
              matrixKK[i,j]<-'lm'
            }
            if (x==7){
              matrixKK[i,j]<-'nn'
            }
            if (x==8){
              matrixKK[i,j]<-'np'
            }
            if (x==9){
              matrixKK[i,j]<-'--'
            }
            if (x==10){
              matrixKK[i,j]<-'--'
            }
            
          }else{matrixKK[i,j]<-'--'}
          
        }
        if (matrixKK$Segregation[i]=='<efxeg>'){
          matrixKK$Segregation[i]<-'<abxcd>'
        }
        if (matrixKK$Segregation[i]=='<hkxhk>'){
          matrixKK$Segregation[i]<-'<abxcd>'
        }
        ###########
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        ############
      }
      return(matrixKK)
    }else{}
  }
  ########################################
  
  ##############################################
  ## CONVERT TO RQTL FORMAT
  ## sup1 argument is the dataframe
  ## p argument is the column where marker data start
  cat("\nConverting to Rqtl format\n")
  dd <- ConvertToRqtl(sup1 = genos,p=6)
  ##############################################
  ##############################################
  ## USE CLOSE MARKERS TO IMPUTE HK MARKERS TO AN
  ## SPECIFIC RQTL VALUE
  cat("\nImputing Aa x Aa markers according to linkage information\n")
  ee <- imputeHKs(SNPdata = dd,start = start)$imputed
  ##############################################
  ##############################################
  ## OBTAIN PARENTAL MAPS
  cat("\nProducing parental maps phased\n")
  ff <- ConvertToJoinMapParentals(rQTLData = ee)
  
  
  ## maternal
  upp <- t(ff$maternal[,c("LG","Position")])
  low <- t(ff$maternal[,c(7:(dim(ff$maternal)[2]))])
  mato <- rbind(upp,low)
  mato2 <- mato
  mato2[which(mato == "lm" | mato=="np", arr.ind = TRUE)] <- "B"
  mato2[which(mato == "ll" | mato=="nn", arr.ind = TRUE)] <- "A"
  
  ## paternal
  upp <- t(ff$paternal[,c("LG","Position")])
  low <- t(ff$paternal[,c(7:(dim(ff$paternal)[2]))])
  pato <- rbind(upp,low)
  pato2 <- pato
  pato2[which(pato == "lm" | pato=="np", arr.ind = TRUE)] <- "B"
  pato2[which(pato == "ll" | pato=="nn", arr.ind = TRUE)] <- "A"
  
  ## paternal
  upp <- t(ee[,c("LG","Position")])
  low <- t(ee[,c(7:(dim(ee)[2]))])
  cato <- rbind(upp,low)
  
  ##############################################
  ## go back to joinmap format but with abxcd formats
  cat("\nProducing joinmap file phased and with abxcd markers from hk's\n")
  hh <- ConvertToJoinMapCP(ee)
  hh$LG <- ee$LG
  
  ff$consen <- NA
  ff$consen <- hh
  
  ff$consenRQTL <- NA
  ff$consenRQTL <- cato
  
  ff$maternalRQTL <- NA
  ff$maternalRQTL <- mato2
  
  ff$paternalRQTL <- NA
  ff$paternalRQTL <- pato2
  
  return(ff)
}