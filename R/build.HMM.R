build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE){
  # build hybrid marker matrix

  if(!is.null(custom.hyb)){
    pheno <- custom.hyb
    found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
    if(found != 3){
      stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
    }
    return.combos.only=FALSE
  }else{
    a <- rownames(M1)
    b <- rownames(M2)
    pheno <- expand.grid(a,b)
    pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
    pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=":")
  }
  
  if(!return.combos.only){
    # check that marker matrices are in -1,0,1 format
    checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
    checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))
    
    checkM1[which(checkM1 > 0)] <- 1
    checkM2[which(checkM2 > 0)] <- 1
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    
    
    ## add markers coming from parents M1
    Z1 <- model.matrix(~Var1-1,pheno);dim(Z1); 
    colnames(Z1) <- gsub("Var1","",colnames(Z1))
    M1 <- M1[colnames(Z1),]
    #M1[1:4,1:4]; Z1[1:4,1:4]; 
    ## add markers coming from parents M2
    Z2 <- model.matrix(~Var2-1,pheno);dim(Z2); 
    colnames(Z2) <- gsub("Var2","",colnames(Z2))
    M2 <- M2[colnames(Z2),]
    #M2[1:4,1:4]; Z2[1:4,1:4];  
    
    ## create the 
    # Z3 <- model.matrix(~hybrid-1,pheno);dim(Z3);
    # colnames(Z3) <- gsub("hybrid","",colnames(Z3))
    # hyb.names <- colnames(Z3)[as.vector(apply(Z3,1,function(x){which(x==1)}))] # names of hybrids
    hyb.names <- pheno$hybrid
    ## marker matrix for hybrids one for each parent
    cat(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))
    
    # M1 <- as(M1, Class="sparseMatrix")
    # M2 <- as(M2, Class="sparseMatrix")
    # Z1 <- as(Z1, Class="sparseMatrix")
    # Z2 <- as(Z2, Class="sparseMatrix")
    
    cat("Extracting M1 contribution\n")
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Md <- Z1 %*% M1;  # was already converted to -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      Md <- 2*Z1 %*% M1 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      Md <- Z1 %*% M1 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    cat("Extracting M2 contribution\n")
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Mf <- Z2 %*% M2;  # was already converted to -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      Mf <- 2*Z2 %*% M2 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      Mf <- Z2 %*% M2 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    ## marker matrix coded as additive -1,0,1
    Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
    rownames(Mdf) <- hyb.names
    #hist(Mdf)
    
    ## dominance matrix for hybrids (0,1 coded)
    Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
    rownames(Delta) <- hyb.names
    #hist(Delta)
    cat("Done!!\n")
    return(list(HMM.add=Mdf, HMM.dom=Delta, data.used=pheno))
    
  }else{
    return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
  }
}