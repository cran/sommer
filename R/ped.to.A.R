# 
# 
# getA <- function(ped, sire="sire", dam="dam", sib="sib"){
#   
#   Dmat<-function (ped){
#     FF <- inbreeding(ped)
#     sire <- ped[,sire]
#     dam <- ped[,dam]
#     Fsire <- ifelse(is.na(sire), -1, FF[sire])
#     Fdam <- ifelse(is.na(dam), -1, FF[dam])
#     ans <- 1 - 0.25 * (2 + Fsire + Fdam)
#     names(ans) <- ped[,sib]
#     ans
#   }
#   relfactor <- function (ped, labs){
#     #stopifnot(is(ped, "pedigree"))
#     if (missing(labs)) 
#       return(Diagonal(x = sqrt(Dmat(ped))) %*% solve(t(as(ped, 
#                                                           "sparseMatrix"))))
#     labs <- factor(labs)
#     stopifnot(all(labs %in% ped[,sib]))
#     rect <- Diagonal(x = sqrt(Dmat(ped))) %*% solve(t(as(ped, 
#                                                          "sparseMatrix")), as(factor(ped[,sib], levels = ped[,sib]), 
#                                                                               "sparseMatrix"))
#     tmpA <- crossprod(rect)
#     tmp <- ped[,sib] %in% labs
#     tmpA <- tmpA[tmp, tmp]
#     orlab <- order(as.numeric(factor(labped <- ped[tmp,sib], 
#                                      levels = labs, ordered = T)))
#     labped <- as.character(labped[orlab])
#     tmpA <- tmpA[orlab, orlab]
#     stopifnot(all.equal(as.character(labped), as.character(labs)))
#     relf <- chol(tmpA)
#     dimnames(relf)[[1]] <- dimnames(relf)[[2]] <- labs
#     relf
#   }
#   
#   #stopifnot(is(ped, "pedigree"))
#   aMx <- crossprod(relfactor(ped))
#   dimnames(aMx)[[1]] <- dimnames(aMx)[[2]] <- ped[,sib]
#   aMx
# }