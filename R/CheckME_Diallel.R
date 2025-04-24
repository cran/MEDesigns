
#' Checking the Properties of a ME-PDC
#'
#' @param design Provide a ME-PDC 
#'
#' @return Generates parameters of the designs along with C matrix, eigenvalues (EVs), degree of fractionations (DF) and canonical efficiency factor (CEF).
#' @export
#'
#' @examples
#' library(MEDesigns)
#' design<-ME_PDC1(10)$ME_PDC
#' CheckME_Diallel(design)
CheckME_Diallel<-function(design){
  Moore_penrose_inverse_mine<-function(matrix){
    N<-as.matrix(matrix)
    NtN<-t(N)%*%N
    NNt<-N%*%t(N)
    eig_val<-eigen(NtN)$values
    positive<-NULL
    for(i in eig_val){
      if(i>10^-7){
        positive<-c(positive,TRUE) 
      }else{
        positive<-c(positive,FALSE) 
      }
    }
    sigma<-1/sqrt(eig_val[eig_val>10^-7])
    u<-eigen(NtN)$vectors[,positive,drop=FALSE]
    v<-eigen(NNt)$vectors[,positive,drop=FALSE]
    return(u%*%diag(sigma,nrow(t(v)))%*%t(v))
  }
diallel_design=as.matrix(design)
lines=max(diallel_design[,2:3])
block_size=length(which(diallel_design[,1]==1))
##crosses vs lines
x1_matrix=matrix(0,nrow(diallel_design),lines)
for(i in 1:nrow(diallel_design)){
  x1_matrix[i,(diallel_design[i,2:3])]<-1
}
#####observation vs block matrix
x2_matrix=matrix(0,nrow(diallel_design),max(diallel_design[,1]))
for(j in 1:nrow(x2_matrix)){
  x2_matrix[j,diallel_design[j,1]]<-1
}
#############3
x_mat=cbind(1,x1_matrix,x2_matrix)
####3
x22=cbind(1,x2_matrix)
###############C matrix
c_mat=(t(x1_matrix)%*%x1_matrix)-(t(x1_matrix)%*%x22)%*%Moore_penrose_inverse_mine(t(x22)%*%x22)%*%(t(x22)%*%x1_matrix)
#######Contrast Matrix
# pos_mat=t(combinat::combn(lines,2))
# p_mat=matrix(0,nrow(pos_mat),lines)
#####
rep=length(which(diallel_design[,2:3]==1))

##########
e1<-(eigen(c_mat)$values)
e1<-e1[e1>(10^-8)]
e2=e1%*%solve(rep)
e3=1/e2
CEF=length(e3)/sum(e3)
#return(CEF)
##Degree of Fractinations
PDC<-diallel_design
v<-max(PDC[,2:3])
N_cdc<-choose(v,2)
s<-t(apply(PDC[,2:3],1,sort))
N_pdc<-nrow(unique(s))
DF<-N_pdc/N_cdc
colnames(PDC)<-c("Blocks","Line 1","Line 2")
list<-list("ME_PDC"=PDC,"Number of lines"=max(c(PDC[,2])),"Number of BLocks"=max(c(PDC[,1])),"Block Size"=length(which(PDC[,1]==1)),"C_matrix"=round(c_mat,3),"EVs"=round(unname(table(e1),3)),"CEF"=round(CEF,3),"DF"=round(DF,3))
return(list)
}

