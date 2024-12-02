
#' ME-CDCs for Even Number of Lines
#'
#' @param lines Number of Lines >=6
#'
#' @return ME-CDCs for an even number of lines along with their parameters, C matrices, eigenvalues (EVs) and canonical efficiency factor (CEF).
#' @export
#'
#' @examples
#' library(MEDesigns)
#' ME_CDC(6)
ME_CDC<-function(lines){
  v=lines
  if(v%%2!=0 || v<6){
    return(message("Please enter an even number of lines."))
  }
initial_vector<-c(v,1:(v-1))
final_matrix<-matrix(NA,v,v)
for(i in 1:length(initial_vector)){
  if(initial_vector[i]%%2==0){
    vec1<-NULL
    for(j in 1:(v)){
      vec1<-c(vec1,(initial_vector[i]-j+1))
    }
    final_matrix[,i]<-vec1
  }
  if(initial_vector[i]%%2!=0){
    vec2<-NULL
    for(k in 1:(v)){
      vec2<-c(vec2,(initial_vector[i]+k-1))
    }
    final_matrix[,i]<-vec2
  }
}
final_matrix<-final_matrix%%v
final_matrix[final_matrix==0]<-v
################################################
half_matrix<-final_matrix[1:(v/2),]
##########
list_regions<-list()
for(l in 1:(v/2)){
  list_regions<-append(list_regions,list(half_matrix[,((2*l)-1):((2*l))]))
}
######################
############
zigzag<-function(matrix){
  seq1<-c()
  seq2<-c()
  for(i in 1:(v/2)){
      if(i%%2!=0){
        seq1<-c(seq1,matrix[i,1])
        seq2<-c(seq2,matrix[i,2])
      }else{
        seq1<-c(seq1,matrix[i,2])
        seq2<-c(seq2,matrix[i,1])
      }
  }
  combine<-c(seq1,seq2)
  return(combine)
}
#zigzag(list_regions[[4]])
#######################################
fixed_region<-function(matrix){
  vec<-rbind(t(t(matrix[,1])),t(t(matrix[,2])))
  return(vec)
}
##################
#fixed_region(list_regions[[4]])
####################################
######
##################################
diallel_crosses_blockwise<-NULL
possible_combinations<-t(combn((v/2),2))
for(i in 1:nrow(possible_combinations)){
diallel_crosses_blockwise<-rbind(diallel_crosses_blockwise,cbind(fixed_region(list_regions[[possible_combinations[i,1]]]),zigzag(list_regions[[possible_combinations[i,2]]])))
}
################
Blocks<-rep(1:(nrow(possible_combinations)),each=v)
###########
PDC<-cbind(Blocks,diallel_crosses_blockwise)
###################CEF code

# PDC_CEF<-function(PDC){
#   diallel_design=as.matrix(PDC)
#   lines=max(diallel_design[,2:3])
#   block_size=length(which(diallel_design[,1]==1))
#   ##crosses vs lines
#   x1_matrix=matrix(0,nrow(diallel_design),lines)
#   for(i in 1:nrow(diallel_design)){
#     x1_matrix[i,(diallel_design[i,2:3])]<-1
#   }
#   #####observation vs block matrix
#   x2_matrix=matrix(0,nrow(diallel_design),max(diallel_design[,1]))
#   for(j in 1:nrow(x2_matrix)){
#     x2_matrix[j,diallel_design[j,1]]<-1
#   }
#   #############3
#   x_mat=cbind(1,x1_matrix,x2_matrix)
#   ####3
#   x22=cbind(1,x2_matrix)
#   ###############C matrix
#   c_mat=(t(x1_matrix)%*%x1_matrix)-(t(x1_matrix)%*%x22)%*%MASS::ginv(t(x22)%*%x22)%*%(t(x22)%*%x1_matrix)
#   #######Contrast Matrix
#   # pos_mat=t(combinat::combn(lines,2))
#   # p_mat=matrix(0,nrow(pos_mat),lines)
#   #####
#   rep=length(which(diallel_design[,2:3]==1))
#   ##########
#   e1<-(eigen(c_mat)$values)
#   e1<-e1[e1>(10^-8)]
#   e2=e1/rep
#   e3=1/e2
#   CEF=length(e3)/sum(e3)
#   return(CEF)
# }

rep=length(which(PDC[,2:3]==1))
##########################################################################
result<-CheckME_Diallel(PDC)
colnames(PDC)<-c("Blocks","Line 1","Line 2")
list<-list("ME_CDC"=PDC,"Number of Lines"=lines,"Number of BLocks"=max(Blocks),"Block Size"=length(which(PDC[,1]==1)),"C_matrix"=round(result$C_matrix,3),"CEF"=round(result$CEF,3))
return(list)

}
