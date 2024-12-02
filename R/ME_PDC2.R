#' ME PDCs for Composite Number of Lines 
#'
#' @param p Any value (p>=3)
#' @param q Any value (q>=3)
#' @return This function will provide ME-PDCs for a composite number, v(= pq) along with basic parameters, C matrix, eigenvalues (EVs), degree of fractionations (DF) and canonical efficiency factor (CEF).
#' @export
#'
#' @examples
#' library(MEDesigns)
#' ME_PDC2(3,3)
ME_PDC2<-function(p,q){
  v=p*q
  lines=v
  if(p<3 || q<3){
    return(message("Please enter a composite number(p*q) of lines where p,q >= 3."))
  }
  #########factors for p and q
  # find_all_factors <- function(v) {
  #   factors <- list()
  # 
  #   # Iterate over possible values of p starting from 3
  #   for (p in 3:floor(sqrt(v))) {
  #     if (v %% p == 0) {
  #       q <- v / p
  #       # Check if both p and q are greater than or equal to 3
  #       if (q >= 3) {
  #        # factors <- append(factors, list(c(p, q)))
  #         factors <- append(factors, list(c(q, p)))  # Add the reverse pair as well
  #       }
  #     }
  #   }
  # 
  #   # If no factors are found
  #   if (length(factors) == 0) {
  #     return(FALSE)
  #   } else {
  #     return(factors)
  #   }
  # }
# ###########
# ##########
#   
# p<-find_all_factors(v)[[1]][1]
# q<-find_all_factors(v)[[1]][2]
#v=p*q

scheme<-matrix(1:v,p,q,byrow=TRUE)
###########matrix rotation rownwise
row_rotation<-function(matrix){
vector<-c(1:nrow(matrix))
mat<-NULL
for(i in 0:(length(vector)-1)){
  mat<-rbind(mat,vector+i)
}
mat<-mat%%max(vector)
mat[mat==0]<-max(vector)
return(mat)
}
#########################
base_rotation<-row_rotation(scheme)
list_of_regions<-list()
for(i in 1:ncol(base_rotation)){
  list_of_regions<-append(list_of_regions,list(scheme[c(base_rotation[,i]),]))
}
############# matrix rotation for each columnwise rotation
col_rotation<-function(matrix){
  vector<-c(1:ncol(matrix))
  mat<-NULL
  for(i in 0:(length(vector)-1)){
    mat<-rbind(mat,vector+i)
  }
  mat<-mat%%max(vector)
  mat[mat==0]<-max(vector)
  return(mat)
}
####################################
base_rotation_for_column<-col_rotation(scheme)
finallist<-list()
for(i in list_of_regions){
  blanklist<-list()
  for(j in 1:nrow(base_rotation_for_column)){
    blanklist<-append(blanklist,list(i[,c(base_rotation_for_column[j,])]))
  }
  finallist<-append(finallist,list(blanklist))
}

#############
################################## final regions of rectangular 3rd AS
rectangular_scheme<-finallist
  for(k in 1:length(finallist)){
    for(l in 1:length(finallist[[k]])){
      rectangular_scheme[[k]][[l]]<-rectangular_scheme[[k]][[l]][-1,-1]
    }
  }
  ################################take row wise regions
rowwise_regions<-list()
for(k in 1:length(rectangular_scheme[[1]])){
  rowwise_regions<-append(rowwise_regions,
                          list(lapply(rectangular_scheme,function(l) l[[k]])))
}

###############
combinations<-t(combn(length(rowwise_regions[[1]]),2))
############
###for 1<- 11,22,33 2<-12,23,31 and 3<-13,21,32 need to generate for zigzag
mat1<-matrix(1:nrow(rowwise_regions[[1]][[1]]),nrow(rowwise_regions[[1]][[1]]),ncol(rowwise_regions[[1]][[1]]))
zigzag<-list()
for(i in 1:ncol(mat1)){
  zigzag<-append(zigzag,list(cbind(mat1[,i],mat1[,i]+(i-1))))
}
#########
for(i in 1:length(zigzag)){
  zigzag[[i]][,2]<- zigzag[[i]][,2]%%ncol(rowwise_regions[[1]][[1]])
  zigzag[[i]][,2][zigzag[[i]][,2]==0]<-ncol(rowwise_regions[[1]][[1]])
}
##############
##################
#crosses<-cbind(c(rowwise_regions[[1]][[1]]),)
# zigzag<-list()
# for(i in 1:ncol(mat1)){
#   zigzag<-append(zigzag,list(cbind(mat1[,i],mat2[,i])))
# }
##############
zigzag_mat<-do.call(rbind,zigzag)
zigzag_elements<-function(matrix1,matrix2){
zigzag_elements<-NULL
  for(x in 1:nrow(matrix2)){
    zigzag_elements<-c(zigzag_elements,matrix1[matrix2[x,1],matrix2[x,2]])
}
return(zigzag_elements)
}

#######################################

####################################
######
##################################
diallel_crosses_blockwise<-NULL
possible_combinations<-combinations
for(r in 1:length(rowwise_regions)){
for(i in 1:nrow(possible_combinations)){
  fixed_mat<-rowwise_regions[[r]][[possible_combinations[i,1]]]
  manupulate_mat<-rowwise_regions[[r]][[possible_combinations[i,2]]]
  diallel_crosses_blockwise<-rbind(diallel_crosses_blockwise,cbind(c(fixed_mat),zigzag_elements(manupulate_mat,zigzag_mat)))
}
  }
################
rep_vec<-NULL
for(i in 1:v){
  rep_vec<-c(rep_vec,length(which(diallel_crosses_blockwise==i)))
}
######################
Blocks<-rep(1:(nrow(diallel_crosses_blockwise)/(length(rowwise_regions[[1]][[1]]))),each=length(rowwise_regions[[1]][[1]]))
###########
PDC<-cbind(Blocks,diallel_crosses_blockwise)
################################################CEF code

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
##########################################################################
##Degree of Fractinations
N_cdc<-choose(v,2)
s<-t(apply(PDC[,2:3],1,sort))
N_pdc<-nrow(unique(s))
DF<-N_pdc/N_cdc
result<-CheckME_Diallel(PDC)
colnames(PDC)<-c("Blocks","Line 1","Line 2")
list<-list("ME_PDC"=PDC,"Number of Lines"=p*q,"Number of BLocks"=max(Blocks),"Block Size"=length(rowwise_regions[[1]][[1]]),"C_matrix"=round(result$C_matrix,3),"CEF"=round(CheckME_Diallel(PDC)$CEF,3),"DF"=round(DF,3))
return(list)
}
