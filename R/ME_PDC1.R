
#' ME-PDCs for Even Number of Lines
#'
#' @param lines Number of Lines >=6
#'
#' @return ME-PDCs for an even number of lines along with their parameters, C matrices, eigenvalues (EVs), degree of fractionations (DF) and canonical efficiency factor (CEF).
#' @export
#'
#' @examples
#' library(MEDesigns)
#' ME_PDC1(6)
ME_PDC1<-function(lines){
  v=lines
  if(v%%2!=0 || v<6){
    return(message("Please enter an even number of lines(>=6)."))
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
#########
v_hlf<-v/2
diallel<-NULL
for(i in 1:(v/2)){
  zigzag1<-c()
  zigzag2<-c()
  select_reg<-half_matrix[,((2*i-1):(2*i))]
  for(j in 1:(v/2)){
  if(j%%2!=0){
  zigzag1<-c(zigzag1,select_reg[j,1])
  zigzag2<-c(zigzag2,select_reg[j,2]) 
  }else{
    zigzag1<-c(zigzag1,select_reg[j,2]) 
    zigzag2<-c(zigzag2,select_reg[j,1])
  }
  }
  # print(zigzag1)
  # print(zigzag2)
  #########make pairs
  pairs1<-NULL
  pairs2<-NULL
  for(k in 1:(length(zigzag1)-1)){
    pairs1<-rbind(pairs1,zigzag1[c(k,(k+1))])
    pairs1<-rbind(pairs1,zigzag2[c(k,(k+1))])
  }
  diallel<-rbind(diallel,rbind(pairs1,pairs2))
#########
# crosses<-rbind(cbind(rep(ele1,v_hlf-1),select_reg[2:v_hlf,2]),
  #                cbind(rep(ele2,v_hlf-1),select_reg[2:v_hlf,1]) )
  # diallel<-rbind(diallel,crosses)
}
 final_diallel<-cbind(rep(1:((v_hlf)),each=2*(v_hlf-1)),diallel)
#return(CheckME_Diallel(final_diallel))
 ##Degree of Fractinations
 PDC=final_diallel
 N_cdc<-choose(lines,2)
 s<-t(apply(PDC[,2:3],1,sort))
 N_pdc<-nrow(unique(s))
 DF<-N_pdc/N_cdc
 result<-CheckME_Diallel(PDC)
 colnames(PDC)<-c("Blocks","Line 1","Line 2")
 list<-list("ME_PDC"=PDC,"Number of Lines"=lines,"Number of BLocks"=max(c(PDC[,1])),"Block Size"=length(which(PDC[,1]==1)),"C_matrix"=round(result$C_matrix,3),"CEF"=round(result$CEF,3),"DF"=round(DF,3))
return(list)
 }
