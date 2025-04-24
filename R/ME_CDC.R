
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
  if((v/2)%%2==0 || v<6){
    return(message("Please enter even number of lines >=6"))
  }
  ##############
  initial_vector<-c(1,(v):2)
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
  #}
  #################3
  blockwise_crosses<-list()
  for(i in 1:(length(list_regions)-1)){
    for(j in (i+1):length(list_regions)){
      b1=rbind(cbind((list_regions[[i]][,1]),(list_regions[[j]][,1])),
               cbind((list_regions[[i]][,2]),(list_regions[[j]][,2])))
      b2=rbind((list_regions[[i]]),(list_regions[[j]]))
      blockwise_crosses<-append(blockwise_crosses,list(rbind(b1,b2)))
    }
  }
  ###############
  final_blockwise_crosses<-blockwise_crosses
  ##########
  sort_mat<-function(mat){
    for(i in 1:nrow(mat)){
      mat[i,]<-sort(mat[i,])
    }
    return(mat)
  }
  #final_blockwise_crosses<-lapply(final_blockwise_crosses,sort_mat)
  #final_blockwise_crosses<-lapply(final_blockwise_crosses,function(mat)mat[order(mat[,1]),])
  design=cbind(rep(1:((v/2)*((v/2)-1)),each=v),do.call(rbind,final_blockwise_crosses))
  lm=(CheckME_Diallel(design))
  names(lm)[1]<-"ME_CDC"
  return(lm)
}

