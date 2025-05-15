##CODE FOR OBTAINING LEFT DUASE EMBEDDINGS FROM LLM-DATA

library(MASSExtra)
library(MatrixExtra)
library(irlba)
library(MVA)
library(smacof)
library(wordspace)

set.seed(1234)

#declaring model parameters
n_model<-30
n_time<-30
n_component<-768
n_replica<-8
seq_nnbr<-c(seq(1,28,3),29)
seq_time<-1:n_time


nnbr<-4

#generic function to store an input as a matrix
store_as_matrix<-function(X)
{
  Y<-as.matrix(X)
  
  rownames(Y)<-NULL
  colnames(Y)<-NULL
  
  return(Y)
}



#creating adjacency matrix from vector of inter-neighbour communication 
create_matrix<-function(x)
{
  M<-matrix(0,nrow=length(x),ncol=length(x))
  for(i in 1:length(x))
  {
    M[i,x[i]]<-1
  }
  return(M)
}





######Regression per replica amalgamating info from all time#####




  
  A_list<-list()
  
  for(replica in 1:n_replica)
  {
    
    file_rd_id<-c("nb",as.character(nnbr),"_",as.character(replica-1),".csv")
    DD<-t(read.csv(paste(file_rd_id,collapse = ''))[,-1])
    
    DD<-store_as_matrix(DD)
    
    
    
    
    DD_col_list<-lapply(seq_len(ncol(DD)), function(j) DD[,j])
    
    A_list[[replica]]<-store_as_matrix(do.call(cbind,
                                               lapply(DD_col_list,create_matrix)))
    
  }
  
  
  
  AA_grand<-store_as_matrix(do.call(rbind,A_list))
  



d<-8

AA_grand_irlba<-irlba(AA_grand,d)
X_left<-AA_grand_irlba$u%*%(diag(AA_grand_irlba$d)^0.5)
Y_right<-AA_grand_irlba$v%*%(diag(AA_grand_irlba$d)^0.5)

print(X_left)


df_X_left<-as.data.frame(X_left)
file_X_out_id<-c("DUXk",".RData")
save(df_X_left,file=paste(file_X_out_id,collapse = ''))


df_Y_right<-as.data.frame(Y_right)
file_Y_out_id<-c("DUYk",".RData")
save(df_Y_right,file=paste(file_Y_out_id,collapse = ''))










