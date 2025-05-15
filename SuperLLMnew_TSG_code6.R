##CODE FOR OBTAINING RESPONSES AND STUDYING VALIDITY OF  
### LINEAR REGRESSION MODEL FROM LLM-DATA 


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
n_replica<-40
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




#finding 2-to-infinity norm of difference between two matrices
f<-function(X,Y)
{
  m<-max(rowNorms(X-Y,method = "euclidean",p=2))
  return(m)
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


### Every replica corresponds to a time series



Ans_grand_list<-list()

  
  Ans<-matrix(,ncol=768)
  
  for(replica in 1:n_replica)
  {

    file_rd_id<-c("LLMdnew",as.character(nnbr),"_",as.character(replica-1),".csv")
    DD<-read.csv(paste(file_rd_id,collapse = ''))[,-1]
    
    DD_mat<-store_as_matrix(DD)
    
    DD_mat<-DD_mat[((n_time-1)*n_model+1):(n_time*n_model),]

    
    #Ans<-rbind(Ans,DD_mat)
    
    Ans_grand_list[[replica]]<-DD_mat
  }
  
  #Ans<-Ans[-1,]
  


M<-store_as_matrix(do.call(rbind,Ans_grand_list))
print(M)

n_series<-n_replica

Dt_list<-lapply(1:n_series, 
       function(series)
         dist(M[((series-1)*n_model+1):(series*n_model),],method="euclidean",diag=TRUE,upper=TRUE,p=2)
       )



#Assigning a response to every TSG
y<-sapply(1:n_series, 
               function(series)
                 mean(store_as_matrix(Dt_list[[series]])[upper.tri(store_as_matrix(Dt_list[[series]]))])
                 )

print(y)










### Obtaining left DUASE embeddings

load("DUXk.RData")

print(df_X_left)
X_left<-store_as_matrix(df_X_left)
print(X_left)


X_left_list<-lapply(1:n_series,
                    function(series)
                      X_left[((series-1)*n_model+1):(series*n_model),])

Dt<-proxy::dist(X_left_list,method = f,diag=TRUE,upper=TRUE)
print(Dt)







#raw-stress minimization
MM<-mds(Dt,ndim = 1,type = "ratio",
        weightmat = NULL,
        init = "torgerson")


#finding the scaling factor
dl<-as.vector(MM$delta)
dh<-as.vector(MM$dhat)
dvec<-dh/dl
dvec<-na.omit(dvec)
fac<-mean(dvec)


#raw-stress embeddings
ZZ<-as.vector(MM$conf)
Z<-ZZ/fac



#regressors
x<-Z

print(x)





#studying validity of simple linear regression model linking responses to regressors
df_lm<-data.frame(y,x)
print(df_lm)

fit<-lm(y~.,data = df_lm)
summary(fit)




Zmm<-matrix(Z,nrow=n_model*n_replica,byrow=FALSE)

print(Zmm)

df_Zmm<-as.data.frame(Zmm)

save(df_Zmm,file="res.RData")

load("res.RData")

print(df_Zmm)



