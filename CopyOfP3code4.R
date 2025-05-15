###CODE TO SHOW THAT SAMPLE MEAN OF $(\hat{y}_{true}-\hat{y}_{sub})^2$
#### CONVERGES TO ZERO AS SAMPLE SIZE GROWS


library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()


library(Matrix)
library(MASS)
library(MVA)
library(irlba)
library(igraph)
library(smacof)
library(pracma)
library(lsbclust)
library(wordspace)


RNGkind("L'Ecuyer-CMRG")
set.seed(1234)




e <- new.env()
e$libs <- c("Matrix","MASS","lsbclust",
            "MVA","irlba","igraph","pracma","wordspace",
            "smacof",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust,.libPaths(libs))











#sequence of common index K
K_vec<-seq(1,20,1)

#embedding dimension of DUASE
d<-2

#number of labeled points
s<-5

#regression parameters
alpha<-2.0
beta<-8.0
sigma_ep<-10^(-2)

#number of Monte Carlo trials
n_trial<-100



#number of nodes in each graph as function of K
n_func<-function(K)
{
  return(as.integer(15+(K-1)^1.5))
}

#number of graphs in each TSG as function of K
M_func<-function(K)
{
  return(as.integer(8+(K-1)))
}

#number of TSGs as function of K
N_func<-function(K)
{
  return(as.integer(10+(K-1)))
}






#function to generate an adjacency matrix from a given probability matrix
gen_adj_from_prob<-function(P)
{
  n<-nrow(P)
  pvec<-c(P)
  avec<-rbinom(length(pvec),1,pvec)
  A<-matrix(avec,nrow=n)
  
  return(A)
}








#function to store an entity as a matrix
store_as_matrix<-function(X)
{
  Y<-as.matrix(X)
  
  rownames(Y)<-NULL
  colnames(Y)<-NULL
  
  return(Y)
}


#performing dual unfolded adjacency spectral embedding
DUASE<-function(A_grand_mat,n,d)
{
  N<-as.integer(nrow(A_grand_mat)/n)
  M<-as.integer(ncol(A_grand_mat)/n)
  
  A_grand_irlba<-irlba(A_grand_mat,d)
  X_left<-A_grand_irlba$u%*%(diag(A_grand_irlba$d)^0.5)
  Y_right<-A_grand_irlba$v%*%(diag(A_grand_irlba$d)^0.5)
  
  return(list(X_left,Y_right))
}


#finding 2-to-infinity norm of difference between two matrices
f<-function(X,Y)
{
  m<-max(rowNorms(X-Y,method = "euclidean",p=2))
  return(m)
}


clusterExport(clust,list("K_vec","d","s","alpha","beta","n_trial"))
clusterExport(clust,list("n_func","M_func","N_func","gen_adj_from_prob",
              "store_as_matrix","DUASE","f"))


val_vec<-vector()



for(K in K_vec)
{
  
  n<-n_func(K)
  N<-N_func(K)
  M<-M_func(K)
  
  clusterExport(clust,list("n","N","M"))
  
  err_vec<-vector()
  
  E<-foreach(trial=1:n_trial,.combine = 'c') %dopar%
  {
    
    
    #generating regressors
    t_lab<-runif(s,min=0,max=1)
    t_aux<-runif(N-s,min=0,max=1)
    t_vec<-c(t_lab,t_aux)
    
    
    #generating responses
    y<-alpha+beta*t_lab+rnorm(s,mean=0,sd=sigma_ep)
    
    
    
    #generating left-DUASE-emb for grand-P-matrix
    XP_list<-lapply(t_vec,function(t_val) matrix(t_val/d^0.5,nrow=n,ncol=d))
    XP<-store_as_matrix(do.call(rbind,XP_list))
    
    
    
    #generating right-DUASE-emb for grand-P-matrix
    YP_list<-lapply(1:M,
                    function(i) matrix(runif(1,min = 0.2,max=0.5),nrow=n,ncol=d)) 
    YP<-store_as_matrix(do.call(rbind,YP_list))
    
    
    #grand probability matrix
    P_grand<-XP%*%t(YP)
    
    
    
    
    #generating the grand adjacency matrix
    A_grand<-gen_adj_from_prob(P_grand)
    
    
    
    
    #finding left and right DUASE embeddings of the grand-A-matrix
    X_left<-DUASE(A_grand,n,d)[[1]]
    
    
    #storing the X-matrices per regressor as elements in a list
    X_list<-lapply(1:N,function(i) X_left[((i-1)*n+1):(i*n),] )
    
    #finding pairwise distances between DUASE-emb-matrices
    D<-proxy::dist(X_list,method = f,diag=TRUE,upper=TRUE)
    Dt<-store_as_matrix(D)
    
    
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
    z_vec<-ZZ/fac
    
    z_lab<-z_vec[1:s]
    
    
    
    #true estimates of regression parameters
    beta_true<-cov(y,t_lab)/var(t_lab)
    alpha_true<-mean(y)-beta_true*mean(t_lab)
    
    
    #substitute estimate of regression parameters
    beta_sub<-cov(y,z_lab)/var(z_lab)
    alpha_sub<-mean(y)-beta_sub*mean(z_lab)
    
    
    
    #true predicted response
    y_hat_true<-alpha_true+beta_true*t_vec[s+1]
    
    
    #substitute predicted response
    y_hat_sub<-alpha_sub+beta_sub*z_vec[s+1]
    
    
    err<-(y_hat_sub-y_hat_true)^2
    
    
    err
  }
  
  val<-mean(E)
  val_vec<-c(val_vec,val)
  
  dec<-c(K,val)

  print(dec)  
}
  

stopCluster(clust)




#KK_vec<-K_vec[1:length(val_vec)]
df<-as.data.frame(cbind(K_vec,val_vec))
save(df,file = "P3dat4_pred_res_const.RData")


df<-get(load("P3dat4_pred_res_const.RData"))

print(df)





library(ggplot2)
library(reshape2)
library(latex2exp)

#df<-as.data.frame(cbind(KK_vec,val_vec))
g<-ggplot(df,aes(x=K_vec,y=val_vec)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0,7.5)) +
  xlab(TeX("K")) +
  ylab(TeX("sample E$(\\hat{y}_{sub}-\\hat{y}_{true})^2$"))

g


ggsave(g,file="P3plot4_pred_res_const.pdf",
     height = 3, width = 5,
     units = "in",dpi = 500)


