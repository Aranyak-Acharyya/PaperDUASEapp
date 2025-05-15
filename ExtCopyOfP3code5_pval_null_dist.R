##CODE FOR TYPE-I ERROR CONTROL UNDER NULL HYPOTHESIS, REFER TO FIG-6##


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
library(stats)


RNGkind("L'Ecuyer-CMRG")
set.seed(1234)




e <- new.env()
e$libs <- c("Matrix","MASS","MVA","irlba","igraph",
            "smacof","pracma","lsbclust","wordspace",
            "stats",.libPaths())
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
beta<-0.0
sigma_ep<-10^(-2)

#number of Monte Carlo trials
n_trial<-100


#level of significance
lvl<-0.05

#threshold for F-statistic
thres<-qf(1-lvl,1,s-2)



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


#experiment for a single $K$
K<-10

n<-n_func(K)
N<-N_func(K)
M<-M_func(K)

clusterExport(clust,list("n","N","M"))

err_vec<-vector()

PS<-foreach(trial=1:n_trial,.combine = 'rbind') %dopar%
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
    
    
    #model with true regressors
    fit_true<-lm(y~t_lab)
    summary(fit_true)
    p_true<-summary(fit_true)$coefficients[2,4] 
    
    
    #model with 1d embeddings
    fit_sub<-lm(y~z_lab)
    summary(fit_sub)
    p_sub<-summary(fit_sub)$coefficients[2,4] 
    
    
    print(p_true)
    
    
    
    
    store<-c(p_true,p_sub)
    store
    
  }


stopCluster(clust)

print(PS)

p_true_vec<-PS[,1]
p_sub_vec<-PS[,2]

#creating a vector of Unif(0,1) values for reference
ref_vec<-runif(100,min=0,max=1)


#computing empirical CDFs of p-values for true regressors and isomap embeddings
ecdf_true<-ecdf(p_true_vec)
ecdf_sub<-ecdf(p_sub_vec)
ecdf_ref<-ecdf(ref_vec)
ecdf_unif<-function(x)
{
  return(x)
}

arg_vec<-seq(0,1,0.001)
ecdf_true_vec<-sapply(seq(0,1,0.001),ecdf_true)
ecdf_sub_vec<-sapply(seq(0,1,0.001),ecdf_sub)
ecdf_ref_vec<-sapply(seq(0,1,0.001),ecdf_unif)

library(ggplot2)
library(reshape2)
library(latex2exp)
library(melt)


df<-data.frame(arg_vec,ecdf_true_vec,ecdf_sub_vec,ecdf_ref_vec)
print(df)
dfm<-melt(df,id.vars='arg_vec')
g<-ggplot(dfm, aes(x=arg_vec, y=value, 
                colour = variable)) +
  #geom_point() +
  geom_line() +
  ylab(TeX("Empirical CDF")) +
  xlab(TeX("values")) +
  scale_colour_manual(values = c("black","blue","orange"),
                      name="Colours",
                      labels=unname(TeX(c(
                        "p-values from true regressors",
                        "p-values from raw-stress embeddings",
                        "CDF of U(0,1)$"
                      )))) 

g

ggsave(g,file="P3rev1_plot_pval_distn.pdf",
       height = 3, width = 6.5,
       units = "in",dpi = 500)








df<-data.frame(x=c(p_true_vec,p_sub_vec,ref_vec),g=gl(3,length(p_true_vec)))
print(df)




ggplot(df, aes(x,colour = g)) +
  stat_ecdf(geom = "step") +
  xlab("regressor values") +
  ylab("Empirical CDF") 






  


df<-as.data.frame(cbind(p_true_vec,ecdf_true))
g<-ggplot(df,aes(x=p_true_vec,y=ecdf_true)) +
  geom_point() +
  geom_line() +
  xlab(TeX("K")) +
  ylab(TeX("$|\\hat{\\pi}-\\pi^*|$"))


g









#repeat experiment by varying $K$
for(K in K_vec)
{
  K<-3
  
  n<-n_func(K)
  N<-N_func(K)
  M<-M_func(K)
  
  clusterExport(clust,list("n","N","M"))
  
  err_vec<-vector()
  
  PS<-foreach(trial=1:n_trial,.combine = 'rbind') %dopar%
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
    
    
    #model with true regressors
    fit_true<-lm(y~t_lab)
    summary(fit_true)
    p_true<-summary(fit_true)$coefficients[2,4] 
    
    
    #model with 1d embeddings
    fit_sub<-lm(y~z_lab)
    summary(fit_sub)
    p_sub<-summary(fit_sub)$coefficients[2,4] 
    
    
    print(p_true)
    
    
    
    
    store<-c(p_true,p_sub)
    store
    
  }
  
  print(PS)
  
  ecdf_sub<-ecdf(PS[,2])
  
  #power of test based on true F-statistic
  power_true<-mean(ifelse(PS[,1]>thres,1,0))
  
  #power of test based on substitute F-statistic
  power_sub<-mean(ifelse(PS[,2]>thres,1,0))
  
  #difference in true power and substitute power, stored in a vector
  val<-abs(power_true-power_sub)
  val_vec<-c(val_vec,val)
  
  dec<-c(K,val)

  print(dec)  
}
  

stopCluster(clust)


df<-as.data.frame(cbind(K_vec,val_vec))
save(df,file = "P3dat5_power_conv.RData")


df<-get(load("P3dat5_power_conv.RData"))

print(df)



library(ggplot2)
library(reshape2)
library(latex2exp)

#df<-as.data.frame(cbind(K_vec,val_vec))
g<-ggplot(df,aes(x=K_vec,y=val_vec)) +
  geom_point() +
  geom_line() +
  xlab(TeX("K")) +
  ylab(TeX("$|\\hat{\\pi}-\\pi^*|$"))

g


ggsave(g,file="P3plot5_power_conv.pdf",
     height = 3, width = 5,
     units = "in",dpi = 500)


