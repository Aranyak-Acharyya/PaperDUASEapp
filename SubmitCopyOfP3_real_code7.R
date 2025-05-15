####CODE__FOR__REAL__DATA__ANALYSIS OF DROSOPHILA########


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




load("missingEdge-tsg-forAranyak-2.RData")
print(data)
df1<-data[,-2]
print(df1)
col.order<-c("model","score","tg","missing_edge")
df1[,col.order]

load("missingEdge-tsg-forAranyak.RData")
print(data)
df2<-data
data<-rbind(df1,df2)
k_vec<-seq(1,nrow(data),1)
n<-140
d<-3

print(data)

data$model




#specifying model-type
mspec<-vector()
for(i in 1:143)
{
  if(i<=130)
  {
    if(i%%10!=0)
    {
      q<-i%/%10
      mspec[i]<-(q+1)
    }
    if(i%%10==0)
    {
      q<-i%/%10
      mspec[i]<-q
    }
  }
  if(i>130)
  {
    mspec[i]<-(i-130)
  }
}





data$model<-mspec


print(data$model)
length(mspec)



#function to store an entity as a matrix
store_as_matrix<-function(X)
{
  Y<-as.matrix(X)
  
  rownames(Y)<-NULL
  colnames(Y)<-NULL
  
  return(Y)
}


#choose threshold for censoring adjacency matrix
p<-0.25



#function to store a graph as its adjacency matrix
store_graph_as_matrix<-function(g)
{
  weight<-E(g)$weight
  
  
  p<-0.25
  
  
  A<-as_adjacency_matrix(
    g,
    type = "both",
    attr = "weight",
    names = FALSE
  )
  
  
  vec_modA<-as.vector(abs(A))
  vnz_modA<-vec_modA[!vec_modA %in% c(0)]
  
  
  qq<-quantile(vnz_modA,probs=p)
  
  
  
  AA<-matrix(ifelse(as.vector(abs(A))>qq,1,0),
             nrow=n,byrow=FALSE)
  
  diag(AA)<-0
  
  return(AA)
}



#function for dual unfolded adjacency spectral embedding
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




#function to estimate sparsity factor from a single adjacency matrix
estimate_sparsity<-function(A)
{
  return(mean(A[upper.tri(A,diag=FALSE)]))
}


n<-gorder(data$tg[[1]][[1]])


N<-length(data$model)

M<-length(data$tg[[1]])


clusterExport(clust,list("d","n","p","k_vec","N","M"))
clusterExport(clust,list("DUASE","f"))








#calculating the grand-A-matrix and global sparsity estimate simultaneously
D_grand<-foreach(k=k_vec) %dopar%
  {
    AA_list<-lapply(1:M,function(i) store_graph_as_matrix(data$tg[[k]][[i]]))
    
    A_TSG<-store_as_matrix(do.call(cbind,AA_list))
    
    m<-mean(sapply(AA_list,estimate_sparsity))
    
    ll<-list(A_TSG,m)
    
    ll
    
  }

#global estimate for sparsity factor
rho_hat<-mean(sapply(k_vec,function(k) D_grand[[k]][[2]]))



#grand-A-matrix
A_grand_list<-lapply(k_vec,function(k) D_grand[[k]][[1]])
A_grand<-store_as_matrix(do.call(rbind,A_grand_list))






print(A_grand)






d<-2

XA<-DUASE(A_grand,n,d)[[1]]

XA<-(1/(rho_hat^0.5))*XA


#storing the X-matrices per regressor as elements in a list
XA_list<-lapply(1:N,function(i) XA[((i-1)*n+1):(i*n),] )

#finding pairwise distances between DUASE-emb-matrices
D<-proxy::dist(XA_list,method = f,diag=TRUE,upper=TRUE)
Dt<-store_as_matrix(D)

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
z_vec<-ZZ/fac


#proxy/substitute regression model
y<-data$score
z<-z_vec
df<-as.data.frame(cbind(k_vec,y,z,mspec))

mod<-lm(y~z,data=df)

summary(mod)





stopCluster(clust)







library(ggplot2)
library(GGally)
library(reshape2)
library(latex2exp)
library(ggpmisc)
library(ggpubr)
library(broom)
library(np)

x<-seq(0,10,length.out=5)
print(x)


mnx<-(-5.0)
mxx<-10.0

mny<-0.50
mxy<-2.50








gr1<-ggplot(df,aes(x=z,y=y)) + 
  geom_point() +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$"))+
  scale_x_continuous(name = TeX("$\\hat{z}_i$"),
                     breaks = seq(mnx,mxx,length.out=11),
                     labels = seq(mnx,mxx,length.out=11),
                     limits = c(mnx,mxx)) +
  scale_y_continuous(name = TeX("$y_i$"),
                     breaks = seq(mny,mxy,length.out=11),
                     labels = seq(mny,mxy,length.out=11),
                     limits = c(mny,mxy)) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_correlation(mapping = use_label(c("R", "n","p","R2")),
                   label.x = "left",
                   label.y = "top") 


gr1<-ggplot(df,aes(x=z,y=y)) + 
  geom_point() +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$"))+
  scale_x_continuous(name = TeX("$\\hat{z}_i$"),
                     breaks = seq(mnx,mxx,length.out=11),
                     labels = seq(mnx,mxx,length.out=11),
                     limits = c(mnx,mxx)) +
  scale_y_continuous(name = TeX("$y_i$"),
                     breaks = seq(mny,mxy,length.out=11),
                     labels = seq(mny,mxy,length.out=11),
                     limits = c(mny,mxy)) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_correlation(mapping = use_label(c("R", "n","p","R2")),
                   label.x = "left",
                   label.y = "top") 


gr1






ggsave(file="P3plot7_real.pdf",plot=gr1,
       height=3,width=5,units = "in",dpi=500)





z_modnp<-z
print(z)
print(z_modnp)
df_modnp<-data.frame(y,z,z_modnp)
print(df_modnp)


modnp<-npreg(y~z_modnp,
             regtype="ll",
             bwmethod="cv.aic",
             gradients=TRUE,
             data=df_modnp)

summary(modnp)

plot(modnp,gradients=TRUE,plot.errors.method="asymptotic")

#yhatnp<-predict(modnp,newdata=df)



z_modnp<-seq(min(z)-0.005,max(z)+0.005,length.out=400)
print(z_modnp)

yhatnp<-predict(modnp,newdata=as.data.frame(z_modnp))

print(yhatnp)

dfnp2<-data.frame(z_modnp,yhatnp)



#dfnp<-data.frame(z,y,yhatnp)

#dfnpm<-reshape::melt(dfnp,id.vars = 'z')

grnp<-ggplot(df,aes(x=z,y=y)) +
  geom_point() +
  geom_line(data=dfnp2,aes(x=z_modnp,y=yhatnp),color="blue") +
  xlab(TeX("$\\hat{z}_i$")) +
  ylab(TeX("$y_i$"))

grnp

ggsave(file="P3res_vs_embedding_np.pdf",
       grnp,
       width=5,
       height=3,
       units="in",
       dpi = 500)


gr_test<-ggplot(df,aes(x=z,y=y)) + 
  geom_point() +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$"))+
  scale_x_continuous(name = TeX("$\\hat{z}_i$"),
                     breaks = seq(mnx,mxx,length.out=11),
                     labels = seq(mnx,mxx,length.out=11),
                     limits = c(mnx,mxx)) +
  scale_y_continuous(name = TeX("$y_i$"),
                     breaks = seq(mny,mxy,length.out=11),
                     labels = seq(mny,mxy,length.out=11),
                     limits = c(mny,mxy)) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth",
              colour="blue") +
  #geom_smooth(method="loess",se=F,colour="red")+
  geom_line(data=dfnp2,aes(x=z_modnp,y=yhatnp),colour="red")+
  stat_correlation(mapping = use_label(c("R", "n","p","R2")),
                   label.x = "left",
                   label.y = "top") 

gr_test

ggsave(file="P3plot7_real_new.pdf",plot=gr_test,
       height=3,width=5,units = "in",dpi=500)






















p<-summary(mod)$coefficients[2,4]

rh<-cor(y,z)

res<-c(tim,p,rh,lambda)

print(res)








library(ggplot2)
library(reshape2)
library(latex2exp)
library(ggpmisc)
library(ggpubr)
library(broom)


lb<-c("one","two","three","four",
      "five","six","seven","eight",
      "nine","ten","eleven","twelve",
      "thirteen")

#gr4<-ggplot(df,aes(x=z,y=y,label=lb)) + 
#  geom_point() +
#  geom_text(label=lb) +
#  ylab(TeX("$y_i$")) +
#  xlab(TeX("$\\hat{z}_i$"))+
#  stat_smooth(method = "lm", 
#              formula = y ~ x, 
#              geom = "smooth") 

length(y)







gr4<-ggplot(df,aes(x=z,y=y)) + 
  geom_point() +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$"))+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_correlation(mapping = use_label(c("R", "n","p","R2")),
                   label.x = "left",
                   label.y = "top") 
  #stat_poly_eq(aes(label =  paste(after_stat(eq.label), "*\" with \"*", 
  #                                after_stat(rr.label), "*\", \"*",
  #                               after_stat(p.value.label), "*\".\"",
  #                                sep = "")),
  #             formula = y~x, size = 3)




gr4









ggsave(file="P2pairs_all_t40.pdf",
       plot=gr1,
       width=5,height=5,
       units="in",dpi=500)


ggsave(file="P2res_vs_embedding_all_t40.pdf",
       plot=gr4,
       width=5,height=4,
       units="in",dpi=500)




data$model<-as.factor(data$model)

gr6<-ggplot(data,aes(x=model,y=score,fill=model)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  labs(x="model",y="learning score") +
  theme(legend.position="none")

gr6


ggsave(file="P2model_vs_score.pdf",
       plot=gr6,
       width=5,height=4,
       units="in",dpi=500)






library(gridExtra)
library(grid)
library(ggplotify)



gr5<-arrangeGrob(gr2,gr4,ncol=2)

plot(gr5)



ggsave(file="P2intro_plot_all_t40.pdf",gr5,
       width=7,height=3.5,
       units="in",dpi=500)



#####THE__END__FOR__REAL__DATA__SECTION__NEXT__COMES__NONPARAMETRIC########



library(ggplot2)
library(reshape2)
library(latex2exp)
library(np)

z_modnp<-z
print(z)
print(z_modnp)
df_modnp<-data.frame(y,z,z_modnp)
print(df_modnp)


modnp<-npreg(y~z_modnp,
             regtype="ll",
             gradients=TRUE,
             data=df_modnp)

summary(modnp)

plot(modnp,gradients=TRUE,plot.errors.method="asymptotic")

#yhatnp<-predict(modnp,newdata=df)



z_modnp<-seq(min(z)-0.005,max(z)+0.005,length.out=400)
print(z_modnp)

yhatnp<-predict(modnp,newdata=as.data.frame(z_modnp))

print(yhatnp)

dfnp2<-data.frame(z_modnp,yhatnp)



#dfnp<-data.frame(z,y,yhatnp)

#dfnpm<-reshape::melt(dfnp,id.vars = 'z')

grnp<-ggplot(df,aes(x=z,y=y)) +
  geom_point() +
  geom_line(data=dfnp2,aes(x=z_modnp,y=yhatnp),color="blue") +
  xlab(TeX("$\\hat{z}_i$")) +
  ylab(TeX("$y_i$"))

grnp

ggsave(file="P2res_vs_embedding_np.pdf",
       grnp,
       width=5,
       height=4,
       units="in")







library(gridExtra)
library(grid)

#grid.arrange(gr2,gr4,ncol=2)









gr5<-arrangeGrob(gr2,gr4,ncol=2)

plot(gr5)



ggsave(file="P2intro_plot.pdf",gr5,
       width=7,height=3.5,
       units="in",dpi=500)







