library(ggplot2)
library(ggpubr)
source("Algroithm.r")
####import datasets
mRNA=read.csv("LUAD_mRNA.csv",row.names=1)
miRNA=read.csv("LUAD_miRNA.csv",row.names=1)
methy=read.csv("LUAD_methy.csv",row.names=1)

dim=c(67,67,66)
#tuning parameters by BIC
lam=seq(0.01,1,0.05)
LUAD=list(mRNA,miRNA,methy)
BIC_list=list()
bic_mat=matrix(0,length(lam),length(lam))
for (i in 1:length(lam)) {
  for (j in 1:length(lam)) {
    for (k in 1:length(lam)) {
      res=Hier(LUAD,lambda1 = lam[i],lambda2 = lam[j],lambda3 = lam[k],dim = dim)
    }
    bic_mat[j,k]=res[[5]]
  }
  BIC_list[[i]]=bic_mat
}

after_tuning=Hier(LUAD,lambda1 = lam[3],lambda2 = lam[3],lambda3 = lam[3],dim = dim)

XZ1_sps=sign(abs(after_tuning[[1]]))
XZ2_sps=sign(abs(after_tuning[[2]]))
XZ3_sps=sign(abs(after_tuning[[3]]))

###graph 1
jpeg(file="Hier_1.jpeg",width=3900,height=3900,res=500)
diag(XZ1_sps)=0
set.seed(2)
Hier_1=graph_from_adjacency_matrix(sign(abs(XZ1_sps)),mode="undirected")
E(Hier_1)$color='blue'
l=layout_with_kk(Hier_1) 
plot(Hier_1,vertex.size=2,vertex.label.cex=1,
     vertex.label.dist=1.2,edge.width=1.6,layout=l)
dev.off()
###graph 2
jpeg(file="Hier_2.jpeg",width=3900,height=3900,res=500)
diag(XZ2_sps)=0
Hier_1=graph_from_adjacency_matrix(sign(abs(XZ2_sps)),mode="undirected")
E(Hier2)$color='green'
plot(Hier2,vertex.size=2,vertex.label.cex=1,
     vertex.label.dist=1.2,edge.width=1.6,layout=l)
dev.off()
###graph 3
jpeg(file="real_3.jpeg",width=3900,height=3900,res=500)
diag(XZ3_sps)=0
Hier3=graph_from_adjacency_matrix(sign(abs(XZ3_sps)),mode="undirected")
E(Hier3)$color='orange'
plot(Hier3,vertex.size=2,vertex.label.cex=1,
     vertex.label.dist=1.2,edge.width=1.6,layout=l)
dev.off()
