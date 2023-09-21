
setwd("E:/singlecall_immune/RCC/myeloid")
dir.create("network")
library(stringr)
a1<-read.csv('pyscenicresult/regulons.txt', sep=',',
             encoding='UTF-8',header = F)
a1<-t(a1)
a2<-c()

for(i in 1:nrow(a1)){
  if(str_detect(a1[i,1],"Regulon")){
    a2<-c(a2,i)
  }
}
a2<-c(a2,(nrow(a1)+1))

a3<-c()
for(i in 1:nrow(a1)){
  if(str_detect(a1[i,1],"gene2weight")){
    a3<-c(a3,i)
  }
}

tflist<-data.frame(TF=NA,Target=NA)

for (i in 1:(length(a2)-1)) {
  t1<-str_split(a1[a2[i],1],"'",simplify = T)[2]
  t2<-str_split(a1[a3[i],1],"'",simplify = T)[2]
  tflist<-rbind(tflist,c(t1,t2))
  f1<-a3[i]+1
  f2<-a2[i+1]-10
  for (j in f1:f2) {
    t3<-str_split(a1[j,1],"'",simplify = T)[2]
    tflist<-rbind(tflist,c(t1,t3))
  }
}
tflist<-tflist[-1,]

tflist$TF<-str_split(tflist$TF,"\\(",simplify = T)[,1]
tflist$TF<-str_split(tflist$TF," ",simplify = T)[,1]

write.csv(tflist,file = "network/myeloid_TF.csv")
save(tflist,file = "myeloid_TF.Rdata")

myeloid_clin<-clinical[which(clinical$Lineage=="Myeloid"),]
myeloid_clin$ICB_Response<-ifelse(myeloid_clin$ICB_Response=="ICB_SD","ICB_NR",
                                  ifelse(myeloid_clin$ICB_Response=="ICB_PD","ICB_NR",
                                         ifelse(myeloid_clin$ICB_Response=="ICB_NE","ICB_NR",myeloid_clin$ICB_Response)))

save(myeloid_clin,file = "髓系细胞临床信息.Rdata")
library(igraph)
g<-graph_from_data_frame(tflist,directed = F)


g1<-degree(g)

g1<-data.frame(g1)

g2<-closeness(g,mode = "all")
g2<-data.frame(g2)

g3<-betweenness(g)
g3<-data.frame(g3)

g4<-evcent(g,scale = F)$vector
g4<-data.frame(g4)

g5<-page.rank(g)$vector
g5<-data.frame(g5)

tuopu<-cbind(g1,g2,g3,g4,g5)

save(g,tuopu,tflist,file = "TF.Rdata")
write.csv(tuopu,file = "network/拓扑系数.csv")
write.table(tuopu,file = "network/拓扑系数.txt",sep = "\t",quote = F)

q1<-data.frame(gene1=rownames(g1),g1=g1$g1)
q1<-q1[order(q1$g1),]
q2<-data.frame(gene2=rownames(g2),g2=g2$g2)
q2<-q2[order(q2$g2),]
q3<-data.frame(gene3=rownames(g3),g3=g3$g3)
q3<-q3[order(q3$g3),]
q4<-data.frame(gene4=rownames(g4),g4=g4$g4)
q4<-q4[order(q4$g4),]
q5<-data.frame(gene5=rownames(g5),g5=g5$g5)
q5<-q5[order(q5$g5),]

qstatistic<-data.frame(gene=q1$gene1,V=NA,Q=NA)
for(i in 1:nrow(q1)){
  r5<-(which(q5$gene5==q1$gene1[i])-1)/nrow(q1)
  r4<-(which(q4$gene4==q1$gene1[i])-1)/nrow(q1)
  r3<-(which(q3$gene3==q1$gene1[i])-1)/nrow(q1)
  r2<-(which(q2$gene2==q1$gene1[i])-1)/nrow(q1)
  r1<-(i-1)/nrow(q1)
  v0=1
  v1=(-1)^0*v0*r5
  v2=(-1)^0*v1*r4+(-1)^1*v0/2*r4^2
  v3=(-1)^0*v2*r3+(-1)^1*v1/2*r3^2+(-1)^2*v0/6*r3^3
  v4=(-1)^0*v3*r2+(-1)^1*v2/2*r2^2+(-1)^2*v1/6*r2^3+
    (-1)^3*v0/24*r2^4
  v5=(-1)^0*v4*r1+(-1)^1*v3/2*r1^2+(-1)^2*v2/6*r1^3+
    (-1)^3*v1/24*r1^4+(-1)^4*v0/120*r1^5
  Q=factorial(5)*v5
  qstatistic[i,2]<-v5
  qstatistic[i,3]<-Q
}

#
library(limma)
annodata2<-myeloid_clin[which(myeloid_clin$ICB_Response!="NoICB"),]
hubdata2<-exprMat[,annodata2$NAME]
design<-model.matrix(~0+factor(annodata2$ICB_Response))
rownames(design)<-annodata2$NAME
colnames(design)<-levels(factor(annodata2$ICB_Response))
contrast.matrix<-makeContrasts(paste0(unique(annodata2$ICB_Response),collapse = "-"),
                               levels = design)

fit<-lmFit(hubdata2,design)
fit1<-contrasts.fit(fit,contrast.matrix)
fit1<-eBayes(fit1)
tempoutput_Responsed<-topTable(fit1,coef = 1,n=Inf)
deg_Responsed<-tempoutput_Responsed[which(tempoutput_Responsed$adj.P.Val<0.05 & 
                                      abs(tempoutput_Responsed$logFC)>0.5),]

qstatistic1<-qstatistic
qstatistic1<-qstatistic1[order(qstatistic1$Q,decreasing =T),]
qstatistic1<-qstatistic1[qstatistic1$Q<1,]
hubgene<-qstatistic1$gene[1:46]
dir.create("hubgene")

save(tempoutput_Responsed,deg_Responsed,
     file = "hubgene/差异分析结果.Rdata")
save(hubgene,file = "关键基因.Rdata")
save(qstatistic,qstatistic1,file = "Q统计量.Rdata")

library(pheatmap)
library(RColorBrewer)

annodata2<-myeloid_clin[which(myeloid_clin$ICB_Response!="NoICB"),]
annodata2<-annodata2[order(annodata2$ICB_Response),]
hubdata2<-exprMat[hubmye,annodata2$NAME]

color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata2)))
n[n>1]=1
n[n< -1]= -1

n1<-n[,annodata2$NAME[which(annodata2$ICB_Response=="ICB_PR")]]
p1<-pheatmap(n1,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
n2<-n[,annodata2$NAME[which(annodata2$ICB_Response=="ICB_NR")]]
p2<-pheatmap(n2,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
pr<-colnames(n1)[p1$tree_col[["order"]]]
nr<-colnames(n2)[p2$tree_col[["order"]]]
annodata2<-annodata2[c(nr,pr),]
annodata<-data.frame(row.names = annodata2$NAME,
                     ICB_Response=annodata2$ICB_Response)
pheatmap(n[,c(nr,pr)],show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,
         legend_breaks=c(-1,0,1),
         legend_labels=c("-1","0","1"),
         main = "",filename = "hubgene/46gene_ICB_Response.pdf",
         width = 6,height = 6)

hubdata3<-exprMat[hubgene,annodata2$NAME]
n=t(scale(t(hubdata3)))
n[n>1]=1
n[n< -1]= -1
n1<-n[,annodata2$NAME[which(annodata2$ICB_Response=="ICB_PR")]]
p1<-pheatmap(n1,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
n2<-n[,annodata2$NAME[which(annodata2$ICB_Response=="ICB_NR")]]
p2<-pheatmap(n2,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
pr<-colnames(n1)[p1$tree_col[["order"]]]
nr<-colnames(n2)[p2$tree_col[["order"]]]
annodata2<-annodata2[c(nr,pr),]
annodata<-data.frame(row.names = annodata2$NAME,
                     ICB_Response=annodata2$ICB_Response)
pheatmap(n[,c(nr,pr)],show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,
         legend_breaks=c(-1,0,1),
         legend_labels=c("-1","0","1"),
         main = "",filename = "hubgene/10gene_ICB_Response.pdf",
         width = 6,height = 6)


################################
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

dir.create("cluster")
scRNAsubmy<-readRDS("scRNAsubmy.rds")
scRNAsubmy<- FindVariableFeatures(scRNAsubmy, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsubmy)
scRNAsubmy <- ScaleData(scRNAsubmy, features = scale.genes)
scRNAsubmy <- RunPCA(scRNAsubmy, features = VariableFeatures(scRNAsubmy))
ElbowPlot(scRNAsubmy, ndims=30, reduction="pca")
pc.num=1:25

scRNAsubmy<- FindNeighbors(scRNAsubmy, dims = pc.num) 
scRNAsubmy <- FindClusters(scRNAsubmy, resolution = 0.9)
table(scRNAsubmy@meta.data$seurat_clusters)
scRNAsubmy = RunTSNE(scRNAsubmy, dims = pc.num)
scRNAsubmy <- RunUMAP(scRNAsubmy, dims = pc.num)
p1 = DimPlot(scRNAsubmy, group.by="celltype", label=F, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsubmy, group.by="celltype", label=F, label.size=5, reduction='umap')
plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("cluster/tSNE_UMAP细胞再聚类.pdf", plot = plotc, width = 10, height = 5)
saveRDS(scRNAsubmy,'scRNAsubmy.rds')

setwd("E:/singlecall_immune/RCC/myeloid")
annodata<-myeloid_clin[which(myeloid_clin$ICB_Response!="NoICB"),]
hubdata<-exprMat[hubgene,annodata$NAME]
library(ConsensusClusterPlus)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 5,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "E:/singlecall_immune/RCC/consensus/myeloid/",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "png")

res_tree<-results[[2]]$consensusClass
res_tree<-data.frame(res_tree)
save(results,res_tree,file = "E:/singlecall_immune/RCC/consensus/myeloid/result.Rdata")

library(pheatmap)
library(RColorBrewer)
res_tree<-data.frame(sample=row.names(res_tree),cluster=res_tree$res_tree)
res_tree<-res_tree[order(res_tree$cluster),]
hubdata<-hubdata[,res_tree$sample]
annodata<-data.frame(row.names = res_tree$sample,
                     cluster=res_tree$cluster)
annodata$cluster<-as.character(annodata$cluster)
annodata1<-myeloid_clin[row.names(annodata),]
annodata$ICB_Responsed<-annodata1$ICB_Response

color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,
         legend_breaks=c(-1,0,1),
         legend_labels=c("-1","0","1"),
         main = "",
         filename = "E:/singlecall_immune/RCC/consensus/myeloid/10hubgene热图.pdf",
         width = 6,height = 6)
