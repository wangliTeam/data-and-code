###############

setwd("E:/singlecall_immune/RCC/lymphoid/")

data<-exprMatlym[which(apply(exprMatlym,1,
                             function(x){return(sum(x!=0))})>ncol(exprMatlym)*0.01),]
dim(data)
save(data,file="data.Rdata")
write.csv(data,file="data.csv")

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
tflist<-na.omit(tflist)
write.csv(tflist,file = "network/lymphoid_TF.csv")
save(tflist,file = "lymphoid_TF.Rdata")


lymphoid_clin<-clinical[which(clinical$Lineage=="Lymphoid"),]
lymphoid_clin$ICB_Response<-ifelse(lymphoid_clin$ICB_Response=="ICB_SD","ICB_NR",
                                   ifelse(lymphoid_clin$ICB_Response=="ICB_PD","ICB_NR",
                                         ifelse(lymphoid_clin$ICB_Response=="ICB_NE","ICB_NR",
                                                lymphoid_clin$ICB_Response)))
save(lymphoid_clin,file = "淋巴细胞临床信息.Rdata")


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


library(limma)
setwd("/pub6/temp/wangli/yxz/RCC/lymphoid")
load("lymclincal.Rdata")
load("exprMatlym.Rdata")

annodata2<-lymphoid_clin[which(lymphoid_clin$ICB_Response!="NoICB"),]
hubdata2<-exprMatlym[,annodata2$NAME]
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
save(tempoutput_Responsed,deg_Responsed,
     file = "hubgene/差异分析结果.Rdata")


load("差异分析结果.Rdata")

qstatistic1<-qstatistic
qstatistic1<-qstatistic1[order(qstatistic1$Q,decreasing =T),]
qstatistic1<-qstatistic1[qstatistic1$Q<1,]
hubgene<-qstatistic1$gene[1:32]
dir.create("hubgene")
save(hubgene,file = "关键基因.Rdata")
save(qstatistic,qstatistic1,file = "Q统计量.Rdata")



library(stringr)
library(limma)
cd8clin<-lymphoid_clin[which(str_detect(lymphoid_clin$FinalCellType,"T cell")),]
cd8clin<-cd8clin[which(cd8clin$FinalCellType!="MitoHigh CD8+ T cell"),]
annodata<-cd8clin[,c(1,3,7)]
annodata$cell<-ifelse(annodata$FinalCellType=="41BB-Lo CD8+ T cell",
                      "progenitor","terminally")

cd8data<-exprMatlym[,annodata$NAME]

design<-model.matrix(~0+factor(annodata$cell))
rownames(design)<-annodata$NAME
colnames(design)<-levels(factor(annodata$cell))
contrast.matrix<-makeContrasts(paste0(unique(annodata$cell),collapse = "-"),
                               levels = design)
fit<-lmFit(cd8data,design)
fit1<-contrasts.fit(fit,contrast.matrix)
fit1<-eBayes(fit1)
tempoutput<-topTable(fit1,coef = 1,n=Inf)
degcd8<-tempoutput[which(tempoutput$adj.P.Val<0.05 & abs(tempoutput$logFC)>0.5),]
save(tempoutput,degcd8,annodata,cd8clin,cd8data,
     file = "CD8Tcell差异分析.Rdata")



library(pheatmap)
library(RColorBrewer)

annodata2<-lymphoid_clin[which(lymphoid_clin$ICB_Response!="NoICB"),]
annodata2<-annodata2[order(annodata2$ICB_Response),]
hubdata2<-data[hublym,annodata2$NAME]

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

hubdata3<-data[hubgene,annodata2$NAME]
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

setwd("/pub6/temp/wangli/yxz/RCC/lymphoid/")
scRNAsublym<-readRDS("scRNAsublym.rds")
scRNAsublym <- FindVariableFeatures(scRNAsublym, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsublym)
scRNAsublym <- ScaleData(scRNAsublym, features = scale.genes)
scRNAsublym <- RunPCA(scRNAsublym, features = VariableFeatures(scRNAsublym))
p1<-ElbowPlot(scRNAsublym, ndims=30, reduction="pca")
dir.create("cluster")
ggsave("cluster/pca.pdf",p1,height=5,width=5)
pc.num=1:25
scRNAsublym<- FindNeighbors(scRNAsublym, dims = pc.num)

scRNAsublym <- FindClusters(scRNAsublym, resolution = 0.9)
table(scRNAsublym@meta.data$seurat_clusters)
scRNAsublym = RunTSNE(scRNAsublym, dims = pc.num)
scRNAsublym <- RunUMAP(scRNAsublym, dims = pc.num)
p1 = DimPlot(scRNAsublym, group.by="celltype", label=F, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsublym, group.by="celltype", label=F, label.size=5, reduction='umap')
plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("cluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
saveRDS(scRNAsublym,'scRNAsublym.rds')

setwd("E:/singlecall_immune/RCC/lymphoid/")

annodata<-lymphoid_clin[which(lymphoid_clin$ICB_Response!="NoICB"),]
hubdata<-data[hubgene,annodata$NAME]
library(ConsensusClusterPlus)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 5,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "E:/singlecall_immune/RCC/consensus/lymphoid/",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "png")
