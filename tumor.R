setwd("E:/singlecall_immune/RCC")

library(data.table)
rcc_data<-fread("data/ccRCC_scRNASeq_NormalizedCounts.txt",sep = "\t",header = TRUE,na.strings = "NA")
row.names(rcc_data)<-rcc_data[,1]
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)

mydata<-readRDS("mydata.rds")

library(celldex)
library(SingleR)
load('testdata.Rdata')
load("clusters.Rdata")
top<-read.csv("top_p.adj.csv")
top_gene<-unique(top$gene)
refdata<-celldex::HumanPrimaryCellAtlasData()
save(refdata,file = "HumanPrimaryCell_refdata.Rdata")
testdata1<-testdata[top_gene,]


cellpred<-SingleR(test = testdata1,ref = refdata,labels = refdata$label.main,
                  method = "cluster",clusters = clusters,
                  assay.type.test = "logcounts",assay.type.ref = "logcounts")
save(cellpred,file = "cellpred.Rdata")



setwd("E:/singlecall_immune/RCC/tumor")

data<-cellexpr[which(apply(cellexpr,1,
                          function(x){return(sum(x!=0))})>ncol(cellexpr)*0.01),]
dim(data)


library(stringr)
a1<-read.table("pyscenic_result/regulons.txt")
a1<-t(a1)
a2<-c()

for(i in 1:nrow(a1)){
  if(str_detect(a1[i,1],"Regulon")){
    a2<-c(a2,i)
  }
}
a2<-c(a2,nrow(a1))

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
  t<-(a2[i+1]-a2[i]-12-2)%/%3
  for (j in 1:t) {
    t3<-a1[a2[i]+1+3*j,1]
    tflist<-rbind(tflist,c(t1,t3))
  }
}
tflist<-tflist[-1,]

#tflist$TF<-str_split(tflist$TF,"\\(",simplify = T)[,1]
write.csv(tflist,file = "network/tumor_TF.csv")

library(igraph)
g<-graph_from_data_frame(tflist,directed = F)

save(g,tflist,file = "tumor/TF.Rdata")


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
colnames(tuopu)<-c("degree","closeness","betweenness","eigenvalue","PageRank")
save(g,tuopu,tflist,file = "tumor/TF.Rdata")
write.csv(tuopu,file = "tumor/拓扑系数.csv")
write.table(tuopu,file = "tumor/拓扑系数.txt",sep = "\t",quote = F)



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

qstatistic1<-qstatistic
qstatistic1<-qstatistic1[order(qstatistic1$Q,decreasing = T),]
qstatistic1<-qstatistic1[qstatistic1$Q<1,]
save(qstatistic,qstatistic1,file = "Q统计量.Rdata")
hubgene<-qstatistic1$gene[1:46]
save(hubgene,file = "关键基因.Rdata")



setwd("E:/singlecall_immune/RCC/tumor")


library(limma)
annodata2<-tumor_clin[which(tumor_clin$ICB_Response!="NoICB"),]
hubdata2<-cellexpr[,annodata2$NAME]
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

################################
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

require(ggplot2)
require(SCENIC)
#require(SCopeLoomR)
require(AUCell)

regulonAUC <- importAUCfromText('pyscenic_result/auc_mtx.tsv')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
scRNAsubcell<-readRDS("tumor/scRNAsub.rds")
cellInfo = scRNAsubcell@meta.data
cellexpr<-scRNAsubcell@assays$RNA@scale.data


regulon<-regulonAUC@assays@data@listData$AUC
save(cellexpr,cellInfo,file = "表达谱及细胞信息.Rdata")
save(regulon,regulonAUC,file="SCENIC_AUC.Rdata")

scRNAsubcell <- FindVariableFeatures(scRNAsubcell, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsubcell)
scRNAsubcell <- ScaleData(scRNAsubcell, features = scale.genes)
scRNAsubcell <- RunPCA(scRNAsubcell, features = VariableFeatures(scRNAsubcell))
ElbowPlot(scRNAsubcell, ndims=30, reduction="pca")
pc.num=1:20
scRNAsubcell <- FindNeighbors(scRNAsubcell, dims = pc.num) 
scRNAsubcell <- FindClusters(scRNAsubcell, resolution = 0.9)
table(scRNAsubcell@meta.data$seurat_clusters)
scRNAsubcell = RunTSNE(scRNAsubcell, dims = pc.num)
scRNAsubcell <- RunUMAP(scRNAsubcell, dims = pc.num)

p1 = DimPlot(scRNAsubcell, group.by="celltype", label=F, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsubcell, group.by="celltype", label=F, label.size=5, reduction='umap')
plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("tumor/tSNE_UMAP细胞再聚类.pdf", plot = plotc, width = 10, height = 5)
saveRDS(scRNAsubcell,'tumor/scRNAsubcell.rds')


annodata1<-data.frame(row.names = rownames(cellInfo),cell=rownames(cellInfo),celltype=cellInfo$celltype)
annodata1<-annodata1[order(annodata1$celltype),]
annodata1<-data.frame(row.names = rownames(annodata1),celltype=annodata1$celltype)
hubdata<-cellexpr[hubtum,]
hubdata<-hubdata[,rownames(annodata1)]
save(hubgene,hubdata,file = "关键基因及表达谱.Rdata")

library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,color=color,
         legend_breaks=c(-1,0,1),
         legend_labels=c("-1","0","1"),
         main = "",filename = "46gene热图.pdf",
         width = 6,height = 6)

hubdata1<-sapply(split(rownames(annodata1), annodata1$celltype),
                 function(cells) rowMeans(hubdata[,cells]))

color<- colorRampPalette(c('#87CEFA','white','red'))(100)

pheatmap(hubdata1,
         show_rownames = T,show_colnames = T, 
         cluster_rows = T,cluster_cols = T,
         color=color,
         main = "",filename = "hubgene/46gene热图_ave.pdf",
         width = 4,height = 8)

#,filename = "46gene热图.pdf",width = 6,height = 6


clinical<-read.csv("data/Clinical dataset.csv")
rownames(clinical)<-clinical$NAME
clinical<-clinical[,-8]
save(clinical,file = "clinical.Rdata")

table(clinical$ICB_Response)
table(clinical$ICB_Exposed)
tumor_clin<-clinical[row.names(cellInfo),]

tumor_clin$ICB_Response<-ifelse(tumor_clin$ICB_Response=="ICB_SD","ICB_NR",
                                ifelse(tumor_clin$ICB_Response=="ICB_PD","ICB_NR",
                                       ifelse(tumor_clin$ICB_Response=="ICB_NE","ICB_NR",tumor_clin$ICB_Response)))

save(tumor_clin,file = "肿瘤细胞临床信息.Rdata")

annodata2<-tumor_clin[which(tumor_clin$ICB_Response!="NoICB"),]
annodata2<-annodata2[order(annodata2$ICB_Response),]
hubdata2<-hubdata[,annodata2$NAME]

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

hubdata3<-cellexpr[hubgene,annodata2$NAME]
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



hubdata1<-sapply(split(rownames(annodata2), annodata2$ICB_Response),
                 function(cells) rowMeans(hubdata2[,cells]))

color<- colorRampPalette(c('#87CEFA','white','red'))(100)

pheatmap(hubdata1,
         show_rownames = T,show_colnames = T, 
         cluster_rows = T,cluster_cols = T,
         color=color,
         main = "",filename = "hubgene/46gene_ICB_Response_ave.pdf",
         width = 4,height = 8)
##########
annodata3<-tumor_clin[order(tumor_clin$ICB_Exposed),]
hubdata3<-hubdata[,annodata3$NAME]

n=t(scale(t(hubdata3)))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = data.frame(row.names = annodata3$NAME,
                                     ICB_Exposed=annodata3$ICB_Exposed),
         color=color,
         legend_breaks=c(-1,0,1),
         legend_labels=c("-1","0","1"),
         main = "",filename = "46gene_ICB_Exposed.pdf",
         width = 6,height = 6)

hubdata1<-sapply(split(rownames(tumor_clin), tumor_clin$ICB_Exposed),
                 function(cells) rowMeans(hubdata[,cells]))

color<- colorRampPalette(c('#87CEFA','white','red'))(100)

pheatmap(hubdata1,
         show_rownames = T,show_colnames = T, 
         cluster_rows = T,cluster_cols = T,
         color=color,
         main = "",filename = "hubgene/46gene_ICB_Exposed_ave.pdf",
         width = 4,height = 8)


dir.create("tumor/enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
ego_ALL <- enrichGO(gene          = hubgene,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'enrich/enrichGO.csv') 

ego_CC <- enrichGO(gene          = hubgene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_cc<-data.frame(ego_CC)
write.csv(ego_cc,'enrich/enrichGOcc.csv') 

ego_MF <- enrichGO(gene          = hubgene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_mf<-data.frame(ego_MF)
write.csv(ego_mf,'enrich/enrichGOmf.csv') 

ego_BP <- enrichGO(gene          = hubgene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp<-data.frame(ego_BP)
write.csv(ego_bp,'enrich/enrichGObp.csv') 

save(ego_all,ego_ALL,ego_bp,ego_BP,ego_cc,
     ego_CC,ego_mf,ego_MF,file = "enrich/enrichGO_result.Rdata")
#ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
#ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
#ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

#p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
#p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
#p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
#plotc <- p_BP/p_CC/p_MF
#ggsave('enrich/enrichGO.pdf', plotc, width = 12,height = 10)

genelist <- bitr(hubgene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
eKEGG<-data.frame(ekegg)
save(ekegg,eKEGG,file = "enrich/enrichKEGG_result.Rdata")

p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave("enrich/enrichKEGG.pdf", plot = plotc, width = 12, height = 10)


setwd("E:/singlecall_immune/RCC/tumor")
library(GSVA)
library(GSVAdata)
geneset <- getGmt("GSVA/h.all.v7.5.1.symbols.gmt")

expr_gsva<-cellexpr[which(rowSums(cellexpr)>10),]

va_hub<-gsva(as.matrix(hubdata),geneset,mx.diff=FALSE, verbose=TRUE)
va_all<-gsva(expr_gsva,geneset,mx.diff=FALSE, verbose=TRUE)

va_hub_score<-sapply(split(rownames(cellInfo), cellInfo$celltype),
                     function(cells) rowMeans(va_hub[,cells]))
library(stringr)
row.names(va_hub_score)<-str_sub(row.names(va_hub_score),10,100)

va_all_score<-sapply(split(rownames(cellInfo), cellInfo$celltype),
                     function(cells) rowMeans(va_all[,cells]))
row.names(va_all_score)<-str_sub(row.names(va_all_score),10,100)
save(va_all,va_all_score,file = "GSVA/GSVA_result_allgene.Rdata")
save(va_hub,va_hub_score,file = "GSVA/GSVA_result_hubgene.Rdata")

va_hub_score_scaled <- t(scale(t(va_hub_score), center = T, scale=T))
va_hub_score_scaled[va_hub_score_scaled>1]<- 1
va_hub_score_scaled[va_hub_score_scaled< -1]<- -1
pheatmap::pheatmap(va_hub_score_scaled, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color = NA) 
                   #treeheight_row=10, treeheight_col=10, border_color=NA)

va_all_score_scaled <- t(scale(t(va_all_score), center = T, scale=T))
va_all_score_scaled[va_all_score_scaled>1]<- 1
va_all_score_scaled[va_all_score_scaled< -1]<- -1
pheatmap::pheatmap(va_all_score_scaled, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color = NA) 

va_hub_ICB_exposed<-sapply(split(rownames(tumor_clin), tumor_clin$ICB_Exposed),
                     function(cells) rowMeans(va_hub[,cells]))
rownames(va_hub_ICB_exposed)<-str_sub(rownames(va_hub_ICB_exposed),10,100)
va_hub_ICB_exposed<-t(scale(t(va_hub_ICB_exposed),center = T,scale = T))
va_hub_ICB_exposed[va_hub_ICB_exposed>1]<- 1
va_hub_ICB_exposed[va_hub_ICB_exposed< -1]<- -1
pheatmap::pheatmap(va_hub_ICB_exposed, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color = NA,filename = "GSVA/GSVA_ICB_Exposed.pdf",
                   height = 6,width = 5) 

va_ICB<-va_hub[,annodata2$NAME]
va_hub_ICB_responsed<-sapply(split(rownames(annodata2), annodata2$ICB_Response),
                           function(cells) rowMeans(va_ICB[,cells]))
rownames(va_hub_ICB_responsed)<-str_sub(rownames(va_hub_ICB_responsed),10,100)
va_hub_ICB_responsed<-t(scale(t(va_hub_ICB_responsed),center = T,scale = T))
va_hub_ICB_responsed[va_hub_ICB_responsed>1]<- 1
va_hub_ICB_responsed[va_hub_ICB_responsed< -1]<- -1
pheatmap::pheatmap(va_hub_ICB_responsed, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color = NA,filename = "GSVA/GSVA_ICB_Responsed.pdf",
                   height = 6,width = 5) 

##################

library(AUCell)
dir.create("TF_AUC")
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                    function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled[regulonActivity_byCellType_Scaled>1]<-1
regulonActivity_byCellType_Scaled[regulonActivity_byCellType_Scaled< -1]<- -1

rownames(regulonActivity_byCellType_Scaled)<-str_split(rownames(regulonActivity_byCellType_Scaled)," ",simplify = T)[,1]
t1<-intersect(row.names(regulonActivity_byCellType_Scaled),hubgene)
n<-regulonActivity_byCellType_Scaled[t1,]
rownames(n)<-paste0(row.names(n),"(+)")
pheatmap::pheatmap(n, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color=NA,filename = "TF_AUC/AUC_cell.pdf",
                   width = 3,height = 6)

regulonActivity_byCellType <- sapply(split(rownames(tumor_clin), tumor_clin$ICB_Exposed),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled[regulonActivity_byCellType_Scaled>1]<-1
regulonActivity_byCellType_Scaled[regulonActivity_byCellType_Scaled< -1]<- -1

rownames(regulonActivity_byCellType_Scaled)<-str_split(rownames(regulonActivity_byCellType_Scaled)," ",simplify = T)[,1]
t1<-intersect(row.names(regulonActivity_byCellType_Scaled),hubgene)
n<-regulonActivity_byCellType_Scaled[t1,]
rownames(n)<-paste0(row.names(n),"(+)")
pheatmap::pheatmap(n, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color=NA,filename = "TF_AUC/AUC_ICB_exposed.pdf",
                   width = 3,height = 6)

regulon_ICB<-regulon[,tumor_clin$NAME[which(tumor_clin$ICB_Response!="NoICB")]]
regulonActivity <- sapply(split(rownames(annodata2), annodata2$ICB_Response),
                                     function(cells) rowMeans(regulon_ICB[,cells]))
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
regulonActivity_Scaled[regulonActivity_Scaled>1]<-1
regulonActivity_Scaled[regulonActivity_Scaled< -1]<- -1

rownames(regulonActivity_Scaled)<-str_split(rownames(regulonActivity_Scaled)," ",simplify = T)[,1]
t1<-intersect(row.names(regulonActivity_Scaled),hubgene)
n<-regulonActivity_Scaled[t1,]
rownames(n)<-paste0(row.names(n),"(+)")
pheatmap::pheatmap(n, #fontsize_row=3, 
                   cluster_rows = T,cluster_cols = F,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   border_color=NA,filename = "TF_AUC/AUC_ICB_responsed.pdf",
                   width = 3,height = 6)

######
library(monocle)
scRNAsubcell<-readRDS("scRNAsubcell.rds")
tumor_clin<-tumor_clin[rownames(scRNAsubcell@meta.data),]
scRNAsubcell@meta.data$ICB_Exposed<-tumor_clin$ICB_Exposed
scRNAsubcell@meta.data$ICB_Reponsed<-tumor_clin$ICB_Response
saveRDS(scRNAsubcell,"pseudotime/scRNAsubcell.rds")

mycds<-readRDS("pseudotime/mycds.rds")

diff_test <- differentialGeneTest(mycds[hubgene,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, pval < 0.05))
p1<-plot_pseudotime_heatmap(mycds[hubgene,], cluster_rows = T,
                        num_clusters = 4,
                        show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap.pdf", plot = p1, width = 5, height = 6)

immune_checkpoint<-c("HAVCR2", "LAG3", "PDCD1", "TIGIT", "CD27", "CTLA4", 
                     "TNFRSF9", "CD226", "CD244", "TMIGD2", "CD7", "KLRG1")
p1<-plot_pseudotime_heatmap(mycds[immune_checkpoint,], cluster_rows = T,
                            num_clusters = 3,
                            show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_immunecheckpoint.pdf", plot = p1, width = 5, height = 6)


s.genes <- c("STAT1","MSC","SOX4","YBX1")
p1 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
ggsave("pseudotime/genes_celltype.pdf", plot = p1, width = 4, height = 6)

p2 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "ICB_Exposed")
ggsave("pseudotime/genes_ICB_Exposed.pdf", plot = p2, width = 4, height = 6)

p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "ICB_Reponsed")
ggsave("pseudotime/genes_ICB_Responsed.pdf", plot = p3, width = 4, height = 6)

expr_var<-cellexpr[rownames(tuopu),]
net_var<-sapply(split(rownames(tumor_clin),tumor_clin$ICB_Exposed),
                function(cells) rowMeans(expr_var[,cells]))

expr_var_response<-expr_var[,tumor_clin$NAME[which(tumor_clin$ICB_Response!="NoICB")]]
tumor_clin_ICB<-tumor_clin[which(tumor_clin$ICB_Response!="NoICB"),]

net_var_response<-sapply(split(rownames(tumor_clin_ICB),tumor_clin_ICB$ICB_Response),
                         function(cells) rowMeans(expr_var_response[,cells]))
net_expr<-cbind(net_var_response,net_var)
save(net_expr,expr_var,expr_var_response,
     file = "network/节点表达值.Rdata")
write.csv(net_expr,file = "network/net_expr1.csv") 
net_expr[net_expr>1.5]<-1.5
net_expr[net_expr< -1.5]<- -1.5
write.csv(net_expr,file = "network/net_expr.csv") 
write.csv(hubgene,file = "network/hubgene.csv")
write.table(hubgene,file = "network/hubgene.txt",row.names = F,col.names = F)



library(limma)
design<-model.matrix(~0+factor(annodata2$ICB_Response))
rownames(design)<-annodata2$NAME
colnames(design)<-levels(factor(annodata2$ICB_Response))
contrast.matrix<-makeContrasts(paste0(unique(annodata2$ICB_Response),collapse = "-"),
                               levels = design)
fit<-lmFit(hubdata2,design)
fit1<-contrasts.fit(fit,contrast.matrix)
fit1<-eBayes(fit1)
tempoutput_Responsed<-topTable(fit1,coef = 1,n=Inf)
deg_Responsed<-tempoutput_Responsed[tempoutput_Responsed$adj.P.Val<0.05 & 
                                      abs(tempoutput_Responsed$logFC)>0.5,]



design<-model.matrix(~0+factor(tumor_clin$ICB_Exposed))
rownames(design)<-tumor_clin$NAME
colnames(design)<-levels(factor(tumor_clin$ICB_Exposed))
contrast.matrix<-makeContrasts(paste0(unique(tumor_clin$ICB_Exposed),collapse = "-"),
                               levels = design)
fit<-lmFit(hubdata,design)
fit1<-contrasts.fit(fit,contrast.matrix)
fit1<-eBayes(fit1)
tempoutput_Exposed<-topTable(fit1,coef = 1,n=Inf)

deg_Exposed<-tempoutput_Exposed[tempoutput_Exposed$adj.P.Val<0.05 & 
                                      abs(tempoutput_Exposed$logFC)>0.5,]


save(tempoutput_Exposed,tempoutput_Responsed,
     deg_Exposed,deg_Responsed,
     file = "hubgene/差异分析结果.Rdata")

write.table(row.names(tempoutput_Responsed)[which(tempoutput_Responsed$adj.P.Val<0.05)],
            file = "hubgene/tempoutput_Responsed.txt",
            sep = "\t",quote = F,row.names = F,col.names = F)

write.table(row.names(tempoutput_Exposed)[which(tempoutput_Exposed$adj.P.Val<0.05)],
            file = "hubgene/tempoutput_Exposed.txt",
            sep = "\t",quote = F,row.names = F,col.names = F)

write.table(row.names(deg_Responsed),
            file = "hubgene/deg_Responsed.txt",
            sep = "\t",quote = F,row.names = F,col.names = F)

write.table(row.names(deg_Exposed),
            file = "hubgene/deg_Exposed.txt",
            sep = "\t",quote = F,row.names = F,col.names = F)


for(i in 1:4667){
  r7<-(which(q7$gene7==q1$gene1[i])-1)/nrow(q1)
  r6<-(which(q6$gene6==q1$gene1[i])-1)/nrow(q1)
  r5<-(which(q5$gene5==q1$gene1[i])-1)/nrow(q1)
  r4<-(which(q4$gene4==q1$gene1[i])-1)/nrow(q1)
  r3<-(which(q3$gene3==q1$gene1[i])-1)/nrow(q1)
  r2<-(which(q2$gene2==q1$gene1[i])-1)/nrow(q1)
  r1<-(i-1)/nrow(q1)
  v0=1
  v1=(-1)^0*v0*r7
  v2=(-1)^0*v1*r6+(-1)^1*v0/2*r6^2
  v3=(-1)^0*v2*r5+(-1)^1*v1/2*r5^2+(-1)^2*v0/6*r5^3
  v4=(-1)^0*v3*r4+(-1)^1*v2/2*r4^2+(-1)^2*v1/6*r4^3+
    (-1)^3*v0/24*r4^4
  v5=(-1)^0*v4*r3+(-1)^1*v3/2*r3^2+(-1)^2*v2/6*r3^3+
    (-1)^3*v1/24*r3^4+(-1)^4*v0/120*r3^5
  v6<-(-1)^0*v5*r2+(-1)^1*v4/2*r2^2+(-1)^2*v3/6*r2^3+
    (-1)^3*v2/24*r2^4+(-1)^4*v1/120*r2^5+(-1)^5*v0/720*r2^6
  v7<-(-1)^0*v6*r1+(-1)^1*v5/2*r1^2+(-1)^2*v4/6*r1^3+
    (-1)^3*v3/24*r1^4+(-1)^4*v2/120*r1^5+
    (-1)^5*v1/720*r1^6+(-1)^6*v0/5040*r1^7
  Q=factorial(7)*v7
  qstatistic[i,2]<-v7
  qstatistic[i,3]<-Q
}
####################################
setwd("E:/singlecall_immune/RCC/lymphoid/")
hublym<-hubgene
deglym<-tempoutput_Responsed[hubgene,]
setwd("E:/singlecall_immune/RCC/myeloid")
hubmye<-hubgene
degmye<-tempoutput_Responsed[hubgene,]
setwd("E:/singlecall_immune/RCC/tumor/")
hubtum<-hubgene
degtum<-tempoutput_Responsed[hubtum,]
#Reduce(intersect,list(hublym,hubmye,hubtum))
#intersect(intersect(hublym,hubmye),hubtum)

hubgene<-c(rownames(deglym)[which(deglym$adj.P.Val<0.05&abs(deglym$logFC)>0.5)],
           rownames(degmye)[which(degmye$adj.P.Val<0.05&abs(degmye$logFC)>0.5)],
           rownames(degtum)[which(degtum$adj.P.Val<0.05&abs(degtum$logFC)>0.5)])
hubgene<-hubgene[!duplicated(hubgene)]
write.table(hubgene,file = "E:/singlecall_immune/RCC/consensus/汇总的基因.txt",
            row.names = F,quote = F,col.names = F)
save(hubgene,hublym,hubmye,hubtum,deglym,degmye,degtum,
     file = "E:/singlecall_immune/RCC/consensus/hubgene汇总.Rdata")



