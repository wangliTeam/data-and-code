setwd("E:/singlecall_immune/RCC/")
immunegene<-read.table("E:/singlecall_immune/RCC/bulkdata/GeneList.txt",header = T,sep = "\t")
save(immunegene,file = "bulkdata/immunegene.Rdata")
load("E:/singlecall_immune/RCC/tumor/Q统计量.Rdata")
a<-intersect(immunegene$Symbol,qstatistic1$gene)
q_tumor<-qstatistic1[which(qstatistic1$gene%in%a),]
load("E:/singlecall_immune/RCC/myeloid/Q统计量.Rdata")
a<-intersect(immunegene$Symbol,qstatistic1$gene)
q_myeloid<-qstatistic1[which(qstatistic1$gene%in%a),]
load("E:/singlecall_immune/RCC/lymphoid/Q统计量.Rdata")
a<-intersect(immunegene$Symbol,qstatistic1$gene)
q_lymphoid<-qstatistic1[which(qstatistic1$gene%in%a),]
three<-Reduce(intersect,list(q_lymphoid$gene,q_myeloid$gene,q_tumor$gene))

hubtum<-q_tumor$gene[1:13]
hubmye<-q_myeloid$gene[1:17]
hublym<-q_lymphoid$gene[1:12]
hubgene<-c(hublym,hubmye,hubtum)
hubgene<-hubgene[!duplicated(hubgene)]
save(hubgene,hublym,hubmye,hubtum,q_lymphoid,q_myeloid,q_tumor,
     file = "E:/singlecall_immune/RCC/consensus/汇总4.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
p1<-ggplot(q_tumor[1:14,],aes(gene,Q,fill=gene))+
  geom_bar(stat="identity")+xlab("")+
  ggtitle("Tumor")+ylab("Q statistic")+
  theme(axis.ticks.length=unit(0.5,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill="none")+
  coord_flip()+theme_bw()
ggsave("tumor/tumor_Q.pdf",p1,width = 4,height = 3)
p2<-ggplot(q_myeloid[1:17,],aes(gene,Q,fill=gene))+
  geom_bar(stat="identity")+xlab("")+
  ggtitle("Myeloid")+ylab("Q statistic")+
  theme(axis.ticks.length=unit(0.5,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill="none")+
  coord_flip()+theme_bw()
ggsave("myeloid/myeloid_Q.pdf",p2,width = 4,height = 4)
p3<-ggplot(q_lymphoid[1:12,],aes(gene,Q,fill=gene))+
  geom_bar(stat="identity")+xlab("")+
  ggtitle("Lymphoid")+ylab("Q statistic")+
  theme(axis.ticks.length=unit(0.5,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill="none")+
  coord_flip()+theme_bw()
ggsave("lymphoid/lymphoid_Q.pdf",p3,width = 4,height = 3)


dir.create("enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
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
ego_BP@result<-ego_BP@result[c(-3),]
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("")
ggsave('enrich/enrichGO.pdf', p_BP, width = 8,height = 4)

genelist <- bitr(hubgene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
eKEGG<-data.frame(ekegg)
write.csv(eKEGG,'enrich/enrichKEGG.csv')
save(ekegg,eKEGG,ego_bp,ego_BP,file = "enrich/GO_KEGG.Rdata")
p1 <- barplot(ekegg, showCategory=10)
ggsave("enrich/enrichKEGG.pdf", plot = p1, width = 8, height = 4)


setwd("E:/singlecall_immune/RCC/consensus/summary4/")
hubdata<-logCPM[hubgene,]
library(ConsensusClusterPlus)
hubdata<-as.matrix(hubdata)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 5,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "E:/singlecall_immune/RCC/consensus/summary4/RCC",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "pdf")
save(results,file = "RCC/result.Rdata")
res_tree<-results[[3]]$consensusClass
res_tree<-data.frame(res_tree)
res_tree<-data.frame(sample=row.names(res_tree),cluster=res_tree$res_tree)
annodata<-data.frame(row.names = res_tree$sample,cluster=res_tree$cluster,ICB=clinical_donor$ICB_Response)
annodata<-annodata[order(annodata$cluster,annodata$ICB),]
annodata$cluster<-as.character(annodata$cluster)
hubdata<-hubdata[,rownames(annodata)]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,main = "",
         filename = "RCC/hubgene.pdf",
         width = 5,height = 5)
p1<-pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color)

gn=rownames(n)[p1$tree_row[["order"]]]
hubdata<-as.data.frame(hubdata)
hubdata_ave<-sapply(split(rownames(annodata)[1:7],annodata$cluster[1:7]),
                    function(cells) rowMeans(hubdata[,cells]))
hubdata_ave<-as.data.frame(hubdata_ave)
hubdata_ave$'3'<-hubdata$P906
n<-t(scale(t(hubdata_ave)))
pheatmap(n[gn,],cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,height = 5,width = 2.5,
         filename = "RCC/hubgene_ave.pdf")

g<-c("STAT1","IRF1","IRF7","IRF9","NR1H3","STAT3","FGR","CSK",
     "NR4A1","NR4A2","JUN","NR3C1","FOS","NR2E1","JUND","RARG","ENG","NFKB1",
     "IL6","SOCS3","TNFAIP3","TCF7L2","RARA","NAMPT","NR2F1")
n<-n[g,]
pheatmap(n,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,height = 5,width = 2,
         filename = "RCC/hubgene_功能排序.pdf")

#############
load("E:/singlecall_immune/RCC/data/donor_clin.Rdata")
clinical2$cluster<-NA
for(i in 1:8){
  clinical2$cluster<-ifelse(clinical2$donor_id==row.names(annodata)[i],
                            annodata$cluster[i],clinical2$cluster)
}
clinical2$cluster<-as.character(clinical2$cluster)

library(stringr)
clinical2$celltype<-ifelse(str_detect(clinical2$FinalCellType,"CD8"),"CD8+ T cell",clinical2$FinalCellType)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"TAM"),"TAM",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"Helper"),"T-Helper",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"MitoHigh NK"),"NK cell",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"FGFBP2"),"NK cell",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"DC"),"Dendritic cell",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"Monocyte"),"Monocyte",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"Macrophage"),"Macrophage",clinical2$celltype)
clinical2$celltype<-ifelse(str_detect(clinical2$celltype,"Myeloid"),"Macrophage",clinical2$celltype)
save(clinical2,clinical_donor,file = "E:/singlecall_immune/RCC/data/donor_clin.Rdata")

load("E:/singlecall_immune/RCC/RCCsingledata.Rdata")
hubdata<-rcc_data[hubtum,]
clinical3<-clinical2[which(clinical2$Lineage=="Putative Tumor"),]
clinical3<-clinical3[order(clinical3$cluster,clinical3$ICB_Response),]
annodata1<-data.frame(row.names = clinical3$NAME,
                      ICB=clinical3$ICB_Response,
                      Cluster=clinical3$cluster)
hubdata<-hubdata[,clinical3$NAME]
hubdata<-log2(hubdata+1)
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubdata,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,
         color=color,main = "",
         filename = "tumor/hubtumor.pdf",
         height = 3,width = 5)

library(ggplot2)
library(reshape2)
library(ggpubr)
hubtumor<-t(hubdata)
hubtumor<-as.data.frame(hubtumor)
hubtumor$Cluster<-annodata1$Cluster
hubtumor_melt<-melt(hubtumor,id.vars = "Cluster")
p<-ggplot(data = hubtumor_melt, aes(x = variable, y = value, fill = Cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=Cluster), label = "p.signif")
ggsave("tumor/hubgene_box.pdf",p,height = 4,width = 6)

dir.create("tumor/enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
ego_BP <- enrichGO(gene          = hubtum,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp<-data.frame(ego_BP)
write.csv(ego_bp,'tumor/enrich/enrichGObp.csv')
genelist <- bitr(hubtum, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
write.csv(genelist,file = "tumor/enrich/gene_entrezID.csv",row.names = F)
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
eKEGG<-data.frame(ekegg)
write.csv(eKEGG,'tumor/enrich/enrichKEGG.csv')
save(ekegg,eKEGG,ego_bp,ego_BP,file = "tumor/enrich/GO_KEGG.Rdata")
ego_BP@result<-ego_BP@result[c(-2,-5),]
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("")
ggsave('tumor/enrich/enrichGO.pdf', p_BP, width = 6,height = 3)
ekegg@result<-ekegg@result[c(3,5,11,13,19,20,28,46),]
p1 <- barplot(ekegg, showCategory=10)
ggsave("tumor/enrich/enrichKEGG.pdf", plot = p1, width = 7, height = 3)

hubdata<-rcc_data[hubmye,]
clinical3<-clinical2[which(clinical2$Lineage=="Myeloid"),]
clinical3<-clinical3[order(clinical3$cluster,clinical3$ICB_Response),]
annodata1<-data.frame(row.names = clinical3$NAME,
                      ICB=clinical3$ICB_Response,
                      Cluster=clinical3$cluster)
hubdata<-hubdata[,clinical3$NAME]
hubdata<-log2(hubdata+1)
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubdata,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,
         color=color,main = "",
         filename = "myeloid/hubmyeloid.pdf",
         width = 5,height = 3)

hubmyeloid<-t(hubdata)
hubmyeloid<-as.data.frame(hubmyeloid)
hubmyeloid$Cluster<-annodata1$Cluster
hubmyeloid_melt<-melt(hubmyeloid,id.vars = "Cluster")
p<-ggplot(data = hubmyeloid_melt, aes(x = variable, y = value, fill = Cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=Cluster), label = "p.signif")
ggsave("myeloid/hubgene_box.pdf",p,height = 4,width = 6)

dir.create("myeloid/enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
ego_BP <- enrichGO(gene          = hubmye,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp<-data.frame(ego_BP)
write.csv(ego_bp,'myeloid/enrich/enrichGObp.csv')
genelist <- bitr(hubmye, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
write.csv(genelist,file = "myeloid/enrich/gene_entrezID.csv",row.names = F)
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
eKEGG<-data.frame(ekegg)
write.csv(eKEGG,'myeloid/enrich/enrichKEGG.csv')
save(ekegg,eKEGG,ego_bp,ego_BP,file = "myeloid/enrich/GO_KEGG.Rdata")
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("")
ggsave('myeloid/enrich/enrichGO.pdf', p_BP, width = 6,height = 3)
ekegg@result<-ekegg@result[c(3,5,9,10,12,14:16,20),]
p1 <- barplot(ekegg, showCategory=10)
ggsave("myeloid/enrich/enrichKEGG.pdf", plot = p1, width = 7, height = 3)

hubdata<-rcc_data[hublym,]
clinical3<-clinical2[which(clinical2$Lineage=="Lymphoid"),]
clinical3<-clinical3[order(clinical3$cluster,clinical3$ICB_Response),]
annodata1<-data.frame(row.names = clinical3$NAME,
                      ICB=clinical3$ICB_Response,
                      Cluster=clinical3$cluster)
hubdata<-hubdata[,clinical3$NAME]
hubdata<-log2(hubdata+1)
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubdata,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,
         color=color,main = "",
         filename = "lymphoid/hublymphoid.pdf",
         width = 5,height = 3)

library(ggplot2)
library(reshape2)
hublymphoid<-t(hubdata)
hublymphoid<-as.data.frame(hublymphoid)
hublymphoid$Cluster<-annodata1$Cluster
hublymphoid_melt<-melt(hublymphoid,id.vars = "Cluster")
p<-ggplot(data = hublymphoid_melt, aes(x = variable, y = value, fill = Cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") + ylab("Expression") + theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=Cluster), label = "p.signif")
ggsave("lymphoid/hubgene_box.pdf",p,height = 4,width = 6)
dir.create("lymphoid/enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
ego_BP <- enrichGO(gene          = hublym,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp<-data.frame(ego_BP)
write.csv(ego_bp,'lymphoid/enrich/enrichGObp.csv')
genelist <- bitr(hublym, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
write.csv(genelist,file = "lymphoid/enrich/gene_entrezID.csv",row.names = F)
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
eKEGG<-data.frame(ekegg)
write.csv(eKEGG,'lymphoid/enrich/enrichKEGG.csv')
save(ekegg,eKEGG,ego_bp,ego_BP,file = "lymphoid/enrich/GO_KEGG.Rdata")
ego_BP@result<-ego_BP@result[c(6:16),]
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("")
ggsave('lymphoid/enrich/enrichGO.pdf', p_BP, width = 6,height = 3)
ekegg@result<-ekegg@result[c(7,8,19,20,27,25,23),]
p1 <- barplot(ekegg, showCategory=10)
ggsave("lymphoid/enrich/enrichKEGG.pdf", plot = p1, width = 7, height = 3)


hubdata<-rcc_data[hubgene,]
annodata1<-data.frame(row.names = clinical2$NAME,cluster=clinical2$cluster,
                      lineage=clinical2$Lineage)
annodata1<-annodata1[order(annodata1$lineage,annodata1$cluster),]
annodata1<-annodata1[which(annodata1$lineage!="Normal Tissue"),]
hubdata<-hubdata[,rownames(annodata1)]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('white','red'))(100)
n<-log2(hubdata+1)
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,
         color=color,
         main = "",filename = "RCC/marker_lineage.pdf",
         height = 5,width = 7)

hubdata<-rcc_data[hubgene,]
annodata1<-data.frame(row.names = clinical2$NAME,cluster=clinical2$cluster,ICB=clinical2$ICB_Response)
annodata1<-annodata1[order(annodata1$cluster,annodata1$ICB),]
hubdata<-hubdata[,rownames(annodata1)]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('white','red'))(100)
n<-log2(hubdata+1)
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,
         color=color,
         main = "",filename = "RCC/marker_hot.pdf",
         height = 5,width = 7)
hubdata_ave<-sapply(split(rownames(annodata1),annodata1$lineage),
                    function(cells) rowMeans(n[,cells]))
n1<-t(scale(t(hubdata_ave)))
pheatmap(n1, #fontsize_row=3, 
         cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,filename = "RCC/marker_lineage_ave.pdf",
         height = 5,width = 3)

clinical2$cluster<-NA
for(i in 1:8){
  clinical2$cluster<-ifelse(clinical2$donor_id==row.names(annodata)[i],
                            annodata$cluster[i],clinical2$cluster)
}
clinical2$cluster<-as.character(clinical2$cluster)
clinical2$number<-1
library(ggplot2)
library(reshape2)
library(ggthemes)
p<-ggplot(clinical2[,c(7,10,11)],aes(cluster,number,fill=Lineage))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("% of total cells")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.5,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  coord_flip()+theme_bw()
ggsave("RCC/cluster_lineage.pdf",p,height = 2,width = 6)

p<-ggplot(clinical2[,c(9:11)],aes(cluster,number,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("% of total cells")+
  theme(axis.ticks.length=unit(0.5,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  coord_flip()+theme_bw()
ggsave("RCC/cluster_cell.pdf",p,height = 2,width = 9)


patient<-as.data.frame(table(clinical2$donor_id))
patient$Var1<-as.character(patient$Var1)
patient$cluster<-NA
patient$ICB<-NA
for (i in 1:8) {
  patient$cluster<-ifelse(patient$Var1==rownames(annodata)[i],annodata$cluster[i],patient$cluster)
  patient$ICB<-ifelse(patient$Var1==rownames(annodata)[i],annodata$ICB[i],patient$ICB)

}
patient<-patient[order(patient$cluster,patient$ICB),]

cellper<-matrix(0,4,3,dimnames = list(unique(clinical2$Lineage),
                                      c("C1/C2","C2/C3","C1/C3")))
cellper<-as.data.frame(cellper)
for (i in 1:4) {
  n1<-clinical2[which(clinical2$cluster=="1"),]
  n2<-clinical2[which(clinical2$cluster=="2"),]
  n3<-clinical2[which(clinical2$cluster=="3"),]
  c1<-sum(n1$Lineage==rownames(cellper)[i])/nrow(n1)
  c2<-sum(n2$Lineage==rownames(cellper)[i])/nrow(n2)
  c3<-sum(n3$Lineage==rownames(cellper)[i])/nrow(n3)
  cellper[i,1]<-c1/c2
  cellper[i,2]<-c2/c3
  cellper[i,3]<-c1/c3
}

cellper_ran<-matrix(0,4000,3,dimnames = list(c(1:4000),
                                      c("C1/C2","C2/C3","C1/C3")))
cellper_ran<-as.data.frame(cellper_ran)
cellper_ran$cluster<-rep(unique(clinical2$Lineage),1000)
for (i in 1:1000) {
  n1<-sample(clinical2$NAME,19090,replace = F)
  m<-clinical2[-which(clinical2$NAME%in%n1),]
  n2<-sample(m$NAME,12766,replace = F)
  m<-m[-which(m$NAME%in%n2),]
  n3<-sample(m$NAME,2470,replace = F)
  b1<-clinical2[clinical2$NAME%in%n1,]
  b2<-clinical2[clinical2$NAME%in%n2,]
  b3<-clinical2[clinical2$NAME%in%n3,]
  for (j in 1:4) {
    c1<-sum(b1$Lineage==rownames(cellper)[j])/nrow(b1)
    c2<-sum(b2$Lineage==rownames(cellper)[j])/nrow(b2)
    c3<-sum(b3$Lineage==rownames(cellper)[j])/nrow(b3)
    cellper_ran[(i-1)*4+j,1]<-c1/c2
    cellper_ran[(i-1)*4+j,2]<-c2/c3
    cellper_ran[(i-1)*4+j,3]<-c1/c3
    }
}

cellper$p1<-NA
cellper$p2<-NA
cellper$p3<-NA
for (i in 1:4) {
  a<-cellper_ran[which(cellper_ran$cluster==rownames(cellper)[i]),]
  for (j in 1:3) {
    if(cellper[i,j]>1){
      cellper[i,j+3]<-sum(a[,j]>cellper[i,j])/1000
    }else{cellper[i,j+3]<-sum(a[,j]<cellper[i,j])/1000}
  }
}
cellper$lineage<-rownames(cellper)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
per_melt<-melt(cellper[,c(1:3,7)],id.vars = "lineage")
p<-ggplot(data = per_melt, aes(x = lineage, y = value, fill = variable))+
  geom_bar(stat = "identity",width = 0.5,position = "dodge")+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("Fold Change") + theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())
ggsave("RCC/cluster_lineage_FC.pdf",p,height = 4,width = 4)
write.csv(cellper,file = "RCC/lineage_FC_P.csv")
save(cellper,cellper_ran,file = "RCC/lineage_FCandP.Rdata")

cellper<-matrix(0,18,3,dimnames = list(unique(clinical2$celltype),
                                      c("C1/C2","C2/C3","C1/C3")))
cellper<-as.data.frame(cellper)
for (i in 1:18) {
  n1<-clinical2[which(clinical2$cluster=="1"),]
  n2<-clinical2[which(clinical2$cluster=="2"),]
  n3<-clinical2[which(clinical2$cluster=="3"),]
  c1<-sum(n1$celltype==rownames(cellper)[i])/nrow(n1)
  c2<-sum(n2$celltype==rownames(cellper)[i])/nrow(n2)
  c3<-sum(n3$celltype==rownames(cellper)[i])/nrow(n3)
  cellper[i,1]<-c1/c2
  cellper[i,2]<-c2/c3
  cellper[i,3]<-c1/c3
}

cellper_ran<-matrix(0,18000,3,dimnames = list(c(1:18000),
                                             c("C1/C2","C2/C3","C1/C3")))
cellper_ran<-as.data.frame(cellper_ran)
cellper_ran$cluster<-rep(unique(clinical2$celltype),1000)
for (i in 1:1000) {
  n1<-sample(clinical2$NAME,19090,replace = F)
  m<-clinical2[-which(clinical2$NAME%in%n1),]
  n2<-sample(m$NAME,12766,replace = F)
  m<-m[-which(m$NAME%in%n2),]
  n3<-sample(m$NAME,2470,replace = F)
  b1<-clinical2[clinical2$NAME%in%n1,]
  b2<-clinical2[clinical2$NAME%in%n2,]
  b3<-clinical2[clinical2$NAME%in%n3,]
  for (j in 1:18) {
    c1<-sum(b1$celltype==rownames(cellper)[j])/nrow(b1)
    c2<-sum(b2$celltype==rownames(cellper)[j])/nrow(b2)
    c3<-sum(b3$celltype==rownames(cellper)[j])/nrow(b3)
    cellper_ran[(i-1)*18+j,1]<-c1/c2
    cellper_ran[(i-1)*18+j,2]<-c2/c3
    cellper_ran[(i-1)*18+j,3]<-c1/c3
  }
}

cellper$p1<-NA
cellper$p2<-NA
cellper$p3<-NA
for (i in 1:18) {
  a<-cellper_ran[which(cellper_ran$cluster==rownames(cellper)[i]),]
  for (j in 1:3) {
    if(cellper[i,j]>1){
      cellper[i,j+3]<-sum(a[,j]>cellper[i,j])/1000
    }else{cellper[i,j+3]<-sum(a[,j]<cellper[i,j])/1000}
  }
}
cellper$celltype<-rownames(cellper)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
per_melt<-melt(cellper[,c(1:3,7)],id.vars = "celltype")
per_melt<-per_melt[c(-33,-51),]
p<-ggplot(data = per_melt, aes(x = celltype, y = value, fill = variable))+
  geom_bar(stat = "identity",width = 0.5,position = "dodge")+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("Fold Change") + theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())
ggsave("RCC/cluster_celltype_FC.pdf",p,height = 4,width = 7)
write.csv(cellper,file = "RCC/celltype_FC_P.csv")
save(cellper,cellper_ran,file = "RCC/celltype_FCandP.Rdata")

#####################################
lim<-clinical2[,c(1,10)]
save(lim,file = "RCC/cluster_limma_group.Rdata")

setwd("/pub6/temp/wangli/yxz/RCC")
library(limma)
load("RCCsingledata.Rdata")
load("RCC/cluster_limma_group.Rdata")
data<-rcc_data[which(apply(rcc_data,1,function(x){return(sum(x!=0))})>ncol(rcc_data)*0.05),]
lim$cluster<-ifelse(lim$cluster=="1","C1","C23")
design<-model.matrix(~0+factor(lim$cluster))
colnames(design)<-levels(factor(lim$cluster))
rownames(design)<-lim$cluster
contrast.matrix<-makeContrasts(paste0(unique(lim$cluster),collapse = "-"),levels = design)
fit<-lmFit(data,design)
fit1<-contrasts.fit(fit,contrast.matrix)
fit1<-eBayes(fit1)
tempoutput<-topTable(fit1,coef = 1,n=Inf)
save(tempoutput,file = "summary4/limma_result.Rdata")
save(data,file="summary4/scRNA_seq%5.Rdata")

setwd("E:/singlecall_immune/RCC/consensus/summary4/")

filename<-list.files(path = 'E:/singlecall_immune/RCC/cancerSEA/',
                     pattern = "*.txt")
a<-read.table(file = paste0("E:/singlecall_immune/RCC/cancerSEA/",filename[1]),header = T)
scRNASCP<-readRDS("E:/singlecall_immune/RCC/scRNASCP.rds")
clinical3<-clinical2[clinical2$NAME%in%rownames(scRNASCP@meta.data),]
rownames(clinical3)<-clinical3$NAME
clinical3<-clinical3[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical3$cluster
scRNASCP@meta.data$celltype<-clinical3$celltype
scRNASCP@meta.data$ICB<-clinical3$ICB_Response
library(Seurat)
library(tidyverse)
library(stringr)
genename<-rownames(scRNASCP)
for(i in 1:14){
  a<-read.table(file = paste0("E:/singlecall_immune/RCC/cancerSEA/",filename[i]),header = T)
  ALPHA<-str_split(filename[i],"\\.",simplify = T)[,1]
  c<-intersect(a$GeneName,genename)
  c<-as.data.frame(c)
  colnames(c)<-ALPHA
  c1<-as.list(c)
  scRNASCP<-AddModuleScore(object = scRNASCP,features = c1,
                           name = ALPHA)
  colnames(scRNASCP@meta.data)[9+i]<-ALPHA
}

geneset_score<-scRNASCP@meta.data
geneset_score<-geneset_score[,-(1:5)]
save(geneset_score,file = "E:/singlecall_immune/RCC/cancerSEA/cancerSEA基因集特征分数.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
sea_melt<-melt(geneset_score[,c(3,5:18)],id.vars = "cluster")
p<-ggplot(data = sea_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Score") +
  theme_classic() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(sea_melt$value)$stats[c(1, 5)])*1.05)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 0.2)
ggsave("RCC/cell_state.pdf",p,height = 4,width = 6)


scRNASCP<-readRDS("E:/singlecall_immune/RCC/tumor/scRNAsubcell.rds")
clinical3<-clinical2[clinical2$NAME%in%rownames(scRNASCP@meta.data),]
rownames(clinical3)<-clinical3$NAME
clinical3<-clinical3[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical3$cluster
scRNASCP@meta.data$celltype<-clinical3$celltype
scRNASCP@meta.data$ICB<-clinical3$ICB_Response
library(Seurat)
library(tidyverse)
library(stringr)
genename<-rownames(scRNASCP)
for(i in 1:14){
  a<-read.table(file = paste0("E:/singlecall_immune/RCC/cancerSEA/",filename[i]),header = T)
  ALPHA<-str_split(filename[i],"\\.",simplify = T)[,1]
  c<-intersect(a$GeneName,genename)
  c<-as.data.frame(c)
  colnames(c)<-ALPHA
  c1<-as.list(c)
  scRNASCP<-AddModuleScore(object = scRNASCP,features = c1,
                           name = ALPHA)
  colnames(scRNASCP@meta.data)[10+i]<-ALPHA
}

geneset_score<-scRNASCP@meta.data
geneset_score<-geneset_score[,-(1:5)]
save(geneset_score,file = "E:/singlecall_immune/RCC/cancerSEA/cancerSEA基因集特征分数_tumor.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
sea_melt<-melt(geneset_score[,c(4,6:19)],id.vars = "cluster")
p<-ggplot(data = sea_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Score") +
  theme_classic() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(sea_melt$value)$stats[c(1, 5)])*2)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 0.4)
ggsave("RCC/cell_state_tumor.pdf",p,height = 4,width = 6)


scRNASCP<-readRDS("E:/singlecall_immune/RCC/myeloid/scRNAsubmy.rds")
clinical3<-clinical2[clinical2$NAME%in%rownames(scRNASCP@meta.data),]
rownames(clinical3)<-clinical3$NAME
clinical3<-clinical3[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical3$cluster
scRNASCP@meta.data$celltype<-clinical3$celltype
scRNASCP@meta.data$ICB<-clinical3$ICB_Response
library(Seurat)
library(tidyverse)
library(stringr)
genename<-rownames(scRNASCP)
for(i in 1:14){
  a<-read.table(file = paste0("E:/singlecall_immune/RCC/cancerSEA/",filename[i]),header = T)
  ALPHA<-str_split(filename[i],"\\.",simplify = T)[,1]
  c<-intersect(a$GeneName,genename)
  c<-as.data.frame(c)
  colnames(c)<-ALPHA
  c1<-as.list(c)
  scRNASCP<-AddModuleScore(object = scRNASCP,features = c1,
                           name = ALPHA)
  colnames(scRNASCP@meta.data)[10+i]<-ALPHA
}

geneset_score<-scRNASCP@meta.data
geneset_score<-geneset_score[,-(1:5)]
save(geneset_score,file = "E:/singlecall_immune/RCC/cancerSEA/cancerSEA基因集特征分数_myeloid.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
sea_melt<-melt(geneset_score[,c(4,6:19)],id.vars = "cluster")
p<-ggplot(data = sea_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Score") +
  theme_classic() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(sea_melt$value)$stats[c(1, 5)])*1.5)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 0.3)
ggsave("RCC/cell_state_myeloid.pdf",p,height = 4,width = 6)


scRNASCP<-readRDS("E:/singlecall_immune/RCC/lymphoid/scRNAsublym.rds")
clinical3<-clinical2[clinical2$NAME%in%rownames(scRNASCP@meta.data),]
rownames(clinical3)<-clinical3$NAME
clinical3<-clinical3[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical3$cluster
scRNASCP@meta.data$celltype<-clinical3$celltype
scRNASCP@meta.data$ICB<-clinical3$ICB_Response
library(Seurat)
library(tidyverse)
library(stringr)
genename<-rownames(scRNASCP)
for(i in 1:14){
  a<-read.table(file = paste0("E:/singlecall_immune/RCC/cancerSEA/",filename[i]),header = T)
  ALPHA<-str_split(filename[i],"\\.",simplify = T)[,1]
  c<-intersect(a$GeneName,genename)
  c<-as.data.frame(c)
  colnames(c)<-ALPHA
  c1<-as.list(c)
  scRNASCP<-AddModuleScore(object = scRNASCP,features = c1,
                           name = ALPHA)
  colnames(scRNASCP@meta.data)[10+i]<-ALPHA
}

geneset_score<-scRNASCP@meta.data
geneset_score<-geneset_score[,-(1:5)]
save(geneset_score,file = "E:/singlecall_immune/RCC/cancerSEA/cancerSEA基因集特征分数_lymphoid.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
sea_melt<-melt(geneset_score[,c(4,6:19)],id.vars = "cluster")
p<-ggplot(data = sea_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +ylab("Score") +
  theme_classic() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(sea_melt$value)$stats[c(1, 5)])*1.05)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 0.2)
ggsave("RCC/cell_state_lymphoid.pdf",p,height = 4,width = 6)

###################################
load("E:/singlecall_immune/RCC/data/pseudo-bulk.Rdata")
donorlog<-donor_expr[rowSums(donor_expr)>10,]
donorlog<-log2(donorlog+1)
#使用pseudo——bulk，C1，C2，C3  GSVA分析
dir.create("RCC/GSVA")
library(GSVA)
library(GSVAdata)
library(stringr)
geneset <- getGmt("E:/singlecall_immune/RCC/GSVA/h.all.v7.5.1.symbols.gmt")  
va_kirc<-gsva(as.matrix(donorlog),geneset,mx.diff=FALSE, verbose=TRUE)
row.names(va_kirc)<-str_sub(rownames(va_kirc),10,100)
va_kirc<-va_kirc[,rownames(annodata)]
save(va_kirc,file = "RCC/GSVA/GSVA_donor.Rdata")

library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n<-t(scale(t(va_kirc)))
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,main = "",
         legend_breaks = c(-1.5,0,1.5),
         legend_labels = c("-1.5","0","1.5"),
         filename = "RCC/GSVA/GSVA.pdf",
         height = 8,width = 6)

library(limma)
lim<-data.frame(row.names = row.names(annodata),
                cluster=c("C1","C1","C1","C2","C2","C2","C2","C3"))
design<-model.matrix(~0+factor(lim$cluster))
colnames(design)=levels(factor(lim$cluster))
rownames(design)<-rownames(lim)
contrast.matrix<-makeContrasts("C2-C3",levels=design)
fit <- lmFit(va_kirc,design)
fit23 <- contrasts.fit(fit, contrast.matrix)
fit23 <- eBayes(fit23)  
tempOutput23 = topTable(fit23, coef=1, n=Inf)
contrast.matrix<-makeContrasts("C1-C2",levels=design)
fit <- lmFit(va_kirc,design)
fit12 <- contrasts.fit(fit, contrast.matrix)
fit12 <- eBayes(fit12)  
tempOutput12 = topTable(fit12, coef=1, n=Inf)
contrast.matrix<-makeContrasts("C1-C3",levels=design)
fit <- lmFit(va_kirc,design)
fit13 <- contrasts.fit(fit, contrast.matrix)
fit13 <- eBayes(fit13)  
tempOutput13 = topTable(fit13, coef=1, n=Inf)

lim<-data.frame(row.names = row.names(annodata),
                cluster=c("C1","C1","C1","C23","C23","C23","C23","C23"))
design<-model.matrix(~0+factor(lim$cluster))
colnames(design)=levels(factor(lim$cluster))
rownames(design)<-rownames(lim)
contrast.matrix<-makeContrasts("C1-C23",levels=design)
fit <- lmFit(va_kirc,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)  
tempOutput = topTable(fit, coef=1, n=Inf)
save(tempOutput,tempOutput12,tempOutput23,tempOutput13,file = "RCC/GSVA/GSVAlimma_result.Rdata")

deggsva<-c(row.names(tempOutput)[which(tempOutput$P.Value<0.05)],
           row.names(tempOutput12)[which(tempOutput12$P.Value<0.05)],
           row.names(tempOutput13)[which(tempOutput13$P.Value<0.05)],
           row.names(tempOutput23)[which(tempOutput23$P.Value<0.05)])
deggsva<-deggsva[!duplicated(deggsva)]
degdata<-va_kirc[deggsva,]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n<-t(scale(t(degdata)))
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,main = "",
         legend_breaks = c(-1.5,0,1.5),
         legend_labels = c("-1.5","0","1.5"),
         filename = "RCC/GSVA/GSVA_DEG.pdf",
         height = 3,width = 6)


require(ggplot2)
require(SCENIC)
require(AUCell)

load("E:/singlecall_immune/RCC/tumor/SCENIC_AUC.Rdata")
clinical2<-clinical2[order(clinical2$cluster),]
clinical3<-clinical2[which(clinical2$Lineage=="Putative Tumor"),]

hubauc<-regulon[hubtum,]
hubauc<-na.omit(hubauc)
library(pheatmap)
library(RColorBrewer)
a1<-data.frame(row.names = clinical3$NAME,cluster=clinical3$cluster,ICB=clinical3$ICB_Response)
a1<-a1[order(a1$cluster,a1$ICB),]
hubauc<-hubauc[,row.names(a1)]
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubauc,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = a1,
         color=color,main = "",
         filename = "scenicAUC/tumorhubauc.pdf",
         height = 3,width = 7)

hubauc1<-t(hubauc)
hubauc1<-as.data.frame(hubauc1)
hubauc1<-hubauc1[clinical3$NAME,]
hubauc1$cluster<-clinical3$cluster
#hubauc1$cluster<-ifelse(hubauc1$cluster!="1","2_3",hubauc1$cluster)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
auc_melt<-melt(hubauc1,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.48)
ggsave("scenicAUC/tumor/AUCtumor差异.pdf",p,height = 4,width = 6)

hubauc2<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response!="NoICB")],]
auc_melt<-melt(hubauc2,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.45)
ggsave("scenicAUC/tumor/AUCtumor_PRNR差异.pdf",p,height = 4,width = 6)

hubauc3<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response=="NoICB")],]
auc_melt<-melt(hubauc3,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.48)
ggsave("scenicAUC/tumor/AUCtumor_NOICB差异.pdf",p,height = 4,width = 6)

load("E:/singlecall_immune/RCC/myeloid/SCENIC_AUC.Rdata")
clinical3<-clinical2[which(clinical2$Lineage=="Myeloid"),]
hubauc<-regulon[hubmye,]
hubauc<-na.omit(hubauc)
a1<-data.frame(row.names = clinical3$NAME,cluster=clinical3$cluster,ICB=clinical3$ICB_Response)
a1<-a1[order(a1$cluster,a1$ICB),]
hubauc<-hubauc[,row.names(a1)]
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubauc,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = a1,
         color=color,main = "",
         filename = "scenicAUC/myehoidhubauc.pdf",
         height = 3,width = 7)

hubauc1<-t(hubauc)
hubauc1<-as.data.frame(hubauc1)
hubauc1<-hubauc1[clinical3$NAME,]
hubauc1$cluster<-clinical3$cluster
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
auc_melt<-melt(hubauc1,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.5)
ggsave("scenicAUC/myeloid/AUCmyeloid差异.pdf",p,height = 4,width = 6)

hubauc2<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response!="NoICB")],]
auc_melt<-melt(hubauc2,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.6)
ggsave("scenicAUC/myeloid/AUCmyeploid_PRNR差异.pdf",p,height = 4,width = 6)

hubauc3<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response=="NoICB")],]
auc_melt<-melt(hubauc3,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.45)
ggsave("scenicAUC/myeloid/AUCmyeloid_NOICB差异.pdf",p,height = 4,width = 6)

load("E:/singlecall_immune/RCC/lymphoid/SCENIC_AUC.Rdata")
clinical3<-clinical2[which(clinical2$Lineage=="Lymphoid"),]
hubauc<-regulon[hublym,]
hubauc<-na.omit(hubauc)
a1<-data.frame(row.names = clinical3$NAME,cluster=clinical3$cluster,ICB=clinical3$ICB_Response)
a1<-a1[order(a1$cluster,a1$ICB),]
hubauc<-hubauc[,rownames(a1)]
color<- colorRampPalette(c('white','red'))(100)
pheatmap(hubauc,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = a1,
         color=color,main = "",
         filename = "scenicAUC/lymphoidhubauc.pdf",
         height = 2.5,width = 7)

hubauc1<-t(hubauc)
hubauc1<-as.data.frame(hubauc1)
hubauc1<-hubauc1[clinical3$NAME,]
hubauc1$cluster<-clinical3$cluster
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
auc_melt<-melt(hubauc1,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.32)
ggsave("scenicAUC/lymphoid/AUClymphoid差异.pdf",p,height = 4,width = 6)

hubauc2<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response!="NoICB")],]
auc_melt<-melt(hubauc2,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.38)
ggsave("scenicAUC/lymphoid/AUClymphoid_PRNR差异.pdf",p,height = 4,width = 6)

hubauc3<-hubauc1[clinical3$NAME[which(clinical3$ICB_Response=="NoICB")],]
auc_melt<-melt(hubauc3,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.25)
ggsave("scenicAUC/lymphoid/AUClymphoid_NOICB差异.pdf",p,height = 4,width = 6)

load("E:/singlecall_immune/RCC/tumor/SCENIC_AUC.Rdata")
auctumor<-regulon[hubgene,]
auctumor[is.na(auctumor)]<-0
row.names(auctumor)<-hubgene
load("E:/singlecall_immune/RCC/myeloid/SCENIC_AUC.Rdata")
aucmyeloid<-regulon[hubgene,]
aucmyeloid[is.na(aucmyeloid)]<-0
row.names(aucmyeloid)<-hubgene
load("E:/singlecall_immune/RCC/lymphoid/SCENIC_AUC.Rdata")
auclymphoid<-regulon[hubgene,]
auclymphoid[is.na(auclymphoid)]<-0
row.names(auclymphoid)<-hubgene
clinical3<-clinical2[which(clinical2$Lineage!="Normal Tissue"),]
a1<-data.frame(row.names = clinical3$NAME,cluster=clinical3$cluster,ICB=clinical3$ICB_Response)
a1<-a1[order(a1$cluster,a1$ICB),]
aucall<-cbind(auctumor,aucmyeloid,auclymphoid)
aucall<-aucall[which(rowSums(aucall)!=0),]
aucall<-aucall[,row.names(a1)]
color<- colorRampPalette(c('white','red'))(100)
pheatmap(aucall,
         show_rownames = T,show_colnames = F,
         cluster_rows = T,cluster_cols = F,
         annotation_col = a1,
         color=color,main = "",
         filename = "scenicAUC/all_hubauc.pdf",
         height = 4,width = 7)
aucall1<-t(aucall)
aucall1<-as.data.frame(aucall1)
aucall1<-aucall1[clinical3$NAME,]
aucall1$cluster<-clinical3$cluster
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggpubr)
auc_melt<-melt(aucall1,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +
  ylab("AUCell") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.35)
ggsave("scenicAUC/AUCall差异.pdf",p,height = 3,width = 6)

aucall2<-aucall1[clinical3$NAME[which(clinical3$ICB_Response!="NoICB")],]
auc_melt<-melt(aucall2,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("AUCell") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.35)
ggsave("scenicAUC/AUCall_PRNR差异.pdf",p,height = 3,width = 6)

aucall3<-aucall1[clinical3$NAME[which(clinical3$ICB_Response=="NoICB")],]
auc_melt<-melt(aucall3,id.vars = "cluster")
p<-ggplot(data = auc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,position = position_dodge(0.5),
               outlier.shape = NA)+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("AUCell") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(auc_melt$value)$stats[c(1, 5)])*1.2)+
  stat_compare_means(aes(group=cluster),label = "p.signif",label.y = 0.32)
ggsave("scenicAUC/AUCall_NOICB差异.pdf",p,height = 3,width = 6)


enrich<-function(target,cluster){
  ego_GO <- enrichGO(gene          = target,
                     #universe     = row.names(dge.celltype),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05) 
  genelist <- bitr(target, fromType="SYMBOL",
                   toType="ENTREZID", OrgDb='org.Hs.eg.db')
  genelist <- pull(genelist,ENTREZID)               
  ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
  save(ekegg,ego_GO,file = paste0("RCC/enrich/",cluster,".Rdata"))
}

load("E:/singlecall_immune/RCC/tumor/TF.Rdata")
tf_tumor<-tflist
load("E:/singlecall_immune/RCC/myeloid/TF.Rdata")
tf_myeloid<-tflist
load("E:/singlecall_immune/RCC/lymphoid/TF.Rdata")
tf_lymphoid<-tflist
dir.create("RCC/enrich")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
target<-tf_tumor[which(tf_tumor$TF=="NR2E1"),2]
enrich(target,"tumor_NR2E1")

z<-c("STAT1","IRF1","NR1H3","STAT3","JUN",
     "IRF7","JUND","FOS","TCF7L2","RARA")
for (i in 1:10) {
  target<-tf_myeloid[which(tf_myeloid$TF==z[i]),2]
  enrich(target,paste0("myeloid_",z[i]))
}
z<-c("STAT1","IRF1","NR3C1","STAT3","JUN",
     "IRF7","JUND","FOS","IRF9")
for (i in 1:9) {
  target<-tf_lymphoid[which(tf_lymphoid$TF==z[i]),2]
  enrich(target,paste0("lymphoid_",z[i]))
}
setwd("E:/singlecall_immune/RCC/consensus/summary4/")

net_expr<-as.data.frame(net_expr)
data<-net_expr[tf_lymphoid[which(tf_lymphoid$TF=="NR3C1"),2],]
target<-rownames(net_expr)[which(data$C1<data$C2)]
enrich(target,"lymphoid_NR3C1_C2")


library(pheatmap)
library(RColorBrewer)
gn<-c("FGR","IRF1","CSK","IRF7","NR1H3",  
      "STAT1","IRF9","TNFAIP3","NR3C1","NAMPT",  
      "NR4A2","NR4A1","FOS","SOCS3","JUND",   
      "JUN","STAT3","RARG","IL6","NR2F1",  
      "NR2E1","NFKB1","RARA","TCF7L2","ENG")
fun<-read.csv("RCC/gene_function.csv",row.names = 1,header = T)
fun[is.na(fun)]<-0
fun<-as.data.frame(t(fun))
fun<-fun[,gn]
annodata2<-data.frame(gene=hubgene,tumor=NA,lymphoid=NA,myeloid=NA)
annodata2[which(annodata2$gene%in%hubtum),2]<-"tumor"
annodata2[which(annodata2$gene%in%hublym),3]<-"lymphoid"
annodata2[which(annodata2$gene%in%hubmye),4]<-"myeloid"
rownames(annodata2)<-annodata2$gene
annodata2<-annodata2[,-1]
annodata2<-annodata2[gn,]
pheatmap(fun,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("white","orange"))(100),
         border_color = 4,height = 4,width = 8,
         annotation_col = annodata2,
         filename = "RCC/gene功能特征.pdf")
pheatmap(fun,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("white","orange"))(100),
         border_color = 4,annotation_col = annodata2)


a<-intersect(tflist$TF,hubtum)
allgene<-c()
for (i in 1:11) {
  tfgene<-tflist[tflist$TF==a[i],2]
  write.table(tfgene,file = paste0("tumor/network/",a[i],".txt"),
              row.names = F,col.names = F)
  allgene<-c(allgene,tfgene)
}
allgene<-allgene[!duplicated(allgene)]
write.table(allgene,file = "tumor/network/allgene.txt",
            row.names = F,col.names = F)

clinical3<-clinical2[which(clinical2$Lineage=="Putative Tumor"),]
expr_var<-rcc_data[rownames(tuopu),clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1","C2","C3")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "tumor/network/expression.csv")
save(net_expr,file = "tumor/network/expression.Rdata")
clinical3<-clinical3[which(clinical3$ICB_Response!="NoICB"),]
expr_var<-expr_var[,clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1_PR","C2_NR","C3_NR")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "tumor/network/expression_ICB.csv")


a<-intersect(tflist$TF,hubmye)
allgene<-c()
for (i in 1:11) {
  tfgene<-tflist[tflist$TF==a[i],2]
  write.table(tfgene,file = paste0("myeloid/network/",a[i],".txt"),
              row.names = F,col.names = F)
  allgene<-c(allgene,tfgene)
}
allgene<-allgene[!duplicated(allgene)]
write.table(allgene,file = "myeloid/network/allgene.txt",
            row.names = F,col.names = F)

clinical3<-clinical2[which(clinical2$Lineage=="Myeloid"),]
expr_var<-rcc_data[rownames(tuopu),clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1","C2","C3")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "myeloid/network/expression.csv")
save(net_expr,file = "myeloid/network/expression.Rdata")
clinical3<-clinical3[which(clinical3$ICB_Response!="NoICB"),]
expr_var<-expr_var[,clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1_PR","C2_NR","C3_NR")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "myeloid/network/expression_ICB.csv")


a<-intersect(tflist$TF,hublym)
allgene<-c()
for (i in 1:9) {
  tfgene<-tflist[tflist$TF==a[i],2]
  write.table(tfgene,file = paste0("lymphoid/network/",a[i],".txt"),
              row.names = F,col.names = F)
  allgene<-c(allgene,tfgene)
}
allgene<-allgene[!duplicated(allgene)]
write.table(allgene,file = "lymphoid/network/allgene.txt",
            row.names = F,col.names = F)

clinical3<-clinical2[which(clinical2$Lineage=="Lymphoid"),]
expr_var<-rcc_data[rownames(tuopu),clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1","C2","C3")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "lymphoid/network/expression.csv")
save(net_expr,file = "lymphoid/network/expression.Rdata")
clinical3<-clinical3[which(clinical3$ICB_Response!="NoICB"),]
expr_var<-expr_var[,clinical3$NAME]
net_expr<-sapply(split(clinical3$NAME,clinical3$cluster),
                 function(cells) rowMeans(expr_var[,cells]))
net_expr<-as.data.frame(net_expr)
colnames(net_expr)<-c("C1_PR","C2_NR","C3_NR")
net_expr<-t(scale(t(net_expr)))
write.csv(net_expr,file = "lymphoid/network/expression_ICB.csv")



setwd("/pub6/temp/wangli/yxz/RCC")
library(monocle)
scRNASCP<-readRDS("scRNASCP.rds")
clinical3<-clinical2[clinical2$NAME%in%rownames(scRNASCP@meta.data),]
rownames(clinical3)<-clinical3$NAME
clinical3<-clinical3[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical3$cluster
scRNASCP@meta.data$celltype<-clinical3$celltype
scRNASCP@meta.data$ICB<-clinical3$ICB_Response
scRNASCP@meta.data$donor<-clinical3$donor_id
dir.create("pseudotime")
data <- as(as.matrix(scRNASCP@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNASCP@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
var.genes <- VariableFeatures(scRNASCP)
mycds <- setOrderingFilter(mycds, var.genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
plot1 <- plot_cell_trajectory(mycds, color_by = "cluster")
ggsave("pseudotime/cluster.pdf", plot = plot1, width = 6, height = 5)
plot1 <- plot_cell_trajectory(mycds, color_by = "ICB")
ggsave("pseudotime/ICB.pdf", plot = plot1, width = 6, height = 5)
plot1 <- plot_cell_trajectory(mycds, color_by ="donor")
ggsave("pseudotime/donor.pdf", plot = plot1, width = 6, height = 5)
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/pseudotime.pdf", plot = plot3, width = 6, height = 5)
save(mycds,file="pseudotime/pesudo_mycds.Rdata")
plot1 <- plot_cell_trajectory(mycds, color_by = "Lineage")
ggsave("pseudotime/lineage.pdf", plot = plot1, width = 6, height = 5)

annogene<-read.csv("pseudotime/hubgene.csv",row.names = 1)
annogene<-annogene[order(annogene$lineage1,annogene$lineage2,annogene$lineage3),]
annogene[annogene==""]<-NA
p1 <- plot_pseudotime_heatmap(mycds[rownames(annogene),], num_clusters=1,show_rownames=T,
                             cluster_rows = F)
ggsave("pseudotime/hotmap.pdf",p1,height=6,width=4)


library(Seurat)
library(tidyverse)
library(stringr)
rownames(clinical2)<-clinical2$NAME
clinical2<-clinical2[rownames(scRNASCP@meta.data),]
scRNASCP@meta.data$cluster<-clinical2$cluster
scRNASCP@meta.data$celltype<-clinical2$celltype
scRNASCP@meta.data$donor<-clinical2$donor_id
p1 = DimPlot(scRNASCP, group.by="cluster", label=T, label.size=3, reduction='tsne')
ggsave("RCC/tsne_cluster.pdf",p1,width = 6,height = 5)
p1 = DimPlot(scRNASCP, group.by="celltype", label=T, label.size=3, reduction='tsne')
ggsave("RCC/tsne_celltype.pdf",p1,width = 7,height = 5)
p1 = DimPlot(scRNASCP, group.by="donor", label=T, label.size=3, reduction='tsne')
ggsave("RCC/tsne_donor.pdf",p1,width = 6,height = 5)

library(ConsensusClusterPlus)
hubdata<-kirc_gene1[hubgene,]
hubdata<-as.matrix(hubdata)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 10,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "E:/singlecall_immune/RCC/consensus/summary4/TCGA/",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "pdf")

res_tree<-results[[3]]$consensusClass
res_tree<-data.frame(res_tree)
res_tree<-data.frame(sample=row.names(res_tree),cluster=res_tree$res_tree)
res_tree<-res_tree[order(res_tree$cluster),]
save(results,res_tree,file = "TCGA/TCGA_result.Rdata")

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
hubclin<-merge(clinical2,res_tree,by.x = "patient_id",by.y = "sample")
hubclin$cluster<-as.character(hubclin$cluster)
fit_im<-survfit(Surv(time,status)~cluster,data = hubclin)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = T,
              risk.table = TRUE,
              risk.table.col = "strata",
              linetype = "strata",
              ggtheme = theme_bw(),
              legend = c(0.8,0.75), 
              legend.title = "", 
              legend.labs = c("1", "2","3"),
              title="Overall Survival")

dir.create("TCGA/clinical")
clin_feature<-clinical1[,c(1,8:11)]
colnames(clin_feature)<-c("sample","M","N","Stage","T")
library(stringr)
clin_feature$T<-str_sub(clin_feature$T,1,2)
clin_feature$sample<-gsub("-",".",clin_feature$sample)
clin_feature<-merge(clin_feature,res_tree,by = "sample")
clin_feature$cluster<-as.character(clin_feature$cluster)
clin_feature$number<-1
save(clin_feature,clinical1,clinical2,clinical,file = "TCGA/clinical/clinical.Rdata")
library(ggplot2)
library(reshape2)
library(ggthemes)
a<-clin_feature[,c(2,6,7)]
a<-a[which(a$M!="'--"),]
p<-ggplot(a,aes(M,number,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("TCGA/clinical/M_cluster.pdf",p,height = 4,width = 4)
b<-table(a$M,a$cluster)
P_M<-matrix(0,3,3,dimnames = list(c("M0","M1","MX"),
                                  c("M0","M1","MX")))
for (i in 1:3) {
  for (j in 1:3) {
    P_M[i,j]<-chisq.test(b[c(i,j),])$p.value
  }
}
P_M<-round(P_M,2)
pheatmap(P_M, cluster_rows = F,cluster_cols = F,
         display_numbers = TRUE,number_color = "blue",
         filename = "TCGA/clinical/M_chisq.pdf",
         height = 1.5,width = 3)

a<-clin_feature[,c(3,6,7)]
a<-a[which(a$N!="'--"),]
p<-ggplot(a,aes(N,number,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("TCGA/clinical/N_cluster.pdf",p,height = 4,width = 4)
b<-table(a$N,a$cluster)
P_N<-matrix(0,3,3,dimnames = list(c("N0","N1","NX"),
                                  c("N0","N1","NX")))
for (i in 1:3) {
  for (j in 1:3) {
    P_N[i,j]<-chisq.test(b[c(i,j),])$p.value
  }
}
P_N<-round(P_N,2)
P_N[2,2]<-1
pheatmap(P_N, cluster_rows = F,cluster_cols = F,
         display_numbers = TRUE,number_color = "blue",
         filename = "TCGA/clinical/N_chisq.pdf",
         height = 1.5,width = 3)

a<-clin_feature[,c(4,6,7)]
a<-a[which(a$Stage!="'--"),]
p<-ggplot(a,aes(Stage,number,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("TCGA/clinical/Stage_cluster.pdf",p,height = 4,width = 5)
b<-table(a$Stage,a$cluster)
P_Stage<-matrix(0,4,4,dimnames = list(c("Stage I","Stage II","Stage III","Stage IV"),
                                  c("Stage I","Stage II","Stage III","Stage IV")))
for (i in 1:4) {
  for (j in 1:4) {
    P_Stage[i,j]<-chisq.test(b[c(i,j),])$p.value
  }
}
P_Stage<-round(P_Stage,2)
pheatmap(P_Stage, cluster_rows = F,cluster_cols = F,
         display_numbers = TRUE,number_color = "blue",
         filename = "TCGA/clinical/Stage_chisq.pdf",
         height = 2,width = 4)

a<-clin_feature[,c(5,6,7)]
a<-a[which(a$T!="'--"),]
p<-ggplot(a,aes(T,number,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("TCGA/clinical/T_cluster.pdf",p,height = 4,width = 5)
b<-table(a$T,a$cluster)
P_T<-matrix(0,4,4,dimnames = list(c("T1","T2","T3","T4"),
                                  c("T1","T2","T3","T4")))
for (i in 1:4) {
  for (j in 1:4) {
    P_T[i,j]<-chisq.test(b[c(i,j),])$p.value
  }
}
P_T<-round(P_T,2)
P_T[4,4]<-1
pheatmap(P_T, cluster_rows = F,cluster_cols = F,
         display_numbers = TRUE,number_color = "blue",
         filename = "TCGA/clinical/T_chisq.pdf",
         height = 2,width = 4)
save(P_T,P_M,P_N,P_Stage,file = "TCGA/clinical/分期卡方检验结果.Rdata")
##############################
tide_score<-merge(kirc_tide,res_tree,by.x = "Sample",by.y = "sample")
tide_score$cluster<-as.character(tide_score$cluster)
library(ggplot2)
library(ggpubr)
my_comparisons<-list(c("1","2"),c("1","3"),c("2","3"))
p<-ggboxplot(tide_score,x="cluster",y="scale",
             title = "",ylab = "TIDE prediction score",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 3)
ggsave("TCGA/TIDE.pdf",p,height = 6,width = 5)


library(pheatmap)
library(RColorBrewer)
hubdata<-kirc_gene1[hubgene,]
hubdata<-hubdata[,res_tree$sample]
annodata3<-data.frame(row.names = res_tree$sample,
                      cluster=res_tree$cluster)
annodata3$cluster<-as.character(annodata3$cluster)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
n[n>2]<-2
n[n< -2]<- -2
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color,
         main = "",filename = "TCGA/hubgene热图.pdf",
         width = 6,height = 6)
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color)

hubkirc<-kirc_gene1[hubgene,]
library(ggplot2)
library(reshape2)
library(ggpubr)
hubkirc<-t(hubkirc)
hubkirc<-data.frame(hubkirc)
hubkirc<-hubkirc[res_tree$sample,]
hubkirc$cluster<-res_tree$cluster
hubkirc_melt<-melt(hubkirc,id.vars = "cluster")
hubkirc_melt$cluster<-as.character(hubkirc_melt$cluster)
p<-ggplot(data = hubkirc_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.2,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(hubkirc_melt$value)$stats[c(1, 5)])*1.05)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 8.5)
ggsave("TCGA/hubgene_箱型图.pdf",p,width=8,height=4)


immune_checkpoint<-c("ADORA2A","BTLA","BTNL2","CD160","CD200","CD200R1","CD7",
                     "CD244","CD226","CD27","CD274","CD276","CD28","CD40","CD40LG",
                     "CD44","CD48","CD70","CD80","CD86","CTLA4","HAVCR2",
                     "HHLA2","ICOS","ICOSLG","IDO1","IDO2","KIR3DL1","KLRG1","LAG3",
                     "LAIR1","LGALS9","NRP1","PDCD1","PDCD1LG2","TIGIT","TMIGD2",
                     "TNFRSF14","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8",
                     "TNFRSF9","TNFSF14","TNFSF15","TNFSF18","TNFSF4",
                     "TNFSF9","VSIR","VTCN1")
checkdata<-kirc_gene1[immune_checkpoint,]
checkdata<-na.omit(checkdata)
checkdata<-checkdata[,res_tree$sample]
annodata3<-data.frame(row.names = res_tree$sample,
                      cluster=res_tree$cluster)
annodata3$cluster<-as.character(annodata3$cluster)
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(checkdata)))
n[n>2]=2
n[n< -2]= -2
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color,
         legend_breaks = c(-2,0,2),
         legend_labels = c("-2","0","2"),
         main = "",filename = "TCGA/immune_checkpoint.pdf",
         width = 6,height = 10)

library(ggplot2)
library(reshape2)
library(ggpubr)
hubcheck<-t(checkdata)
hubcheck<-as.data.frame(hubcheck)
hubcheck$cluster<-annodata3$cluster
hub_melt<-melt(hubcheck,id.vars = "cluster")
p<-ggplot(data = hub_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +
  ylab("Expression") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=cluster), label = "p.signif")
ggsave("TCGA/immunecheckpoint_box.pdf",p,width = 12,height = 4)


ICD<-c("CALR","CXCL10","FPR1","HGP","IFNAR2",'IFNB1','IFNE',
       'IFNK','IFNW1','MET','P2RX7','TLR4','LRP1','EIF2A',
       'EIF2AK3','EIF2AK2','EIF2AK4',"EIF2AK1",'HMGB1',
       'ANXA1','PANX1','P2RY2','IFNA1','IFNA2','TLR3')
ICDdata<-kirc_gene1[ICD,]
ICDdata<-na.omit(ICDdata)
ICDdata<-ICDdata[,res_tree$sample]
annodata3<-data.frame(row.names = res_tree$sample,
                      cluster=res_tree$cluster)
annodata3$cluster<-as.character(annodata3$cluster)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(ICDdata)))
n[n>2]=2
n[n< -2]= -2
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color,
         legend_breaks = c(-2,0,2),
         legend_labels = c("-2","0","2"),
         main = "",filename = "TCGA/ICD.pdf",
         width = 6,height = 6)

library(ggplot2)
library(reshape2)
library(ggpubr)
hubICD<-t(ICDdata)
hubICD<-as.data.frame(hubICD)
hubICD$cluster<-annodata3$cluster
hub_melt<-melt(hubICD,id.vars = "cluster")
p<-ggplot(data = hub_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.3,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","red"))+
  xlab("") +
  ylab("Expression") +
  theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=cluster), label = "p.signif")
ggsave("TCGA/ICD_box.pdf",p,width = 7,height = 4)


load("E:/singlecall_immune/RCC/GSVA/TCGA/GSVA_result.Rdata")
###########
va_kirc<-va_kirc[,res_tree$sample]
res_tree$subtype<-ifelse(res_tree$cluster==1,"C1","C23")
library(limma)
design <- model.matrix(~0+factor(res_tree$subtype))
colnames(design)=levels(factor(res_tree$subtype))
rownames(design)=res_tree$sample
contrast.matrix<-makeContrasts("C1-C23",levels=design)
fit <- lmFit(va_kirc,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit) 
tempOutput = topTable(fit, coef=1, n=Inf)
#############
res_tree$cluster<-as.character(res_tree$cluster)
library(pheatmap)
library(RColorBrewer)
annodata3<-data.frame(row.names = res_tree$sample,cluster=res_tree$cluster)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n<-t(scale(t(va_kirc)))
n<-n[,res_tree$sample]
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color,
         main = "",
         filename = "TCGA/GSVA.pdf",
         width = 8,height = 8)

va_kirc_ave<-sapply(split(rownames(annodata3),annodata3$cluster),
                    function(cells) rowMeans(n[,cells]))
n<-t(scale(t(va_kirc_ave)))
pheatmap(n, #fontsize_row=3, 
         cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,filename = "TCGA/GSVA_ave.pdf",
         height = 8,width = 5)

g<-egmt@result$core_enrichment[8]
g<-str_split(g,"/",simplify = T)
g<-c(g)
expressionMat<-kirc_gene1[g,]
expressionMat<-na.omit(expressionMat)
g1<-colMeans(expressionMat)
g1<-as.data.frame(g1)
g1$cluster<-NA
g1<-g1[res_tree$sample,]
g1$cluster<-res_tree$cluster

library(ggplot2)
library(ggpubr)
p<-ggboxplot(g1,x="cluster",y="g1",
             title = "",ylab = "Expression",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  stat_compare_means()
ggsave("TCGA/糖酵解.pdf",p,width = 5,height = 5)


load("E:/singlecall_immune/RCC/cibersort/TCGA_cell_result.Rdata")
result_tree<-result[,1:22]
result_tree<-t(result_tree)
result_tree<-result_tree[,rownames(annodata3)]
library(pheatmap)
library(RColorBrewer)
n=t(scale(t(result_tree)))
n[n>2]=2
n[n< -2]= -2
n1<-n[,row.names(annodata3)[which(annodata3$cluster=="1")]]
p1<-pheatmap(n1,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
n2<-n[,row.names(annodata3)[which(annodata3$cluster=="2")]]
p2<-pheatmap(n2,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
n3<-n[,row.names(annodata3)[which(annodata3$cluster=="3")]]
p3<-pheatmap(n3,
             cluster_rows = F,cluster_cols = T,
             show_rownames = T,show_colnames = F)
r1<-colnames(n1)[p1$tree_col[["order"]]]
r2<-colnames(n2)[p2$tree_col[["order"]]]
r3<-colnames(n3)[p3$tree_col[["order"]]]
annodata2<-data.frame(row.names = c(r1,r2,r3),cluster=annodata3$cluster)
pheatmap(n[,c(r1,r2,r3)],show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation = annodata2,
         main = "",filename = "TCGA/immunecell.pdf",
         height = 4,width = 7)

immune_cluster<-result
immune_cluster<-as.data.frame(immune_cluster)
immune_cluster<-immune_cluster[rownames(annodata3),]
immune_cluster$cluster<-annodata3$cluster
immune_melt<-melt(immune_cluster[,c(1:22,26)],id.vars = "cluster")
p<-ggplot(data = immune_melt, aes(x = variable, y = value, fill = cluster))+
  geom_boxplot(width=0.2,outlier.shape = NA,
               position = position_dodge(0.5))+
  scale_fill_manual(values = c("#1c9e77","#d95f02","blue"))+
  xlab("") +ylab("score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  coord_cartesian(ylim = (boxplot.stats(immune_melt$value)$stats[c(1, 5)])*3)+
  stat_compare_means(aes(group=cluster), label = "p.signif",label.y = 0.37)
ggsave("TCGA/immunecell_箱型图.pdf",p,width=8,height=4)

####TMB
dir.create("TCGA/TMB")
library(readr)
library(stringr)
kirc<-read_tsv("E:/singlecall_immune/RCC/bulkdata/TCGA/TCGA-KIRC.mutect2_snv.tsv/TCGA-KIRC.mutect2_snv.tsv")
kirc<-as.data.frame(kirc)
colnames(kirc) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
kirc$Entrez_Gene_Id =1
kirc$Center ='ucsc'
kirc$NCBI_Build ='GRCh38'
kirc$NCBI_Build ='GRCh38'
kirc$Strand ='+'
kirc$Variant_Classification = kirc$effect
tail(sort(table(kirc$Variant_Classification )))
kirc$Tumor_Seq_Allele1 = kirc$Reference_Allele
kirc$Variant_Type = ifelse(
  kirc$Reference_Allele %in% c('A','C','T','G') & kirc$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(kirc$Variant_Type )
library(maftools)
tcga.kirc = read.maf(maf = kirc,
                     vc_nonSyn=names(tail(sort(table(kirc$Variant_Classification )))))
save(tcga.kirc,kirc,file = "TCGA/TMB/mutation.Rdata")
yb<-tcga.kirc@clinical.data$Tumor_Sample_Barcode
yb<-str_sub(yb,1,12)
yb<-gsub("-",".",yb)
yb_feature<-clin_feature[clin_feature$sample%in%yb,]
tcga.kirc@clinical.data$cluster<-NA
tcga.kirc@clinical.data$sample<-gsub("-",".",str_sub(tcga.kirc@clinical.data$Tumor_Sample_Barcode,1,12))
for (i in 1:332) {
  tcga.kirc@clinical.data$cluster<-ifelse(tcga.kirc@clinical.data$sample==yb_feature$sample[i],
                                          yb_feature$cluster[i],tcga.kirc@clinical.data$cluster)
}
pdf("TCGA/TMB/All.oncoplot.clin.pdf", width=8, height=5)
oncoplot(maf = tcga.kirc, top = 20,clinicalFeatures = "cluster",sortByAnnotation = T)
dev.off()

pdf("TCGA/TMB/All.oncoplot.hubgene.pdf", width=8, height=6)
oncoplot(maf = tcga.kirc, top = 20,clinicalFeatures = "cluster",
         sortByAnnotation = T,genes = hubgene)
dev.off()

library(dplyr)
total<-table(kirc$Tumor_Sample_Barcode)
total<-as.data.frame(total)
colnames(total)<-c("sample","mutation ")
total$sample<-str_sub(total$sample,1,12)
total$sample<-gsub("-",".",total$sample)
total<-merge(total,res_tree,id.x="sample",id.y="sample")
library(ggplot2)
library(ggpubr)
p<-ggboxplot(total,x="cluster",y="mutation ",
             title = "",ylab = "Number of Mut Genes",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(total$`mutation `)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = 150)
ggsave("TCGA/TMB/gene_mutation.pdf",p,height = 5,width = 4)


use_cols <- c("Tumor_Sample_Barcode", "effect", "Hugo_Symbol","vaf")

mut_type <- c("5'UTR", "3'UTR", "3'Flank", "5'Flank", 
              "Intron", "IGR","Silent", "RNA", "Splice_Region")
df<-kirc[,use_cols]
data1 <- df %>% subset(!effect %in% mut_type)

word<-c("deletion","insertion","missense","nonsense","frameshift")
num<-c()
for(i in 1:5){
  a<-which(str_detect(data1$effect,word[i]))
  num<-c(num,a)
}
num<-num[!duplicated(num)]
data2<-data1[num,]
data <- data2 %>%
  mutate(vaf = vaf) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))

data$Tumor_Sample_Barcode<-str_sub(data$Tumor_Sample_Barcode,1,12)
data$Tumor_Sample_Barcode<-gsub("-",".",data$Tumor_Sample_Barcode)
colnames(data)[1]<-"sample"
data<-merge(data,res_tree,id.x="sample",id.y="sample")
library(ggplot2)
library(ggpubr)
p<-ggboxplot(data,x="cluster",y="TMB",
             title = "",ylab = "TMB",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(data$TMB)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = 3)
ggsave("TCGA/TMB/TMB.pdf",p,height = 5,width = 4)
save(kirc,data,total,file = "E:/singlecall_immune/RCC/TMB/TCGA_TMB.Rdata")


#RECA
library(ConsensusClusterPlus)
hubdata<-reca_cancer[hubgene,]
hubdata<-as.matrix(hubdata)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 10,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "E:/singlecall_immune/RCC/consensus/summary4/RECA",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "pdf")

res_tree<-results[[3]]$consensusClass
res_tree<-data.frame(res_tree)
res_tree<-data.frame(sample=row.names(res_tree),cluster=res_tree$res_tree)
res_tree<-res_tree[order(res_tree$cluster),]
save(results,file = "RECA/result.Rdata")

recaclin<-merge(clinical[,c(1,6,17)],res_tree,by.x = "icgc_donor_id",by.y = "sample")
colnames(recaclin)[1:3]<-c("patient_id","vital_status","time")
recaclin$status<-ifelse(recaclin$vital_status=="alive",0,1)
fit_im<-survfit(Surv(time,status)~cluster,data = recaclin)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              ggtheme = theme_bw(), 
              title="Overall Survival")

tide_score<-merge(reca_tide,res_tree,by.x = "Sample",by.y = "sample")
tide_score$cluster<-as.character(tide_score$cluster)
library(ggplot2)
library(ggpubr)
my_comparisons<-list(c("1","2"),c("1","3"),c("2","3"))
p<-ggboxplot(tide_score,x="cluster",y="scale",
             title = "",ylab = "TIDE prediction score",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = -2)
ggsave("RECA/TIDE.pdf",p,height = 6,width = 5)

#GSE167573
setwd("E:/singlecall_immune/RCC/consensus/summary4")
hubgene<-read.table("hubgene.txt",header = F)
hubgene<-hubgene$V1
source("E:/singlecall_immune/RCC/TIDE/TIDE.R")
library(openxlsx)
gse167573_clin<-read.xlsx("GSE167573/GSE167573os.xlsx")
clinical<-gse167573_clin
rownames(clinical)<-gse167573_clin$sample
clinical<-clinical[,-1]
gse167573_expr<-log2(gse167573_expr+1)
gse167573_expr1<-gse167573_expr[,gse167573_clin$sample]
exclusion<-read.table("E:/ketizu/毕设plus/TIDE/exclusion1.txt", sep='\t', header=T, check.names=F, quote=NULL)
gse167573_tide<-TIDE(gse167573_expr1,clinical,exclusion)

library(ConsensusClusterPlus)
hubdata<-gse167573_expr1[hubgene,]
hubdata<-as.matrix(hubdata)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 10,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "GSE167573/",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "pdf")

res_tree<-results[[3]]$consensusClass
res_tree<-data.frame(res_tree)
res_tree<-data.frame(sample=row.names(res_tree),cluster=res_tree$res_tree)
res_tree<-res_tree[order(res_tree$cluster),]
save(results,file = "GSE167573/result.Rdata")

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
gseclin<-merge(gse167573_clin,res_tree,by = "sample")
gseclin$cluster<-as.character(gseclin$cluster)
fit_im<-survfit(Surv(time,status)~cluster,data = gseclin)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = T,
              risk.table = TRUE,
              risk.table.col = "strata",
              linetype = "strata",
              ggtheme = theme_bw(),
              legend = c(0.8,0.75), 
              legend.title = "", 
              legend.labs = c("1", "2","3"),
              title="Overall Survival")
library(pheatmap)
library(RColorBrewer)
hubdata<-hubdata[,res_tree$sample]
annodata3<-data.frame(row.names = res_tree$sample,
                      cluster=res_tree$cluster)
annodata3$cluster<-as.character(annodata3$cluster)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
n[n>4]<-4
n[n< -4]<- -4
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata3,
         color=color,
         main = "",filename = "GSE167573/hubgene热图.pdf",
         width = 6,height = 6)

gsetide<-gse167573_tide[[5]]
tide_score<-merge(gsetide,res_tree,by.x = "Sample",by.y = "sample")
tide_score$cluster<-as.character(tide_score$cluster)
library(ggplot2)
library(ggpubr)
my_comparisons<-list(c("1","2"),c("1","3"),c("2","3"))
p<-ggboxplot(tide_score,x="cluster",y="scale",
             title = "",ylab = "TIDE prediction score",
             color = "cluster",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = -4)
ggsave("GSE167573/TIDE.pdf",p,height = 6,width = 5)

hubdata_ave<-sapply(split(rownames(annodata3),annodata3$cluster),
                    function(cells) rowMeans(hubdata[,cells]))
hubdata_ave<-as.data.frame(hubdata_ave)
n<-t(scale(t(hubdata_ave)))
pheatmap(n,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,height = 5,width = 2.5,
         filename = "GSE167573/hubgene_ave.pdf")
save(gse167573_tide,file = "GSE167573/tide.Rdata")
#PMID
hubgene<-read.table("hubgene.txt",header = F)
hubgene<-hubgene$V1
pmidos<-read.xlsx("PMID32472114/PMID32472114os.xlsx")
pmidIR<-read.xlsx("PMID32472114/PMID32472114RNAclinic.xlsx")
pmidIR<-na.omit(pmidIR)
pmidos<-na.omit(pmidos)
library(ConsensusClusterPlus)
hubdata<-pmid[hubgene,]
hubdata<-as.matrix(hubdata)
d = sweep(hubdata,1, apply(hubdata,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK = 10,reps = 100,
                              pItem = 0.8,pFeature = 1,
                              title = "PMID32472114/",
                              clusterAlg = "hc",distance = 'pearson',
                              seed = 1262118388.71279,plot = "pdf")
pmid_tree<-results[[3]]$consensusClass
pmid_tree<-data.frame(pmid_tree)
pmid_tree<-data.frame(sample=row.names(pmid_tree),cluster=pmid_tree$pmid_tree)
pmid_tree<-pmid_tree[order(pmid_tree$cluster),]
save(results,file = "PMID32472114/result.Rdata")

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
pmidclin<-merge(pmidos,pmid_tree,by.x = "RNA_ID",by.y = "sample")
pmidclin$cluster<-as.character(pmidclin$cluster)
fit_im<-survfit(Surv(OS,OS_CNSR)~cluster,data = pmidclin)
fit_im<-survfit(Surv(PFS,PFS_CNSR)~cluster,data = pmidclin)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = T,
              risk.table = TRUE,
              risk.table.col = "strata",
              linetype = "strata",
              ggtheme = theme_bw(),
              legend = c(0.8,0.75), 
              legend.title = "", 
              legend.labs = c("1", "2","3"),
              title="Overall Survival")
pmidIR[grep("PR",pmidIR$ORR),"ORR"]<-"CR/PR"
pmidIR[grep("CR",pmidIR$ORR),"ORR"]<-"CR/PR"
pmidIR[grep("D",pmidIR$ORR),"ORR"]<-"PD/SD"
pmidIR<-merge(pmidIR,pmid_tree,by.x = "RNA_ID",by.y = "sample")

library(pheatmap)
library(RColorBrewer)
hubdata<-hubdata[,pmid_tree$sample]
annodata<-data.frame(row.names = pmid_tree$sample,
                      cluster=pmid_tree$cluster)
annodata$cluster<-as.character(annodata$cluster)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n=t(scale(t(hubdata)))
n[n>4]<-4
n[n< -4]<- -4
pheatmap(n,
         show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,
         color=color,
         main = "",filename = "PMID32472114/hubgene热图.pdf",
         width = 6,height = 6)
hubdata_ave<-sapply(split(rownames(annodata),annodata$cluster),
                    function(cells) rowMeans(hubdata[,cells]))
hubdata_ave<-as.data.frame(hubdata_ave)
n<-t(scale(t(hubdata_ave)))
pheatmap(n,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("blue","white","red"))(100),
         border_color = NA,height = 5,width = 2.5,
         filename = "PMID32472114/hubgene_ave.pdf")
pmidIR<-pmidIR[which(pmidIR$ORR!="NE"),]
pmidIR$num<-1
pmidIR$cluster<-as.character(pmidIR$cluster)
library(ggplot2)
library(reshape2)
library(ggthemes)
p<-ggplot(pmidIR,aes(cluster,num,fill=ORR))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("PMID32472114/ORR_cluster1.pdf",p,height = 4,width = 4)
p<-ggplot(pmidIR,aes(ORR,num,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab("percentage")+xlab("")+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.2,'cm'),
        legend.text = element_text(size = 8))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()
ggsave("PMID32472114/ORR_cluster2.pdf",p,height = 4,width = 3)

table(pmidIR$ORR,pmidIR$cluster)
chisq.test(table(pmidIR$ORR,pmidIR$cluster))
#####
tcga_tree<-res_tree
hubtcga<-kirc_gene1[hubgene,]
hubtcga_ave<-sapply(split(tcga_tree$sample,tcga_tree$cluster),
                    function(cells) rowMeans(hubtcga[,cells]))

hubgse<-gse167573_expr1[hubgene,]
hubgse_ave<-sapply(split(rownames(annodata3),annodata3$cluster),
                    function(cells) rowMeans(hubgse[,cells]))

hubpmid<-pmid[hubgene,]
hubpmid_ave<-sapply(split(pmid_tree$sample,pmid_tree$cluster),
                    function(cells) rowMeans(hubpmid[,cells]))

tcga_gse<-data.frame(matrix(0,3,3,dimnames = list(c("T1","T2","T3"),
                                                  c("G1","G2","G3"))))
for(i in 1:3){
  for(j in 1:3){
    tcga_gse[i,j]<-paste(cor(hubtcga_ave[,i],hubgse_ave[,j],method = "spearman"),";",
                         cor.test(hubtcga_ave[,i],hubgse_ave[,j],method = "spearman")$p.value)
  }
}
tcga_pmid<-data.frame(matrix(0,3,3,dimnames = list(c("T1","T2","T3"),
                                                   c("P1","P2","P3"))))
for(i in 1:3){
  for(j in 1:3){
    tcga_pmid[i,j]<-paste(cor(hubtcga_ave[,i],hubpmid_ave[,j],method = "spearman"),";",
                          cor.test(hubtcga_ave[,i],hubpmid_ave[,j],method = "spearman")$p.value)
  }
}
