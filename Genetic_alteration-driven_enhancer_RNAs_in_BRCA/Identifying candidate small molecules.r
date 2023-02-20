################################################################
#   Prediction of drug-eRNAs pairs based on GSEA               #
################################################################
##############################################################################
#     Recognition of enhancer target genes based on Spearman correlation     #
##############################################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt",sep="\t",header=T)
load("/pub5/xiaoyun/Jobs/J4/7.TCGA count/1.lncRNApcg/2.tpm/OSF_BRCA_All_tpm.RData")
ls()
all_TCGA <- data.frame(pcg)
t_1 <- all_TCGA[,grep("11A",colnames(all_TCGA))]
t_2 <- all_TCGA[,grep("11B",colnames(all_TCGA))]
t_3 <- all_TCGA[,grep("01A",colnames(all_TCGA))]
t_4 <- all_TCGA[,grep("01B",colnames(all_TCGA))]
t_5 <- t_4[,setdiff(substring(colnames(t_4),1,12),substring(colnames(t_3),1,12))]
t_p <- cbind(t_1,t_2,t_3,t_5)
colnames(t_p)<- substring(colnames(t_p),1,16)
colnames(t_p) <- gsub(".11A","_normal",colnames(t_p))
colnames(t_p) <- gsub(".11B","_normal",colnames(t_p))
colnames(t_p) <- gsub(".01A","_tumor",colnames(t_p))
colnames(t_p) <- gsub(".01B","_tumor",colnames(t_p))
setdiff(colnames(enhancer_profile),colnames(t_p))
"TCGA.AR.A0U1" 
setdiff(colnames(t_p),colnames(enhancer_profile))
"TCGA.E2.A14U"
#enhancer_p was the expression profile of 1207 samples containing 5441 enhancers, and tcga_p was the expression profile of 1207 samples containing 19655 genes
enhancer_p <- enhancer_profile[,-which(colnames(enhancer_profile)=="TCGA.AR.A0U1_tumor")]
tcga_p <- t_p[,-which(colnames(t_p)=="TCGA.E2.A14U_tumor")]
e_p <- enhancer_p[,intersect(colnames(tcga_p),colnames(enhancer_p))]
t_p <- tcga_p[,intersect(colnames(tcga_p),colnames(enhancer_p))]
save(e_p,file="e_p.RData")
save(t_p,file="t_p.RData")
t_p<- as.matrix(t_p)
m <- t(t_p)
enhancer_names <- rownames(e_p)
library(psych)
n <- t(e_p)
#Correlation analysis
for(i in 1:5441){
   cor_res_1 <- corr.test(m,n[,i],method="spearman")
   cor_res <- cbind(cor_res_1$r,cor_res_1$p)
   colnames(cor_res) <- c(enhancer_names[i],"p")
   filenamecsv <- paste("cor",i,".csv",sep="")
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/all_cor")
   write.csv(cor_res,file=filenamecsv)
}
dir_Dgene<-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/all_cor"
CORenhgene <-list.files(dir_Dgene, full.names=T, recursive =T)
r_p <- data.frame(ENSG=NA,cor=NA,cor_p=NA,enhancer=NA)
for(p in 1:length(CORenhgene))
{
  Dsample <- read.table(CORenhgene[p], header=T, sep=",", as.is=T, quote="")
  Dsample$enhancer <- colnames(Dsample)[2]
  colnames(Dsample) <- c("ENSG","cor","cor_p","enhancer")
  r_p <- rbind(r_p,Dsample)
}
r_p <- r_p[-1,]
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对")
save(r_p,file="r_p.RData")
write.table(r_p,"r_P.txt",sep="\t",row.names=T,col.names=T)
#####################################################################
#  genetic alteration-driven eRNA-gene pairs（cor>0.3,fdr<0.0001）  #
#####################################################################
library(tidyverse)
all_enhancer_gene_cor <- separate(data=r_p,col=enhancer,into=c("no","chr","start","end",sep="[.]"))
all_enhancer_gene_cor <- all_enhancer_gene_cor[,-c(4,8)]
library(tidyr)
all_enhancer_gene_cor <- unite(all_enhancer_gene_cor, "enhancer", start, end, sep = "-")
all_enhancer_gene_cor <- unite(all_enhancer_gene_cor, "enhancer", chr, enhancer, sep = ":")
save(all_enhancer_gene_cor,file="all_enhancer_gene_cor.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/analysis_enhancer.RData")
enhancer_gene_cor <- all_enhancer_gene_cor[which(all_enhancer_gene_cor$enhancer %in% analysis_enhancer$enhancer),]
#p-value adjust，and cor>0.3,fdr<0.0001
enhancer_gene_cor <- enhancer_gene_cor[order(enhancer_gene_cor$cor_p),]  
enhancer_gene_cor$cor_fdr <-  p.adjust(as.numeric(enhancer_gene_cor$cor_p),method = "fdr",n=length(enhancer_gene_cor$cor_p))
enhancer_gene_cor_0.3 <- enhancer_gene_cor[which(as.numeric(enhancer_gene_cor$cor)>0.3& as.numeric(enhancer_gene_cor$cor_fdr<0.0001)),]
save(enhancer_gene_cor_0.3,file="enhancer_gene_cor_0.3.RData")
#########################################
#     Differentially expressed gene     #
#           Gene ID conversion          #
#########################################
#Gene ID conversion
library(org.Hs.eg.db)
t_p1 <- read.table("/pub6/temp/wangli/zcy/result/tcga_p.txt",sep="\t",header=T)
gene_ensgid <- rownames(t_p1)
gene_symple <- bitr(gene_ensgid, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
write.table(gene_symple,"/pub6/temp/wangli/zcy/result/gene_symple.txt",sep="\t")   #19293个
load("/pub5/xiaoyun/Jobs/J4/7.TCGA count/1.lncRNApcg/3.DEseq2/OSF_BRCA_DEseq2_Allres_.RData")
deseq <- DEseq2.res
library(DESeq2)
t_p1 <- read.table("/pub6/temp/wangli/zcy/result/tcga_p.txt",sep="\t",header=T)
d_gene <- data.frame(intersect(rownames(t_p1),rownames(deseq)))
deseq_gene <- deseq[intersect(rownames(t_p1),rownames(deseq)),]
#Differentially expressed gene
diff_gene_deseq2 <- subset(deseq_gene, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))\
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
diff_gene_deseq2_symple <- merge(diff_gene_deseq2,gene_symple,by.x="row.names",by.y="ENSEMBL")\
colnames(diff_gene_deseq2_symple)[1] <- "ENSEMBL"
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对")
save(diff_gene_deseq2_symple,file="diff_gene_deseq2_symple.RData")
#####################################################################
#   genetic alteration-driven eRNA-DEG pairs（cor>0.3,fdr<0.0001）  #
#####################################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对")
load("diff_gene_deseq2_symple.RData")
load("analysis_enhancer.RData")
load("all_enhancer_gene_cor.RData")
all_enhancer_gene_cor$ENSG <- gsub("[\"]","",all_enhancer_gene_cor$ENSG)
diffenhancer_diffgene_cor <- all_enhancer_gene_cor[which((all_enhancer_gene_cor$enhancer %in% analysis_enhancer$enhancer) & (all_enhancer_gene_cor$ENSG %in% diff_gene_deseq2_symple$ENSEMBL)),]
diffenhancer_diffgene_cor <- diffenhancer_diffgene_cor[order(diffenhancer_diffgene_cor$cor_p),]
diffenhancer_diffgene_cor$cor_fdr <-  p.adjust(as.numeric(diffenhancer_diffgene_cor$cor_p),method = "fdr",n=length(diffenhancer_diffgene_cor$cor_p))
diffenhancer_diffgene_cor_0.3 <- diffenhancer_diffgene_cor[which(as.numeric(diffenhancer_diffgene_cor$cor)>0.3& as.numeric(diffenhancer_diffgene_cor$cor_fdr<0.0001)),]  
save(diffenhancer_diffgene_cor_0.3,file="diffenhancer_diffgene_cor_0.3.RData")
##################################################
#              Process Cmap raw data             #
##################################################
library(R.oo)
library(R.methodsS3)
library(R.utils) 
library(BiocGenerics)  
library(parallel)  
library(Biobase)
library(affy)  
library(u133aaofav2cdf)
setwd("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据")
unzip("cmap_build02.volume1of7.zip",exdir = "解压zip/")   
unzip("cmap_build02.volume2of7.zip",exdir = "解压zip/")
unzip("cmap_build02.volume3of7.zip",exdir = "解压zip/")
unzip("cmap_build02.volume4of7.zip",exdir = "解压zip/")
unzip("cmap_build02.volume5of7.zip",exdir = "解压zip/")
unzip("cmap_build02.volume6of7.zip",exdir = "解压zip/")
unzip("cmap_build02.volume7of7.zip",exdir = "解压zip/")
bz2_file <- list.files(path="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压zip",pattern="*.bz2")
setwd("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压zip")
for(i in 1:7056){
    bunzip2(bz2_file[i])   
}
dir.files<-list.files(path = "/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压zip" ,pattern="*.CEL")
count<-length(dir.files)
dir.names <-dir.files
for(i in 1:7056){
        print(i);
        setwd("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压zip")
        AffyData<-ReadAffy(filenames = dir.files[i])
        txtnames <- paste(gsub(".CEL","",dir.names[i]),".txt",sep="")
        exprsSet.MAS5 <- mas5(AffyData)
		setwd("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压txt")
        write.exprs(exprsSet.MAS5, file=txtnames)
}
dir_cmap <-"/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/解压txt"
primary_cmap <-list.files(dir_cmap, full.names=T, recursive =T)
data_1 <- read.table(primary_cmap[1], header=T, sep="\t")
for(p in 2:length(primary_cmap))
{
  print(p);
  data_p <- read.table(primary_cmap[p], header=T, sep="\t")
  data_1 <- merge(data_1,data_p,by="X")
}
all_cmap <- data_1
save(all_cmap,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/all_cmap.RData")
write.table (colnames(all_cmap), file ="aacmap.txt",sep="\t")
instance <- read.table("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/cmap_instances.csv",sep=",",quote="\"",header=T)
colnames(all_cmap) <- gsub(".CEL","",colnames(all_cmap))
colnames(all_cmap) <- gsub("^X","",colnames(all_cmap))
write.table(all_cmap,"/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/file 1.csv",sep=",",row.names=T,col.names=T)
###########################################
#   Processing drug sample data--Cmap     #
###########################################
perturbation_sample <- data.frame(perturbation_scan_id=instance[,"perturbation_scan_id"],cmap_name=instance[,"cmap_name"],class="perturbation")
perturbation_sample$perturbation_scan_id <- sub("^'","",perturbation_sample$perturbation_scan_id)
perturbation_cmap <- all_cmap[,which(colnames(all_cmap) %in% perturbation_sample$perturbation_scan_id)]
rownames(perturbation_cmap) <- all_cmap[,1]
perturbation_pro_1 <- merge(t(perturbation_cmap),perturbation_sample,by.x="row.names",by.y="perturbation_scan_id") 
save(perturbation_pro_1,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/perturbation_pro_1.RData")
#The expression values of the same drug at different concentrations were averaged
my_drug <- factor(perturbation_pro_1$cmap_name)
perturbation_pro <- apply(perturbation_pro_1[,2:22269],2,function(x){tapply(as.numeric(x),my_drug,mean)})  
save(perturbation_pro,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/perturbation_pro.RData")
###########################################
#    Process normal sample data--Cmap     # 
###########################################
nonper_cmap <- all_cmap[,setdiff(colnames(all_cmap)[-1],perturbation_sample$perturbation_scan_id)] 
rownames(nonper_cmap) <- all_cmap[,1]
vehicle_sample <- data.frame(vehicle_scan_id=instance[,"vehicle_scan_id4"],perturbation_scan_id=instance[,"perturbation_scan_id"],cmap_name=instance[,"cmap_name"])
vehicle_sample$vehicle_scan_id <- sub("^[.]","",vehicle_sample$vehicle_scan_id)
vehicle_sample$vehicle_scan_id <- sub("^'","",vehicle_sample$vehicle_scan_id)
vehicle_sample$perturbation_scan_id <- sub("^'","",vehicle_sample$perturbation_scan_id)
vehicle_pro_1 <- merge(t(nonper_cmap),vehicle_sample,by.x="row.names",by.y="vehicle_scan_id") 
save(vehicle_pro_1,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro_1.RData")
vehicle_sample$per_scan_id <- sub("\\d$","",vehicle_sample$perturbation_scan_id) 
vehicle_sample2 <- vehicle_sample[-which(vehicle_sample$vehicle_scan_id %in% vehicle_pro_1$Row.names),] 
nonper_nonveh <- nonper_cmap[,setdiff(colnames(nonper_cmap),vehicle_pro_1[,1])]
nonper_nonveh <- as.data.frame(t(nonper_nonveh))
nonper_nonveh$pri_sample <- rownames(nonper_nonveh)
nonper_nonveh$sample <- sub("\\d$","",rownames(nonper_nonveh))
vehicle_pro_2 <- merge(nonper_nonveh,vehicle_sample2,by.x="sample",by.y="per_scan_id",all.y=T)  
save(vehicle_pro_2,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro_2.RData")
colnames(vehicle_pro_1)[1] <- "vehicle_scan_id"
vehicle_pro_3 <- unique(rbind(vehicle_pro_1[,c(2:22269,22271,1)],vehicle_pro_2[,c(2:22269,22273,22271)])) 
save(vehicle_pro_3,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro_3.RData")
vehicle_pro_4 <- vehicle_pro_3[complete.cases(vehicle_pro_3),]   #4986
vehicle_pro_5 <- vehicle_pro_3[-which(vehicle_pro_3[,1]!="NA"),] #1281
vehicle_scan_id <- factor(vehicle_pro_4$vehicle_scan_id)
vehicle_pro_6 <- apply(vehicle_pro_4[,1:22268],2,function(x){tapply(as.numeric(x),vehicle_scan_id,mean)})
save(vehicle_pro_6,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro_6.RData")
vehicle_pro_samplemean <- merge(vehicle_pro_6,vehicle_sample,by.x="row.names",by.y="vehicle_scan_id")
save(vehicle_pro_samplemean,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro_samplemean.RData")
#####Take an average of different samples of the same drug
my_drug <- factor(vehicle_pro_samplemean$cmap_name)
vehicle_pro <- apply(vehicle_pro_samplemean[,2:22269],2,function(x){tapply(as.numeric(x),my_drug,mean)})  
save(vehicle_pro,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro.RData")
############################################################
#    Order of difference before and after disturbance      # 
#              vehicle_pro and perturbation_pro            #
#                         log2(per/veh)                    #
############################################################ 
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/vehicle_pro.RData")
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/perturbation_pro.RData")
probe_1 <- read.table("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/GPL96-15653.txt",sep="\t",skip=16,header=T,fill=T,quote="")
probe_2 <- read.table("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/GPL3921-25447.txt",sep="\t",skip=16,header=T,fill=T,quote="")
probe_1 <- probe_1[,c("ID","Gene.Symbol")]  
probe_2 <- probe_2[,c("ID","Gene.Symbol")]   
probe_anno <- unique(rbind(probe_1,probe_2)) 
save(probe_anno,file="/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/probe_anno.RData")
#The probe ID of cmap was matched with the gene ID, and the probe without gene ID was deleted
probe_anno_na <- subset(probe_anno,Gene.Symbol !="")  
#Replace the probe in the expression profile with the Gene.Symbol
vehicle_pro <- t(vehicle_pro)
perturbation_pro <- t(perturbation_pro)
#Take the average of multiple probes corresponding to a gene
vehicle_pro_gene <- merge(vehicle_pro,probe_anno_na,by.x="row.names",by.y="ID")   
perturbation_pro_gene <- merge(perturbation_pro,probe_anno_na,by.x="row.names",by.y="ID")  #24447*1311
my_gene <- factor(vehicle_pro_gene$Gene.Symbol)
result_vehicle_pro <- apply(vehicle_pro_gene[,2:1310],2,function(x){tapply(as.numeric(x),my_gene,mean)}) #15334*1309
result_perturbation_pro <- apply(perturbation_pro_gene[,2:1310],2,function(x){tapply(as.numeric(x),my_gene,mean)})
setwd("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/")
save(result_vehicle_pro,file="result_vehicle_pro.RData")
save(result_perturbation_pro,file="result_perturbation_pro.RData")
############################################################ 
#                       log2(per/veh)                      #
############################################################
result_perturbation_pro <- result_perturbation_pro[order(rownames(result_perturbation_pro)),]
result_perturbation_pro <- result_perturbation_pro[,order(colnames(result_perturbation_pro))]
result_vehicle_pro <- result_vehicle_pro[order(rownames(result_vehicle_pro)),]
result_vehicle_pro <- result_vehicle_pro[,order(colnames(result_vehicle_pro))]
log2_per_veh <- matrix(data = 0,15334,1309)
for(i in 1:1309){
 print(i);
 log2_per_veh[,i]<-log2(as.numeric(result_perturbation_pro[,i])/as.numeric(result_vehicle_pro[,i]))
} 
rownames(log2_per_veh) <- rownames(result_perturbation_pro)
colnames(log2_per_veh) <- colnames(result_perturbation_pro)
save(log2_per_veh,file="log2_per_veh.RData")
foldchange_per_veh <- matrix(data = 0,15334,1309)
for(i in 1:1309){
 print(i);
 foldchange_per_veh[,i]<-as.numeric(result_perturbation_pro[,i])/as.numeric(result_vehicle_pro[,i])
} 
rownames(foldchange_per_veh) <- rownames(result_perturbation_pro)
colnames(foldchange_per_veh) <- colnames(result_perturbation_pro)
save(foldchange_per_veh,file="foldchange_per_veh.RData")
drug <- data.frame(colnames(foldchange_per_veh))
gene <- data.frame(rownames(foldchange_per_veh))
foldchange_drug_gene <- data.frame(drug_name=NA,gene=NA,foldchange=NA)
for(i in 1:1309){
  print(i);
  new_1 <- data.frame(drug[i,1],gene[,1],foldchange_per_veh[,i])
  colnames(new_1) <- c("drug_name","gene","foldchange")
  foldchange_drug_gene <- rbind(foldchange_drug_gene,new_1)
}
foldchange_drug_gene <- foldchange_drug_gene[-1,]
foldchange_drug_gene$foldchange=round(as.numeric(foldchange_drug_gene$foldchange),digits = 3)
save(foldchange_drug_gene,file="foldchange_drug_gene.RData")
###########################################################
#                  GSEA analysis，FDR<0.05                #
###########################################################
install.packages("stringi")
library(stringi)
install.packages("devtools") 
library(devtools)
setwd("D:/")
devtools::install_local("clusterProfiler-master.zip")
library(clusterProfiler)
install.packages("tidyverse") 
library(tidyverse)
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
diff_r_p_0.3 <- unique(diffenhancer_diffgene_cor[,c(4,6)]) 
colnames(diff_r_p_0.3) <- c("enhancer","gene")
save(diff_r_p_0.3,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diff_r_p_0.3.RData")
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/log2_per_veh.RData")
intersect(rownames(log2_per_veh),diff_r_p_0.3$gene)  
intersect(probe_1[,2],diff_r_p_0.3$gene) 
library(clusterProfiler)
library(tidyverse)
for(i in 1:1309){
  one <- log2_per_veh[,i]
  new <- sort(one,decreasing=T)	  
  gsea <- GSEA(new,TERM2GENE=diff_r_p_0.3,verbose=F,minGSSize=1,maxGSSize=1000,pvalueCutoff=1,pAdjustMethod="fdr",nPermSimple = 10000,eps=0)
  new_gsea <- gsea %>% as_tibble() %>% arrange(desc(NES))
  setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/result_log2_per_veh")
  filenamecsv=paste("gsea",i,".txt",sep="")
  while (nrow(new_gsea) > 0){
    drug <- colnames(log2_per_veh)[i]
    new_gsea1 <- cbind(drug,new_gsea)
    write.table (new_gsea1, file =filenamecsv,sep="\t")
    break;
  }
}
dir_Dgene<-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/result_log2_per_veh"
GSEAgene <-list.files(dir_Dgene, full.names=T, recursive =T)
Dset_0.3 <-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
p=1
Dsample<- read.table(GSEAgene[p], header=T, sep="\t")
t <- colnames(Dsample)
colnames(Dset_0.3) <- t
for(p in 1:length(GSEAgene))
{
  Dsample<- read.table(GSEAgene[p], header=T, sep="\t")
  Dset_0.3<-rbind(Dset_0.3,Dsample)
}
Dset_log2_per_veh <- Dset_0.3[-1,]
save(Dset_log2_per_veh,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/Dset_log2_per_veh.Rdata")
GSEA_enh_drug_final <- Dset_log2_per_veh[which(Dset_log2_per_veh$p.adjust<0.05),]
dim(GSEA_enh_drug_final)                 #100856
length(unique(GSEA_enh_drug_final$drug)) #1252
length(unique(GSEA_enh_drug_final$ID))   #238
range(GSEA_enh_drug_final$NES)           #-3.705169 3.265990
length(which(GSEA_enh_drug_final$NES>0)) #32582 
length(which(GSEA_enh_drug_final$NES<0)) #68274  
save(GSEA_enh_drug_final,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/GSEA_enh_drug_final.Rdata")
write.table(GSEA_enh_drug_final, file ="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/GSEA_enh_drug_final.txt",sep="\t")
######################################################################
#   Prediction of drug-enhancer relationship pairs based on network  #
#                       eRNA-drug-gene network                       #
######################################################################
############################################
#               drug-target network        #
#       DRUGBANK and DGIdb database        #
############################################
#drugbank data
DB_dt <- read.csv("/pub6/temp/wangli/zcy/network/drug_target/drugbank/all/all.csv",sep=",",header = T,fill=T)
DB_dt1 <- DB_dt[,c(3,12,13)]
DB_dt2 <- data.frame(DB_dt1[which(DB_dt1[,2]=="Humans"),])
DB_dt3 <- data.frame(DB_dt2[which(DB_dt2[,1]!=""),])
bb <- data.frame(NA,NA)
colnames(bb) <- c("gene","drug_id")
for(i in 1:3219){
  q <- as.character(DB_dt3[i,3])
  q1 <- c(unlist(strsplit(q,"[; ]")))
  qq <- data.frame(gene=DB_dt3[i,1],drug_id=q1)
  bb <- rbind(bb,qq)
}
bb1 <- bb[-1,]
bb2 <- data.frame(bb1[which(bb1[,2]!=""),])
write.table(bb2, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/drugbank_dt_human.csv",sep=",",row.names=F,col.names=T)
bb3 <- data.frame("DrugBank.ID"=bb2[,2],"gene"=bb2[,1])
drug <- unique(bb3[,1])   
DB_all_d <- read.csv("/pub6/temp/wangli/zcy/network/drug_target/drugbank/all/drug_links.csv",sep=",",header = T,fill=T)
DB_all_d1 <- DB_all_d[,1:2]
DB_all_dt <- merge(bb3,DB_all_d1,by="DrugBank.ID")
colnames(DB_all_dt) <-c("DrugBank_ID","gene","drug")
write.table(DB_all_dt, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_all_dt.csv",sep=",",row.names=F,col.names=T)
DB_approved_d <- read.csv("/pub6/temp/wangli/zcy/network/drug_target/drugbank/approve/all.csv",sep=",",header = T,fill=T)
DB_approved_dt1 <- DB_approved_d[,c(3,12,13)]
DB_approved_dt2 <- data.frame(DB_approved_dt1[which(DB_approved_dt1[,2]=="Humans"),])
DB_approved_dt3 <- data.frame(DB_approved_dt2[which(DB_approved_dt2[,1]!=""),])  ###2500对
aa <- data.frame(NA,NA)
colnames(aa) <- c("gene","drug_id")
for(i in 1:2500){
  q <- as.character(DB_approved_dt3[i,3])
  q1 <- c(unlist(strsplit(q,"[; ]")))
  qq <- data.frame(gene=DB_approved_dt3[i,1],drug_id=q1)
  aa <- rbind(aa,qq)
}
aa1 <- aa[-1,]
aa2 <- data.frame(aa1[which(aa1[,2]!=""),])  ###9942
write.table(aa2, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/drugbank_approved_dt_human.csv",sep=",",row.names=F,col.names=T)
aa3 <- data.frame("DrugBank.ID"=aa2[,2],"gene"=aa2[,1])
drug <- unique(aa3[,1])   ###1900
gene <- unique(aa3[,2])   ###2155
DB_approved_all_dt <- merge(aa3,DB_all_d1,by="DrugBank.ID")
colnames(DB_approved_all_dt) <- c("DrugBank_ID","gene","drug")
write.table(DB_approved_all_dt, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_approved_all_dt.csv",sep=",",row.names=F,col.names=T)
DB_all_dt <- as.matrix(read.table("/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_all_dt.csv",sep=",",header = T))
DB_approved_all_dt <- as.matrix(read.table("/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_approved_all_dt.csv",sep=",",header = T))
c_drug <- as.matrix(read.table("/pub6/temp/wangli/zcy/network/drug-cmap.txt",sep="\t",quote="",fill=T,header = F))
colnames(c_drug) <- c("error_drug","drug")
DB_all_dt[,3] <- tolower(DB_all_dt[,3])
DB_approved_all_dt[,3] <- tolower(DB_approved_all_dt[,3])
c_drug[,2] <- tolower(c_drug[,2])
write.table(c_drug, file ="/pub6/temp/wangli/zcy/network/c_drug.txt",sep="\t",row.names=F,col.names=T)
save(c_drug,file="/pub6/temp/wangli/zcy/network/c_drug.Rdata")
write.table(DB_all_dt, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_all_dt.csv",sep=",",row.names=F,col.names=T)
save(DB_all_dt,file="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_all_dt.Rdata")
write.table(DB_approved_all_dt, file ="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_approved_all_dt.csv",sep=",",row.names=F,col.names=T)
save(DB_approved_all_dt,file="/pub6/temp/wangli/zcy/network/drug_target/drugbank/result/DB_approved_all_dt.Rdata")
a <- intersect(c_drug[,2],DB_all_dt[,3])  
match_drug_all <- merge(c_drug,DB_all_dt,by="drug")
length(unique(match_drug_all[,1]))
a <- intersect(c_drug[,2],DB_approved_all_dt[,3]) 
#Drug target relationship pairs in DGIdb database were read, and relationship pairs without gene symple were removed
DGIdb <- read.table("/pub6/temp/wangli/zcy/network/drug_target/DGIdb/DGIdb_interaction.txt",sep="\t",quote="",fill=T,header=T)
DGIdb_1 <- DGIdb[,c(1,2,5,7,8)]         ##85832
DG <- DGIdb_1[which(DGIdb_1[,1]!=""),]  ##82706
DG1 <- DG[which(DG[,5]!=""),]           ##58995
DG2 <- as.matrix(DG1[,c(1,5)])
DG2[,2] <- tolower(DG2[,2]) 
DG2[,2] <-iconv(DG2[,2],"WINDOWS-1252","UTF-8")  
colnames(DG2) <- c("gene","drug")
length(unique(DG2[,2]))  ###10691
length(unique(DG2[,1]))  ###3231
write.table(DG2,"/pub6/temp/wangli/zcy/network/drug_target/DGIdb/result/DGIdb_drug_gene.txt",sep="\t",row.names=F,col.names=T)
save(DG2,file="/pub6/temp/wangli/zcy/network/drug_target/DGIdb/result/DGIdb_drug_gene.Rdata")
two_database <- union(DG2[,2],DB_all_dt[,3]) ##14239
length(unique(two_database))
a <- intersect(c_drug[,2],two_database)  ###850
#DGIdb--drug_claim_primary_name
DG_1<- DG[which(DG[,4]!=""),]           ##82706
DG_2 <- as.matrix(DG_1[,c(1,4)])
DG_2[,2] <- tolower(DG_2[,2]) 
DG_2[,2] <-iconv(DG_2[,2],"WINDOWS-1252","UTF-8")  
DG_2[,2] <- tolower(DG_2[,2]) 
colnames(DG_2) <- c("gene","drug")
length(unique(DG_2[,2]))  ###25387
length(unique(DG_2[,1]))  ###3953
write.table(DG_2,"/pub6/temp/wangli/zcy/network/drug_target/DGIdb/result/DGIdb_nndrug_gene.txt",sep="\t",row.names=F,col.names=T)
save(DG_2,file="/pub6/temp/wangli/zcy/network/drug_target/DGIdb/result/DGIdb_nndrug_gene.Rdata")
#Finally, the drug name in the drug claim primary name column is selected
a <- intersect(DG2[,2],DG_2[,2]) 
D <- union(DG2[,2],DG_2[,2])
two_database <- union(D,DB_all_dt[,3])
b <- intersect(c_drug[,2],two_database)  
#DGIdb and drugbank databases were merged
two_dd <- rbind(DB_all_dt[,2:3],DG_2)
two_dd_1 <- unique(two_dd)
length(unique(two_dd_1[,1]))   ##4903
length(unique(two_dd_1[,2]))   ##28712
write.table(two_dd_1,"/pub6/temp/wangli/zcy/network/drug_target/two_dd.txt",sep="\t",row.names=F,col.names=T)
save(two_dd_1,file="/pub6/temp/wangli/zcy/network/drug_target/two_dd.Rdata")
################################
#       PPI----HINT database   #
################################
#HINT(http://hint.yulab.org/)
##high-quality interactions
HomoSapiens_binary_hq <- read.table("/pub6/temp/wangli/zcy/network/ppi/text/HomoSapiens_binary_hq.txt",sep="\t",header = T)
HomoSapiens_cocomp_hq <- read.table("/pub6/temp/wangli/zcy/network/ppi/text/HomoSapiens_cocomp_hq.txt",sep="\t",header = T)
#After evaluation, these two data are added and de-duplicated to be HomoSapiens_binary_hq
htb_hq <- read.table("/pub6/temp/wangli/zcy/network/ppi/text/HomoSapiens_htb_hq.txt",sep="\t",header = T)
lcb_hq <- read.table("/pub6/temp/wangli/zcy/network/ppi/text/HomoSapiens_lcb_hq.txt",sep="\t",header = T)
binary <- unique(HomoSapiens_binary_hq[,3:4])
cocomp <- unique(HomoSapiens_cocomp_hq[,3:4])
ppi_HINT <- rbind(binary,cocomp)
ppi_HINT <- unique(ppi_HINT)
write.table(binary,"/pub6/temp/wangli/zcy/network/ppi/result/binary.csv",sep=",",row.names=F,col.names=T)
write.table(cocomp,"/pub6/temp/wangli/zcy/network/ppi/result/cocomp.csv",sep=",",row.names=F,col.names=T)
write.table(ppi_HINT,"/pub6/temp/wangli/zcy/network/ppi/result/ppi_HINT.csv",sep=",",row.names=F,col.names=T)
#Remove data from the loop
ppi_HINT_nloop <- subset(ppi_HINT,ppi_HINT[,1]!=ppi_HINT[,2])
write.table(ppi_HINT_nloop,"/pub6/temp/wangli/zcy/network/ppi/result/ppi_HINT_nloop.txt",sep="\t",row.names=F,col.names=T)
double_gene_1 <- ppi_HINT[grep(pattern = "[|]",ppi_HINT[,1]),]
double_gene_1 <- data.frame(Gene_A=double_gene_1[,2],Gene_B=double_gene_1[,1])
double_gene_2 <- ppi_HINT[grep(pattern = "[|]",ppi_HINT[,2]),]
double_g <- rbind(double_gene_1,double_gene_2)
write.table(double_g,"/pub6/temp/wangli/zcy/network/ppi/result/double_g.txt",sep="\t",row.names =F)
double_dg <- double_g[grep(pattern = "[|]",double_g[,1]),]
write.table(double_dg,"/pub6/temp/wangli/zcy/network/ppi/result/double_dg.txt",sep="\t",row.names =F)
#The data is converted to GMT format(clusterProfiler--read.gmt)
 a <- read.gmt("/pub6/temp/wangli/zcy/network/ppi/text/dg.gmt")
 a1 <- a[-which(a[,2]==""),]  
 a2 <- data.frame(a1[,2],a1[,1])  
 colnames(a2)=c("Gene_A","Gene_B")
write.table(a2,"/pub6/temp/wangli/zcy/network/ppi/result/dg_2.txt",sep="\t",row.names =F)
double_ndg <- double_g[-grep(pattern = "[|]",double_g[,1]),]
d2_g <- unique(rbind(a2,double_ndg))
write.table(d2_g,"/pub6/temp/wangli/zcy/network/ppi/result/d2_g.txt",sep="\t",row.names =F)
b <- read.gmt("/pub6/temp/wangli/zcy/network/ppi/text/d2_g.gmt")
b1 <- b[-which(b[,2]==""),]  
colnames(b1)=c("Gene_A","Gene_B")
s_ppi_nloop <- subset(b1,b1[,1]!=b1[,2])
colnames(s_ppi_nloop)=c("Gene_A","Gene_B") 
write.table(s_ppi_nloop,"/pub6/temp/wangli/zcy/network/ppi/result/s_ppi_nloop.txt",sep="\t",row.names =F)
ss_1 <- ppi_HINT_nloop[-grep(pattern = "[|]",ppi_HINT_nloop[,1]),]
ss_2 <- ss_1[-grep(pattern = "[|]",ss_1[,2]),]
ss_2 <- unique(ss_2)
write.table(ss_2,"/pub6/temp/wangli/zcy/network/ppi/result/ss_2.txt",sep="\t",row.names =F)
final_HINT <- rbind(ss_2,s_ppi_nloop)
final_HINT_1 <- unique(final_HINT)
write.table(final_HINT_1,"/pub6/temp/wangli/zcy/network/ppi/result/final_HINT.txt",sep="\t",row.names =F)

#########################################
#         eRNA-drug-gene network        #
#########################################
ppi_HINT <- read.table("/pub6/temp/wangli/zcy/network/ppi/result/final_HINT.txt",sep="\t",header=T)
colnames(ppi_HINT) <- c("node1","node2")
load("/pub6/temp/wangli/zcy/network/drug_target/two_dd.Rdata")
dg <- two_dd_1
colnames(ppi_HINT) <- c("node1","node2")
dg_1 <- data.frame(dg[,2],dg[,1])
colnames(dg_1) <- c("node1","node2")
net <- rbind(ppi_HINT,dg_1) 
net_1 <- data.frame(node1=dg_1[,1],class="+",node2=dg_1[,2])
net_2 <- data.frame(node1=ppi_HINT[,1],class="-",node2=ppi_HINT[,2])
dgg_network <- rbind(net_1,net_2)
write.table(dgg_network,"/pub6/temp/wangli/zcy/network/result/dgg_network.txt",sep="\t",row.names =F,col.names=T)
save(dgg_network,file="/pub6/temp/wangli/zcy/network/result/dgg_network.Rdata")
gggg <- data.frame(a=unique(dgg_network[,3]))
gggg1 <- data.frame(a=unique(ppi_HINT[,1]))
ll <- rbind(gggg,gggg1)
length(unique(ll[,1]))
length(unique(dg[,2]))
########################################################################
#                       network analysis----SPL                        #
########################################################################
library(BiRewire)
library(igraph)
#drug-gene-gene network
load("/pub6/temp/wangli/zcy/network/result/dgg_network.Rdata")
dgg=birewire.induced.bipartite(dgg_network) 
#eRNA-gene
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
diff_r_p_0.3 <- unique(diffenhancer_diffgene_cor[,c(4,6)]) 
colnames(diff_r_p_0.3) <- c("enhancer","gene")
node1 <- data.frame(from=diff_r_p_0.3[,1],to=diff_r_p_0.3[,2])  
all_diff_enhancer <- data.frame(unique(diff_r_p_0.3[,1]))
load("/pub6/temp/wangli/zcy/network/drug_target/two_dd.Rdata")
dg <- two_dd_1
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/随机网络")
for(i in 1:500){
   print(i);
   new_deg <- birewire.rewire.dsg(dgg,verbose=FALSE)  
   #random network
   r_net <- birewire.build.dsg(new_deg,delimitators=list(negative='-',positive='+'))
   node <- data.frame(from=r_net[,1],to=r_net[,3])
   n <- rbind(node,node1)    
   r_g <- graph_from_data_frame(n,directed=FALSE)
   r_all_path <- shortest.paths(r_g,mode="all")
   r_drug_path <- r_all_path[which(rownames(r_all_path) %in% dg[,2]),]
   diff_r_drug_enhancer_path <- r_drug_path[,which(colnames(r_drug_path) %in% all_diff_enhancer[,1])] 
   filenametxt=paste("diff_r_path",i,".txt",sep="")
   write.table (diff_r_drug_enhancer_path, file =filenametxt,sep="\t")   
}
####################################################
# Shortest path of eRNA- drug pair in real network #
####################################################
library(BiRewire)
library(igraph)
load("/pub6/temp/wangli/zcy/network/result/dgg_network.Rdata")
dgg_network <- dgg_network
node <- data.frame(from=dgg_network[,1],to=dgg_network[,3])   
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
diff_r_p_0.3 <- unique(diffenhancer_diffgene_cor[,c(4,6)]) 
colnames(diff_r_p_0.3) <- c("enhancer","gene")
node1 <- data.frame(from=diff_r_p_0.3[,1],to=diff_r_p_0.3[,2])
diff_n <- rbind(node,node1)        ######362272
diff_nn <- diff_n[which(diff_n[,1]!="<NA>"),]    ###### 362270
diff_g <- graph_from_data_frame(diff_nn,directed=FALSE)
all_diff_enhancer <- data.frame(unique(diff_r_p_0.3[,1]))
load("/pub6/temp/wangli/zcy/network/drug_target/two_dd.Rdata")
dg <- two_dd_1
diff_all_path <- shortest.paths(diff_g,mode="all")   
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/")
save(diff_all_path,file="diff_all_path.Rdata")  ######47217*47217
diff_drug_path <- diff_all_path[which(rownames(diff_all_path) %in% dg[,2]),]
diff_drug_enhancer_path <- diff_drug_path[,which(colnames(diff_drug_path) %in% all_diff_enhancer[,1])]    #########28711*259
save(diff_drug_enhancer_path,file="diff_drug_enhancer_path.Rdata")
write.table (diff_drug_enhancer_path, file ="diff_drug_enhancer_path.txt",sep="\t")
load("diff_drug_enhancer_path.Rdata") 
load("/pub6/temp/wangli/zcy/network/c_drug.Rdata")
colnames(diff_drug_enhancer_path) <- gsub("X.","",colnames(diff_drug_enhancer_path))
diff_my_drug_enhancer_path <- diff_drug_enhancer_path[which(rownames(diff_drug_enhancer_path) %in% c_drug[,2]),]
save(diff_my_drug_enhancer_path,file="diff_my_drug_enhancer_path.Rdata")
###################################################################
#      standardized shortest path length of eRNA-drug pair        #
###################################################################
diff_drug_enhancer_path <- read.table("diff_drug_enhancer_path.txt",header=T, sep="\t",quote="")
colnames(diff_drug_enhancer_path) <- gsub("X.","",colnames(diff_drug_enhancer_path))
diff_drug_enhancer_path_0.3 <- diff_drug_enhancer_path[order(rownames(diff_drug_enhancer_path)),]
diff_drug_enhancer_path_0.3 <- diff_drug_enhancer_path_0.3[,order(colnames(diff_drug_enhancer_path_0.3))]
save(diff_drug_enhancer_path_0.3,file="diff_drug_enhancer_path_0.3.Rdata")
ran_net <- "/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/随机网络"
net <- list.files(ran_net, full.names=T, recursive =T)
diff_equl_d_0.3 <- matrix(data = 0,28711,259) 
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/equl_npl")
for(p in 1:500)
{
  print(p);
  random_net <- read.table(net[p], header=T, sep="\t",quote="")
  colnames(random_net) <- gsub("X.","",colnames(random_net))
  r_row <- random_net[order(rownames(random_net)),]
  r_random <- r_row[,order(colnames(r_row))]
  diff_net <- diff_drug_enhancer_path_0.3-r_random 
  for(i in 1:259) {
  diff_equl_d_0.3[which(diff_net[,i]>=0),i] <- diff_equl_d_0.3[which(diff_net[,i]>=0),i]+1 
  }
  filenametxt=paste("diff_equl_r_r",p,".txt",sep="")
  write.table (diff_equl_d_0.3, file =filenametxt,sep="\t")
}
save(diff_equl_d_0.3,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_equl_d_0.3.Rdata")  
diff_equl_NPL_0.3 <- diff_equl_d_0.3/500
p=1
random_net <- read.table(net[p], header=T, sep="\t",quote="")
colnames(random_net) <- gsub("X.","",colnames(random_net))
r_row <- random_net[order(rownames(random_net)),]
r_random <- r_row[,order(colnames(r_row))]
rownames(diff_equl_NPL_0.3) <- rownames(r_random)
colnames(diff_equl_NPL_0.3) <- colnames(r_random)
save(diff_equl_NPL_0.3,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_equl_NPL_0.3.Rdata")
dd <- data.frame(rownames(diff_equl_NPL_0.3))
ee <- data.frame(colnames(diff_equl_NPL_0.3))
diff_equl_final_all_NPL_0.3 <- data.frame(drug_name=NA,enhancer_name=NA,SPL=NA)
for(i in 1:259){
  new_1 <- data.frame(dd[,1],ee[i,1],diff_equl_NPL_0.3[,i])
  colnames(new_1) <- c("drug_name","enhancer_name","SPL")
  diff_equl_final_all_NPL_0.3 <- rbind(diff_equl_final_all_NPL_0.3,new_1)
}
diff_equl_final_all_NPL_0.3 <- diff_equl_final_all_NPL_0.3[-1,]
save(diff_equl_final_all_NPL_0.3,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_equl_final_all_NPL_0.3.Rdata")
load("/pub6/temp/wangli/zcy/network/c_drug.Rdata")
diff_equl_final_all_NPL_0.3[,1] <- gsub("[\"]","",diff_equl_final_all_NPL_0.3[,1])
diff_equl_my_final_NPL_0.3 <- diff_equl_final_all_NPL_0.3[which(diff_equl_final_all_NPL_0.3[,1] %in% c_drug[,2]),]
dim(diff_equl_my_final_NPL_0.3)                      
length(unique(diff_equl_my_final_NPL_0.3[,1]))      
diff_equl_my_final_NPL_0.3_2 <- diff_equl_final_all_NPL_0.3[grep("1,5-isoquinolinediol",diff_equl_final_all_NPL_0.3[,1]),]
diff_equl_my_final_NPL_0.3_3 <- rbind(diff_equl_my_final_NPL_0.3,diff_equl_my_final_NPL_0.3_2)    
save(diff_equl_my_final_NPL_0.3_3,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_equl_my_final_NPL_0.3_3.Rdata")
####################################################
# Merge the result of PL with the result of SPL    #
####################################################
#Read the PL value of the real network
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_PL_0.3.Rdata")
#Read the result of the SPL
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/diff_equl_my_final_NPL_0.3_3.Rdata")
library(tidyverse)
diff_equl_NPL <- separate(data=diff_equl_my_final_NPL_0.3_3,col=enhancer_name,into=c("chr","start","end",sep="[.]"))
diff_equl_NPL <- diff_equl_NPL[,-5]
library(tidyr)
diff_equl_NPL <- unite(diff_equl_NPL, "enhancer_name", start, end, sep = "-")
diff_equl_NPL <- unite(diff_equl_NPL, "enhancer_name", chr, enhancer_name, sep = ":")
diff_equl_NPL$drug_name <- gsub("[\\]","",diff_equl_NPL$drug_name)
diff_PL_0.3$drug_name <- gsub("[\"]","",diff_PL_0.3$drug_name)
equl_NPL_PL <- merge(diff_PL_0.3,diff_equl_NPL,by=c("drug_name","enhancer_name"))
save(equl_NPL_PL,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/equl_NPL_PL.RData")
write.table (equl_NPL_PL,"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/equl_NPL_PL.txt",sep="\t",row.names=FALSE)
#######################################################
#    Merge the result of NES with the result of SPL   #
#######################################################
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/GSEA_enh_drug_final.Rdata")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/equl_NPL_PL.RData")
colnames(GSEA_enh_drug_final)[1:2] <- c("drug_name","enhancer_name")
GSEA_enh_drug_final[,1] <- tolower(GSEA_enh_drug_final[,1])
equl_NES_NPL_PL <- merge(GSEA_enh_drug_final,equl_NPL_PL,by=c("drug_name","enhancer_name")) 
save(equl_NES_NPL_PL,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/equl_NES_NPL_PL.RData") 
equl_NES_NPL_PL_mi <- equl_NES_NPL_PL
equl_NES_NPL_PL_mi$MI <- equl_NES_NPL_PL_mi$NES/equl_NES_NPL_PL_mi$SPL
equl_NES_NPL_PL_mi_inf <- equl_NES_NPL_PL_mi[-which(equl_NES_NPL_PL_mi$PL=="Inf"),]
save(equl_NES_NPL_PL_mi,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/equl_NES_NPL_PL_mi.RData") 
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/NPL_PL.RData")
NES_NPL_PL <- merge(GSEA_enh_drug_final,equl_NPL_PL,by=c("drug_name","enhancer_name")) 
save(NES_NPL_PL,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/NES_NPL_PL.RData") 
NES_NPL_PL_mi <- NES_NPL_PL
NES_NPL_PL_mi$MI <- NES_NPL_PL_mi$NES/NES_NPL_PL_mi$SPL
NES_NPL_PL_mi_inf <- NES_NPL_PL_mi[-which(NES_NPL_PL_mi$PL=="Inf"),]
save(NES_NPL_PL_mi,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/NES_NPL_PL_mi.RData") 
##############################################################################################
#  Construct the properties of drugs, enhancers, and genes in the network--Cytoscape         #
##############################################################################################
#Drug attribute
load("/pub6/temp/wangli/zcy/network/drug_target/two_dd.Rdata")
drug_attribute <- data.frame(node=unique(two_dd_1[,2]),attribute="drug")
#eRNA attribute
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
diff_r_p_0.3 <- unique(diffenhancer_diffgene_cor[,c(4,6)]) 
colnames(diff_r_p_0.3) <- c("enhancer","gene")
enhancer_attribute <- data.frame(node=unique(diff_r_p_0.3[,1]),attribute="enhancer")
#Gene attribute
load("/pub6/temp/wangli/zcy/network/result/dgg_network.Rdata")
dgg_network <- dgg_network
node <- data.frame(from=dgg_network[,1],to=dgg_network[,3])   
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
diff_r_p_0.3 <- unique(diffenhancer_diffgene_cor[,c(4,6)]) 
colnames(diff_r_p_0.3) <- c("enhancer","gene")
node1 <- data.frame(from=diff_r_p_0.3[,1],to=diff_r_p_0.3[,2])
diff_n <- rbind(node,node1)        ######362272
diff_nn <- diff_n[which(diff_n[,1]!="<NA>"),] 
save(diff_nn,file="/pub6/temp/wangli/zcy/network/result/diff_nn.RData")
library(igraph)
diff_gg <- graph_from_data_frame(diff_nn,directed = FALSE)
n <- data.frame(node=get.vertex.attribute(diff_gg)[[1]])
non_gene <- c(drug_attribute[,1],enhancer_attribute[,1])
gene <- setdiff(n[,1],non_gene)
gene_attribute <- data.frame(node=unique(gene),attribute="gene")  #17553
#all nodes attribute
node_attribute <- rbind(drug_attribute,enhancer_attribute,gene_attribute)
a<-data.frame(node="1,5-isoquinolinediol",attribute="drug")
node_attribute <- rbind(node_attribute,a)
write.table(node_attribute, file ="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/node_attribute.txt",sep="\t",row.names=F,col.names=T,quote=F)
node_attribute <- read.table("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/node_attribute.txt",sep="\t",header=T)
#################################
#          eRNA attribute       #
#################################
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/analysis_enhancer.RData")
analysis_enhancer$Mutation_1[analysis_enhancer$Mutation.Status=="Mutation"] <- 1
analysis_enhancer$Mutation_1[analysis_enhancer$Mutation.Status=="Non_mutation"] <- 0
analysis_enhancer$Methylation_1[analysis_enhancer$Methylation.Status=="Hypomethylation"] <- 1
analysis_enhancer$Methylation_2[analysis_enhancer$Methylation.Status=="Hypermethylation"] <- 1
analysis_enhancer$Methylation_1[analysis_enhancer$Methylation.Status!="Hypomethylation"] <- 0
analysis_enhancer$Methylation_2[analysis_enhancer$Methylation.Status!="Hypermethylation"] <- 0
analysis_enhancer$CNV_1[analysis_enhancer$CNV.Status=="Del"] <- 1
analysis_enhancer$CNV_2[analysis_enhancer$CNV.Status=="Amp"] <- 1
analysis_enhancer$CNV_1[analysis_enhancer$CNV.Status!="Del"] <- 0
analysis_enhancer$CNV_2[analysis_enhancer$CNV.Status!="Amp"] <- 0
analysis_enhancer$regualation_1[analysis_enhancer$regualation_class=="Downregulated"] <- 1
analysis_enhancer$regualation_2[analysis_enhancer$regualation_class=="Upregulated"] <- 1
analysis_enhancer$regualation_1[analysis_enhancer$regualation_class!="Downregulated"] <- 0
analysis_enhancer$regualation_2[analysis_enhancer$regualation_class!="Upregulated"] <- 0
analysis_enhancer_network <- analysis_enhancer
write.table(analysis_enhancer_network, file ="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/analysis_enhancer_network.txt",sep="\t",row.names=F,col.names=T,quote=F)
analysis_enhancer_network <- read.table("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/analysis_enhancer_network.txt",sep="\t",header=T)
###################################################
#  Properties of differentially expressed genes   #
###################################################
diff_gene_deseq2 <- read.table("/pub6/temp/wangli/zcy/diff_result/diff_gene_deseq2.txt",sep="\t") 
diff_gene_deseq2 <- diff_gene_deseq2[,c(2,6)]
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
gene_symple2 <- gene_symple[!duplicated(gene_symple[,1]),]
colnames(gene_symple2) <- c("gene","symple")
rownames(gene_symple2) <- gene_symple2[,1]
gene_FC <- merge(gene_symple2,diff_gene_deseq2,by="row.names",sort=FALSE)
gene_FC <- gene_FC[,3:4]
colnames(gene_FC) <- c("gene_name","log2FC")
write.table(gene_FC, file ="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/gene_FC.txt",sep="\t",row.names=F,col.names=T,quote=F)
gene_FC <- read.table("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/gene_FC.txt",sep="\t",header=T)
##########################################################################
#        Attributes of all edges in the network                          #
##########################################################################
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diffenhancer_diffgene_cor_0.3.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/diff_result/gene_symple.txt",sep="\t")
colnames(gene_symple) <- c("ENSG","symple")
diffenhancer_diffgene_cor <- merge(diffenhancer_diffgene_cor_0.3,gene_symple,by="ENSG")
enhancer_gene_attribute_1 <- data.frame(node1=diffenhancer_diffgene_cor[,4],node2=diffenhancer_diffgene_cor[,6],cor=round(as.numeric(diffenhancer_diffgene_cor[,2]),digits = 4))
enhancer_gene_attribute_2 <- data.frame(node1=diffenhancer_diffgene_cor[,6],node2=diffenhancer_diffgene_cor[,4],cor=round(as.numeric(diffenhancer_diffgene_cor[,2]),digits = 4))
enhancer_gene_label <- rbind(enhancer_gene_attribute_1,enhancer_gene_attribute_2)
enhancer_gene_label$cor <- as.numeric(enhancer_gene_label$cor)
write.table(enhancer_gene_label, file ="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_gene_label.txt",sep="\t",row.names=F,col.names=T,quote=F)     
enhancer_gene_label <- read.table("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_gene_label.txt",sep="\t",header=T)
#####################################################
#       Analysis--Final eRNA-drug pairs             #
#####################################################
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/NES_NPL_PL_mi.RData") 
enhancer_drug_network <- NES_NPL_PL_mi[,c("drug_name","enhancer_name","NES","SPL","PL","MI")]
enhancer_drug_network$regulation[which(enhancer_drug_network$NES<0)] <- "Inhibition"
enhancer_drug_network$regulation[which(enhancer_drug_network$NES>0)] <- "Promotion"
write.table(enhancer_drug_network,"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_drug_network.txt",sep="\t",row.names=FALSE)
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/analysis_enhancer.RData")
enhancer_drug_network_analysis <- merge(enhancer_drug_network,analysis_enhancer,by.x="enhancer_name",by.y="enhancer")
enhancer_drug_network_promoter <- enhancer_drug_network_analysis[which(enhancer_drug_network_analysis$MI>0 & enhancer_drug_network_analysis$regualation_class=="Downregulated"),]    #26011
enhancer_drug_network_inhibite <- enhancer_drug_network_analysis[which(enhancer_drug_network_analysis$MI<0 & enhancer_drug_network_analysis$regualation_class=="Upregulated"),]       #11776
enhancer_drug_network_promoter_inhibite <- rbind(enhancer_drug_network_promoter,enhancer_drug_network_inhibite)
save(enhancer_drug_network_promoter,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_drug_network_promoter.RData")
save(enhancer_drug_network_inhibite,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_drug_network_inhibite.RData")
save(enhancer_drug_network_promoter_inhibite,file="/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_drug_network_promoter_inhibite.RData")
#top5 SI-value of eRNA-drug
enhancer <- data.frame(enhancer=unique(enhancer_drug_network_promoter_inhibite$enhancer_name))
predict_enh_top5 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
colnames(predict_enh_top5) <- colnames(enhancer_drug_network_promoter_inhibite)
for(i in 1:dim(enhancer)[1]){
    a <- enhancer_drug_network_promoter_inhibite[which(enhancer_drug_network_promoter_inhibite$enhancer_name==enhancer[i,1]),]
    a_1 <- a[order(abs(a$MI),decreasing=TRUE)[1:5],]
	predict_enh_top5 <- rbind(predict_enh_top5,a_1)
}
predict_enh_top5 <- predict_enh_top5[-1,]    
predict_enh_top5 <- predict_enh_top5[which(predict_enh_top5$drug_name!="NA"),]
write.table(predict_enh_top5,"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/predict_enh_top5.txt",sep="\t",row.names=FALSE,,quote=F)
#########################
#         Venn diagram  #
#########################
q=6    
m=343     
n=1309-m    
k=8       
phyper(q-1,m,n,k,lower.tail=F) #0.005322767

 install.packages("VennDiagram")
 library(grid)
 library(futile.logger)
 library(VennDiagram)
venn.plot <- draw.pairwise.venn(343,8,6,c("Specific_Predicted_drug", "FDA_approved_brca_drug"),  
  fill = c("blue", "yellow"),  
  scaled = FALSE,
  main = "Venn Diagram of Drug",
  lty = "blank",  
  cex = 2, 
  cat.cex = 1.5, 
  cat.pos = c(0, 0),
  cat.dist = c(0.03,0.03), 
  cat.just = list(c(0, 0), c(0, 0)), 
  ext.pos = 0,  
  ext.dist = -0.05, ext.length = 0.85, 
  ext.line.lwd = 2, ext.line.lty = "dashed", 
  alpha=0.3, euler.d=T )
grid.draw(venn.plot)
######################################################
#          Density distribution map                  #
######################################################
setwd("C:/Users/Dell/Desktop/result3_method/数据")
load("NPL_PL.RData")
predict_enh_top5 <- read.table("predict_enh_top5.txt",sep="\t",header=T)
load("Dset_log2_per_veh.Rdata")
Dset_log2_per_veh[,1] <- tolower(Dset_log2_per_veh[,1])
pre_all_gsea <- merge(Dset_log2_per_veh[,c(1,2,8)],predict_enh_top5[,1:2],by.x=c("drug","ID"),by.y=c("drug_name","enhancer_name"))
non_pre_all_gsea <- Dset_log2_per_veh[which(Dset_log2_per_veh$p.adjust==1),c(1,2,8)]
pre_all_gsea$class <- "Specific Predictions"
non_pre_all_gsea$class <- "Non-Predictions"
non_pre_all_gsea_2 <- non_pre_all_gsea[sample(1:6709,1139,replace=FALSE),]
data_gsea <- rbind(pre_all_gsea,non_pre_all_gsea)
data_gsea$drug <- tolower(data_gsea$drug)
colnames(data_gsea)[1:2] <- c("drug_name","enhancer_name")
unequl_NPL <- merge(NPL_PL,data_gsea,by=c("drug_name","enhancer_name"))
ggplot(unequl_NPL,aes(x = SPL))+geom_density(aes(fill=class),alpha=0.4)+labs(x="Standardized path length",title="Predicted vs. Non-Predicted Drug-Enhancer Path lengths")+xlim(0,0.25)+ylim(0,75) 
wilcox.test(unequl_NPL[which(unequl_NPL$class=="Specific Predictions"),]$SPL,unequl_NPL[which(unequl_NPL$class=="Non-Predictions"),]$SPL,paired = FALSE,alternative = "two.sided")
p-value < 2.2e-16
ggviolin(unequl_NPL,x ="class",y="PL",fill="class",color="class")
#######################################
#       PL distribution map           #
#######################################
setwd("C:/Users/Dell/Desktop/result3_method/数据")
load("NPL_PL.RData")
PL <- merge(NPL_PL,data_gsea,by=c("drug_name","enhancer_name"))
PL$PL <- factor(PL$PL)
ggdensity(unequl_NPL, x="PL",fill="class",color="class"),
         palette = "npg", #杂志nature的配色
         sort.val = "desc", #下降排序
         sort.by.groups=FALSE, #不按组排序
         x.text.angle=60)
#####################################################################################
#                    enhancer_subnetwork---function                                 #
#    Extraction of enhancers by the shortest path subnetwork of drug candidates     #
#####################################################################################
load("/pub6/temp/wangli/zcy/network/result/diff_nn.RData")
library(igraph)
diff_gg <- graph_from_data_frame(diff_nn,directed = FALSE)
n <- data.frame(node=get.vertex.attribute(diff_gg)[[1]])
node_1 <- data.frame(node=n,vids=seq(1,length(n[,1]),by=1))
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_drug_network_promoter_inhibite.RData")
enhancer_gene_label <- read.table("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-药物网络/enhancer_gene_label.txt",sep="\t",header = T)
colnames(enhancer_gene_label)[1:2] <- c("node1","node2")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/hi-c/all_enhancer_symple_hic_cor.RData")
enhancer_symple_hic_1 <- all_enhancer_symple_hic_cor[,c(2,8,4,5,3)]
colnames(enhancer_symple_hic_1)[1:4] <- c("node1","node2","target1","target2")
enhancer_symple_hic_2 <- all_enhancer_symple_hic_cor[,c(8,2,5,4,3)]
colnames(enhancer_symple_hic_2)[1:4] <- c("node1","node2","target1","target2")
enhancer_symple_hic <- rbind(enhancer_symple_hic_1,enhancer_symple_hic_2)
load("/pub6/temp/wangli/zcy/network/drug_target/drug_gene_annotation.RData")
drug_gene_annotation_1 <- drug_gene_annotation
colnames(drug_gene_annotation_1)[1:2] <- c("node1","node2") 
drug_gene_annotation_2 <- drug_gene_annotation[,c(2,1,3:10)]
colnames(drug_gene_annotation_2)[1:2] <- c("node1","node2") 
drug_gene_annotation <- rbind(drug_gene_annotation_1,drug_gene_annotation_2)
drug_gene_annotation <- unique(drug_gene_annotation[,c(1,2,5,6,10)])
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/foldchange_drug_gene_1.RData") 
foldchange_drug_gene$type[which(foldchange_drug_gene$foldchange<1)] <- "Inhibition"
foldchange_drug_gene$type[which(foldchange_drug_gene$foldchange>1)] <- "Promoter"
foldchange_drug_gene_1 <- foldchange_drug_gene
colnames(foldchange_drug_gene_1)[1:2] <- c("node1","node2") 
foldchange_drug_gene_2 <- foldchange_drug_gene[,c(2,1,3,4)]
colnames(foldchange_drug_gene_2)[1:2] <- c("node1","node2") 
drug_gene_FC_annotation <- rbind(foldchange_drug_gene_1,foldchange_drug_gene_2)
library(dplyr)
enhancer_subnetwork <- function(x,y,z,m){
      input_enhancer <- x
	  target_gene <- z
      a <- enhancer_drug_network_promoter_inhibite[which(enhancer_drug_network_promoter_inhibite$enhancer_name==input_enhancer),]
	  if(m=="high"){
      a_1 <- a[order(a$MI,decreasing = F),] 
	  }else
	  {
	  a_1 <- a[order(a$MI,decreasing = T),]  
	  }
	  a_1 <- a_1[which(a_1$PL!="NA"),]
      enh_enh <- data.frame(NA)
	  if(length(unique(a_1$drug_name))>=y){
	  num_drug=y
	  }else 
	  {
	  num_drug=length(unique(a_1$drug_name))
	  }
      for(i in 1:num_drug){
          enh <- all_shortest_paths(diff_gg,from=input_enhancer,to=a_1[i,2])
          enh_1 <- as.data.frame(t(sapply(enh$res, as_ids)))
          enh_enh<- plyr::rbind.fill(enh_enh,enh_1)   
      }
      enh_enh_1 <- enh_enh[-1,]
      enh_enh_1 <- enh_enh_1[,-1]
      v <- data.frame(node=c(unique(unlist(as.list(enh_enh_1))),target_gene))
      vids <- merge(v,node_1,by="node")
      vids_enh <- as.vector(vids[,2])
      subgraph_enh <- induced.subgraph(diff_gg,vids_enh)
      edge_enh <- get.edgelist(subgraph_enh)
      colnames(edge_enh) <- c("node1","node2")
      edge_enh_1 <- merge(enhancer_gene_label,edge_enh,by=c("node1","node2"),all.y=T)
	  edge_enh_2 <- merge(edge_enh_1,drug_gene_FC_annotation,by=c("node1","node2"),all.x=T)
	  print(edge_enh_2)
      filename <- paste("edge_",x,".txt",sep="")
      write.table(edge_enh_2, file =filename,sep="\t",row.names=F,col.names=T,quote=F)
	  if(m=="high"){
      a_2 <- a[order(a$MI,decreasing = F),] 
	  }else
	  {
	  a_2 <- a[order(a$MI,decreasing = T),]  
	  }
	  enhancer_drug_table <- a_2
	  enhancer_drug_table$NL <- enhancer_drug_table$NES/enhancer_drug_table$PL
	  filename <- paste("enhancer_drug_table_",x,".txt",sep="")
	  write.table(enhancer_drug_table, file =filename,sep="\t",row.names=F,col.names=T,quote=F)
}
##############################################################
#                   GSEA_function---function                 #
#                     GSEA plot---eRNA-drug                  #
##############################################################
library(enrichplot)
library(tidyverse)
library(clusterProfiler)
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/log2_per_veh.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diff_r_p_0.3.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/GSEA_enh_drug_final.Rdata")
GSEA_function <- function(m,n){
     filenamepdf <- paste("GSEA.plot_",n,".pdf",sep="")
     pdf(filenamepdf,width=30,height=20)
     drug_name <- m
	 enhancer_name <- n
	 splots <- list()
	 GSEA_enh_drug_final$drug <- tolower(GSEA_enh_drug_final$drug)
	 colnames(log2_per_veh) <- tolower(colnames(log2_per_veh))
  for(i in 1:length(m))
  {
     a_1 <-GSEA_enh_drug_final[which(GSEA_enh_drug_final$ID==n & GSEA_enh_drug_final$drug==m[i]),]
     one <- log2_per_veh[,which(colnames(log2_per_veh)==m[i])]
     new <- sort(one,decreasing=T)
     gene <- diff_r_p_0.3[which(diff_r_p_0.3$enhancer==n),]
     gsea <- GSEA(new,TERM2GENE=gene,verbose=F,minGSSize=1,maxGSSize=1000,pvalueCutoff=0.05,pAdjustMethod="fdr",nPermSimple = 10000) 
	 title_name <- paste(n,"Targets after",m[i],"Perturbation",sep=" ")
     splots[[i]] <- gseaplot2(gsea,1,color="red",pvalue_table = F,ES_geom = "line",subplots = 1:3,title=title_name)
   }
	p <- plot_grid(plotlist = splots,labels = c('A','B','C','D','E'),align="hv",ncol=3,nrow =2)
	print(p)
    dev.off()
}
########################################################################
#           enhancer_diff_analysis_function------Function              #
#                   plot---Differential expression eRNA                #
########################################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA")
eRNA_profile <- read.table("TCGA_BRCA_enhancer_profile.txt",sep="\t",header=T)
normal_exp <- eRNA_profile[,grep("normal",colnames(eRNA_profile))]  
tumor_exp <- eRNA_profile[,grep("tumor",colnames(eRNA_profile))]    
normal_sample <- gsub("_normal","",colnames(normal_exp))
tumor_sample <- gsub("_tumor","",colnames(tumor_exp))
match_sample <- intersect(tumor_sample,normal_sample)
match_tumor_exp <- tumor_exp[,pmatch(match_sample,colnames(tumor_exp))]
match_normal_exp <- normal_exp[,pmatch(match_sample,colnames(normal_exp))] 
paired_expr_data <- cbind(match_normal_exp,match_tumor_exp)
save(paired_expr_data,file="paired_expr_data.RData")
library(ggpubr)
library(ggplot2)
enhancer_diff_analysis_function <- function(x){
     filenamepdf <- paste("enhdiff.plot_",x,".pdf",sep="")
     pdf(filenamepdf,width=10,height=20) 
     enhancer_name <- x
     expr_data <- paired_expr_data[which(rownames(paired_expr_data)==enhancer_name),]
     expr_data <- data.frame(t(expr_data))
     expr_data$class <- "normal_sample"
     expr_data$class[grepl("tumor",rownames(expr_data))] <- "tumor_sample"
     colnames(expr_data)[1] <- "exp"
     expr_data$exp <- as.numeric(expr_data$exp) 
	 #FC value
     a <- length(paired_expr_data)/2
     log2_FC <- log2((mean(as.numeric(expr_data[(a+1):(a*2),1])))/(mean(as.numeric(expr_data[1:a,1]))))
	 log2_FC <- signif(log2_FC,3)
	 log2_FC <- paste("Log2Foldchange = ",log2_FC,sep="")
     title_name <- x		 
     plot1 <- ggplot(expr_data, aes(y=exp,x=class,fill=class)) + geom_boxplot(outlier.shape = NA) + labs(y="FPKM value",x="",title=title_name) + stat_compare_means(method="t.test",paired=T)+annotate("text",y = max(expr_data$exp)*0.8,x=1, label = log2_FC)
	 plot2 <- ggpaired(expr_data, x = "class", y = "exp",
         color = "class", line.color = "gray", line.size = 0.4, palette = "jco")+ labs(y="FPKM value",x="",title=title_name) + stat_compare_means(method="t.test",paired=T)+annotate("text",y = max(expr_data$exp)*0.8,x=1,label = log2_FC)
     plot3 <- ggarrange(plot1,plot2,ncol=1,nrow=2)
	 print(plot3)
     dev.off()
}
#######################################################################
#               cmap_drug_gene_function------Function                 #
#     plot---Differential expression of genes after drug action       #
#######################################################################
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/result_vehicle_pro.RData")
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/result_perturbation_pro.RData")
result_perturbation_pro <- result_perturbation_pro[order(rownames(result_perturbation_pro)),]
result_perturbation_pro <- result_perturbation_pro[,order(colnames(result_perturbation_pro))]
result_vehicle_pro <- result_vehicle_pro[order(rownames(result_vehicle_pro)),]
result_vehicle_pro <- result_vehicle_pro[,order(colnames(result_vehicle_pro))]
cmap_drug_gene_function <- function(x,y){
        filenamepdf <- paste("cmap_drug_gene_",x,"_",y,".pdf",sep="")
        pdf(filenamepdf,width=10,height=10) #保存图为pdf
        drug <- x
		gene <- y
		expr_data <- data.frame(drug=rbind(result_perturbation_pro[y,x],result_vehicle_pro[y,x]))
        expr_data$class <- c("perturbation","vehicle")
		FC <- expr_data[1,1]/expr_data[2,1]
		FC <- signif(FC,3)
	    FC <- paste("Foldchange = ",FC,sep="")
		title_name <- paste(x,"-",y)
		y_name <- paste(y,"value",sep=" ")
		plot1 <- ggplot(expr_data) + geom_col(aes(y=drug,x=class,fill=class))+labs(y=y_name,x="",title=title_name) +annotate("text",y = max(expr_data$drug)*0.8, x=1, label = FC)
		print(plot1)
     dev.off()
}

##########################################################
#       enhancer_gene_cor_function------Function         #
#    plot---Correlation analysis of eRNA and target gene #
##########################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA")
load("e_p.RData")
load("t_p.RData")
gene_symple <- read.table("/pub6/temp/wangli/zcy/result/gene_symple.txt",sep="\t",header=T)
t_p_symple <- merge(gene_symple,t_p,by.y="row.names",by.x="ENSEMBL")  #19293个
t_p_symple <- t_p_symple[!duplicated(t_p_symple$SYMBOL),] 
rownames(t_p_symple) <- t_p_symple$SYMBOL
t_p_symple <- t_p_symple[,-2]
save(t_p_symple,file="t_p_symple.RData")
library(psych)
library(corrplot)
enhancer_gene_cor_function <- function(x,y){
    filenamepdf <- paste("allcor.plot_",x,".pdf",sep="")
    pdf(filenamepdf,width=10,height=10) 
    enhancer <- x
    target <- y
    enhancer_exp <- data.frame(t(e_p[which(rownames(e_p)==enhancer),]))
    enhancer_exp2 <- data.frame(enhancer_exp[,order(colnames(enhancer_exp))])
    gene_exp <- data.frame(t(t_p_symple[which(rownames(t_p_symple) %in% target),-1]))
    gene_exp2 <- data.frame(gene_exp[,order(colnames(gene_exp))])
    enhancer_gene <- merge(enhancer_exp,gene_exp,by="row.names",all=TRUE)
    cor_result <- corr.test(enhancer_gene[,-1],method="spearman")
    sig.A <- cor.mtest(enhancer_gene[,-1], conf.level = .95)
    corrplot(corr =cor_result$r,type="upper",tl.pos="tp",tl.col="black")
    corrplot(corr = cor_result$r,add=TRUE, type="lower",method="color", col="white",tl.cex=1.2,diag=FALSE,tl.pos="n", cl.pos="n",p.mat = sig.A$p, insig = "label_sig", sig.level = c(.00001, .0001,.001,.05),pch.cex=1, pch.col = "black",tl.col="black",addgrid.col = "black")
    dev.off()
}
antagonist_enhancer_gene_cor_function <- function(x,y){
    filenamepdf <- paste("antagonist_allcor.plot_",x,".pdf",sep="")
    pdf(filenamepdf,width=10,height=10) 
    enhancer <- x
    target <- y
    enhancer_exp <- data.frame(t(e_p[which(rownames(e_p)==enhancer),]))
    enhancer_exp2 <- data.frame(enhancer_exp[,order(colnames(enhancer_exp))])
    gene_exp <- data.frame(t(t_p_symple[which(rownames(t_p_symple) %in% target),-1]))
    gene_exp2 <- data.frame(gene_exp[,order(colnames(gene_exp))])
    enhancer_gene <- merge(enhancer_exp,gene_exp,by="row.names",all=TRUE)
    cor_result <- corr.test(enhancer_gene[,-1],method="spearman")
    sig.A <- cor.mtest(enhancer_gene[,-1], conf.level = .95)
    corrplot(corr =cor_result$r,type="upper",tl.pos="tp",tl.col="black") 
    corrplot(corr = cor_result$r,add=TRUE, type="lower",method="color", col="white",tl.cex=1.2,diag=FALSE,tl.pos="n", cl.pos="n",p.mat = sig.A$p, insig = "label_sig", sig.level = c(.00001, .0001,.001,.05),pch.cex=1, pch.col = "black",tl.col="black",addgrid.col = "black")
    dev.off()
}
###############################################################
#       enhancer_gene_cor_function2------Function             #
#  one plot---Correlation analysis of eRNA and target gene    #
###############################################################
library(RColorBrewer)
enhancer_gene_cor_function2 <- function(x,y){
    filenamepdf <- paste("singlecor.plot_",x,".pdf",sep="")
    pdf(filenamepdf,width=30,height=10*ceiling(length(y)/3)) 
    enhancer <- x
    target <- y
    enhancer_exp <- data.frame(t(e_p[which(rownames(e_p)==enhancer),]))
    enhancer_exp2 <- data.frame(enhancer_exp[,order(colnames(enhancer_exp))])
    gene_exp <- data.frame(t(t_p_symple[which(rownames(t_p_symple) %in% target),-1]))
    gene_exp2 <- data.frame(gene_exp[,order(colnames(gene_exp))])
    enhancer_gene <- merge(enhancer_exp,gene_exp,by="row.names",all=TRUE)
   splots <- list()
   for(i in 1:(dim(enhancer_gene)[2]-2))
   {
     cor_value <- cor(enhancer_gene[,2],enhancer_gene[,i+2],method = c("spearman"))
     eq <- substitute(""~cor~"="~a~"", list(a = signif(cor_value,3))) 
     exp <- as.expression(eq)
	 title_name <- paste("The scatter chart of",enhancer,"and",colnames(enhancer_gene)[i+2],sep=" ")
	 y_name <- paste(colnames(enhancer_gene)[i+2],"expression",sep=" ")
     print(exp)
     palette <- brewer.pal(11,"Spectral")
     enhancer_gene_exp <- data.frame(enhancer=enhancer_gene[,2],gene=enhancer_gene[,i+2])
     splots[[i]] <- ggplot(data = enhancer_gene_exp,aes(x = enhancer,y = gene)) + geom_point()+geom_smooth(method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+annotate("text",x=max(enhancer_gene_exp$enhancer)*0.5,y=max(enhancer_gene_exp$gene)*0.8, label = exp,vjust = 0,nudge_y = 0.1,size = 6,parse = T,color = c("blue4"))+
	   labs(x ='enhancer FPKM value',y= y_name,
        title = title_name)+
		theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey",linetype = 2),
        axis.line = element_line(colour = "black",size = rel(2),arrow = arrow(angle = 30,length = unit(0.1,"inches"))),
        axis.title.y = element_text(size = rel(2),hjust = 0.5),
        axis.title.x = element_text(size = rel(2),hjust = 0.5),
        axis.text.x = element_text(size = rel(2),hjust = 1),
        axis.text.y = element_text(hjust = 1,size = rel(2)),
        axis.ticks = element_line(size = rel(1.3)),
        plot.title = element_text(size = rel(1.8)),
        plot.margin = margin(15,9,9,30))
	}
	library(cowplot)
	p <- plot_grid(plotlist = splots,align="hv",ncol=3,nrow =ceiling(length(y)/3))
	print(p)
	dev.off()
}
###############################################################
#       enhancer_gene_cor_function3------Function             #
# one plot---Correlation analysis of single eRNA and one gene #
###############################################################
enhancer_gene_cor_function3 <- function(x,y){
    filenamepdf <- paste("single.plot_",x,".pdf",sep="")
    pdf(filenamepdf,width=10,height=10) 
    enhancer <- x
    target <- y
    enhancer_exp <- data.frame(t(e_p[which(rownames(e_p)==enhancer),]))
    enhancer_exp2 <- data.frame(enhancer_exp[,order(colnames(enhancer_exp))])
    gene_exp <- data.frame(t(t_p_symple[which(rownames(t_p_symple) %in% target),-1]))
    gene_exp2 <- data.frame(gene_exp[,order(colnames(gene_exp))])
    enhancer_gene <- merge(enhancer_exp,gene_exp,by="row.names",all=TRUE)
     cor_value <- cor(enhancer_gene[,2],enhancer_gene[,3],method = c("spearman"))
     eq <- substitute(""~cor~"="~a~"", list(a = signif(cor_value,3))) 
     exp <- as.expression(eq)
	 title_name <- paste("The scatter chart of",enhancer,"and",colnames(enhancer_gene)[1+2],sep=" ")
	 y_name <- paste(colnames(enhancer_gene)[1+2],"expression",sep=" ")
     print(exp)
     palette <- brewer.pal(11,"Spectral")
     enhancer_gene_exp <- data.frame(enhancer=enhancer_gene[,2],gene=enhancer_gene[,1+2])
     splots <- ggplot(data = enhancer_gene_exp,aes(x = enhancer,y = gene)) + geom_point()+geom_smooth(method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+annotate("text",x=max(enhancer_gene_exp$enhancer)*0.5,y=max(enhancer_gene_exp$gene)*0.8, label = exp,vjust = 0,size = 6,parse = T,color = c("blue4"))+
	   labs(x ='enhancer FPKM value',y= y_name,
        title = title_name)+
		theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey",linetype = 2),
        axis.line = element_line(colour = "black",size = rel(2),arrow = arrow(angle = 30,length = unit(0.1,"inches"))),
        axis.title.y = element_text(size = rel(2),hjust = 0.5),
        axis.title.x = element_text(size = rel(2),hjust = 0.5),
        axis.text.x = element_text(size = rel(2),hjust = 1),
        axis.text.y = element_text(hjust = 1,size = rel(2)),
        axis.ticks = element_line(size = rel(1.3)),
        plot.title = element_text(size = rel(1.8)),
        plot.margin = margin(15,9,9,30))
	p <- splots
	print(p)
	dev.off()
}
###################################
#          eRNA plot              #
###################################
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/paired_expr_data.RData")
library(ggpubr)
library(ggplot2)
library(psych)
library(corrplot)
library(cowplot)  
library(RColorBrewer) 
library(enrichplot)
library(tidyverse)
library(clusterProfiler)
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/NES_NPL_PL_mi.RData") 
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/e_p.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/t_p_symple.RData")
load("/pub6/temp/wangli/zcy/0000——yuanshishuju/cmap数据库数据/log2_per_veh.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/diff_r_p_0.3.RData")
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/GSEA/GSEA_enh_drug_final.Rdata")
#####################################################
#                chr10:9115735-9116358              #
#####################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/结果图/增强子——chr10:9115735-9116358")
#Differential expression analysis
enhancer_diff_analysis_function("chr10:9115735-9116358")
#eRNA-drug-gene network
enhancer_subnetwork("chr10:9115735-9116358",5,c("GATA3","ESR1","FOXA1","AGR3"),"high")
#The correlation diagram of eRNA-gene
enhancer_gene_cor_function("chr10:9115735-9116358",c("CYP2B6","FOXA1","GATA3","ESR1","NR1I3","HSD17B10","EHMT2"))
enhancer_gene_cor_function2("chr10:9115735-9116358",c("CYP2B6","FOXA1","GATA3","ESR1","NR1I3","HSD17B10","EHMT2"))
antagonist_enhancer_gene_cor_function("chr10:9115735-9116358",c("CYP2B6","FOXA1","GATA3","ESR1","AGR3","APC","TOP2A"))
#GSEA plot of drug-eRNA
GSEA_function(c("isoniazid","chlorogenic acid","hydrastinine","sulmazole","pheniramine"),"chr10:9115735-9116358")
#####################################################
#              chr10:3848037-3849536                #
#####################################################
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/结果图/增强子——chr10:3848037-3849536")
#Differential expression analysis
enhancer_diff_analysis_function("chr10:3848037-3849536")
#eRNA-drug-gene network
enhancer_subnetwork("chr10:3848037-3849536",5,c("COL11A1","PLAU","COL10A1","INHBA","FN1","SULF1","KLF6"),"high")
#The correlation diagram of eRNA-gene
enhancer_gene_cor_function("chr10:3848037-3849536",c("LCN1","GPRC5A","COL1A1","FN1","RGS4","SPP1"))
#GSEA plot of drug-eRNA
GSEA_function(c("tretinoin","valproic acid","rottlerin","risperidone","wortmannin"),"chr10:3848037-3849536")
#cmap histogram of changes in gene expression before and after drug action
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/cmap-药物-基因图")
cmap_drug_gene_function("isoniazid","CYP2B6")
cmap_drug_gene_function("chlorogenic acid","GATA3")
cmap_drug_gene_function("chlorogenic acid","FOXA1")
cmap_drug_gene_function("chlorogenic acid","ESR1")
cmap_drug_gene_function("tretinoin","LCN1")
cmap_drug_gene_function("risperidone","RGS4")
cmap_drug_gene_function("rottlerin","RGS4")
cmap_drug_gene_function("valproic acid","COL1A1")
cmap_drug_gene_function("valproic acid","FN1")

