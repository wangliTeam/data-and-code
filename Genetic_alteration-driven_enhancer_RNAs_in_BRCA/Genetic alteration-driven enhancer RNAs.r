#############################################################
#     Treatment of eRNA expression profile in 32 cancer     #
#      Identify active enhancers for each cancer study      #
#############################################################
cancer_32 <- c("TCGA_ACC","TCGA_BLCA","TCGA_BRCA","TCGA_CESC","TCGA_CHOL","TCGA_COADREAD","TCGA_DLBC","TCGA_ESCA","TCGA_GBM","TCGA_HNSC","TCGA_KICH","TCGA_KIRC","TCGA_KIRP","TCGA_LAML","TCGA_LGG","TCGA_LIHC","TCGA_LUAD","TCGA_LUSC","TCGA_MESO","TCGA_OV","TCGA_PAAD","TCGA_PCPG","TCGA_PRAD","TCGA_SARC","TCGA_SKCM","TCGA_STAD","TCGA_TGCT","TCGA_THCA","TCGA_THYM","TCGA_UCEC","TCGA_UCS","TCGA_UVM")
cance_exp <-"/pub6/temp/wangli/zcy/0000——yuanshishuju/TCGA_ALL_cancer_enhancer_60k"
cance_exp_name <-list.files(cance_exp, full.names=T, recursive =T)
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/brca为例")
enhancer_15808_BRCA_1095 <- read.table("enhancer_15808_BRCA_1095_RPKM.txt",sep="\t")
cancer_name <- as.data.frame(cancer_32)
for(p in 1:length(cance_exp_name))
 {
   all_eRNA <- read.table(cance_exp_name[p], header=T, sep="\t")      
   rownames(all_eRNA) <- gsub("FANTOM5_","",rownames(all_eRNA))      
   eRNA <- intersect(rownames(all_eRNA),rownames(enhancer_15808_BRCA_1095))
   enhancer_15436 <- all_eRNA[eRNA,]
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_ALL_CANCER_enhancer_15436")
   filenamecsv=paste(cancer_name[p,],"_15436_RPKM",".txt",sep="") 
   write.table(enhancer_15436, file =filenamecsv,sep="\t")
   num_tumor <- enhancer_15436[,grep("tumor",colnames(enhancer_15436))]
   enh_proportion <- as.data.frame(cbind(rowSums(num_tumor > 0)/ncol(num_tumor)))
   tumor_eRNA <- num_tumor[which(enh_proportion > 0.1),]   
   tumor_eRNA_name <- intersect(rownames(enhancer_15436),rownames(tumor_eRNA))
   enhancer_0.1 <- enhancer_15436[tumor_eRNA_name,]
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_ALL_CANCER_enhancer_15436_10%")
   filenamecsv=paste(cancer_name[p,],"_enhancer_profile",".txt",sep="") 
   write.table (enhancer_0.1, file =filenamecsv,sep="\t")
 }

############################################################################
#    Identify differences, mutations, methylation, CNVS in 12 cancers      #
############################################################################
############################################################
#    Function--Identify differential expression enhancer   #
############################################################
diff_eRNA_function <- function(x)
{
   eRNA_profile <- x
   normal_exp <- eRNA_profile[,grep("normal",colnames(eRNA_profile))]  
   tumor_exp <- eRNA_profile[,grep("tumor",colnames(eRNA_profile))]    
   normal_sample <- gsub("_normal","",colnames(normal_exp))
   tumor_sample <- gsub("_tumor","",colnames(tumor_exp))
   match_sample <- intersect(tumor_sample,normal_sample)
   match_tumor_exp <- tumor_exp[,pmatch(match_sample,colnames(tumor_exp))]
   match_normal_exp <- normal_exp[,pmatch(match_sample,colnames(normal_exp))] 
   paired_expr_data <- cbind(match_normal_exp,match_tumor_exp)
   Pvalue <- rep(0,nrow(paired_expr_data))
   FC <- rep(0,nrow(paired_expr_data))
   a <- length(match_normal_exp)
   for(i in 1:nrow(paired_expr_data))
   {
      if(sd(paired_expr_data[i,1:a])==0 && sd(paired_expr_data[i,(a+1):(a*2)])==0 )
      {
         Pvalue[i] <- "NA"
         log2_FC[i] <- "NA"
      }
      else
	  {
         y=t.test(as.numeric(paired_expr_data[i,(a+1):(a*2)]),as.numeric(paired_expr_data[i,1:a]),paired=T)
         Pvalue[i]<-y$p.value
         FC[i]<-(mean(as.numeric(paired_expr_data[i,(a+1):(a*2)])))/(mean(as.numeric(paired_expr_data[i,1:a])))
       }
    }
   res_paired_enhancer <- data.frame(enhancer=rownames(paired_expr_data),FC=FC,p=Pvalue)
   return(res_paired_enhancer)
}     
#####################################################################
#   Identify differential expression enhancers for each cancer      #
#####################################################################
cance_eRNA_profile <-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_ALL_CANCER_enhancer_15436_10%"
file <- list.files(cance_eRNA_profile, full.names=T, recursive =T)
cancer_name <- as.data.frame(cancer_12)
for(p in 1:length(file))
{
   all_eRNA_profile <- read.table(file[p], header=T, sep="\t") 
   normal_exp <- all_eRNA_profile[,grep("normal",colnames(all_eRNA_profile))]  
   tumor_exp <- all_eRNA_profile[,grep("tumor",colnames(all_eRNA_profile))]   
   normal_sample <- gsub("_normal","",colnames(normal_exp))
   tumor_sample <- gsub("_tumor","",colnames(tumor_exp))
   match_sample <- intersect(tumor_sample,normal_sample)
   if(length(match_sample) > 20){
   res_paired_enhancer <- diff_eRNA_function(all_eRNA_profile)
   #fdr<0.05 and |log2fc|>1
   res_paired_enhancer <- res_paired_enhancer[order(res_paired_enhancer$p),]
   res_paired_enhancer$diff_fdr <- p.adjust(res_paired_enhancer$p,method = "fdr")
   res_paired_enhancer_diff <- subset(res_paired_enhancer,diff_fdr<0.05 & (log2_FC > 1 | log2_FC < -1)) 
   res_paired_enhancer_diff <- res_paired_enhancer_diff[-which(res_paired_enhancer_diff$log2_FC=="Inf"),]   
   print(dim(res_paired_enhancer_diff))
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/TCGA_12_CANCER_diff_enhancer")
   filenamecsv=paste(cancer_name[p,],"_diff_enhancer",".txt",sep="") 
   write.table(res_paired_enhancer_diff, file =filenamecsv,sep="\t")
   }
}
cancer_12 <-  c("TCGA_BRCA","TCGA_COADREAD","TCGA_HNSC","TCGA_KICH","TCGA_KIRC","TCGA_KIRP","TCGA_LIHC","TCGA_LUAD","TCGA_LUSC","TCGA_PRAD","TCGA_STAD","TCGA_THCA")
cance_name_12 <- as.data.frame(cancer_12)
diff_enhancer <-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/TCGA_12_CANCER_diff_enhancer"
file <- list.files(diff_enhancer, full.names=T, recursive =T)
for(i in 1:12){
   diff_12_cancer <- read.table(file[i], header=T, sep="\t")
   library(tidyverse)
   diff_bedtools <- separate(data=diff_12_cancer,col=enhancer,into=c("chr","start","end",sep="[.]"))
   diff_12_bedtools <- diff_bedtools[,1:3]
   diff_12_bedtools$chr <- gsub("chr","",diff_12_bedtools$chr)
   diff_12_bedtools$enhancer <- diff_12_cancer$enhancer
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_data")
   filenamecsv=paste(cance_name_12[i,],"_diff_bedtools",".txt",sep="") 
   write.table(diff_12_bedtools, file =filenamecsv,sep="\t",row.names=F,col.names=F,quote=F)
}
###################################################################
#            Processing of copy number data----q<0.25             #
###################################################################
cancer_12 <-  c("TCGA_BRCA","TCGA_COADREAD","TCGA_HNSC","TCGA_KICH","TCGA_KIRC","TCGA_KIRP","TCGA_LIHC","TCGA_LUAD","TCGA_LUSC","TCGA_PRAD","TCGA_STAD","TCGA_THCA")
cancer_name_12 <- as.data.frame(cancer_12)
GISTIC2 <-"/pub6/temp/wangli/zcy/0000——yuanshishuju/GISTIC2_12_cancer"
file <- list.files(GISTIC2, full.names=T, recursive =T)
for(i in 1:12){
   GISTIC2_12_cancer <- read.table(file[i], header=T, sep="\t")
   GISTIC2_12_cancer_q <- GISTIC2_12_cancer[which(GISTIC2_12_cancer$q> (-log10(0.25))),]       
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/TCGA_12_CANCER_GISTIC2")
   filenamecsv=paste(cancer_name_12[i,],"_GISTIC2_qvalue_0.25",".txt",sep="") 
   write.table(GISTIC2_12_cancer_q, file =filenamecsv,sep="\t",row.names=F,col.names=T,quote=F)
   GISTIC2_12_bedtools <- GISTIC2_12_cancer_q[,c(2:4,1)]
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_data")
   filenamecsv=paste(cancer_name_12[i,],"_GISTIC2_bedtools",".bed",sep="") 
   write.table(GISTIC2_12_bedtools, file =filenamecsv,sep="\t",row.names=F,col.names=F,quote=F)
}
############################################################
#                 Processing of mutation data              #
############################################################
cancer_12 <-  c("TCGA_BRCA","TCGA_COADREAD","TCGA_HNSC","TCGA_KICH","TCGA_KIRC","TCGA_KIRP","TCGA_LIHC","TCGA_LUAD","TCGA_LUSC","TCGA_PRAD","TCGA_STAD","TCGA_THCA")
cancer_name_12 <- as.data.frame(cancer_12)
nutation <-"/pub6/temp/wangli/zcy/0000——yuanshishuju/mutation_12_cancer"
file <- list.files(nutation, full.names=T, recursive =T)
for(i in 1:12){
   mutation_12_cancer <- read.table(file[i], header=T, sep="\t",quote="")
   mutation_12_cancer_data <- mutation_12_cancer[,c(1,2,4:7,9:15,17:18)]
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/TCGA_12_CANCER_mutation")
   filenamecsv=paste(cancer_name_12[i,],"_mutation_data",".txt",sep="") 
   write.table(mutation_12_cancer_data, file =filenamecsv,sep="\t",row.names=F,quote=F)
   mutation_12_bedtools <- mutation_12_cancer_data[,c(4:6,8)]
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_data")
   filenamecsv=paste(cancer_name_12[i,],"_mutation_bedtools",".bed",sep="") 
   write.table(mutation_12_bedtools, file =filenamecsv,sep="\t",row.names=F,col.names=F,quote=F)
}
####################################################
#              Processing of methylation data      #
#            hg38 is converted to hg19             #
####################################################
methyl <-"/pub6/temp/wangli/zcy/0000——yuanshishuju/methyl_12_cancer"
file <- list.files(methyl, full.names=T, recursive =T)
m <- strsplit(file,split="/")
for(i in 1:length(file)){
   methyl_11_cancer <- read.table(file[i], header=T, sep="\t",quote="")
   if(i<12){
         methyl_11_cancer_data <- methyl_11_cancer[,c(3,5,4,14,17:19,6,9)]
   }
   else{
         methyl_11_cancer_data <- methyl_11_cancer[,c(3,5,4,14:17,6,9)]
   }
   methyl_11_cancer_data$name <- paste(methyl_11_cancer_data[,1],methyl_11_cancer_data[,2],sep=":")
   methyl_crossmap <- data.frame(methyl_11_cancer_data[,1],methyl_11_cancer_data[,2],methyl_11_cancer_data[,2],methyl_11_cancer_data$name)
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data")
   cancer_name <- gsub(".txt","",m[[i]][8])
   filenamecsv=paste(cancer_name,"_methyl_crossmap",".bed",sep="")
   write.table(methyl_crossmap, file =filenamecsv,sep="\t",row.names=F,col.names=F,quote=F)
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/methyl_hg38")
   filenamecsv=paste(cancer_name,"_methyl_hg38",".txt",sep="")
   write.table(methyl_11_cancer_data, file =filenamecsv,sep="\t",row.names=F,col.names=T,quote=F)
}
#Coordinate transformation of hg38-hg19 was performed for all methylation data
#####################
#        BRCA       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_BRCA_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_BRCA_methyl.bed
#WGBS
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/WGBS_BRCA_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/WGBS_BRCA_methyl.bed
#####################
#     COADREAD      #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_COAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_COADREAD_methyl.bed
#WGBS
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/WGBS_COAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/WGBS_COADREAD_methyl.bed
#####################
#         HNSC      #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_HNSC_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_HNSC_methyl.bed
#####################
#        KIRC       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_KIRC_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_KIRC_methyl.bed
#####################
#        KIRP       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_KIRP_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_KIRP_methyl.bed
#####################
#        LIHC       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_LIHC_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_LIHC_methyl.bed
#####################
#        LUAD       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_LUAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_LUAD_methyl.bed
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/WGBS_LUAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/WGBS_LUAD_methyl.bed
#####################
#        LUSC       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_LUSC_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_LUSC_methyl.bed
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/WGBS_LUSC_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/WGBS_LUSC_methyl.bed
#####################
#        PRAD       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_PRAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_PRAD_methyl.bed
#####################
#        STAD       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_STAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_STAD_methyl.bed
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/WGBS_STAD_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/WGBS_STAD_methyl.bed
#####################
#        THCA       #
#####################
crossmap bed /pub6/temp/wangli/zcy/0000——yuanshishuju/hg38ToHg19.over.chain.gz /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_data/TCGA_THCA_dms_methyl_crossmap.bed /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_result/TCGA_THCA_methyl.bed
#######################################################################
#                   Integrated processing of methylation data         #
#######################################################################
methyl_hg19 <-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/crossmap_hg19"
file_hg19 <- list.files(methyl_hg19, full.names=T, recursive=T)
methyl_hg38 <-"/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/methyl_hg38"
file_hg38 <- list.files(methyl_hg38, full.names=T, recursive=T)
m <- strsplit(file_hg19,split="/")
for(i in 1:16){
   methyl_11_hg19 <- read.table(file_hg19[i], header=F, sep="\t",quote="")
   colnames(methyl_11_hg19) <- c("methyl_chr","methyl_start","methyl_end","name")
   methyl_11_hg38 <- read.table(file_hg38[i], header=T, sep="\t",quote="")
   methyl_11_hg19_data <- merge(methyl_11_hg19,methyl_11_hg38,by="name") 
   methyl_bedtools <- data.frame(methyl_11_hg19_data[,c(2:4,1)])
   methyl_bedtools$methyl_chr <- gsub("chr","",methyl_bedtools$methyl_chr)
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_data")
   cancer_name <- gsub(".txt","",m[[i]][10])
   filenamecsv=paste(cancer_name,"_bedtools",".bed",sep="")
   write.table(methyl_bedtools, file =filenamecsv,sep="\t",row.names=F,col.names=F,quote=F)
   setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/TCGA_11_CANCER_methyl")
   filenamecsv=paste(cancer_name,"_hg19",".txt",sep="")
   write.table(methyl_11_hg19_data, file =filenamecsv,sep="\t",row.names=F,col.names=T,quote=F)
}

####################################################################
#         Identify dysregulated eRNA----use bedtools               #
####################################################################
#bedtools
cd /pub6/temp/wangli/bedtools2/bin/
#####################
#        BRCA       #
#####################
#Identification of CNV in differential eRNA regions
./bedtools intersect -a bedtools_data/TCGA_BRCA_diff_bedtools.bed -b bedtools_data/TCGA_BRCA_GISTIC2_bedtools.bed -wo > /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_result/diff_cnv/TCGA_BRCA_diff_CNV.bed
#Identification of mutation sites in differential eRNA regions
./bedtools intersect -a bedtools_data/TCGA_BRCA_diff_bedtools.bed -b bedtools_data/TCGA_BRCA_mutation_bedtools.bed -wo > /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_result/diff_mutation/TCGA_BRCA_diff_mutation.bed
#Identification of methylation in differential eRNA regions
./bedtools intersect -a bedtools_data/TCGA_BRCA_diff_bedtools.bed -b bedtools_data/TCGA_BRCA_methyl_bedtools.bed -wo > /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_result/diff_methyl/TCGA_BRCA_diff_methyl.bed
./bedtools intersect -a bedtools_data/TCGA_BRCA_diff_bedtools.bed -b bedtools_data/WGBS_BRCA_methyl_bedtools.bed -wo > /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/TCGA_12_CANCER/bedtools_result/diff_methyl/WGBS_BRCA_diff_methyl.bed
###############################################################
#           第一部分：增强子（遗传驱动的失调增强子）          #
#      1159个差异增强子，2个突变，819个甲基化，2492个CNV      #
###############################################################
#TCGA_BRCA_FANTOM5_60k_eRNA_v2.tsv(eRNA expression profiles of normal and cancer samples from breast cancer)
#TCGA_BRCA_enhancer_profile.txt(active eRNA expression profiles)
#TCGA_BRCA_eRNA_mutation_data.txt(mutation-driven eRNAs)
#TCGA_BRCA_eRNA_methyl_data.txt(methylation-driven eRNAs)
#TCGA_BRCA_eRNA_CNV_data.txt(CNV-driven eRNAs)
#TCGA_BRCA_diff_enhancer.txt(differentially expressed eRNAs)
#############################################################
#    Genetic alterations-driven eRNAs in breast cancer      #
#############################################################
setwd("C:/Users/Dell/Desktop/result1/数据")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile_normal <- enhancer_profile[,grep("normal",colnames(enhancer_profile))]
enhancer_profile_tumor <- enhancer_profile[,grep("tumor",colnames(enhancer_profile))]
enhancer_profile_paired <- enhancer_profile_tumor[,substring(colnames(enhancer_profile_tumor),1,12) %in% substring(colnames(enhancer_profile_normal),1,12)]
enhancer_profile_unpaired_tumor <- enhancer_profile_tumor[,-which(substring(colnames(enhancer_profile_tumor),1,12) %in% substring(colnames(enhancer_profile_normal),1,12))]
enh_profile <- cbind(enhancer_profile_unpaired_tumor,enhancer_profile_paired,enhancer_profile_normal)
enh_profile <- as.matrix(enh_profile)
enh_profile[,1:1208] <- log2(enh_profile[,1:1208])
#-Inf--NA
enh_profile <-apply(enh_profile,2, function(x) {replace(x, is.infinite(x),NA)})
mutation_data <- read.table("TCGA_BRCA_eRNA_mutation_data.txt", header=T, sep="\t")
mutation_heatmap <- data.frame(enhancer=rownames(enh_profile))
mutation_heatmap$Mutation.Status[mutation_heatmap$enhancer %in% mutation_data$enhancer] <- "Mutation"
mutation_heatmap$Mutation.Status[which(is.na(mutation_heatmap$Mutation.Status))] <- "Non_mutation"
methyl_data <- read.table("TCGA_BRCA_eRNA_methyl_data.txt", header=T, sep="\t",quote="")
multiple_enhancer <- factor(methyl_data$enhancer)
methyl_data <- methyl_data[,c("enhancer","Methylation.Difference","Meancase","Meancontrol")]
sum_methyl_data <- data.frame(apply(methyl_data[,-1],2,function(x){tapply(as.numeric(x),multiple_enhancer,sum)}))
methyl_heatmap <- merge(mutation_heatmap,sum_methyl_data,by.x="enhancer",by.y="row.names",all.x=T)
methyl_heatmap$Methylation.Status[methyl_heatmap$Methylation.Difference > 0] <- "Hypermethylation"
methyl_heatmap$Methylation.Status[methyl_heatmap$Methylation.Difference < 0] <- "Hypomethylation"
methyl_heatmap$Methylation.Status[which(is.na(methyl_heatmap$Methylation.Difference))] <- "Non_methylation"
CNV_data <- read.table("TCGA_BRCA_eRNA_CNV_data.txt", header=T, sep="\t")
CNV_data <- unique(CNV_data[,c(8,4)])
m <- data.frame(table(CNV_data$enhancer))
m_28 <- m[which(m[,2]==2),]
CNV_unique <- CNV_data[-which(CNV_data$enhancer %in% m_28[,1]),]
CNV_28 <- data.frame(enhancer=m_28[,1],Type="Amp/Del")
CNV_heatmap <- rbind(CNV_unique,CNV_28) 
colnames(CNV_heatmap)[2] <- "CNV.Status"
CNV_heatmap <- merge(methyl_heatmap[,-c(3:5)],CNV_heatmap,by="enhancer",all.x=T)
CNV_heatmap$CNV.Status[which(is.na(CNV_heatmap$CNV.Status))] <- "Non_CNV"
diff_data <- read.table("TCGA_BRCA_diff_enhancer.txt", header=T, sep="\t")
diff_data$regualation_class[which(diff_data$log2_FC > 0)] <-  "Upregulated"
diff_data$regualation_class[which(diff_data$log2_FC < 0)] <-  "Downregulated"
diff_heatmap <- merge(CNV_heatmap,diff_data[,c(1,5)],by="enhancer",all.x=T)
diff_heatmap$regualation_class[which(is.na(diff_heatmap$regualation_class))] <-  "Non_differential_expression"
type_heatmap <- diff_heatmap[order(diff_heatmap$regualation_class),]
enh_heatmap <- merge(type_heatmap,enh_profile,by.x="enhancer",by.y="row.names",sort=F)
rownames(enh_heatmap) <- enh_heatmap[,1]
enh_heatmap <- as.matrix(enh_heatmap[,-c(1,2,3,4,5)])
save(type_heatmap,file="type_heatmap.RData")
save(enh_heatmap,file="enh_heatmap.RData")
#Genetic alterations-driven eRNAs and differentially expressed eRNAs
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对")
load("type_heatmap.RData")
analysis_enhancer <- type_heatmap[which(type_heatmap$regualation_class!="Non_differential_expression"),]
analysis_enhancer_1 <- analysis_enhancer[which(analysis_enhancer$regualation_class=="Downregulated" & (analysis_enhancer$Methylation.Status=="Hypermethylation" |analysis_enhancer$CNV.Status=="Del" | analysis_enhancer$Mutation.Status=="Mutation")),]  
analysis_enhancer_2 <- analysis_enhancer[which(analysis_enhancer$regualation_class=="Upregulated" & (analysis_enhancer$Methylation.Status=="Hypomethylation" |analysis_enhancer$CNV.Status=="Amp" | analysis_enhancer$Mutation.Status=="Mutation")),]
analysis_enhancer <- rbind(analysis_enhancer_1,analysis_enhancer_2)  
save(analysis_enhancer,file="analysis_enhancer.RData")
##############################################################################################
#    Heatmap of mutation eRNA, methylation eRNA, CNV eRNA and differential eRNA              #
##############################################################################################
#############################################
#      Supplementary Figure S1B--Heatmap    #
#############################################
library(grid)
library(ComplexHeatmap)  
colors <- structure(c("red","gray"), names=c("Mutation","Non_mutation"))
heatmap_mutation <- t(type_heatmap[,2])
colnames(heatmap_mutation) <- type_heatmap[,1]
heatmap_mutation <- as.matrix(heatmap_mutation)
heatmap1 <- Heatmap(heatmap_mutation,
            name="Mutation.Status",    
			col=colors,                         
			row_title="Mutation_enhancer",    
			row_names_side = "left",             
			row_names_gp = gpar(fontsize = 10),
			cluster_rows = F, cluster_columns = F,    
			show_row_names = F,show_column_names = F,use_raster=F)       
colors <- structure(c("#cc340c","#9ec417","gray"),names=c("Upregulated","Downregulated","Non_differential_expression"))
heatmap_regualation <- t(type_heatmap[,5])
colnames(heatmap_regualation) <- type_heatmap[,1]
#anno_histogram，anno_density（密度图），anno_block（盒子），anno_boxplot，anno_barplot
allFC_data <- read.table("TCGA_BRCA_enhancer_allFC.txt", header=T, sep="\t")
allFC <- merge(type_heatmap,allFC_data,by="enhancer",sort=F)
allFC$log2_FC <- log2(allFC$FC)
allFC$log2_FC[is.infinite(allFC$log2_FC)] <- 12
allFC_annotation_data <- as.vector(allFC[,9])
allFC_annotation = HeatmapAnnotation(log2FC = anno_lines(allFC_annotation_data,height=unit(30,"mm")))
heatmap2 <- Heatmap(heatmap_regualation,
            name="regualation_class",    
			col=colors,                         
			row_title="Deregualation_enhancer",    
			row_names_side = "left",             
			row_names_gp = gpar(fontsize = 10),
			cluster_rows = F, cluster_columns = F,    
			show_row_names = F,show_column_names = F,       
			top_annotation = allFC_annotation)
colors <- structure(c("#e20612","#ffd401","gray"),names=c("Hypomethylation", "Hypermethylation","Non_methylation"))
heatmap_Methylation <- t(type_heatmap[,3])
colnames(heatmap_Methylation) <- type_heatmap[,1]
heatmap3 <- Heatmap(heatmap_Methylation,
            name="Methylation.Status",    
			col=colors,                         
			row_title="Methylation_enhancer",    
			row_names_side = "left",             
			row_names_gp = gpar(fontsize = 10),
			cluster_rows = F, cluster_columns = F,    
			show_row_names = F,show_column_names = F)       
colors <- structure(c("#9b3a74","#4da0a0","#f46f20","gray"),names=c("Amp", "Del","Amp/Del","Non_CNV"))
heatmap_CNV <- t(type_heatmap[,4])
colnames(heatmap_CNV) <- type_heatmap[,1]
heatmap4 <- Heatmap(heatmap_CNV,
            name="CNV.Status",    
			col=colors,                         
			row_title="CNV_enhancer",    
			row_names_side = "left",             
			row_names_gp = gpar(fontsize = 10),
			cluster_rows = F, cluster_columns = F,    
			show_row_names = F,show_column_names = F)       			
ht_list <- heatmap2 %v% heatmap3 %v% heatmap4 %v% heatmap1
draw(ht_list, ht_gap = unit(1, "mm")) 
#############################################
#       Supplementary Figure S1D--Heatmap   #
#############################################
library(grid)
library(ComplexHeatmap)  
library(circlize)
library(RColorBrewer)
setwd("C:/Users/Dell/Desktop/result1/数据")
load("type_heatmap.RData") 
load("enh_heatmap.RData")  
sample_type <- gsub("^.*[_]", "", colnames(enh_heatmap))  
sample_type[983:1095] <- "paired_tumor"
sample_type[1:982] <- "unpaired_tumor"
col_fun <- list(sample_type=c("unpaired_tumor"="#0000FF","paired_tumor"="#00FF00","normal"="#FF0000"))
top_annotation <- HeatmapAnnotation(df=data.frame(sample_type = sample_type),col=col_fun) 
type_col = list(regualation_class = c("Upregulated" = "#cc340c", "Downregulated" = "#9ec417",  "Non_differential_expression" = "#44c1f0"),
            Methylation.Status = c("Hypomethylation" = "#e20612", "Hypermethylation" = "#ffd401","Non_methylation"="#00b0eb"),
            CNV.Status = c("Amp" = "#9b3a74", "Del" = "#f46f20","Amp/Del" = "#4da0a0","Non_CNV" = "#44c1f0"),
			Mutation.Status = c("Mutation" = "darkred", "Non_mutation" = "gray")
			)
right_annotation <- rowAnnotation(df=data.frame(regualation_class=type_heatmap$regualation_class,Methylation.Status=type_heatmap$Methylation.Status,CNV.Status=type_heatmap$CNV.Status,Mutation.Status=type_heatmap$Mutation.Status),col=type_col)
enh_heatmap_paired	<- enh_heatmap[,983:1208] 
diffenh_heatmap <- type_heatmap[-which(type_heatmap$regualation_class=="Non_differential_expression"),]
sample_type <- gsub("^.*[_]", "", colnames(enh_heatmap_paired))  
sample_type[1:113] <- "paired_tumor"
sample_type[114:226] <- "paired_normal"
col_fun <- list(sample_type=c("paired_tumor"="#0000FF","paired_normal"="#00FF00"))
top_annotation <- HeatmapAnnotation(df=data.frame(sample_type = sample_type),col=col_fun) 
type_col = list(regualation_class = c("Upregulated" = "#cc340c", "Downregulated" = "#9ec417"),
            Methylation.Status = c("Hypomethylation" = "#e20612", "Hypermethylation" = "#ffd401","Non_methylation"="gray"),
            CNV.Status = c("Amp" = "#9b3a74", "Del" = "#f46f20","Amp/Del" = "#4da0a0","Non_CNV" = "gray"),
			Mutation.Status = c("Mutation" = "darkred", "Non_mutation" = "gray")
			)	
right_annotation <- rowAnnotation(df=data.frame(regualation_class=diffenh_heatmap$regualation_class,Methylation.Status=diffenh_heatmap$Methylation.Status,CNV.Status=diffenh_heatmap$CNV.Status,Mutation.Status=diffenh_heatmap$Mutation.Status),col=type_col)
diffenh_heatmap_paired <- merge(diffenh_heatmap,enh_heatmap_paired,by.x="enhancer",by.y="row.names",sort=F)
rownames(diffenh_heatmap_paired) <- diffenh_heatmap_paired[,1]
diffenh_heatmap_paired <- as.matrix(diffenh_heatmap_paired[,-c(1,2,3,4,5)])	
Heatmap(diffenh_heatmap_paired,name="log2(expression)",na_col = "gray",top_annotation=top_annotation,right_annotation=right_annotation,show_row_names = FALSE, show_column_names = FALSE,cluster_columns = F,cluster_rows = F,use_raster=F)  
###########################################
#            Volcano Plot                 #
###########################################
setwd("C:/Users/Dell/Desktop/result1/数据")
enhancer_allFC <- read.table("TCGA_BRCA_enhancer_allFC.txt", header=T, sep="\t")
enhancer_allFC$FC[which(enhancer_allFC$FC=="Inf")] <- 2000
enhancer_allFC$log10FDR <- -log10(enhancer_allFC$diff_fdr)
enhancer_allFC$log2FC <- log2(enhancer_allFC$FC)
enhancer_allFC$class[enhancer_allFC$diff_fdr<0.05 & enhancer_allFC$log2FC>1] <- "Up"
enhancer_allFC$class[enhancer_allFC$diff_fdr<0.05 & (enhancer_allFC$log2FC< -1)] <- "Down"
enhancer_allFC$class[which(is.na(enhancer_allFC$class))] <- "Normal"
data_volcano <- enhancer_allFC[,c(1,6,5,7)]
colnames(data_volcano)[3] <- '-log10FDR'
data_volcano$class <- factor(data_volcano$class, levels=c('Up','Down','Normal'), order=T)
library(ggplot2)
library(ggrepel)
P_volcano=ggplot(data_volcano,aes(x=log2FC,y=data_volcano[,'-log10FDR']))+
            geom_point(aes(color=class))+
            scale_color_manual(values =c('Up' = 'red', 'Down' = 'blue', 'Normal' = 'grey'))+
            labs(x='log2FC',y='-log10FDR',title="Differential enhancer")+
            geom_hline(yintercept=-log10(0.05),linetype=4)+
            geom_vline(xintercept=c(-1,1),linetype=4)+
            ylim(0,20)+
            theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 18, face = 'bold'),
                    legend.position = 'right',
                    legend.key.size=unit(0.8,'cm'),
                    axis.ticks.x=element_blank(),
                    axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
                    axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
                    axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
                    axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
                    panel.background = element_rect(fill = 'transparent',colour = 'black'),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    plot.background = element_rect(fill = 'transparent',colour = 'black'))
label <- data_volcano[-which(data_volcano$class=="Normal"),]
volcano_label_1 <- label[order(label$log2FC,decreasing=TRUE)[1:10],] 
volcano_label_2 <- label[order(label$log2FC)[1:10],]
volcano_label <- rbind(volcano_label_1,volcano_label_2)
data_volcano$label <- ifelse(data_volcano$enhancer %in% volcano_label$enhancer,data_volcano$enhancer,"")
P_volcano+geom_text_repel(data = data_volcano, aes(x = log2FC, y = data_volcano[,'-log10FDR'],label = data_volcano$label), size = 2,box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"),  segment.color = "black", show.legend = FALSE)
############################################
#                    Circos plot           #
############################################
load("C:/Users/Dell/Desktop/result1/数据/analysis_enhancer.RData")
library(tidyverse)
enhancer <- data.frame(enhancer=unique(analysis_enhancer[,1]))
library(tidyverse)
a <- separate(data=enhancer,col=enhancer,into=c("Chromosome","chromStart","chromEnd",sep="[.]"))
a[,1] <- factor(a[,1])
a[,2] <- as.integer(a[,2])
a[,3] <- as.integer(a[,3])
enhancer[,1] <- factor(enhancer[,1])
enhancer_data <- cbind(a[,-4],enhancer)
a <- merge(enhancer_data,analysis_enhancer,by="enhancer")
#color of eRNA expression values
a$enh_col[a$regualation_class=="Upregulated"]="#cc340c"
a$enh_col[a$regualation_class=="Downregulated"]="#9ec417"
#color of CNV---c("#9b3a74","#4da0a0","#f46f20","gray"),names=c("Amp", "Del","Amp/Del","Non_CNV"))
a$CNV_col="white"
a$CNV_col[a$CNV.Status=="Amp"]= "#9b3a74"
a$CNV_col[a$CNV.Status=="Del"]= "#4da0a0"
a$CNV_col[a$CNV.Status=="Amp/Del"]= "#f46f20"
#color of methylation---c("#e20612","#ffd401","gray"),names=c("Hypomethylation", "Hypermethylation","Non_methylation")
a$methy_col="white"
a$methy_col[a$Methylation.Status=="Hypomethylation"]= "#e20612"
a$methy_col[a$Methylation.Status=="Hypermethylation"]= "#ffd401"
#color of mutation
a$mut_col="white"
a$mut_col[grepl("variant",a$Mutation_effect)]="red"
setwd("C:/Users/Dell/Desktop/result1/数据")
enhancer_allFC <- read.table("TCGA_BRCA_enhancer_allFC.txt", header=T, sep="\t")
enhancer_allFC$FC[which(enhancer_allFC$FC=="Inf")] <- 2000
enhancer_allFC$log10FDR <- -log10(enhancer_allFC$diff_fdr)
enhancer_allFC$log2FC <- log2(enhancer_allFC$FC)
a1 <- merge(a,enhancer_allFC,by="enhancer",all.x=T)
circos_enhancer=a1[,c("enhancer","Chromosome","chromStart","chromEnd","enh_col","CNV_col","methy_col","mut_col","log2FC")]
###########
#  plot   #
###########
library(circlize)
pdf("C:/Users/Dell/Desktop/result1/结果图/1.4_analysis_circos.pdf",height=12, width=12, compress=TRUE)  #from outside to inside:DMR, enhancer,upregulation or downregulation
circos.clear()
par(mar=c(10,10,10,10),mai = c(0,1,0,1))
circos.par(start.degree = 90,cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0),points.overflow.warning = FALSE)
circos.initializeWithIdeogram(plotType = c('ideogram'))    
#Add chromosome
circos.track(ylim = c(0, 1),track.height=0.1, bg.border = 'NA', panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.text(mean(xlim), 0.5, chr, cex = 0.7, col = "black",
              facing = "inside", niceFacing = TRUE)
})
#show linked line,eRNA name  
for(i in seq_len(nrow(circos_enhancer))){ 
  circos.lines(x=c(mean(c(circos_enhancer$chromStart[i],circos_enhancer$chromEnd[i])),mean(c(circos_enhancer$chromStart[i],circos_enhancer$chromEnd[i]))),
               y=c(1,1.5),col = "black",
               sector.index = circos_enhancer$Chromosome[i],track.index = 1)
  circos.text(x=mean(c(circos_enhancer$chromStart[i],circos_enhancer$chromEnd[i])),
                y=c(2,4),col = circos_enhancer$enh_col[i],labels = circos_enhancer$enhancer[i],
                facing = "clockwise", niceFacing = TRUE,cex=0.5,adj = c(0,0.5),
                sector.index = circos_enhancer$Chromosome[i],track.index = 1)
  }  
#show methylation-driven eRNA  
circos.track(ylim=c(0,1),track.height=0.15,bg.col='white',bg.border="black")
for(i in seq_len(nrow(circos_enhancer))){
  circos.rect(xleft = circos_enhancer$chromStart[i],ybottom = 0,xright = circos_enhancer$chromEnd[i],ytop = 1,col=circos_enhancer$methy_col[i],border=circos_enhancer$methy_col[i],sector.index = circos_enhancer$Chromosome[i]) 
  }
#show CNV-driven eRNA 
circos.track(ylim=c(0,1),track.height=0.15,bg.col='white',bg.border="black")
for(i in seq_len(nrow(circos_enhancer))){
  circos.rect(xleft = circos_enhancer$chromStart[i],ybottom = 0,xright = circos_enhancer$chromEnd[i],ytop = 1,col=circos_enhancer$CNV_col[i],border=circos_enhancer$CNV_col[i],sector.index = circos_enhancer$Chromosome[i])  
  }
#show mutation-driven eRNA  
circos.track(ylim=c(0,1),track.height=0.15,bg.col='white',bg.border="black")
for(i in seq_len(nrow(circos_enhancer))){
  circos.rect(xleft = circos_enhancer$chromStart[i],ybottom = 0,xright = circos_enhancer$chromEnd[i],ytop = 1,col=circos_enhancer$mut_col[i],border=circos_enhancer$mut_col[i],sector.index = circos_enhancer$Chromosome[i]) 
  }
 
circos.track(ylim = c(0, 1),track.height=0.3, bg.border = 'NA', panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.lines(xlim,c(0.5,0.5),col = "grey")
})  
#show differentially expressed enhancer
for(i in seq_len(nrow(circos_enhancer))){
  y=circos_enhancer$log2FC[i]/max(abs(circos_enhancer$log2FC))
    circos.lines(x=c(mean(c(circos_enhancer$chromStart[i],circos_enhancer$chromEnd[i])),mean(c(circos_enhancer$chromStart[i],circos_enhancer$chromEnd[i]))),
                 y=c(0.5,0.5+y),type="h",baseline = 0.5,
                 col = ifelse(circos_enhancer$log2FC[i]>0,"red","green"),
                 sector.index = circos_enhancer$Chromosome[i]) 
}

dev.off()
#####################################################
#              DiseaseEnhancer and eRNAs            #
#####################################################
#DiseaseEnhancer database--BrcaEnhancer_loc.csv
setwd("C:/Users/Dell/Desktop/result1/数据/DiseaseEnhancer")
setwd("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/DiseaseEnhancer")
DiseaseEnhancer <- read.table("DiseaseEnhancer.txt",sep="\t")
BrcaEnhancer <- DiseaseEnhancer[which(DiseaseEnhancer[,6]=="Breast cancer"),]
write.table(BrcaEnhancer,"BrcaEnhancer.txt",sep="\t",row.names=F,col.names=F)
BrcaEnhancer_location <- BrcaEnhancer[,c(2:4,1)]
write.table(BrcaEnhancer_location,"BrcaEnhancer_loc.txt",sep="\t",row.names=F,col.names=F,quote=F)
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/type_heatmap.RData")
library(tidyverse)
enhancer_5441 <- separate(data=type_heatmap,col=enhancer,into=c("chr","start","end",sep="[.]"))
enhancer_5441$enhancer <- type_heatmap$enhancer
enhancer_5441 <- enhancer_5441[,c(1:3,9)]
write.table(enhancer_5441,"enhancer_5441.txt",sep="\t",row.names=F,col.names=F,quote=F)
cd /pub6/temp/wangli/bedtools2/bin/
./bedtools intersect -a /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/DiseaseEnhancer/enhancer_5441.bed -b /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/DiseaseEnhancer/BrcaEnhancer_loc.bed -wo > /pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/DiseaseEnhancer/DiseaseEnhancer_brca_5441.bed
DiseaseEnhancer_brca_5441 <- read.table("DiseaseEnhancer_Brca_5441.bed",sep="\t")
colnames(DiseaseEnhancer_brca_5441) <- c("chr","start","end","enhancer","DE_chr","DE_start","DE_end","DE_id")
#genetic alteration-driven eRNAs----BRCA-eRNAs 
load("/pub6/temp/wangli/zcy/0000——BISHE/1.识别遗传改变驱动的增强子活性图谱/BRCA/增强子-基因对/analysis_enhancer.RData")
DiseaseEnhancer_analysis <- merge(DiseaseEnhancer_brca_5441,analysis_enhancer,by="enhancer")
####################################################################
#   Venn diagram--genetic alteration-driven eRNAs and BRCA-eRNAs   #
####################################################################
q=6     
m=431     
n=5441-m    
k=27      
phyper(q-1,m,n,k,lower.tail=F) #0.01686754
venn.plot <- draw.pairwise.venn(431,27,6,c("Specific_enhancer", "DiseaseEnhancer_BRCA_enhancer"), 
  fill = c("red", "green"), 
  scaled = FALSE,
  main = "Venn Diagram of Enhancer",
  lty = "blank", 
  cex = 2, 
  cat.cex = 1.5, 
  cat.pos = c(340, 0), 
  cat.dist = c(0.03,0.03), 
  cat.just = list(c(0, 0), c(0, 0)), 
  ext.pos = 0, 
  ext.dist = -0.05, ext.length = 0.85, 
  ext.line.lwd = 2, ext.line.lty = "dashed", 
  alpha=0.3, euler.d=T  
)
