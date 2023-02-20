#######################################################################################
#          Univariate COX analysis of eRNAs--------cbioportal                           #
#######################################################################################
library(cgdsr)
library(DT) 
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds)
all_TCGA_studies <- getCancerStudies(mycgds)
DT::datatable(all_TCGA_studies)
hnscc2018 <- "brca_tcga"
mycaselist <- getCaseLists(mycgds,hnscc2018)[1,1]
myclinicaldata <- getClinicalData(mycgds,mycaselist) 
save(myclinicaldata,file="E:/生存分析/数据和结果/myclinicaldata.Rdata")
#Retain the required survival data---"AGE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_TUMOR_PATHOLOGIC_PT","DFS_MONTHS","DFS_STATUS","ER_STATUS_BY_IHC","IHC_HER2","HER2_FISH_STATUS","OS_MONTHS","OS_STATUS","PRIMARY_SITE_PATIENT","PR_STATUS_BY_IHC","RACE","SEX"
select <- c("OS_MONTHS","OS_STATUS","DFS_MONTHS","DFS_STATUS","AGE","SEX","ER_STATUS_BY_IHC","IHC_HER2","PR_STATUS_BY_IHC","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT","AJCC_PATHOLOGIC_TUMOR_STAGE","MENOPAUSE_STATUS","HISTOLOGICAL_SUBTYPE","LYMPH_NODES_EXAMINED","MARGIN_STATUS_REEXCISION","RADIATION_TREATMENT_ADJUVANT","HISTORY_NEOADJUVANT_TRTYN")
brca_clinicaldata <- myclinicaldata[,select]
brca_clinicaldata[which(brca_clinicaldata$OS_MONTHS<0),c(1,3)] <- 0
rownames(brca_clinicaldata) <- gsub("[.]01$","",rownames(brca_clinicaldata))
#The status indicator, normally 0=alive, 1=dead
brca_clinicaldata$OS_STATUS <- gsub(":DECEASED","",brca_clinicaldata$OS_STATUS)
brca_clinicaldata$OS_STATUS <- gsub(":LIVING","",brca_clinicaldata$OS_STATUS)
brca_clinicaldata$DFS_STATUS <- gsub(":DiseaseFree","",brca_clinicaldata$DFS_STATUS)
brca_clinicaldata$DFS_STATUS <- gsub(":Recurred/Progressed","",brca_clinicaldata$DFS_STATUS)
colnames(brca_clinicaldata) <- c("OS_MONTHS","OS_STATUS","DFS_MONTHS","DFS_STATUS","AGE","SEX","IHC_ER","IHC_HER2","IHC_PR","AJCC_M","AJCC_N","AJCC_T","AJCC_STAGE","menopause status","histological subtype","lymph_nodes","margin_status","radiation_treatment","neoadjuvant_treatment")
##Delete the data whose survival state is NA
brca_clinicaldata <- brca_clinicaldata[-which(brca_clinicaldata$OS_STATUS==""),]
brca_clinicaldata <- brca_clinicaldata[-grep("[.]06",rownames(brca_clinicaldata)),]
save(brca_clinicaldata,file="C:/Users/Dell/Desktop/result2——生存分析/数据/brca_clinicaldata.Rdata") 
#####################################################################################
#            Count the number of samples for each category of survival data         #
#####################################################################################
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
load("brca_clinicaldata.Rdata")
unique(brca_clinicaldata$IHC_HER2)
unique(brca_clinicaldata$IHC_ER)
unique(brca_clinicaldata$IHC_PR)
unique(brca_clinicaldata$AJCC_M)
unique(brca_clinicaldata$AJCC_N)
unique(brca_clinicaldata$AJCC_T)
unique(brca_clinicaldata$AJCC_STAGE)
##################################################################################
med <- median(brca_clinicaldata[-which(is.na(brca_clinicaldata$AGE)),]$AGE)   #58
medage <- vector(mode="logical",length=length(brca_clinicaldata$AGE))
brca_clinicaldata$medage[which(brca_clinicaldata$AGE >= med)] <- ">=58"
brca_clinicaldata$medage[which(brca_clinicaldata$AGE < med)] <- "<58"
brca_clinicaldata$medage[which(medage=="")] <- NA
brca_clinicaldata$AJCC_STAGE[grep("Stage III[A,B,C]|Stage IV|Stage III$",brca_clinicaldata$AJCC_STAGE)] <- "Stage 3|Stage 4"
brca_clinicaldata$AJCC_STAGE[grep("Stage I[A,B,C]|Stage II[A,B,C]|Stage I$|Stage II$",brca_clinicaldata$AJCC_STAGE)] <- "Stage 1|Stage 2"
brca_clinicaldata$AJCC_STAGE[grep("Stage X$",brca_clinicaldata$AJCC_STAGE)] <- NA
brca_clinicaldata$AJCC_STAGE[which(brca_clinicaldata$AJCC_STAGE=="")] <- NA
brca_clinicaldata$AJCC_M[grep("MX",brca_clinicaldata$AJCC_M)] <- NA
brca_clinicaldata$AJCC_M[grep("cM0",brca_clinicaldata$AJCC_M)] <- "M0"
brca_clinicaldata$AJCC_N[grep("N0|N1",brca_clinicaldata$AJCC_N)] <- "N0|N1"
brca_clinicaldata$AJCC_N[grep("N2|N3",brca_clinicaldata$AJCC_N)] <- "N2|N3"
brca_clinicaldata$AJCC_N[grep("NX",brca_clinicaldata$AJCC_N)] <- NA
brca_clinicaldata$AJCC_T[grep("T1|T2",brca_clinicaldata$AJCC_T)] <- "T1|T2"
brca_clinicaldata$AJCC_T[grep("T3|T4",brca_clinicaldata$AJCC_T)] <- "T3|T4"
brca_clinicaldata$AJCC_T[grep("TX",brca_clinicaldata$AJCC_T)] <- NA
brca_clinicaldata$IHC_HER2[grep("Equivocal|Indeterminate",brca_clinicaldata$IHC_HER2)] <- NA
brca_clinicaldata$IHC_ER[grep("Indeterminate",brca_clinicaldata$IHC_ER)] <- NA
brca_clinicaldata$IHC_PR[grep("Indeterminate",brca_clinicaldata$IHC_PR)] <- NA
###############################################################################################################
#            factor conversion is performed on the covariable, and the first term controlling levels is the reference term                  #
###############################################################################################################
brca_clinicaldata$SEX <- factor(brca_clinicaldata$SEX,levels=c("Male","Female"))
brca_clinicaldata$medage <- factor(brca_clinicaldata$medage,levels=c("<58",">=58"))
brca_clinicaldata$AJCC_M <- factor(brca_clinicaldata$AJCC_M,levels=c("M0","M1"))
brca_clinicaldata$AJCC_N <- factor(brca_clinicaldata$AJCC_N,levels=c("N0|N1","N2|N3"))
brca_clinicaldata$AJCC_T <- factor(brca_clinicaldata$AJCC_T,levels=c("T1|T2","T3|T4"))
brca_clinicaldata$AJCC_STAGE <- factor(brca_clinicaldata$AJCC_STAGE,levels=c("Stage 1|Stage 2","Stage 3|Stage 4"))
brca_clinicaldata$IHC_PR <- factor(brca_clinicaldata$IHC_PR,levels=c("Positive","Negative"))
brca_clinicaldata$IHC_HER2 <- factor(brca_clinicaldata$IHC_HER2,levels=c("Positive","Negative"))
brca_clinicaldata$IHC_ER <- factor(brca_clinicaldata$IHC_ER,levels=c("Positive","Negative"))
brca_clinicaldata$OS_STATUS <- as.numeric(brca_clinicaldata$OS_STATUS)
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
save(brca_clinicaldata,file="brca_clinicaldata_new.Rdata")
################################################################################
#     Univariate cox regression was performed for factors of survival data     #
################################################################################
library(survival)
library("ggplot2")
library("ggpubr")
library("survminer")
library(survival)
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
load("brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile) <- gsub("_tumor","",colnames(enhancer_profile))
brca_clinicaldata_1094 <-  brca_clinicaldata[intersect(rownames(brca_clinicaldata),colnames(enhancer_profile)),]
enhancer_profile_1094 <- enhancer_profile[,intersect(rownames(brca_clinicaldata),colnames(enhancer_profile))]
setwd("C:/Users/Dell/Desktop/result2——生存分析/结果")
write.table(enhancer_profile_1094,"enhancer_profile_1094.txt",sep="\t",row.names=T,col.names=T,quote=F)
write.table(brca_clinicaldata_1094,"brca_clinicaldata_1094.txt",sep="\t",row.names=T,col.names=T,quote=F)
covariates <- c("medage","SEX","IHC_ER","IHC_HER2","IHC_PR","AJCC_M","AJCC_N","AJCC_T","AJCC_STAGE")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS_MONTHS,OS_STATUS)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = brca_clinicaldata_1094)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-round(x$coefficients[ ,5], digits = 4)
                         HR <-round(x$coef[2], digits=4);
                         HR.confint.lower <- round(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- round(x$conf.int[,"upper .95"],4)
                         CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                         res_1094 <- c(p.value,HR,CI)
                         names(res_1094) <- c("p.value","HR", "95% CI for HR")
                         return(res_1094)
                       })
res_1094 <- t(as.data.frame(univ_results, check.names = FALSE))

################################################################################
#           Univariate cox regression was performed for eRNAs                  #
################################################################################
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
library(survival)
load("brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile) <- gsub("_tumor","",colnames(enhancer_profile))
load("C:/Users/Dell/Desktop/result2——生存分析/数据/analysis_enhancer.RData")
enhancer_p <- enhancer_profile[which(rownames(enhancer_profile) %in% analysis_enhancer$enhancer),]
enhancer_p <- scale(enhancer_p) 
os_cox_scale <- data.frame(NA,NA,NA,NA)
colnames(os_cox_scale) <- c("enhancer", "HR", "95% CI for HR", "p.value")
for(i in 1:431){
enh_exp <- t(enhancer_p[i,])
enh_exp_1 <- data.frame(enhancer_name=enh_exp[1,])
follow_enh <- merge(enh_exp_1,brca_clinicaldata,by="row.names")
follow_enh$enhancer_name <- as.numeric(follow_enh$enhancer_name)
cox_enh <- coxph(Surv(OS_MONTHS,OS_STATUS)~ enhancer_name,data=follow_enh)
tmp <- summary(cox_enh)
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
    HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"], digits = 4);
    CI <- paste0(HR.confint.lower, "-", HR.confint.upper);	
	variate <- rownames(enhancer_p)[i];
    os_cox.result <- as.data.frame(cbind(variate, HR, CI, p.value));
    colnames(os_cox.result) <- c("enhancer", "HR", "95% CI for HR", "p.value");
os_cox_scale <- rbind(os_cox_scale,os_cox.result)
}
os_cox_scale <- os_cox_scale[-1,]
os_cox_scale_0.05 <- os_cox_scale[which(os_cox_scale$p.value<0.05),]   
setwd("C:/Users/Dell/Desktop/result2——生存分析/结果")
save(os_cox_scale_0.05,file="os_cox_scale_0.05.RData")
###################################
#           log-rank test         #
###################################
#Log-rank
library(survminer)
library(survival)
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
load("brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile) <- gsub("_tumor","",colnames(enhancer_profile))
setdiff(colnames(enhancer_profile),rownames(brca_clinicaldata))  #"TCGA.BH.A0B2"
enhancer_p <- enhancer_profile[which(rownames(enhancer_profile) %in% analysis_enhancer$enhancer),]
log_rank_431_enh <- data.frame(enhancer=NA,log_rank_p=NA)
for(i in 1:431){
   enh_exp <- t(enhancer_p[i,])
   enh_exp_1 <- data.frame(enhancer_name=enh_exp[,1])
   follow_enh <- merge(brca_clinicaldata,enh_exp_1,by="row.names")
   median_exp <- median(as.numeric(follow_enh[,22]))
   #1 indicates that the eRNA expression value is higher than the median expression value
   #0 indicates that the eRNA expression value is lower than the median expression value
   follow_enh$median_exp_label[follow_enh[,22] > median_exp] <- 1
   follow_enh$median_exp_label[follow_enh[,22] <= median_exp] <- 0
   #Log-rank 
   log_rank_enh <- survdiff(Surv(OS_MONTHS,OS_STATUS)~ median_exp_label,data=follow_enh)
   log_rank <- 1 - pchisq(log_rank_enh$chisq, length(log_rank_enh$n) - 1)
   log_rank_1 <- data.frame(enhancer=rownames(enhancer_p)[i],log_rank_p=log_rank)
   log_rank_431_enh <- rbind(log_rank_431_enh,log_rank_1)
}
log_rank_431_enh <- log_rank_431_enh[-1,]
log_rank_431_0.05 <- log_rank_431_enh[which(log_rank_431_enh$log_rank_p<0.05),] 
setwd("C:/Users/Dell/Desktop/result2——生存分析/结果")
save(log_rank_431_0.05,file="log_rank_431_0.05.RData")
#################################################################################
#      Intersection eRNAs of univariate cox regression and log-rank test        #
#################################################################################
cox_logrank <- intersect(os_cox_scale_0.05$enhancer,log_rank_431_0.05$enhancer)
[1] "chr10:3848037-3849536"    "chr10:5530637-5531333"    "chr6:118759105-118759218"
[4] "chr8:116990613-116990894" "chr8:117175146-117175377" "chr8:117235182-117235423"
setwd("C:/Users/Dell/Desktop/result2——生存分析/结果")
save(cox_logrank,file="cox_logrank.RData")
##################################
#           Survival curve       #
##################################
load("C:/Users/Dell/Desktop/result2——生存分析/结果/log_rank_431_0.05.Rdata")
load("C:/Users/Dell/Desktop/result2——生存分析/数据/brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("C:/Users/Dell/Desktop/result2——生存分析/结果/TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile) <- gsub("_tumor","",colnames(enhancer_profile))
for(i in 1:length(unique(cox_logrank))){
       logrank_p_analysis_enh <- log_rank_p_enh[which(log_rank_p_enh$enhancer %in% cox_logrank[i]),]
       enh_exp <- t(enhancer_profile[which(rownames(enhancer_profile)==cox_logrank[i]),])
       enh_exp_1 <- data.frame(enhancer_name=enh_exp[,1])
       follow_enh <- merge(brca_clinicaldata,enh_exp_1,by="row.names")
	   follow_enh$OS_Year <- follow_enh$OS_MONTHS/12
	   analysis_os <- follow_enh %>% dplyr::select(OS_Year,OS_STATUS,enhancer_name)
	   analysis_os <- na.omit(analysis_os)
	   analysis_os["enhancer_name"] <- ifelse(analysis_os["enhancer_name"] > median(analysis_os[,"enhancer_name"]), 'high', 'low' )
       km <- survfit(Surv(as.numeric(analysis_os$OS_Year),analysis_os$OS_STATUS)~ enhancer_name,data=analysis_os)
	   #setwd("C:/Users/Dell/Desktop/result2——生存分析/结果图")
	   #pdfname <- paste("new_KM_",i,".pdf",sep="")
	   #pdf(pdfname,height=10,width=10)
       ggsurvplot(km, data = analysis_os, pval=TRUE,legend.title = cox_logrank[i], font.main = c(16, "bold", "darkblue"),xlab="Time in years",ylab="Overall survival",xlim = c(0,20),break.time.by = 4,pval.size = 8,risk.table = "absolute", risk.table.y.text.col = T,risk.table.y.text = FALSE,risk.table.fontsize = 5,legend.labs =  c("high", "Low"),  palette = c("#E41A1C","#377EB8"),font.y  = c(16, "bold"),font.x  = c(16, "bold"),.legend = "top",font.legend = c(16, "bold"),
	   surv.median.line = "hv", 
       conf.int = TRUE) 
       #dev.off()
}
###################################################################
#     Correlation analysis between eRNAs and clinical features    #
###################################################################
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
load("brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile_tumor <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile_tumor) <- gsub("_tumor","",colnames(enhancer_profile_tumor))
#six eRNAS
load("C:/Users/Dell/Desktop/result2——生存分析/结果/cox_logrank.RData")
analysis_enhancer_cox_logrank_profile <- enhancer_profile[which(rownames(enhancer_profile) %in% cox_logrank),]
library(ggplot2)
library(ggpubr)
###########################
#       箱式图+点图       #
###########################
for(i in 1:length(cox_logrank)){
   #tumor and normal
   enh <- rownames(analysis_enhancer_cox_logrank_profile)[i]
   eRNA_profile <- enhancer_profile
   normal_exp <- eRNA_profile[,grep("normal",colnames(eRNA_profile))]  
   tumor_exp <- eRNA_profile[,grep("tumor",colnames(eRNA_profile))]    
   normal_sample <- gsub("_normal","",colnames(analysis_enh))
   tumor_sample <- gsub("_tumor","",colnames(analysis_enh))
   match_sample <- intersect(tumor_sample,normal_sample)
   match_tumor_exp <- tumor_exp[,pmatch(match_sample,colnames(tumor_exp))]
   match_normal_exp <- normal_exp[,pmatch(match_sample,colnames(normal_exp))]
   paired_expr_data <- cbind(match_normal_exp,match_tumor_exp)
   analysis_enh <- paired_expr_data[which(rownames(paired_expr_data)==enh),]   
   expr_data <- data.frame(t(analysis_enh))
   expr_data$sample[grepl("tumor",rownames(expr_data))] <- "tumor"
   expr_data$sample[grepl("normal",rownames(expr_data))] <- "normal"
   colnames(expr_data)[1] <- "log2(FPKM+1)"
   expr_data[,1] <- log2(expr_data[,1]+1)
   p_sample <- ggboxplot(expr_data,y="log2(FPKM+1)",x="sample",color="sample",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means(method = "t.test",paired=T)
   #age
   analysis_enh2 <- enhancer_profile_tumor[which(rownames(enhancer_profile_tumor)==enh),]
   expr_data2 <- data.frame(t(analysis_enh2))
   expr_data3 <- merge(expr_data2,brca_clinicaldata,by="row.names")
   colnames(expr_data3)[2] <- "log2(FPKM+1)"
   expr_data3[,2] <- log2(as.numeric(expr_data3[,2])+1)
   p_age <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="medage",color="medage",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #sex
   p_sex <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="SEX",color="SEX",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #ER
   p_ER <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="IHC_ER",color="IHC_ER",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #HER2
   p_HER2 <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="IHC_HER2",color="IHC_HER2",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #PR
   p_PR <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="IHC_PR",color="IHC_PR",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #M
   p_M <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="AJCC_M",color="AJCC_M",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #N
   p_N <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="AJCC_N",color="AJCC_N",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #T
   p_T <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="AJCC_T",color="AJCC_T",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   #STAGE
   p_STAGE <- ggboxplot(expr_data3,y="log2(FPKM+1)",x="AJCC_STAGE",color="AJCC_STAGE",add="jitter",ylab=paste(enh,"\n","log2(FPKM+1)"),xlab="")+stat_compare_means()
   p_result <- ggarrange(p_sample,p_age,p_sex,p_ER,p_HER2,p_PR,p_T,p_N,p_M,p_STAGE,ncol = 5, nrow = 2)
   #ggpar(a,font.title = c(1),font.x = c(1),font.y = c(1))
   setwd("C:/Users/Dell/Desktop/result2——生存分析/结果图2/表达值和临床因素——箱式图+点图_最新版")
   pdfname <- paste("表达值和临床因素_",i,".pdf",sep="")
   ggsave(pdfname,p_result,height=10,width=25)
}
######################################################################
#       Multivariate cox proportional risk regression---6 eRNAs      #
######################################################################
setwd("C:/Users/Dell/Desktop/result2——生存分析/数据")
load("brca_clinicaldata_new.Rdata")
enhancer_profile <- read.table("TCGA_BRCA_enhancer_profile.txt", header=T, sep="\t")
enhancer_profile <- enhancer_profile[,-grep("normal",colnames(enhancer_profile))]
colnames(enhancer_profile) <- gsub("_tumor","",colnames(enhancer_profile))
enhancer_p <- enhancer_profile[which(rownames(enhancer_profile) %in% analysis_enhancer$enhancer),]
enhancer_p <- scale(enhancer_p)
load("C:/Users/Dell/Desktop/result2——生存分析/结果/cox_logrank.RData")
load("C:/Users/Dell/Desktop/result2——生存分析/结果/os_cox_scale_0.05.RData")
os_cox_scale_6 <- os_cox_scale_0.05[which(os_cox_scale_0.05$enhancer %in% cox_logrank),]
enhancer_p <- enhancer_profile[which(rownames(enhancer_profile) %in% analysis_enhancer$enhancer),]
enhancer_p <- scale(enhancer_p)
enhancer_6 <- enhancer_p[which(rownames(enhancer_p) %in% cox_logrank),] 
follow_enh <- merge(t(enhancer_6),brca_clinicaldata,by="row.names")
#Multivariate cox proportional risk regression----eRNAs and Clinical factor
mul_cox <- coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS))~SEX+medage+IHC_ER+IHC_PR+IHC_HER2+AJCC_M+AJCC_N+AJCC_T+AJCC_STAGE+`chr10:3848037-3849536`+`chr10:5530637-5531333`+`chr6:118759105-118759218`+`chr8:116990613-116990894`+`chr8:117175146-117175377`+`chr8:117235182-117235423`,data=follow_enh)
mul_cox1 <- summary(mul_cox)
HR <- round(mul_cox1$coefficients[ ,2], digits = 4);
p.value <- round(mul_cox1$coefficients[ ,5], digits = 4);
HR.confint.lower <- round(mul_cox1$conf.int[ ,"lower .95"], digits = 4);
HR.confint.upper <- round(mul_cox1$conf.int[ ,"upper .95"], digits = 4);
CI <- paste0(HR," (",HR.confint.lower, "-", HR.confint.upper,")");
eRNA6_clinical_result <- data.frame(cbind(p.value,HR,HR.confint.lower,HR.confint.upper,CI))	
write.table(eRNA6_clinical_result,"C:/Users/Dell/Desktop/result2——生存分析/结果/eRNA6_clinical_result.txt",sep="\t",row.names=T,,quote=F)
#Multivariate cox proportional risk regression----eRNAs
a <- follow_enh[,c(8,9,2:7)]
write.table(a, file ="C:/Users/Dell/Desktop/result2——生存分析/结果/six_enhancer_surdata_2.csv",sep=",",row.names=F)
mul_cox <- coxph(Surv(OS_MONTHS,as.numeric(OS_STATUS))~`chr10:3848037-3849536`+`chr6:118759105-118759218`+`chr8:116990613-116990894`,data=a)
ggforest(mul_cox,data=a,main="hazard ratio",
         cpositions=c(0.02,0.22,0.4),
         fontsize=0.8,refLabel="reference")	 
write.table(follow_enh,"C:/Users/Dell/Desktop/result2——生存分析/结果/eRNA6_clinical_new_data.txt",sep="\t",row.names=F,,quote=F)
