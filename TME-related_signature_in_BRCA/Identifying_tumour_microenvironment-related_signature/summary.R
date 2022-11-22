#estimate was used to calculate the corresponding immune score, stromal score, and estimated score
library(estimate)
filterCommonGenes(input.f = "symbol.txt",
                  output.f = "Estimate_gene.gct",
                  id = "GeneSymbol")

estimateScore(input.ds = "Estimate_gene.gct",
              output.ds = "Estimate_score.gct",
              platform = "illumina")
scores=read.table("Estimate_score.gct",skip=2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
#Distribution and cumulative distribution of immune scores, stromal scores, and estimated scores between tumor and normal samples
library(ggplot2)
a<-scores
a<-data.frame(a)
library(stringr)
a$group<-ifelse(as.numeric(str_sub(row.names(a),14,15))<10,"tumor","normal")

a_tumor<-a[a$group=="tumor",]
a_normal<-a[a$group=="normal",]
library(reshape2)
library(ggpubr)
ameta<-melt(a,id.vars = "group")
wilcoxscore<-c(wilcox.test(a[a$group=="tumor","StromalScore"],
                           a[a$group=="normal","StromalScore"])$p.value,
              wilcox.test(a[a$group=="tumor","ImmuneScore"],
                          a[a$group=="normal","ImmuneScore"])$p.value,
              wilcox.test(a[a$group=="tumor","ESTIMATEScore"],
                          a[a$group=="normal","ESTIMATEScore"])$p.value)

anno_a<-data.frame(class=unique(ameta$variable),p=wilcoxscore)
anno_a$h<-c(max(a$StromalScore),max(a$ImmuneScore),max(a$ESTIMATEScore))
anno_a$p=ifelse(anno_a$p<0.001,paste("p<0.001"),round(anno_a$p,3))
p1<-ggplot(data = ameta, aes(x = variable, y = value, fill = group)) + 
  geom_violin(position = position_dodge(width = 1), scale = 'width',trim = FALSE) +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.7, width = 0.2, show.legend = FALSE) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), 
        legend.title = element_blank(), legend.key = element_blank()) +
  labs(x = '', y = 'Score')+
  geom_text(inherit.aes = F,data = anno_a,aes(x=class,y=h,label=p),
            size=2.5,fontface="bold")+
  annotate("text",x=anno_a$class,y=(anno_a$h-300),label="————",vjust=-0.3)+
  scale_fill_manual(values = c("tumor"="orange","normal"="blue"))

ggsave("difference scores.png",p1,width = 8,height = 6)
ggsave("difference scores.pdf",p1,width = 8,height = 6)

p1<-ggplot(ameta, aes(x=value, col=variable,linetype=group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2) + 
  theme_bw() +
  theme(legend.position=c(.75,.25),
        panel.grid = element_blank()) +
  labs(x="Score",
       y="Accumulative proportion",
       col="",
       title = "")+
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"))
  
ggsave("cumulative distribution.png",p1,width = 8,height = 6)
ggsave("cumulative distribution.pdf",p1,width = 8,height = 6)

#Breast cancer samples were grouped based on the median immune score and stromal score, followed by DEseq2 differential analysis
rnaseq<-as.data.frame(scores)
immune_high<-subset(rnaseq,rnaseq$ImmuneScore>median(rnaseq$ImmuneScore))
immune_low<-subset(rnaseq,rnaseq$ImmuneScore<median(rnaseq$ImmuneScore))
stromal_high<-subset(rnaseq,rnaseq$StromalScore>median(rnaseq$StromalScore))
stromal_low<-subset(rnaseq,rnaseq$StromalScore<median(rnaseq$StromalScore))
im_high_symbol<-data_symbol[,rownames(immune_high)]
im_low_symbol<-data_symbol[,rownames(immune_low)]
st_high_symbol<-data_symbol[,rownames(stromal_high)]
st_low_symbol<-data_symbol[,rownames(stromal_low)]

library(stringr)
groop_list1<- ifelse(as.numeric(str_sub(colnames(im_high_symbol),14,15))<10,"tumor","normal")
groop_list1 <- factor(groop_list1,levels = c("normal","tumor"))
table(groop_list1)
colData1<-data.frame(row.names = colnames(im_high_symbol),
                     condition=groop_list1)
groop_list2<- ifelse(as.numeric(str_sub(colnames(im_low_symbol),14,15))<10,"tumor","normal")
groop_list2<- factor(groop_list2,levels = c("normal","tumor"))
table(groop_list2)
colData2<-data.frame(row.names = colnames(im_low_symbol),
                     condition=groop_list2)
groop_list3<- ifelse(as.numeric(str_sub(colnames(st_high_symbol),14,15))<10,"tumor","normal")
groop_list3<- factor(groop_list3,levels = c("normal","tumor"))
table(groop_list3)
colData3<-data.frame(row.names = colnames(st_high_symbol),
                     condition=groop_list3)
groop_list4<- ifelse(as.numeric(str_sub(colnames(st_low_symbol),14,15))<10,"tumor","normal")
groop_list4<- factor(groop_list4,levels = c("normal","tumor"))
table(groop_list4)
colData4<-data.frame(row.names = colnames(st_low_symbol),
                     condition=groop_list4)
save(im_high_symbol,im_low_symbol,st_high_symbol,st_low_symbol,groop_list1,groop_list2,groop_list3,groop_list4,colData1,colData2,colData3,colData4,file = "DESeq.Rdata")


library(stringr)
groop_list<- ifelse(as.numeric(str_sub(row.names(rnaseq),14,15))<10,"tumor","normal")
group_list <- factor(groop_list,levels = c("normal","tumor"))
table(group_list)
colData<-data.frame(row.names = row.names(rnaseq),
                    condition=group_list)
cg1<-row.names(colData)[colData$condition=="tumor"]
score_tumor<-rnaseq[cg1,]
data_tumor<-data_symbol[,cg1]

score_tumor$Immune<- ifelse(score_tumor$ImmuneScore>median(score_tumor$ImmuneScore),"high","low")
table(score_tumor$Immune)
colData1<-data.frame(row.names = row.names(score_tumor),
                     condition=score_tumor$Immune)

score_tumor$Stromal<- ifelse(score_tumor$StromalScore>median(score_tumor$StromalScore),"high","low")
table(score_tumor$Stromal)
colData2<-data.frame(row.names = row.names(score_tumor),
                     condition=score_tumor$Stromal)

suppressPackageStartupMessages(library(DESeq2))
dds1<-DESeqDataSetFromMatrix(countData = data_tumor,colData = colData1,design = ~condition)
dds1_2<-DESeq(dds1)
res1<-results(dds1_2)
DE_im<-as.data.frame(res1)
DE_im$change=as.factor(ifelse(DE_im$padj < 0.05 & abs(DE_im$log2FoldChange) > 1 ,ifelse(DE_im$log2FoldChange > 1  ,'UP','DOWN'),'NOT'))
table(DE_im$change)

dds2<-DESeqDataSetFromMatrix(countData = data_tumor,colData = colData2,design = ~condition)
dds2_2<-DESeq(dds2)
res2<-results(dds2_2)
DE_st<-as.data.frame(res2)
DE_st$change=as.factor(ifelse(DE_st$padj < 0.05 & abs(DE_st$log2FoldChange) > 1 ,ifelse(DE_st$log2FoldChange > 1  ,'UP','DOWN'),'NOT'))
table(DE_st$change)
#Volcano maps, Venn maps and heat maps were created based on differential genes from immune scores and stromal scores
im_up<-row.names(DE_im)[DE_im$change=="UP"]
st_up<-row.names(DE_st)[DE_st$change=="UP"]
im_down<-row.names(DE_im)[DE_im$change=="DOWN"]
st_down<-row.names(DE_st)[DE_st$change=="DOWN"]
im_up<-na.omit(im_up)
st_up<-na.omit(st_up)
im_down<-na.omit(im_down)
st_down<-na.omit(st_down)
up_inter<-intersect(im_up,st_up)
down_inter<-intersect(im_down,st_down)
up_inter<-data.frame(gene=up_inter)
down_inter<-data.frame(gene=down_inter)
write.csv(up_inter,"up.csv")
write.csv(down_inter,"down.csv")
DEG<-rbind(up_inter,down_inter)
write.csv(DEG,file = "DEG.csv")

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(airway)
EnhancedVolcano(DE_im,
                lab = rownames(DE_im),
                x = 'log2FoldChange',
                y = 'padj',
                ylab = 'pvalue',
                title = "Immune Score",
                subtitle = "",
                xlim = c(-8, 8))

EnhancedVolcano(DE_st,
                lab = rownames(DE_st),
                x = 'log2FoldChange',
                y = 'padj',
                ylab = 'pvalue',
                title = "Stromal Score",
                subtitle = "",
                xlim = c(-8, 8))

library(grid)
library(futile.logger)
library(VennDiagram)
vene_list1<-list(Immune = im_up,Stromal = st_up)
venn.diagram(vene_list1, filename = 'up-regulated gene.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 3, fontfamily = 'serif', 
             cat.cex = 0, cat.fontfamily = 'serif')

vene_list2<-list(Immune = im_down,Stromal = st_down)
venn.diagram(vene_list2, filename = 'down-regulated gene.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 3, fontfamily = 'serif', 
             cat.cex = 0, cat.fontfamily = 'serif')

im_symbol<-data_tumor[c(im_up,im_down),]
library(stringr)
score_tumor$Immune<- ifelse(score_tumor$ImmuneScore>median(score_tumor$ImmuneScore),"high","low")
table(score_tumor$Immune)
colData1<-data.frame(row.names = row.names(score_tumor),
                     condition=score_tumor$Immune)
colData1_1<-data.frame(row.names = c(row.names(colData1)[colData1$condition=="high"],
                                     row.names(colData1)[colData1$condition=="low"]),
                       condition=factor(c(rep(c("high"),554),rep(c("low"),555))))

im_symbol1<-cbind(im_symbol[,row.names(colData1)[colData1$condition=="high"]],
                  im_symbol[,row.names(colData1)[colData1$condition=="low"]])
#The top20 differentially expressed genes in up-regulated and down-regulated genes were extracted as examples to draw expression heat map
im_symbol2<-im_symbol1[which(apply(im_symbol1,1,
                                   function(x){return(sum(x>10))})>ncol(im_symbol1)*0.25),]                              
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#436eee','white','#EE0000'))(100)
n1=t(scale(t(im_symbol2)))
n1[n1>1]=1
n1[n1< -1]= -1
ht1 <- pheatmap(n1,show_rownames = F,show_colnames = F, 
                cluster_rows = T,cluster_cols = F,
                annotation_col = colData1_1,color=color,
                main = "Immune Score")
a<-intersect(row.names(im_symbol2),row.names(DE_im_up))
DE_im_up<-DE_im[which(DE_im$change=="UP"),]
DE_im_up<-DE_im_up[order(DE_im_up$baseMean,decreasing = T),]
DE_im_down<-DE_im[which(DE_im$change=="DOWN"),]
DE_im_down<-DE_im_down[order(DE_im_down$pvalue),]
top_im<-im_symbol1[c(rownames(DE_im_up)[1:20],
                     row.names(DE_im_down)[1:20]),]
n1=t(scale(t(top_im)))
n1[n1>1]=1
n1[n1< -1]= -1
ht1 <- pheatmap(n1,show_rownames = T,show_colnames = F,
                cluster_rows = F,cluster_cols = F,
                annotation_col = colData1_1,color=color,
                main = "Immune Score")

st_symbol<-data_tumor[c(st_up,st_down),]
st_symbol1<-cbind(st_symbol[,row.names(colData2)[colData2$condition=="high"]],
                  st_symbol[,row.names(colData2)[colData2$condition=="low"]])

st_symbol2<-st_symbol1[which(apply(st_symbol1,1,
                                   function(x){return(sum(x>10))})>ncol(st_symbol1)*0.25),]
library(stringr)
score_tumor$Stromal<- ifelse(score_tumor$StromalScore>median(score_tumor$StromalScore),"high","low")
table(score_tumor$Stromal)
colData2<-data.frame(row.names = row.names(score_tumor),
                     condition=score_tumor$Stromal)
colData2_1<-data.frame(row.names = c(row.names(colData2)[colData2$condition=="high"],
                                     row.names(colData2)[colData2$condition=="low"]),
                       condition=factor(c(rep(c("high"),554),rep(c("low"),555))))

library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#436eee','white','#EE0000'))(100)
n2=t(scale(t(st_symbol2)))
n2[n2>1]=1
n2[n2< -1]= -1
ht2 <- pheatmap(n2,show_rownames = F,show_colnames = F, 
                cluster_rows = T,cluster_cols = F,
                annotation_col = colData2_1,color=color,
                main = "Stromal Score")

DE_st_up<-DE_st[which(DE_st$change=="UP"),]
DE_st_up<-DE_st_up[order(DE_st_up$baseMean,decreasing = T),]

DE_st_down<-DE_st[which(DE_st$change=="DOWN"),]
DE_st_down<-DE_st_down[order(DE_st_down$pvalue),]
top_st<-st_symbol1[c(rownames(DE_st_up)[1:20],
                     row.names(DE_st_down)[1:20]),]
n2=t(scale(t(top_st)))
n2[n2>1]=1
n2[n2< -1]= -1
ht2 <- pheatmap(n2,show_rownames = T,show_colnames = F,
                cluster_rows = F,cluster_cols = F,
                annotation_col = colData2_1,color=color,
                main = "Stromal Score")


#Tumor samples were grouped based on the median immune score and stromal score to explore their impact on prognosis
score_surv<-merge(clinical2[,c(1,3,4)],score_class[,1:4],
                  by.x = "patient_id",by.y = "Group.1")
library(survival)
library(survminer)
score_surv$stromal<-ifelse(score_surv$StromalScore>median(score_surv$StromalScore),"high","low")
fit_im<-survfit(Surv(time,status)~stromal,data = score_surv)

p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              ggtheme = theme_bw(), 
              palette = c("#E7B800", "#2E9FDF"),
              title="Stromal Survival")


score_surv$immune<-ifelse(score_surv$ImmuneScore>median(score_surv$ImmuneScore),"high","low")
fit_im<-survfit(Surv(time,status)~immune,data = score_surv)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              ggtheme = theme_bw(), 
              palette = c("#E7B800", "#2E9FDF"),
              title="Immune Survival")



score_surv$estimate<-ifelse(score_surv$ESTIMATEScore>median(score_surv$ESTIMATEScore),"high","low")
fit_im<-survfit(Surv(time,status)~estimate,data = score_surv)

p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              ggtheme = theme_bw(), 
              palette = c("#E7B800", "#2E9FDF"),
              title="Estimate Survival")

#Distribution of immune scores, stromal scores and estimated scores among different clinicopathological features
dir.create("estatime_TNM")
score_class<-merge(score_tumor,TNMstage,by.x = "sample",by.y = "patient_id")
score_class<-aggregate(score_class[,2:4],by=list(score_class$sample),FUN = mean)
score_class<-merge(score_class,TNMstage,by.x = "Group.1",by.y = "patient_id")
score_class<-score_class[score_class$M!="MX",]
score_class<-score_class[score_class$N!="NX",]
score_class<-score_class[score_class$Stage!="Stage X",]
score_class<-score_class[score_class$Stage!="'--",]

library(ggplot2)
library(ggsci)
library(ggpubr)

immune_M<-score_class[,c(3,5)]
my_comparisons = list(c("M0","M1"))
p<-ggboxplot(immune_M,x="M",y="ImmuneScore",
          title = "",ylab = "Immune score",
          xlab = "",
          color = "M",palette = "jco")+
  stat_compare_means(label.y = 3000)+
  guides(colour="none")
ggsave("immune_M.png",p,width = 3,height = 3)
ggsave("immune_M.pdf",p,width = 3,height = 3)

immune_T<-score_class[,c(3,8)]
immune_T<-immune_T[order(immune_T$T),]
my_comparisons_T = list(c("T1","T2"),c("T2","T3"),c("T3","T4"))
p<-ggboxplot(immune_T,x="T",y="ImmuneScore",
          title = "",ylab = "Immune score",xlab = "",
          color = "T",palette = "jco")+
  stat_compare_means(label.y = 3000,label.x = 1.2)+
  guides(colour="none")
ggsave("immune_T.png",p,width = 3,height = 3)
ggsave("immune_T.pdf",p,width = 3,height = 3)

immune_N<-score_class[,c(3,6)]
immune_N<-immune_N[order(immune_N$N),]
my_comparisons_N= list(c("N0","N1"),c("N1","N2"),c("N2","N3"))
p<-ggboxplot(immune_N,x="N",y="ImmuneScore",
          title = "",ylab = "Immune score",xlab = "",
          color = "N",palette = "jco")+
  stat_compare_means(label.y = 3000,label.x = 1.2)+
  guides(colour="none")
ggsave("immune_N.png",p,width = 3,height = 3)
ggsave("immune_N.pdf",p,width = 3,height = 3)

immune_stage<-score_class[,c(3,7)]
immune_stage<-immune_stage[order(immune_stage$Stage),]
p<-ggboxplot(immune_stage,x="Stage",y="ImmuneScore",
          title = "",ylab = "Immune score",xlab = "",
          color = "Stage",palette = "jco")+
  stat_compare_means(label.y = -1500,label.x = 1.2)+
  guides(colour="none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
ggsave("immune_stage.png",p,width = 3,height = 3)
ggsave("immune_stage.pdf",p,width = 3,height = 3)


stromal_M<-score_class[,c(2,5)]
my_comparisons = list(c("M0","M1"))
p<-ggboxplot(stromal_M,x="M",y="StromalScore",
          title = "",ylab = "Stromal score",xlab = "",
          color = "M",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 2000)
ggsave("stromal_M.png",p,width = 3,height = 3)
ggsave("stromal_M.pdf",p,width = 3,height = 3)
#

stromal_T<-score_class[,c(2,8)]
stromal_T<-stromal_T[order(stromal_T$T),]
my_comparisons_T = list(c("T1","T2"),c("T2","T3"),c("T3","T4"))
p<-ggboxplot(stromal_T,x="T",y="StromalScore",
          title = "",ylab = "Stromal score",xlab = "",
          color = "T",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 2000,label.x = 1.2)
ggsave("stromal_T.png",p,width = 3,height = 3)
ggsave("stromal_T.pdf",p,width = 3,height = 3)
#
stromal_N<-score_class[,c(2,6)]
stromal_N<-stromal_N[order(stromal_N$N),]
my_comparisons_N= list(c("N0","N1"),c("N1","N2"),c("N2","N3"))
p<-ggboxplot(stromal_N,x="N",y="StromalScore",
          title = "",ylab = "Stromal score",xlab = "",
          color = "N",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 2500,label.x = 1.2)
ggsave("stromal_N.png",p,width = 3,height = 3)
ggsave("stromal_N.pdf",p,width = 3,height = 3)

stromal_stage<-score_class[,c(2,7)]
stromal_stage<-stromal_stage[order(stromal_stage$Stage),]
my_comparisons_stage= list(c("Stage I","Stage II"),c("Stage II","Stage III"),c("Stage III","Stage IV"))
p<-ggboxplot(stromal_stage,x="Stage",y="StromalScore",
          title = "",ylab = "Stromal score",xlab = "",
          color = "Stage",palette = "jco")+
  guides(colour="none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))+
  stat_compare_means(label.y = -2000,label.x = 1.2)
ggsave("stromal_stage.png",p,width = 3,height = 3)
ggsave("stromal_stage.pdf",p,width = 3,height = 3)


estimate_M<-score_class[,c(4,5)]
my_comparisons = list(c("M0","M1"))
p<-ggboxplot(estimate_M,x="M",y="ESTIMATEScore",
          title = "",ylab = "Estimate score",xlab = "",
          color = "M",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 4000)
ggsave("estimate_M.png",p,width = 3,height = 3)
ggsave("estimate_M.pdf",p,width = 3,height = 3)

estimate_T<-score_class[,c(4,8)]
estimate_T<-estimate_T[order(estimate_T$T),]
my_comparisons_T = list(c("T1","T2"),c("T2","T3"),c("T3","T4"))
p<-ggboxplot(estimate_T,x="T",y="ESTIMATEScore",
          title = "",ylab = "Estimate score",xlab = "",
          color = "T",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 4000,label.x = 1.2)
ggsave("estimate_T.png",p,width = 3,height = 3)
ggsave("estimate_T.pdf",p,width = 3,height = 3)


estimate_N<-score_class[,c(4,6)]
estimate_N<-estimate_N[order(estimate_N$N),]
my_comparisons_N= list(c("N0","N1"),c("N1","N2"),c("N2","N3"))
p<-ggboxplot(estimate_N,x="N",y="ESTIMATEScore",
          title = "",ylab = "Estimate score",xlab = "",
          color = "N",palette = "jco")+
  guides(colour="none")+
  stat_compare_means(label.y = 4000,label.x = 1.2)
ggsave("estimate_N.png",p,width = 3,height = 3)
ggsave("estimate_N.pdf",p,width = 3,height = 3)


estimate_stage<-score_class[,c(4,7)]
estimate_stage<-estimate_stage[order(estimate_stage$Stage),]
my_comparisons_stage= list(c("Stage I","Stage II"),c("Stage II","Stage III"),c("Stage III","Stage IV"))
p<-ggboxplot(estimate_stage,x="Stage",y="ESTIMATEScore",
          title = "",ylab = "Estimate score",xlab = "",
          color = "Stage",palette = "jco")+
  guides(colour="none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))+
  stat_compare_means(label.y = -2000,label.x = 1.2)
ggsave("estimate_stage.png",p,width = 3,height = 3)
ggsave("estimate_stage.pdf",p,width = 3,height = 3)

save(immune_M,immune_stage,immune_N,immune_T,stromal_M,
     stromal_N,stromal_stage,stromal_T,estimate_M,
     estimate_N,estimate_stage,estimate_T,file = "临床分期得分.Rdata")



#Organize clinical data
dir.create('survival')
library(readr)
clinical <- read_tsv("clinical/clinical.tsv") 
clinical1<-clinical[,c("case_submitter_id","age_at_index","days_to_death","days_to_last_follow_up","gender","vital_status","race","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage","ajcc_pathologic_t")]
clinical1<-clinical1[!duplicated(clinical1$case_submitter_id),]
clinical2<-clinical1[,c(1,3,4,6)]
clinical2$days_to_death<-ifelse(clinical2$days_to_death=="'--",'0',clinical2$days_to_death)
clinical2$days_to_last_follow_up<-ifelse(clinical2$days_to_last_follow_up=="'--",'0',clinical2$days_to_last_follow_up)

clinical2[,2:3]<-lapply(clinical2[,2:3],as.numeric)
clinical2$time<-clinical2$days_to_death+clinical2$days_to_last_follow_up
clinical2<-clinical2[-which(clinical2$time<=0),]
clinical2<-clinical2[,-(2:3)]
clinical2$status<-ifelse(clinical2$vital_status=='Alive',0,1)
colnames(clinical2)<-c("patient_id","vital_status","time","status")

library(stringr)
surv_symbol<-surv_tumor
surv_symbol<-log2(surv_symbol+1)

#The differentially expressed genes were analyzed by univariate regression analysis, multivariate regression analysis and LASSO regression analysis
dir.create('LASSO')
genediff<-genediff[,which(apply(genediff,2,
                                function(x){return(sum(x>0))})>nrow(genediff)*0.25)]  
genediff$id<-row.names(genediff)
surgene<-merge(clinical2,genediff,by.x = "patient_id",by.y = "id")
library(survival)
result<-matrix(,648,7)
for (i in 5:652){
  group=ifelse(surgene[,i]>median(surgene[,i]),'high','low') 
  c <- coxph(Surv(time,status)~group,surgene)
  result[i-4,1] <- colnames(surgene)[i]
  result[i-4,2] <- summary(c)[[7]][[1]][1]
  result[i-4,3] <- summary(c)[[7]][[2]][1]
  result[i-4,4] <- summary(c)[[8]][3]
  result[i-4,5] <- summary(c)[[8]][4]
  result[i-4,6] <- summary(c)[[7]][[5]][1]
  result[i-4,7] <- summary(c)[[10]][[3]][1]
}
result<-data.frame(result)
colnames(result) <- c("Gene","coef","HR","lower.95","upper.95","cox.p","km.p")
result[,2:7]<-lapply(result[,2:7],as.numeric)
diff_res<-result[which(result$cox.p<0.05),]

genecox<-surgene[,diff_res$Gene]
genecox<-cbind(surgene[,1:4],genecox)
for(i in 5:172){
  genecox[,i]<-ifelse(genecox[,i]>median(genecox[,i]),1,0)
}
mult_models <- coxph(Surv(time, status)~., data = genecox[,3:172])

a<-summary(mult_models)
result_mult<-data.frame(coef=round(as.matrix(a$coefficients[,1]),3),
                        HR=round(as.matrix(a$coefficients[,2]),3),
                        lower.95=round(as.matrix(a$conf.int[,3]),3),
                        upper.95=round(as.matrix(a$conf.int[,4]),3),
                        p.value=as.matrix(a$coefficients[,5]))
library(stringr)
row.names(result_mult)<-ifelse(str_detect(row.names(result_mult),'-'),str_split(row.names(result_mult),'`',simplify = T)[,2],row.names(result_mult))

diff_mult<-result_mult[which(result_mult$p.value<0.05),]

library(glmnet)
y<-surgene[,3:4]
x<-surgene[,diff_res$Gene]
x1<-as.matrix(x)
rownames(y)<-surgene$patient_id
rownames(x1)<-surgene$patient_id
y$time<-as.double(y$time)
y$status<-as.double(y$status)
y1<-Surv(y$time,y$status)
fitcv<-cv.glmnet(x1,y1,family='cox',alpha = 1,nfolds = 5,type.measure = "deviance")
cofficient<-coef(fitcv,s=fitcv$lambda.min)
active.coff<-which(as.numeric(cofficient)!=0)
active.gene<-as.numeric(cofficient)[active.coff]
sig_genecox<-row.names(cofficient)[active.coff]
sig_coef<-data.frame(gene=sig_genecox,coef=active.gene)
#TMERS of each sample were calculated based on LASSO regression analysis results
dir.create("LASSO/10gene")
gene10<-surv_tumor[sig_coef$gene,]
TMEscore<-data.frame(sample=colnames(gene10),score=NA)
for(i in 1:1069){
  c=0
  for(j in 1:6){
    c=c+(gene10[j,i]*sig_coef[j,2])
  }
  TMEscore[i,2]<-c
}
score_surv<-merge(clinical2,TMEscore,by.x = "patient_id",by.y = "sample")
#Based on LASSO regression analysis results, KM survival curve, distribution waterfall plot and survival time distribution scatter plot of TMERS were drawn in TCGA, metabric, GSE21653 and GSE58812, respectively
{
output = "GSE58812"
load(paste0(output,"/TIDE.Rdata"))
load(paste0(output,"expression.Rdata"))
load(paste0(output,"clincal.Rdata"))
mat<-metabric_all
sig_coef<-sig_coef

gene<-mat[sig_coef$gene,]
gene<-na.omit(gene)
sig_coef1<-sig_coef[sig_coef$gene%in%row.names(gene),]

TMEscore<-data.frame(sample=colnames(gene),score=NA)
for(i in 1:ncol(gene)){
  c=sum(gene[,i]*sig_coef1[,2])
  TMEscore[i,2]<-c
}

score_surv<-merge(clinical,TMEscore,by.x = "sample",by.y = "sample")

library(survival)
library(ggplot2)
library(survminer)
score_surv$group<-ifelse(score_surv$score>median(score_surv$score),"High","Low")
fit_im<-survfit(Surv(time,status)~group,data = score_surv)
jpeg(paste0(output,"-TMERS_OS.jpg"),width = 500,height = 600)
p<-ggsurvplot(fit_im,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              ggtheme = theme_bw(), 
              palette = c("#E7B800", "#2E9FDF"),
              title= paste("Overall survival of ",output))
print(p)
dev.off()

tide_score<-merge(score_surv[,c(1,4)],tumor_TIDE,by.x = "sample",by.y = "Sample")
tide_score<-tide_score[,c(-3,-4)]
library(ggplot2)
library(ggpubr)
tide_score$group<-ifelse(tide_score$score>median(tide_score$score),'High Risk','Low Risk')

p<-ggboxplot(tide_score,x="group",y="scale",
             title = "TMERS",ylab = "TIDE prediction score",
             color = "group",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(tide_score$scale)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1,label.x = 1.2)

ggsave(paste0(output,"-TMERS_TIDE.png"),p,width = 4,height = 6)
ggsave(paste0(output,"-TMERS_TIDE.pdf"),p,width = 4,height = 6)

rm(list = ls())

scoreplot<-score_surv
scoreplot$group<-ifelse(scoreplot$status==0,"Alive","Dead")
scoreplot<-scoreplot[order(scoreplot$score),]
scoreplot$num<-1:nrow(scoreplot)
library(ggplot2)
p<-ggplot(scoreplot,aes(x = num,y = score,fill = group))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(legend.position = c(.8,.2))+
  scale_fill_manual(breaks = c("Alive","Dead"),
                    values = c("#EFC000","#0073C2"))+
  labs(x="",y="TMERS")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5))


p<-ggplot()+
  geom_point(data=scoreplot,aes(x=num, y=time,colour=group),size=1)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x="",y="Survival time(month)")+
  scale_color_manual(breaks = c("Alive","Dead"),
                     values = c("#EFC000","#0073C2"))

}

#Distribution of TMERS between clinicopathological features and breast cancer subtypes
TNMplot<-merge(TNMstage,score_surv,by.x = "patient_id",by.y = "sample")
library(ggplot2)
library(ggpubr)
TNMRS_M<-TNMplot[,c(2,8)]
TNMRS_M<-TNMRS_M[TNMRS_M$M!="MX",]

p<-ggboxplot(TNMRS_M,x="M",y="score",
             title = "M classfication",ylab = "TMERS",
             color = "M",palette = "jama",xlab = "",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(TNMRS_M$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(method = "wilcox.test",label.y = -1,label.x = 1.2)

ggsave("TNM/TNMRS_M.png",p,width = 4,height = 4)
ggsave("TNM/TNMRS_M.pdf",p,width = 4,height = 4)

my_comparisons = list(c("M0","M1"))

p<-ggboxplot(TNMRS_M,x="M",y="score",
             title = "M classfication",ylab = "TMERS",
             color = "M",palette = "jama",xlab = "",
             add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
ggsave("TNM/TNMRS_M_2.png",p,width = 5,height = 5)

TNMRS_T<-TNMplot[,c(5,8)]
TNMRS_T<-TNMRS_T[TNMRS_T$T!="TX",]
TNMRS_T<-TNMRS_T[order(TNMRS_T$T),]

p<-ggboxplot(TNMRS_T,x="T",y="score",
             title = "T classfication",ylab = "TMERS",
             color = "T",palette = "jama",xlab = "",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(TNMRS_T$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1.2)
ggsave("TNM/TNMRS_T.png",p,width = 4,height = 4)
ggsave("TNM/TNMRS_T.pdf",p,width = 4,height = 4)

my_comparisons_T = list(c("T1","T2"),c("T2","T3"),c("T3","T4"),c("T1","T4"))
p<-ggboxplot(TNMRS_T,x="T",y="score",
             title = "T classfication",ylab = "TMERS",
             color = "T",palette = "jama",xlab = "",
             add = "jitter")+
  stat_compare_means(comparisons = my_comparisons_T)

ggsave("TNM/TNMRS_T_2.png",p,width = 5,height = 5)
TNMRS_N<-TNMplot[,c(3,8)]
TNMRS_N<-TNMRS_N[TNMRS_N$N!="NX",]
TNMRS_N<-TNMRS_N[order(TNMRS_N$N),]
p<-ggboxplot(TNMRS_N,x="N",y="score",
             title = "N classfication",ylab = "TMERS",
             color = "N",palette = "jama",xlab = "",
             add = "jitter")+
  coord_cartesian(ylim = (boxplot.stats(TNMRS_N$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1.2)
ggsave("TNM/TNMRS_N.png",p,width = 4,height = 4)
ggsave("TNM/TNMRS_N.pdf",p,width = 4,height = 4)

my_comparisons_N= list(c("N0","N1"),c("N1","N2"),c("N2","N3"),c("N0","N3"))
p<-ggboxplot(TNMRS_N,x="N",y="score",
             title = "N classfication",ylab = "TMERS",
             color = "N",palette = "jama",xlab = "",
             add = "jitter")+
  stat_compare_means(comparisons = my_comparisons_N)

ggsave("TNM/TNMRS_N_2.png",p,width = 5,height = 5)

TNMRS_Stage<-TNMplot[,c(4,8)]
TNMRS_Stage<-TNMRS_Stage[TNMRS_Stage$Stage!="Stage X",]
TNMRS_Stage<-TNMRS_Stage[TNMRS_Stage$Stage!="'--",]
TNMRS_Stage<-TNMRS_Stage[order(TNMRS_Stage$Stage),]
p<-ggboxplot(TNMRS_Stage,x="Stage",y="score",
             title = "Stage",ylab = "TMERS",xlab = "",
             color = "Stage",palette = "jama",
             add = "jitter")+
  coord_cartesian(ylim = (boxplot.stats(TNMRS_Stage$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1.2)
ggsave("TNM/TNMRS_Stage.png",p,width = 5,height = 5)
ggsave("TNM/TNMRS_Stage.pdf",p,width = 5,height = 5)

my_comparisons_stage= list(c("Stage I","Stage II"),c("Stage II","Stage III"),c("Stage III","Stage IV"),c("Stage I","Stage IV"))
p<-ggboxplot(TNMRS_Stage,x="Stage",y="score",
             title = "Stage",ylab = "TMERS",xlab = "",
             color = "Stage",palette = "jama",
             add = "jitter")+
  stat_compare_means(comparisons = my_comparisons_stage)
ggsave("TNM/TNMRS_Stage_2.png",p,width = 5,height = 5)

subtype<-read.csv("TNM/BCRA-subtype.csv")
typescore<-merge(score_surv,subtype[,c(1,4,12)],by.x = "sample",by.y = "Sample.ID")

typescore<-typescore[typescore$BRCA_Subtype_PAM50!="Normal",]
colnames(typescore)[6]<-"Subtype"

p<-ggboxplot(typescore,x="Subtype",y="score",
             title = "BRCA Subtype",ylab = "TMERS",xlab = "",
             color = "Subtype",palette = "jama",
             add = "jitter")+
  coord_cartesian(ylim = (boxplot.stats(TNMRS_Stage$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1.2)
ggsave("TNM/TNMRS_subtype.png",p,width = 5,height = 5)
ggsave("TNM/TNMRS_subtype.pdf",p,width = 5,height = 5)

save(TMEscore,TNMRS_M,TNMRS_N,TNMRS_T,TNMstage,TNMplot,
     TNMRS_Stage,typescore,subtype,
     file = "TNM/TNM分期绘图.Rdata")

TNMstage<-TNMstage[TNMstage$M!="MX",]
TNMstage<-TNMstage[TNMstage$N!="NX",]
TNMstage<-TNMstage[TNMstage$Stage!="Stage X",]
TNMstage<-TNMstage[TNMstage$Stage!="'--",]

class_map<-merge(TNMstage,typescore[,c(1,5,6)],by.x = "patient_id",by.y = "sample")

pic_map<-class_map

pic_map$M<-ifelse(pic_map$M=="M0",1,2)
pic_map$N<-ifelse(pic_map$N=="N0",1,
                  ifelse(pic_map$N=="N1",2,
                         ifelse(pic_map$N=="N2",3,4)))

pic_map$T<-ifelse(pic_map$T=="T1",1,
                  ifelse(pic_map$T=="T2",2,
                         ifelse(pic_map$T=="T3",3,4)))

pic_map$Stage<-ifelse(pic_map$Stage=="Stage I",1,
                  ifelse(pic_map$Stage=="Stage II",2,
                         ifelse(pic_map$Stage=="Stage III",3,4)))

pic_map$Subtype<-ifelse(pic_map$Subtype=="LumA",1,
                        ifelse(pic_map$Subtype=="Her2",2,
                               ifelse(pic_map$Subtype=="LumB",3,4)))


anno<-merge(pic_map,TMEscore,by.x="patient_id",by.y = "sample")
anno$group<-ifelse(anno$score>median(anno$score),"High","Low")
anno<-anno[order(anno$group),]
annotation<-data.frame(row.names = anno$patient_id,group=anno$group)

rownames(pic_map)<-pic_map[,1]
pic_map<-data.frame(t(pic_map))
colnames(pic_map)<-pic_map[1,]
pic_map<-pic_map[-1,]
pic_map[,1:850]<-lapply(pic_map[,1:850],as.numeric)
pic_map<-pic_map[,row.names(annotation)]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#F2CD33','#338FCE','#E64B35','#7E6148'))(100)
pheatmap(pic_map,show_rownames = T,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         annotation_col = annotation,color=color,
         legend_breaks = c(1:4), 
         legend_labels = c("1","2","3","4"),
         main = "",filename = "TNM/map.pdf",height = 3,width = 7)

#To confirmed whether the prognostic signature are independent prognostic factors for BRCA, univariate and multivariate Cox regression analyses were performed on BRCA patients.
dir.create("prognostic ")
prognostic<-merge(clinical2,TNMstage,by.x = "patient_id",by.y = "patient_id")

prognostic<-merge(prognostic,clinical1[,c(1,2,5)],by.x = "patient_id",by.y = "case_submitter_id")
colnames(prognostic)[9]<-"age"

prognostic<-merge(prognostic,score_surv[,c(1,4)],by.x = "patient_id",by.y = "sample")

colnames(prognostic)[11]<-"RiskScore"
colnames(prognostic)[5:11]<-c("Pathologic_M_stage","Pathologic_N_stage",
  "Pathologic_tumor_stage","Pathologic_T_stage","Age",
  "Gender","RiskScore")
save(prognostic,file = "prognostic/prognostic_data.Rdata")

library(survival)
prog_cox<-prognostic
prog_cox$RiskScore<-ifelse(prog_cox$RiskScore>median(prog_cox$RiskScore),2,1)
prog_cox$Pathologic_M_stage<-ifelse(prog_cox$Pathologic_M_stage=="M0",1,ifelse(prog_cox$Pathologic_M_stage=="M1",2,0))
prog_cox$Pathologic_N_stage<-ifelse(prog_cox$Pathologic_N_stage=="N0",1,
                   ifelse(prog_cox$Pathologic_N_stage=="N1",2,
                          ifelse(prog_cox$Pathologic_N_stage=="N2",3,
                                 ifelse(prog_cox$Pathologic_N_stage=="N3",4,0))))
prog_cox$Pathologic_T_stage<-ifelse(prog_cox$Pathologic_T_stage=="T1",1,
                   ifelse(prog_cox$Pathologic_T_stage=="T2",2,
                          ifelse(prog_cox$Pathologic_T_stage=="T3",3,
                                 ifelse(prog_cox$Pathologic_T_stage=="T4",4,0))))
prog_cox$Pathologic_tumor_stage<-ifelse(prog_cox$Pathologic_tumor_stage=="Stage I",1,
                       ifelse(prog_cox$Pathologic_tumor_stage=="Stage II",2,
                              ifelse(prog_cox$Pathologic_tumor_stage=="Stage III",3,
                                     ifelse(prog_cox$Pathologic_tumor_stage=="Stage IV",4,0))))

prog_cox$Gender<-ifelse(prog_cox$Gender=="female",2,1)


covariates<-colnames(prog_cox)[5:11]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = prog_cox)})
result<-matrix(0,7,7)

for(i in 1:7){
  x<-summary(univ_models[[i]])
  result[i,1] <- names(univ_models)[i]
  result[i,2] <- round(x[[7]][1],3)
  result[i,3] <- round(x[[7]][2],3)
  result[i,4] <- round(x[[8]][3],3)
  result[i,5] <- round(x[[8]][4],3)
  result[i,6] <- x[[7]][5]
  result[i,7]<-paste0(round(x[[7]][2],3),"(",round(x[[8]][3],3),"-",round(x[[8]][4],3),")")
}
result<-data.frame(result)
colnames(result) <- c("Variable","coef","HR","lower.95","upper.95","p.value","HR(95% CI)")
result[2:6]<-lapply(result[,2:6],as.numeric)

library(forestplot)
result1<-result
result1$p.value<-ifelse(result1$p.value<0.001,"<0.001",round(result1$p.value,3))
result1<-rbind(c("Univariate analysis",NA,NA,NA,NA,"p.value","HR (95% CI)"),result1)
result1[,2:5]<-lapply(result1[,2:5],as.numeric)
#Map the forest
forestplot(labeltext=as.matrix(result1[,c(1,6,7)]),
           mean=result1$HR,
           lower=result1$lower.95,
           upper=result1$upper.95,
           zero=1,
           boxsize=0.2,
           graph.pos=2,
           lwd.zero=1.5,
           lwd.ci=2,
           lineheight = unit(10,'mm'),
           colgap=unit(2,'mm'),
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="<---Favors Beneficial  Favors Harmful--->",
           lwd.xaxis =1,
           lty.ci = "solid",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           line.margin = 0.08,xticks = c(0,1,2,4,6))


prog_cox<-prognostic
prog_cox$RiskScore<-ifelse(prog_cox$RiskScore>median(prog_cox$RiskScore),2,1)
prog_cox<-prog_cox[-which(prog_cox$Pathologic_tumor_stage=="Stage X"),]
prog_cox<-prog_cox[-which(prog_cox$Pathologic_tumor_stage=="'--"),]
prog_cox<-prog_cox[-which(prog_cox$Pathologic_M_stage=="MX"),]
prog_cox<-prog_cox[-which(prog_cox$Pathologic_N_stage=="NX"),]
prog_cox$Pathologic_M_stage<-ifelse(prog_cox$Pathologic_M_stage=="M0",0,1)
prog_cox$Pathologic_N_stage<-ifelse(prog_cox$Pathologic_N_stage=="N0",0,
                   ifelse(prog_cox$Pathologic_N_stage=="N1",1,
                          ifelse(prog_cox$Pathologic_N_stage=="N2",2,3)))
prog_cox$Pathologic_T_stage<-ifelse(prog_cox$Pathologic_T_stage=="T1",1,
                   ifelse(prog_cox$Pathologic_T_stage=="T2",2,
                          ifelse(prog_cox$Pathologic_T_stage=="T3",3,4)))
prog_cox$Pathologic_tumor_stage<-ifelse(prog_cox$Pathologic_tumor_stage=="Stage I",1,
                       ifelse(prog_cox$Pathologic_tumor_stage=="Stage II",2,
                              ifelse(prog_cox$Pathologic_tumor_stage=="Stage III",3,4)))

prog_cox$Gender<-ifelse(prog_cox$Gender=="female",2,1)

mult_model<-coxph(Surv(time,status)~Pathologic_M_stage+Pathologic_tumor_stage
                  +Pathologic_T_stage+Pathologic_N_stage+Age+
                    Gender+RiskScore,data=prog_cox)
a<-summary(mult_model)

result_mult<-data.frame(coef=round(as.matrix(a$coefficients[,1]),3),
                        HR=round(as.matrix(a$coefficients[,2]),3),
                        lower.95=round(as.matrix(a$conf.int[,3]),3),
                        upper.95=round(as.matrix(a$conf.int[,4]),3),
                        p.value=as.matrix(a$coefficients[,5]))

result_mult$variable<-row.names(result_mult)

result_mult$CI<-paste0(result_mult$HR,"(",result_mult$lower.95,"-",result_mult$upper.95,")")
result_mult1<-result_mult

result_mult1$p.value<-ifelse(result_mult1$p.value<0.001,"<0.001",round(result_mult1$p.value,3))
result_mult1<-rbind(c(NA,NA,NA,NA,"p.value","Multivariate analysis","HR (95% CI)"),result_mult1)
result_mult1[,1:4]<-lapply(result_mult1[,1:4],as.numeric) 


forestplot(labeltext=as.matrix(result_mult1[,c(6,5,7)]),
           mean=result_mult1$HR,
           lower=result_mult1$lower.95,
           upper=result_mult1$upper.95,
           zero=1,
           boxsize=0.2,
           graph.pos=2,
           lwd.zero=1.5,
           lwd.ci=2,
           lineheight = unit(10,'mm'),
           colgap=unit(2,'mm'),
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="<---Favors Beneficial  Favors Harmful--->",
           lwd.xaxis =1,
           lty.ci = "solid",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           line.margin = 0.08,xticks = c(0,1,2,3,4))

#combined with prognostic characteristics, age and tumor pathological stage, a comprehensive Nomogram information map was further constructed
library(rms)
library(regplot)
library(survival)
library(survminer)

nom<-prognostic[,c(3,4,7,9,11)]

nom<-nom[-which(nom$Pathologic_tumor_stage=="Stage X"),]
nom<-nom[-which(nom$Pathologic_tumor_stage=="'--"),]

nom$Pathologic_tumor_stage<-ifelse(nom$Pathologic_tumor_stage=="Stage I",1,
              ifelse(nom$Pathologic_tumor_stage=="Stage II",2,
                     ifelse(nom$Pathologic_tumor_stage=="Stage III",3,4)))
res.cox<-coxph(Surv(time,status)~Pathologic_tumor_stage+Age+
                 RiskScore,data = nom)
nom1<-regplot(res.cox, clickable=TRUE, 
              points=TRUE, rank="sd",failtime = c(3650,7300),
              prfail = T,title = "coxph regression")

nom2<-regplot(res.cox,observation=nom[3,], clickable=TRUE, 
              points=TRUE, rank="sd",failtime = c(3650,7300),droplines=T,prfail = T,
              other=(list(bvcol="red",sq="green",obscol="blue")))

#Decision curve analyses (DCA) curves evaluate clinical predictive potentialof nomogram
res.cox2 <- cph(Surv(time,status) ~ Pathologic_tumor_stage+Age+
                  RiskScore, data =nom, x=T, y=T, surv=TRUE, time.inc = 3650)

cal1<-calibrate(res.cox2,cmethod = "KM",method = "boot",u=3650,B=1047,m=250)
par(mar=c(8,5,3,2),cex = 1.0)
plot(cal1,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted OS(%)",
     ylab="Observed OS(%)",
     col=c(rgb(192,98,83,maxColorValue=255)))
title(main="10-year calibration curve")
dev.off()

res.cox3 <- cph(Surv(time,status) ~ Pathologic_tumor_stage+Age+
                  RiskScore, data =nom, x=T, y=T, surv=TRUE, time.inc = 7300)

cal2<-calibrate(res.cox3,cmethod = "KM",method = "boot",u=7300,B=1047,m=250)

par(mar=c(8,5,3,2),cex = 1.0)

plot(cal2,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted OS(%)",
     ylab="Observed OS(%)",
     col=c(rgb(192,98,83,maxColorValue=255)))
title(main="20-year calibration curve")
dev.off()

library(rmda)
nom1<-nom[which(nom$time<3650),]

simple1<- decision_curve(status~Age,data= nom1,
                        family = binomial(link ='logit'),
                        thresholds= seq(0,.4, by = 0.01),
                        confidence.intervals = 0.95,
                        bootstraps = 10)

simple2<- decision_curve(status~RiskScore,data= nom1,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,.4, by = 0.01),
                         confidence.intervals = 0.95,
                         bootstraps = 10)

simple3<- decision_curve(status~Pathologic_tumor_stage,data= nom1,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,.4, by = 0.01),
                         confidence.intervals = 0.95,
                         bootstraps = 10)
List<-list(simple1,simple2,simple3)

plot_decision_curve(List,
                    curve.names=c('Age','RiskScore','Stage'),
                    cost.benefit.axis =FALSE,
                    col= c('red','blue','green',"black",'black'),
                    confidence.intervals=FALSE,
                    lwd = c(4,4,4,2,2),
                    standardize = FALSE,
                    xlab = "Threshold probability")+title("10-years DCA curve")

nom1<-nom[which(nom$time<7300),]
simple1<- decision_curve(status~Age,data= nom1,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,.4, by = 0.01),
                         confidence.intervals = 0.95,
                         bootstraps = 10)
simple2<- decision_curve(status~RiskScore,data= nom1,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,.4, by = 0.01),
                         confidence.intervals = 0.95,
                         bootstraps = 10)

simple3<- decision_curve(status~Pathologic_tumor_stage,data= nom1,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,.4, by = 0.01),
                         confidence.intervals = 0.95,
                         bootstraps = 10)
List<-list(simple1,simple2,simple3)

plot_decision_curve(List,
                    curve.names=c('Age','RiskScore','Stage'),
                    cost.benefit.axis =FALSE,
                    col= c('red','blue','green',"black",'black'),
                    confidence.intervals=FALSE,
                    lwd = c(4,4,4,2,2),
                    standardize = FALSE,
                    xlab = "Threshold probability")+title("20-years DCA curve")


##The prognostic signature correlated with immune cell infiltration in BRCA
dir.create("immune")
library(xCell)
xCell  = xCellAnalysis(as.matrix(data_symbol), rnaseq = TRUE)
xCell<-as.data.frame(t(xCell))
write.table(xCell,file="immune/xCell_fpkm_Charoentong_0600110421_RAW.txt",sep="\t",quote = "")
xcell<-read.table("immune/xCell_fpkm_Charoentong_0600110421_RAW.txt",header = T,sep = "\t",check.names = F,quote = NULL)
rownames(xcell)<-xcell[,1]
xcell<-xcell[,-1]
library(stringr)
tumor_xcell<-xcell
colnames(tumor_xcell)<-gsub("\\.","-",colnames(tumor_xcell))

TMEscore$group<-ifelse(TMEscore$score>median(TMEscore$score),"High Risk","Low Risk")
TMEscore<-TMEscore[order(TMEscore$group),]

annodata<-data.frame(row.names = TMEscore$sample,
                     condition=TMEscore$group)

pic_data<-tumor_xcell[,rownames(annodata)]
library(pheatmap)
library(RColorBrewer)
n=t(scale(t(pic_data)))
n[n>1]=1
n[n < -1]= -1
color<- colorRampPalette(c('#87CEFA','white','#FF4500'))(100)
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata,color=color,
         main = "")
pic_data<-data.frame(t(pic_data))
colnames(pic_data)<-rownames(xcell)
annodata$group<-ifelse(annodata$condition=="High","High Risk","Low Risk")
pic_data<-cbind(pic_data,group=annodata$group)
#Differences in immune cell infiltration between the high and low TMERS groups
library(ggplot2)
library(ggpubr)
for(i in 1:22){
  p<-ggplot(data = pic_data,aes(x=group,y=pic_data[,i],fill=group))+
    geom_violin(position = position_dodge(width = 1), scale = 'width',trim = FALSE)+
    geom_boxplot(position = position_dodge(width = 1), 
                 outlier.size = 0.7, width = 0.2, show.legend = FALSE)+
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), 
          legend.title = element_blank(), legend.key = element_blank(),legend.position = "top") +
    guides(fill="none")+
    labs(x="",y=paste("Relative",colnames(pic_data)[i],"amount"),
         title = paste(colnames(pic_data)[i]))+
    stat_compare_means(method = "t.test",label.y = min(pic_data[,i]),label.x = 1)
  ggsave(paste("immune/cell/",colnames(pic_data)[i],".png",sep = ""),p,width = 2.5,height = 2.5)
  ggsave(paste("immune/cell/",colnames(pic_data)[i],".pdf",sep = ""),p,width = 2.5,height = 2.5)
}
#Correlation between the degree of immune cell infiltration and TMERS, prognostic genes, immune score, stromal score, and estimated score
score_tumor<-score_tumor[,-(4:5)]
row.names(score_tumor)<-gsub("\\.","-",row.names(score_tumor))
score_tumor$sample<-row.names(score_tumor)
score_tumor$sample<-str_sub(score_tumor$sample,1,12)

gene6<-deg[sig_coef$gene,]
gene6<-data.frame(t(gene6))
colnames(gene6)<-sig_coef$gene
TME_cor<-cbind(TMEscore,gene6)
TME_cor<-merge(TME_cor,score_tumor,by.x = "sample",by.y = "sample")
pic_data$sample<-row.names(pic_data)
pic_data<-pic_data[,-23]
TME_cor<-merge(pic_data,TME_cor,by.x = "sample",by.y = "sample")
colnames(TME_cor)[24]<-"RiskScore"
TME_cor<-aggregate(x=TME_cor[,2:33],by=list(TME_cor$sample),FUN = mean)
colnames(TME_cor)[1]<-"sample"
p_matrix<-data.frame(matrix(0,10,22,dimnames = list(colnames(TME_cor)[24:33],colnames(TME_cor)[2:23])))
colnames(p_matrix)<-colnames(TME_cor)[2:23]

cor_matrix<-data.frame(matrix(0,10,22,dimnames = list(colnames(TME_cor)[24:33],colnames(TME_cor)[2:23])))
colnames(cor_matrix)<-colnames(TME_cor)[2:23]
for (i in 24:33) {
  for (j in 2:23) {
    p<-cor.test(TME_cor[,i],
                  TME_cor[,j],method = "spearman")$p.value
    p_matrix[i-23,j-1]<-p
    
    cor<-cor.test(TME_cor[,i],
                  TME_cor[,j],method = "spearman")$estimate
    cor_matrix[i-23,j-1]<-cor
  }
}
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
pheatmap(cor_matrix,
         display_numbers = matrix(ifelse(p_matrix < 0.01,"**",
                                         ifelse(p_matrix<0.05,"*","")), 
                                  nrow(p_matrix)),
         fontsize_number = 10,
         cluster_rows = F,cluster_cols = F,
         show_rownames = T,show_colnames = T,angle_col = 45,
         cellwidth = 20, cellheight = 12,color=color)
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
picdata<-data.frame(t(TME_cor))
colnames(picdata)<-picdata[1,]
picdata<-picdata[-1,]
picdata<-picdata[-33,]
picdata1<-picdata[,c(rownames(annodata)[which(annodata$condition=="High")],
                     rownames(annodata)[which(annodata$condition=="Low")])]

annodata1<-data.frame(row.names = row.names(annodata),condition=annodata[,2])
picdata1[,1:1069]<-lapply(picdata1[,1:1069],as.numeric)

n=t(scale(t(picdata1[24:29,])))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,color=color,
         main = "",cellwidth = 0.3, cellheight = 10)

n=t(scale(t(picdata1[1:22,])))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,color=color,
         main = "",cellwidth = 0.3, cellheight = 10)

n=t(scale(t(picdata1[30:32,])))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1,color=color,
         main = "",cellwidth = 0.3, cellheight = 10)
###############
pic_data<-TME_cor[,c(1,25:33)]
annodata$sample<-row.names(annodata)
annodata1<-merge(annodata[,2:3],class_map[,-6],by.x = "sample",by.y = "patient_id")
annodata1<-annodata1[order(annodata1$group),]
row.names(annodata1)<-annodata1$sample

annodata_high<-annodata1[which(annodata1$group=="High Risk"),]
annodata_low<-annodata1[which(annodata1$group=="Low Risk"),]

annodata_high$N<-annodata_high$N[order(annodata_high$N)]
annodata_high$M<-annodata_high$M[order(annodata_high$M)]
annodata_high$T<-annodata_high$T[order(annodata_high$T)]
annodata_high$Stage<-annodata_high$Stage[order(annodata_high$Stage)]
annodata_high$Subtype<-annodata_high$Subtype[order(annodata_high$Subtype)]

annodata_low$M<-annodata_low$M[order(annodata_low$M)]
annodata_low$N<-annodata_low$N[order(annodata_low$N)]
annodata_low$Stage<-annodata_low$Stage[order(annodata_low$Stage)]
annodata_low$T<-annodata_low$T[order(annodata_low$T)]
annodata_low$Subtype<-annodata_low$Subtype[order(annodata_low$Subtype)]

annodata2<-rbind(annodata_high,annodata_low)

pic_data<-data.frame(t(pic_data))
colnames(pic_data)<-pic_data[1,]
pic_data<-pic_data[-1,]
pic_data<-pic_data[,annodata2$sample]
pic_data[,1:850]<-lapply(pic_data[,1:850],as.numeric)
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)
n1=t(scale(t(pic_data[1:6,])))
n1[n1>1]=1
n1[n1< -1]= -1
n2=t(scale(t(pic_data[7:9,])))
n2[n2>1]=1
n2[n2< -1]= -1
pheatmap(n1,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata2[,2:7],color=color,
         main = "",cellwidth = 0.3, cellheight = 10,
         filename = "immune/gene_pheatmap.pdf",width=8,heigh=4)

pheatmap(n2,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,color=color,
         main = "",cellwidth = 0.3, cellheight = 10,
         filename = "immune/score_pheatmap.pdf",width=8,heigh=4)

pheatmap(n1,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         annotation_col = annodata1[,2:7],color=color,
         main = "",cellwidth = 0.3, cellheight = 10)

library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#87CEFA','white','red'))(100)

pdl<-surv_symbol[c("CD274","CTLA4","PDCD1"),]
pdl<-pdl[,row.names(annodata)]
row.names(pdl)<-c("PD-L1","CTLA-4","PD-1")

n=t(scale(t(pdl)))
n[n>1]=1
n[n< -1]= -1
pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = F,
         color=color,
         main = "",cellwidth = 0.3, cellheight = 10,
         filename = "immune/Immune checkpoint distribution.pdf",width=8,heigh=5)
		 
#The prognostic signature positively correlated with tumor mutation burden
dir.create("TMB")
library(maftools)
laml<-read.maf("TMB/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf")

plotmafSummary(maf = laml,rmOutlier = T,
               addStat = "median",dashboard = T,titvRaw = F)
sample_variant<-laml@variants.per.sample
library(stringr)
sample_variant$Tumor_Sample_Barcode<-str_sub(sample_variant$Tumor_Sample_Barcode,1,12)
sample_variant<-merge(score_surv,sample_variant,by.x = "sample",by.y = "Tumor_Sample_Barcode")
#Distribution of TMB between high and low risk scores groups
sample_variant$group<-ifelse(sample_variant$score>median(sample_variant$score),"High TMERS","Low TMERS")
p<-ggboxplot(sample_variant,x="group",y="Variants",
             title = "",ylab = "TMB",xlab = "",
             color = "group",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(sample_variant$Variants)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = 1)
#Distribution of risk scores between high and low TMB groups
sample_variant$group<-ifelse(sample_variant$Variants>median(sample_variant$Variants),"High TMB","Low TMB")
p<-ggboxplot(sample_variant,x="group",y="score",
             title = "",ylab = "TMERS",xlab = "",
             color = "group",palette = "jama",
             add = "jitter",outlier.shape = NA)+
  coord_cartesian(ylim = (boxplot.stats(sample_variant$score)$stats[c(1, 5)])*1.05)+
  stat_compare_means(label.y = -1.3)
#Correlation between risk score and TMB
p<-ggplot(data=sample_variant, aes(x=score, y=Variants))+geom_point(color="blue")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  stat_cor(data=sample_variant, method = "spearman")+
  labs(x="TMERS",y="TMB",title = "")
ggsave("TMB/TMB_TMERSspearman.pdf",p,width = 4,height = 4)

#Plot the cascade of mutated genes grouped into high and low risk groups
sample_variant<-laml@variants.per.sample
library(stringr)
sample_variant$Sample<-str_sub(sample_variant$Tumor_Sample_Barcode,1,12)
sample_variant<-merge(score_surv,sample_variant,by.x = "sample",by.y = "Sample")

sample_variant$group<-ifelse(sample_variant$score>median(sample_variant$score),"High TMERS","Low TMERS")
sample_variant$Tumor_Sample_Barcode<-as.character(sample_variant$Tumor_Sample_Barcode)

sampeOrder<-sample_variant$Tumor_Sample_Barcode[which(sample_variant$group=="High TMERS")]

#high TMERS
oncoplot(maf = laml,top = 20,borderCol = NULL,keepGeneOrder = T,
         sampleOrder = sampeOrder,titleText = "High TMERS")


sampeOrder<-sample_variant$Tumor_Sample_Barcode[which(sample_variant$group=="Low TMERS")]
#low TMERS
oncoplot(maf = laml,top = 20,borderCol = NULL,
         sampleOrder = sampeOrder,titleText = "Low TMERS")


#The TIDE algorithm was used to estimate the TIDE score for each sample
{
dir.create("metabric")
library(survival)
CTL_genes = c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1')
readmat = function(mat) as.matrix(read.table(mat, sep='\t', header=T, check.names=F, quote=NULL))

writemat = function(result, output)
{
  result = result[!is.na(result[,"p"]),,drop=F]
  FDR = p.adjust(result[,"p"], method="fdr")
  result = cbind(result, FDR)
  write.table(result, file=output, sep='\t', quote=F)
}
#commands = commandArgs(trailingOnly=T)
mat = t(as.matrix(metabric_all))
pivot<-t(metabric_all[CTL_genes,])
survival = data.frame(clinical)
rownames(survival)<-survival$sample
survival<-survival[,-1]
survival = survival[survival[,1] > 0,]
output = "metabric"
common = Reduce(intersect, list(rownames(mat),rownames(survival),rownames(pivot)))
print(paste(length(common), "samples"))
if(length(common) < 20) stop("two few samples")

not_included = setdiff(CTL_genes, colnames(pivot))

if(length(not_included) > 0) stop(paste(c("pivot genes", not_included, "are not included"), ' '))

pivot = rowMeans(pivot[common,CTL_genes])
mat = mat[common,,drop=F]
survival = survival[common,,drop=F]
death_rate = sum(survival[,2])/dim(survival)[1]

if(death_rate < 0.1) stop(paste("death rate", death_rate, "is too low"))

surv = Surv(survival[,1], survival[,2])

if(dim(survival)[2] > 2){
  B = survival[,3:dim(survival)[2], drop=F]
}else{
  B = survival[,c(), drop=F]
}

B = cbind(B, pivot, pivot, pivot)
B = as.data.frame(B)
N_B = dim(B)[2]

coxph.pivot = summary(coxph(surv~., data=B[,c(-N_B, -(N_B-1)), drop=F]))$coef
write.csv(coxph.pivot,file=paste0(output, "_coxph.csv"))

colnames(B)[N_B-1] = "partner"
colnames(B)[N_B] = "Interaction"

features = colnames(mat)
N = length(features)

result_interaction = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_main = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_partial = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_base = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))

for (i in 1:N){
  title = features[i]
  partner = mat[,i]
  
  B[,N_B-1] = partner
  B[,N_B] = partner * pivot
  
  # region 1: model with interaction
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_interaction[i,] = reg.summary["Interaction", c("z", "Pr(>|z|)")]
    result_main[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  
  # region 2: model without interaction
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,-N_B]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_partial[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  
  # region 3: model with only base effect
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,c(-N_B, -(N_B-2)), drop=F]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_base[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
}
writemat(result_interaction, paste0(output, "_interaction.txt"))
writemat(result_main, paste0(output, "_main.txt"))
writemat(result_partial, paste0(output, "_partial.txt"))
writemat(result_base, paste0(output, "_base.txt"))
}
##########The relevance of T Cell Dysfunction was calculated through interaction
result_interaction<-data.frame(result_interaction)
result_diff<-result_interaction[which(result_interaction$p<0.05),]
dys_expr<-surv_symbol[row.names(result_diff),]


BRCA_TIDE<-data.frame(matrix(0,1069,2))
BRCA_TIDE$sample<-colnames(dys_expr)
colnames(BRCA_TIDE)<-c("Dysfunction","Exclusion","Sample")
for(i in 1:1069){
  BRCA_TIDE[i,1]<-cor(dys_expr[,i],result_diff[,1],method="pearson")
}
#############Calculate the correlation of T Cell Exclusion
for(i in 1:1069){
  BRCA_TIDE[i,2]<-cor(excl_expr[,i],excl_ave[,4],method="pearson")
}
tumor_TIDE<-merge(BRCA_TIDE,CTL_surv,by.x = "Sample",by.y = "sample")
tumor_TIDE$TIDE<-ifelse(tumor_TIDE$CTL>mean(tumor_TIDE$CTL),tumor_TIDE$Dysfunction,tumor_TIDE$Exclusion)
tumor_TIDE$group<-ifelse(tumor_TIDE$CTL>mean(tumor_TIDE$CTL),"High","Low")
tumor_TIDE$scale<-scale(tumor_TIDE$TIDE)

#############Correlation between risk score and PD-1,PD-L1,CTLA4
#PD-L1
pdl<-surv_symbol[c("CD274","CTLA4","PDCD1"),]
pdl<-data.frame(t(pdl))
colnames(pdl)<-c("PD-L1","CTLA4","PD-1")
pdl<-cbind(TMEscore,pdl)
library(ggplot2)
library(ggpubr)
p<-ggplot(data=pdl, aes(x=score, y=pdl[,5]))+
  geom_point(color="#EACD50")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  stat_cor(data=pdl, method = "spearman")+
  labs(x="TMERS",y="PD-L1",title = "PD-L1")

ggsave("gene/TMERS_PD-L1.png",p,width = 4,height = 4)
ggsave("gene/TMERS_PD-L1.pdf",p,width = 4,height = 4)

p<-ggplot(data=pdl, aes(x=score, y=pdl[,6]))+
  geom_point(color="#EACD50")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  stat_cor(data=pdl, method = "spearman")+
  labs(x="TMERS",y="CTLA4",title = "CTLA4")

ggsave("gene/TMERS_CTLA4.png",p,width = 4,height = 4)
ggsave("gene/TMERS_CTLA4.pdf",p,width = 4,height = 4)

p<-ggplot(data=pdl, aes(x=score, y=pdl[,7]))+
  geom_point(color="#EACD50")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  stat_cor(data=pdl, method = "spearman")+
  labs(x="TMERS",y="PD-1",title = "PD-1")

ggsave("gene/TMERS_PD-1.png",p,width = 4,height = 4)
ggsave("gene/TMERS_PD-1.pdf",p,width = 4,height = 4)