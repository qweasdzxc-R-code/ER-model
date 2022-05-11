

# ER prognostic model R transcripts

# part 1: outline
#(1)TCGA exp download
#(2) clini download
#(3) Target gene-exp+clin combi
#(4) GSE14520 download-sort
#(5) TCGA-GSE14520 batch match
#(6) TCGA-Normal-Tumor-diff
#(7) GEO-Normal-Tumor-diff
#(8) TCGA-GEO-uni cox
#(9) Lasso modeling
#(10) uni/multi cox survival -KM
#(11) ROC
#(12) Nomogram
#(13) risk-clinical-differ
#(14) GO-KEGG-GSEA
# part 2: outline
#(15) clusting analysis
#(16) immune cell infiltration analysis
#(17) clinical information sorting-exp-combi




rm(list=ls()) 
work_dir <- "#/#/#" 
setwd(work_dir)

#(15) clusting analysis#######################################################################################
library("org.Hs.eg.db")
library("clusterProfiler")
library(SummarizedExperiment)  
library(TCGAbiolinks)  

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")

## try http:// if https:// URLs are not supported
#source("http://bioconductor.org/biocLite.R")
#biocLite("ConsensusClusterPlus")
#browseVignettes("ConsensusClusterPlus")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

rm(list=ls())  

rt_1=read.table("0_TCGA_training-riskscore.txt",sep="\t",header=T,check.names=F)              
fix(rt_1)
dim(rt_1)
rt_2<-rt_1[,c("id", "KPNA2", "FMO3", "SPP1", "KIF2C", "LPCAT1" )]
fix(rt_2)
dim(rt_2)   
rt_3<-t(rt_2)
fix(rt_3)
dim(rt_3)  
colnames(rt_3)<- rt_3[1,] 
rt_4<-rt_3[-1,] 
fix(rt_4)
dim(rt_4)  
library(limma)
rt_4=as.matrix(rt_4)
rownames(rt_4)=rt_4[,1]
exp=rt_4[,2:ncol(rt_4)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
fix(data)
dim(data)  
data=as.matrix(data)
maxK=9  
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="png")
clusterNum=2                 
cluster_2=results[[clusterNum]][["consensusClass"]]
cluster_2<-as.data.frame(cluster_2)
fix(cluster_2)
dim(cluster_2)
write.table(cluster_2,file="1_TCGA-two-group_cluster(PCA).txt",sep="\t",quote=F,col.names = F)
write.table(cluster_2,file="2_TCGA-cluster-surv-trait-enrichment.txt",sep="\t",quote=F)

clusterNum=3                 
cluster_3=results[[clusterNum]][["consensusClass"]]
cluster_3<-as.data.frame(cluster_3)
fix(cluster_3)
dim(cluster_3)
write.table(cluster_3,file="3_TCGA-three-group_cluster(PCA).txt",sep="\t",quote=F,col.names = F)
write.table(cluster_3,file="4_TCGA-three-group_cluster(surv-trait-enrichment).txt",sep="\t",quote=F)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("ggplot2")
#pca analysis
library(limma)
rt=read.table("0_TCGA-LIHC_batch-2.txt",sep="\t",header=T,check.names=F)    
fix(rt)
dim(rt)  
TCGA_1=avereps(rt)
TCGA_1=TCGA_1[rowMeans(TCGA_1)>0.5,]   
fix(TCGA_1)
dim(TCGA_1) 
dim(data)  
HCC_sample<-colnames(data) 
TCGA_2<-TCGA_1[,HCC_sample] 
dim(TCGA_2) 
TCGA_3<-t(TCGA_2) 
fix(TCGA_3)
dim(TCGA_3) 
TCGA_3.class <- rownames(TCGA_3)
TCGA_3.pca <- prcomp(TCGA_3, scale. = TRUE)                             
library(ggplot2)
cluster_2a=read.table("1_TCGA-two-group_cluster(PCA).txt",sep="\t",header=F,check.names=F)    
group=paste0("cluster",as.vector(cluster_2a[,2]))
pcaPredict_2a=predict(TCGA_3.pca)
fix(pcaPredict_2a)
dim(pcaPredict_2a) 
PCA_2a = data.frame(PCA1 = pcaPredict_2a[,1], PCA2 = pcaPredict_2a[,2],group=group)
tiff(file="6_TCGA-PCA_two-group.tiff",width =20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data = PCA_2a, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
library(ggplot2)
cluster_3a=read.table("3_TCGA-three-group_cluster(PCA).txt",sep="\t",header=F,check.names=F)     
group=paste0("cluster",as.vector(cluster_3a[,2]))
pcaPredict_3a=predict(TCGA_3.pca)
fix(pcaPredict_3a)
dim(pcaPredict_3a) 
PCA_3a = data.frame(PCA1 = pcaPredict_3a[,1], PCA2 = pcaPredict_3a[,2],group=group)
tiff(file="7_TCGA-PCA_three-group.tiff",width =20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data = PCA_3a, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
rm(list=ls())  
library(survival)
cluster_1=read.table("0_TCGA-cluster-surv-trait-enrichment.txt",header=T,sep="\t",check.names=F,row.names=1)
cli_1=read.table("0_TCGA_training-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
fix(cli_1)
dim(cli_1)
colnames(cli_1)
cli_2<-cli_1[,c("os_time", "os_stat")] 
fix(cli_2)
dim(cli_2) 
cli_3<-cli_2[row.names(cluster_1),] 
dim(cli_3)
fix(cli_3) 
aaa <- data.frame(cli_3,cluster_1[match(rownames(cluster_1),rownames(cli_3)),])
dim(aaa)
fix(aaa)
colnames(aaa)[3]="cluster"
dim(aaa)
fix(aaa) 
write.table(aaa,file="1_TCGA-two-cluster-OS-comb.txt",sep="\t",quote=F)
clusterNum=2    
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }

fit_2 <- survfit(Surv(os_time, os_stat) ~ cluster, data = aaa)
#summary(fit)
tiff(file="2_TCGA_cluster-OS-differ(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_2 , 
     lwd=4, 
     col=c("#0066FF","#FF1493"),
     xlab="Time (year)",  
     mark.time=T,
     bg="yellow",
     ylab="Survival rate",
     
     main=paste("TCGA-OS (p=", pValue ,")",sep=""))
legend("topright", 
       paste0("cluster",1:clusterNum), 
       lwd=4, 
       bty = "n", 
       col=c("#0066FF","#FF1493"))
dev.off()
library(survival)
library("survminer") 
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit2a <- survfit(Surv(os_time, os_stat) ~ cluster, data = aaa)
tiff(file="2_TCGA_cluster-OS-differ(ggsurvplot).tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) #设置图的大小
surPlot=ggsurvplot(fit2a, 
                   data=aaa,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("cluster1", "cluster2"),
                   legend.title="TCGA_OS",
                   xlab="os_time(years)",
                   break.os_time.by = 1,
                   risk.table.title="",
                   palette=c("#0066FF","#2E8B57"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
cluster_2=read.table("0_TCGA-three-group_cluster(surv-trait-enrichment).txt",header=T,sep="\t",check.names=F, row.names=1)
fix(cluster_2)
dim(cluster_2) 
cli_1=read.table("0_TCGA_training-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
fix(cli_1)
dim(cli_1) 
colnames(cli_1)
cli_2<-cli_1[,c("os_time", "os_stat")] 
fix(cli_2)
dim(cli_2) 
cli_3<-cli_2[row.names(cluster_2),] 
dim(cli_3)
fix(cli_3) 
bbb <- data.frame(cli_3,cluster_2[match(rownames(cluster_2),rownames(cli_3)),])
dim(bbb)
fix(bbb)
colnames(bbb)[3]="cluster"
dim(bbb)
fix(bbb) 
write.table(aaa,file="3_TCGA-three-group_OS-surv-comb.txt",sep="\t",quote=F)
clusterNum=3    
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = bbb) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }
fit_3 <- survfit(Surv(os_time, os_stat) ~ cluster, data = bbb)
#summary(fit)
tiff(file="4_TCGA_three-group-OS-surv-diff(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_3 , 
     lwd=4, 
     col=c("#0066FF","#FF1493","#006600"),
     xlab="Time (year)",  
     mark.time=T,
     ylab="Survival rate",
     main=paste("TCGA-OS (p=", pValue ,")",sep=""))
legend( "bottomleft",    #"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
        paste0("cluster",1:clusterNum), 
        lwd=4, 
        bty = "n", 
        col=c("#0066FF","#FF1493","#006600"))
dev.off()
rm(list=ls())  
cluster_1=read.table("0_TCGA-cluster-surv-trait-enrichment.txt",sep="\t",header=T,check.names=F, row.names=1) 
fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
fix(cluster_2)
dim(cluster_2) 
cli_1=read.table("0_TCGA_OS_model-gene_riskscore_chr.txt",sep="\t",header=T,check.names=F, row.names=1)
fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
fix(cli_2)
dim(cli_2) 

cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
fix(cluster_cli_2) 
write.table(cluster_cli_2, file="1_TCGA_OS_cluster（twogroup）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_2"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
TCGA_OS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_2","gender", "age","stage")])  
fix(TCGA_OS_clin_rt)
dim(TCGA_OS_clin_rt)  
colnames(TCGA_OS_clin_rt)
write.table(TCGA_OS_clin_rt, file="2_TCGA_OS_cluster（twogroup）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
TCGA_OS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
fix(TCGA_OS_cluster_rt)
dim(TCGA_OS_cluster_rt) 
colnames(TCGA_OS_cluster_rt)
write.table(TCGA_OS_cluster_rt, file="2_TCGA_OS_cluster（twogroup）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
field="cluster_2"  
flag1="1"
flag2="2"
trainFlag=TCGA_OS_clin_rt[TCGA_OS_clin_rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=TCGA_OS_clin_rt[TCGA_OS_clin_rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)
fix(newTable) 
newLabels=c("cluster_2")
for(i in 2:(ncol(TCGA_OS_clin_rt)) ){ 
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=fisher.test(tableStat)  
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels 
#fix(newTable)
colnames(TCGA_OS_clin_rt)=newLabels
#fix(TCGA_OS_clin_rt) 
dim(TCGA_OS_clin_rt)
write.table(TCGA_OS_clin_rt,file="3_TCGA_OS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",quote=F)
#install.packages("pheatmap")
library(pheatmap)
rt=read.table("2_TCGA_OS_cluster（twogroup）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("3_TCGA_OS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)

Type=Type[order(Type$cluster_2),] 
fix(Type)
rt_1=rt[,row.names(Type)] 
fix(rt_1)
TCGA_OS_cluster_rt_1=read.table("2_TCGA_OS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
fix(TCGA_OS_cluster_rt_1)
gxj_table=table(TCGA_OS_cluster_rt_1[,c("stage","cluster_2")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="4_clustertwogroup_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("TCGA-OS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
TCGA_OS_cluster_rt_2=read.table("2_TCGA_OS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
fix(TCGA_OS_cluster_rt_2)
TCGA_OS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(TCGA_OS_cluster_rt_2)){
  nameStat=colnames(TCGA_OS_cluster_rt_2)[i]
  tableStat=table(TCGA_OS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  TCGA_OS_cluster_rt_2StatOut=rbind(TCGA_OS_cluster_rt_2StatOut,tableStatOut)}
gxj_StatOut=TCGA_OS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="5_TCGA_OS_cluster（twogroup）-num.xls",sep="\t",quote=F,row.names=F)
rm(list=ls())  
cluster_1=read.table("0_TCGA-three-group_cluster(surv-trait-enrichment).txt",sep="\t",header=T,check.names=F, row.names=1) 
fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
fix(cluster_2)
dim(cluster_2) 
cli_1=read.table("0_TCGA_OS_model-gene_riskscore_chr.txt",sep="\t",header=T,check.names=F, row.names=1)
fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
fix(cli_2)
dim(cli_2) 
cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
fix(cluster_cli_2) 
write.table(cluster_cli_2, file="6_TCGA_OS__cluster（three-group）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_3"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
TCGA_OS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_3","gender", "age","stage")])  
fix(TCGA_OS_clin_rt)
dim(TCGA_OS_clin_rt)  
colnames(TCGA_OS_clin_rt)
write.table(TCGA_OS_clin_rt, file="7_TCGA_OS__cluster（three-group）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
TCGA_OS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
fix(TCGA_OS_cluster_rt)
dim(TCGA_OS_cluster_rt) 
colnames(TCGA_OS_cluster_rt)
write.table(TCGA_OS_cluster_rt, file="8_TCGA_OS__cluster（three-group）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
#install.packages("pheatmap")
library(pheatmap)
rt=read.table("8_TCGA_OS__cluster（three-group）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("7_TCGA_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)
Type=Type[order(Type$cluster_3),] 
fix(Type)
rt_1=rt[,row.names(Type)] 
fix(rt_1)
TCGA_OS_cluster_rt_1=read.table("7_TCGA_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
fix(TCGA_OS_cluster_rt_1)
gxj_table=table(TCGA_OS_cluster_rt_1[,c("stage","cluster_3")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="9_cluster（three-group）_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("TCGA-OS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
TCGA_OS_cluster_rt_2=read.table("7_TCGA_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
fix(TCGA_OS_cluster_rt_2)
TCGA_OS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(TCGA_OS_cluster_rt_2)){
  nameStat=colnames(TCGA_OS_cluster_rt_2)[i]
  tableStat=table(TCGA_OS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  TCGA_OS_cluster_rt_2StatOut=rbind(TCGA_OS_cluster_rt_2StatOut,tableStatOut)}
gxj_StatOut=TCGA_OS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="10_TCGA_OS__cluster（three-group）-num.xls",sep="\t",quote=F,row.names=F)
rm(list=ls())  
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0_TCGA-cluster-surv-trait-enrichment.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1)  
com = subset(group_1, cluster_2 =="1"|cluster_2  =="2",  select=cluster_2 )  
dim(com)   
#fix(com)   
colnames(com)[1]<-c("group")
head(com,5)
gene_1<-read.table("0_TCGA-LIHC_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
#fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
#fix(group_list) 
rownames(group_list) = rownames(com)        
colnames(group_list)[1] = "description"   
dim(group_list) 
#fix(group_list) 
design = model.matrix(~ description + 0, group_list)   
#fix(design)
colnames(design) = levels(as.factor(c("cluster1","cluster2")))   
factor(design) 
head(design)  
cont.matrix = makeContrasts(cluster2-cluster1, levels=design)   
fdr = 0.05
logFC = 1
fit = lmFit(eSet, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="BH", sort.by="logFC", resort.by = "P" ,n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "Pos" 
tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Neg"  
tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "NO"
colnames(tT)[ncol(tT)] = "expression"
dim(tT)
write.table(tT,file="1_TCGA_cluster2-cluster1_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
model_gene_1<-read.table("0_lasso_(geneCoef).txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(model_gene_1)  
#fix(model_gene_1)  
model_gene<-row.names(model_gene_1) 
model_gene
diff_gene<-row.names(deg) 
diff_gene[1:10]  
aa_1<-setdiff(diff_gene, model_gene)  
dim(deg)  
deg_1<-deg[aa_1,] 
dim(deg_1)    
write.table(deg_1,file="2_TCGA_cluster2-cluster1_diff-gene.xls",sep="\t",quote=F)
write.table(deg_1,file="2_TCGA_cluster2-cluster1_diff-gene.txt",sep="\t",quote=F)
group_reg <- factor(deg_1$expression)
labels=levels(group_reg)
gxj<-summary(group_reg)
gxj
fix(tT)
dim(tT)
aa_2<-setdiff(row.names(tT), model_gene) 
tT_1<-tT[aa_2,]
fix(tT_1)
dim(tT_1)
temp1 = tT_1[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("Pos","Neg","NO"), order=T)
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+ 
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080" ))+  
  labs(x=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]),y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  #xlim(-5,5)+  
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = "#FFE1FF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "#FFE1FF"))
tiff(file="3_TCGA_cluster-volcano_style1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_TCGA_cluster-volcano_style2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+  
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080"))
dev.off()
temp2 = data.frame(log10(rowMeans(eSet))) 
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT_1),rownames(temp2)),] 
temp2 = data.frame(temp2,tT_1[,c("logFC","expression")])
colnames(temp2)=c("log10Exp","logFC","expression")
temp2$expression=factor(temp2$expression, levels=c("Pos","Neg","NO"), order=T)
P_MA=ggplot(temp2,aes(x=temp2$log10Exp,y=temp2$logFC))+
  geom_point(alpha=0.8, size=6, 
             aes(color=temp2$expression))+
  scale_color_manual(values =c("Pos" = "#FF4500", "NO" = "gray","Neg" = "#483D8B"))+
  labs(x="log10Exp",y=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+
  # xlim(0.75,1.15)+  
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "black"))
tiff(file="4_TCGA_cluster-MA.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)  
print(P_MA)
dev.off()
dim(eSet) 
DEG_list = rownames(deg_1)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),] 
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_TCGA_cluster-all-gene-heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset, color = colorRampPalette(c('blue2',  'white',  'red'))(10000),
                   cluster_rows = T, cluster_cols=T,                  
                   scale ="row" , 
                   show_rownames=F,show_colnames=F, 
                   treeheight_row = 0,  
                   clustering_distance_rows = "correlation",  
                   clustering_distance_cols = "euclidean",   
                   border_color = "grey",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()
target_gene<-row.names(deg_1)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_TCGA_cluster-target-gene-heatmap_(TOP10-remove-hub).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset_2, color = colorRampPalette(c('#000066',  'white',  '#FF0000'))(60), 
                   cluster_rows = T, #clustering_distance_rows = "correlation", 
                   cluster_cols=T,   clustering_distance_cols = "euclidean",          
                   scale ="row" , 
                   show_rownames=T,  
                   show_colnames= F, 
                   treeheight_row = 0,  
                   border_color = "White",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()
#KPNA2; FMO3; SPP1; KIF2C; LPCAT1
model_gene_1<-read.table("0_lasso_(geneCoef).txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(model_gene_1)  
#fix(model_gene_1)  
deg_eset_2<-eSet[row.names(model_gene_1),]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_TCGA_cluster-model-heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset_2, color = colorRampPalette(c('#000066',  'white',  '#FF0000'))(60), 
                   cluster_rows = T, #clustering_distance_rows = "correlation", 
                   cluster_cols=T,   clustering_distance_cols = "euclidean",          
                   scale ="row" , 
                   show_rownames=T,  
                   show_colnames= F, 
                   treeheight_row = 0,  
                   border_color = "White",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()
rm(list=ls())  
library("colorspace")
library("stringi")
library("rlang")
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
library("pathview")
library("GOplot")
clx_kegg1="9_FC_map"   
dir.create(paste(workdir,clx_kegg1,sep="/"))  
clx_kegg2="10_Pvalue_map"   
dir.create(paste(workdir,clx_kegg2,sep="/"))  
gxj <- read.table("0_TCGA_cluster2-cluster1_diff-gene-remove-hub-gene.txt",header=T,comment.char = "",check.names=F,sep="\t", row.names=1 ) 
dim(gxj)
fix(gxj)
colnames(gxj)
gxj$GeneSymbol<-row.names(gxj)  
transfered_ID <- bitr(gxj$GeneSymbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
transfered_ID<-data.frame(transfered_ID) 
fix(transfered_ID)
Final_ID_acorss<- merge(gxj, transfered_ID,  by.x="GeneSymbol", by.y = "SYMBOL") 
dim(Final_ID_acorss)  
fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"ENTREZID"])==F,]   
dim(Final_ID_acorss)  
aaa<-duplicated(Final_ID_acorss$ENTREZID)   
Final_ID_acorss<-Final_ID_acorss[!aaa,]
dim(Final_ID_acorss) 
fix(Final_ID_acorss)
write.table(Final_ID_acorss, file = "0_ID_change.txt",sep="\t", row.names =F, quote = F)
gxj <- read.table("0_ID_change.txt",header=T,comment.char = "",check.names=F,sep="\t")
fix(gxj)
DEG_list <- gxj$ENTREZID
kk <- enrichGO(gene = DEG_list,
               keyType ="ENTREZID",  
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="0_all_GO.txt",sep="\t",quote=F,row.names = F)
pdf(file="0_GO_barplot.pdf",width = 30,height = 15)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="0_GO_bubble.pdf",width =30,height = 15)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
ego_MF <- enrichGO(gene          = DEG_list ,
                   keyType       ="ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "MF",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,  
                   readable      = TRUE)  
write.table(ego_MF, file = "1_GO_MF_enrichment.txt",sep="\t", row.names =F, quote = F)
#ego_results<-summary(ego)
tiff(file="2_GO_MF_barplot_top10.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GO_MF_barplot_top20.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_MF_dotplot_top10.tiff",width = 100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=10)
dev.off()
tiff(file="3_GO_MF_dotplot_top20.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=20)
dev.off()
tiff(file="4_GO_MF_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF)
dev.off()
tiff(file="4_GO_MF_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_MF_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_MF_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_MF_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GO_MF_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_MF_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GO_MF_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_MF,showCategory = 20, foldChange =NULL)
dev.off()
tiff(file="7_GO_MF_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_MF, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
ego_BP <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "BP",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.table(ego_BP, file = "1_GO_BP_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GO_BP_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GO_BP_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_BP_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=10)
dev.off()
tiff(file="3_GO_BP_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=20)
dev.off()
tiff(file="4_GO_BP_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP)
dev.off()
tiff(file="4_GO_BP_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_BP_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_BP_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_BP_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GO_BP_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_BP_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GO_BP_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_BP,showCategory = 20, foldChange =NULL)
dev.off()
tiff(file="7_GO_BP_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_BP, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
ego_CC <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "CC",  
                   pAdjustMethod = "BH",
                   #universe = names(geneList), 
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE) 
write.table(ego_CC, file = "1_GO_CC_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GO_CC_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_CC_CC_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_CC_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=10)
dev.off()
tiff(file="3_GO_CC_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=20)
dev.off()
tiff(file="4_GO_CC_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC)
dev.off()
tiff(file="4_GO_CC_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_CC_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_CC_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
#cnetplot(ego_CC, showCategory = 5, categorySize="pvalue", foldChange=geneList)
cnetplot(ego_CC, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_CC_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GO_CC_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GO_CC_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GO_CC_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_CC,showCategory = 20, foldChange = NULL)
dev.off()
tiff(file="7_GO_CC_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_CC, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
library("pathview")
kegg <- enrichKEGG(gene         = gxj$ENTREZID,
                   organism     = 'hsa',   
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
write.table(kegg, file = "8_KEGG_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="9_KEGG_barplot.tiff",width = 35,height = 22,units ="cm",compression="lzw",bg="white",res=300)
barplot(kegg, drop = T, showCategory = 20,color = "pvalue")
dev.off()
tiff(file="10_KEGG_dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)  
dotplot(kegg,title="Enrichment KEGG_dot")
dev.off()
setwd(paste(workdir,clx_kegg1,sep="/"))
map_ids = as.matrix(gxj['logFC'])
for (pathway_id in kegg$ID){
  map = pathview(gene.data  = map_ids,
                 pathway.id = pathway_id,
                 species    = "hsa", kegg.native = TRUE)}
setwd(paste(workdir,clx_kegg2,sep="/"))
grd<-gxj$P.Value   
names(grd)<-gxj$ENTREZID 
for(i in kegg$ID){
  pv.out <- pathview(gene.data = -log10(grd), 
                     pathway.id = i, 
                     species = "hsa", 
                     out.suffix = "pathview",
                     limit=list(gene=5, cpd=5)  
  )}
rm(list=ls())  
library("colorspace")
library("stringi")
library("rlang")
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOplot")
library("stringr")
gxj<- read.table("0_ER-surv_cluster2_vs._cluster1_diff-gene --EXCEL-sort-Symbol_logFC.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(gxj)
fix(gxj)
colnames(gxj)[c(1)]=c("ID")  
colnames(gxj)
transfered_ID <- bitr(gxj$ID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL" ), OrgDb="org.Hs.eg.db")
transfered_ID<-data.frame(transfered_ID) 
#fix(transfered_ID)
Final_ID_acorss<- merge(gxj, transfered_ID,  by.x="ID", by.y = "SYMBOL") 
dim(Final_ID_acorss)  
#fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"ENSEMBL"])==F,]   
aaa<-duplicated(Final_ID_acorss$ENSEMBL)   
Final_ID_acorss<-Final_ID_acorss[!aaa,]
dim(Final_ID_acorss) 
#fix(Final_ID_acorss)
write.table(Final_ID_acorss, file = "0_ID_change.txt",sep="\t", row.names =F, quote = F)
gxj_GO <- read.table("0_ID_change.txt",header=T,comment.char = "",check.names=F,sep="\t")
#fix(gxj_GO)
ego_MF <- enrichGO(gene          = gxj_GO$ENSEMBL,
                   keyType       ="ENSEMBL", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "MF",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,  
                   readable      = TRUE)  
write.table(ego_MF, file = "1_GO_MF_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_BP  <- enrichGO(gene          = gxj_GO$ENSEMBL,
                    keyType       ="ENSEMBL", 
                    OrgDb         = org.Hs.eg.db, 
                    ont           = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,  
                    readable      = TRUE)  
write.table(ego_BP , file = "1_GO_BP_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_CC  <- enrichGO(gene          = gxj_GO$ENSEMBL,
                    keyType       ="ENSEMBL", 
                    OrgDb         = org.Hs.eg.db, 
                    ont           = "CC",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,  
                    readable      = TRUE)  

write.table(ego_CC, file = "1_GO_CC_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_MF<-data.frame(ego_MF)
dim(ego_MF)
#fix(ego_MF)
colnames(ego_MF)
GO_MF<-ego_MF[, c(1,2,8,6)]  
dim(GO_MF)
#fix(GO_MF)
GO_MF$geneID<-str_replace_all (GO_MF$geneID, "/", ",")  
names(GO_MF)<-c("ID", "Term", "Genes", "adj_pval")  
GO_MF$Category<-"MF"  
dim(GO_MF)
#fix(GO_MF)

ego_BP<-data.frame(ego_BP)
dim(ego_BP)
#fix(ego_BP)
colnames(ego_BP)
GO_BP<-ego_BP[, c(1,2,8,6)]  
dim(GO_BP)
#fix(GO_BP)
GO_BP $geneID<-str_replace_all (GO_BP$geneID, "/", ",")  
names(GO_BP)<-c("ID", "Term", "Genes", "adj_pval")  
GO_BP$Category<-"BP"  
dim(GO_BP)
#fix(GO_BP)

ego_CC<-data.frame(ego_CC)
dim(ego_CC)
#fix(ego_CC)
colnames(ego_CC)
GO_CC<-ego_CC[ c(1,2,8,6)]  
dim(GO_CC)
#fix(GO_CC)
GO_CC$geneID<-str_replace_all (GO_CC$geneID, "/", ",")  
names(GO_CC)<-c("ID", "Term", "Genes", "adj_pval")  
GO_CC$Category<-"CC"  
dim(GO_CC)
#fix(GO_CC)
enrichment_MF_BP_CC<-rbind(GO_MF,GO_BP, GO_CC)   
dim(enrichment_MF_BP_CC)
fix(enrichment_MF_BP_CC)
write.table(enrichment_MF_BP_CC, file = "2_enrichment_MF_BP_CC.txt",sep="\t", row.names =F, quote = F)
  enrichment_MF_BP_CC_top10<- read.table("2_enrichment_MF_BP_CC - top10.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(enrichment_MF_BP_CC_top10)
fix(enrichment_MF_BP_CC_top10)
circ<- circle_dat(enrichment_MF_BP_CC_top10, gxj) 
dim(circ)
fix(circ)  
tiff(file="1_GOplot_1.tiff",width =70,height =50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, labels = 2)  
dev.off()
reduced_circ <- reduce_overlap(circ, overlap = 0.75) # Reduce redundant terms with a gene overlap >= 0.75...
tiff(file="1_GOplot_1_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, labels = 3) 
dev.off()
tiff(file="1_GOplot_2.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3) 
dev.off()
tiff(file="1_GOplot_2_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3) 
dev.off()
tiff(file="1_GOplot_3.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)   
dev.off()
tiff(file="1_GOplot_3_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)   
dev.off()
tiff(file="2_GOplot_BP_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'BP'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)   #主要由第一颜色决定,order.by.zscore=F表示P值展示.
dev.off()  
tiff(file="2_GOplot_BP_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'BP'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)  #可以改变颜色
dev.off()
tiff(file="2_GOplot_CC_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'CC'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)   #主要由第一颜色决定,order.by.zscore=F表示P值展示.
dev.off()  
tiff(file="2_GOplot_CC_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'CC'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)  #可以改变颜色
dev.off()
tiff(file="2_GOplot_MF_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'MF'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)   #主要由第一颜色决定,order.by.zscore=F表示P值展示.
dev.off()  
tiff(file="2_GOplot_MF_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'MF'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)  #可以改变颜色
dev.off()
tiff(file="2_GOplot_all_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)
dev.off()
tiff(file="2_GOplot_all_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)+
  theme(axis.text.x = element_text(size=20))
dev.off()
tiff(file="3_GOplot_.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
IDs <- c("GO:0008017", "GO:0008392", "GO:0044282", "GO:0046364","GO:0062023","GO:0060205")  #选择感兴趣的-result进行展示
GOCircle(circ, nsub = IDs, zsc.col = c('yellow', 'black', 'cyan'),lfc.col = c('purple', 'orange'))
dev.off()
tiff(file="3_GOplot_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOCircle(circ, nsub = 6, zsc.col = c('yellow', 'red', 'cyan'),lfc.col = c('blue', 'green'))   
dev.off() 
#fix(gxj)
target_gene<-gxj 
target_process<-unique(enrichment_MF_BP_CC_top10$Term)    
tiff(file="4-chord_all-more.tiff",width =45,height = 50,units ="cm",compression="lzw",bg="white",res=300)
chord <- chord_dat(circ, target_gene, target_process)
GOChord(chord, gene.order = 'logFC', gene.space = 0.15, gene.size = 10)  
dev.off()
tiff(file="4-chord_all-less.tiff",width =45,height =50,units ="cm",compression="lzw",bg="white",res=300)
chord <- chord_dat(circ, target_gene, target_process)
GOChord(chord, limit = c(2, 5), gene.order = 'logFC', gene.space = 0.20, 
        lfc.col=c('red','gray','cyan'),
        gene.size = 6)  
dev.off()    
target_process_2<-c( "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",
                     "small molecule catabolic process",
                     "collagen-containing extracellular matrix")  
tiff(file="4-chord_GO_Okokok.tiff",width =45,height =50,units ="cm",compression="lzw",bg="white",res=300)
chord_1 <- chord_dat(circ, target_gene, target_process_2)
GOChord(chord_1, space = 0.02, 
        lfc.col=c('#FF0088', 'white','green'),
        #ribbon.col=brewer.pal(length(target_process_2), "Set2"),
        ribbon.col=colorRampPalette(c("orange", "#3399FF"))(length(target_process_2)),
        gene.order = 'logFC', lfc.min=-3, lfc.max=3,  
        gene.space = 0.20, gene.size = 6,border.size=1,
        process.label=15)

dev.off()
tiff(file="5-heatmap-GOHeat_Count_all.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
#fix(chord)
#colnames(chord) 
abc<-data.frame (dim(chord))  
chord_no_logFC<-chord[,-abc[2,]] 
GOHeat(chord_no_logFC, nlfc = 0)  #the chord object without logFC column to create the heatmap
dev.off()
tiff(file="5-heatmap-GOHeat_logFC_all.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
GOHeat(chord, nlfc = 1, fill.col = c('red', 'white', 'green'))   #create the heatmap with logFC values and user-defined colour scale
dev.off()
tiff(file="5-heatmap-GOHeat_Count_term.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
abc_1<-data.frame (dim(chord_1))  
chord_no_logFC_1<-chord_1[,-abc_1[2,]] 
GOHeat(chord_no_logFC_1, nlfc = 0)  #the chord object without logFC column to create the heatmap
dev.off()
tiff(file="5-heatmap-GOHeat_logFC_term.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
GOHeat(chord_1, nlfc = 1, fill.col = c('red', 'white', 'green'))   #create the heatmap with logFC values and user-defined colour scale
dev.off()
tiff(file="6-cluster-GOCluster_logFCcluster_all.tiff",width =120,height = 120,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process, clust.by = 'logFC', term.width = 0.5)  cluster
dev.off()
tiff(file="6-cluster-GOCluster_logFCcluster_all.tiff",width =120,height = 120,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process)  cluster
dev.off()
target_process_2
tiff(file="6-cluster-GOCluster_logFCcluster_.tiff",width =70,height = 80,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process_2, clust.by = 'logFC', term.width = 2)   
dev.off()
tiff(file="6-cluster-GOCluster_termcluster_.tiff",width =70,height = 80,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process_2, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
dev.off()
tiff(file="7-veen-GOVenn.tiff",width =40,height = 30,units ="cm",compression="lzw",bg="white",res=300)
l1 <- subset(circ, term == "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor", c(genes,logFC))
l2 <- subset(circ, term == "small molecule catabolic process",  c(genes,logFC))
l3 <- subset(circ, term ==  "collagen-containing extracellular matrix", c(genes,logFC))
GOVenn(l1,l2,l3, label = c('oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor', 'small molecule catabolic process', 'collagen-containing extracellular matrix'))
dev.off()
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")
## try http:// if https:// URLs are not supported
#source("http://bioconductor.org/biocLite.R")
#biocLite("ConsensusClusterPlus")
#browseVignettes("ConsensusClusterPlus")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
rm(list=ls())  
rt_1=read.table("0_GSE14520_OS_valid-riskscore.txt",sep="\t",header=T,check.names=F)              
#fix(rt_1)
dim(rt_1)
rt_2<-rt_1[,c("id", "KPNA2", "FMO3", "SPP1", "KIF2C", "LPCAT1" )]
#fix(rt_2)
dim(rt_2)  
rt_3<-t(rt_2)
#fix(rt_3)
dim(rt_3) 
colnames(rt_3)<- rt_3[1,]
rt_4<-rt_3[-1,]  
#fix(rt_4)
dim(rt_4) 
library(limma)
rt_4=as.matrix(rt_4)
rownames(rt_4)=rt_4[,1]
exp=rt_4[,2:ncol(rt_4)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
#fix(data)
dim(data)  
data=as.matrix(data)
maxK=9  
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="png")

clusterNum=2                 
cluster_2=results[[clusterNum]][["consensusClass"]]
cluster_2<-as.data.frame(cluster_2)
#fix(cluster_2)
dim(cluster_2)
write.table(cluster_2,file="1_GEO-two-group_cluster(PCA).txt",sep="\t",quote=F,col.names = F)
write.table(cluster_2,file="2_GEO-cluster-surv-trait-enrichment.txt",sep="\t",quote=F)
clusterNum=3                 
cluster_3=results[[clusterNum]][["consensusClass"]]
cluster_3<-as.data.frame(cluster_3)
#fix(cluster_3)
dim(cluster_3)
write.table(cluster_3,file="3_GEO-three-group_cluster(PCA).txt",sep="\t",quote=F,col.names = F)
write.table(cluster_3,file="4_GEO-three-group_cluster(surv-trait-enrichment).txt",sep="\t",quote=F)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("ggplot2")
#pca analysis
library(limma)
rt=read.table("0_GSE14520_batch-2.txt",sep="\t",header=T,check.names=F)    
#fix(rt)
dim(rt)  
GEO_1=avereps(rt)
GEO_1=GEO_1[rowMeans(GEO_1)>0.5,]   
#fix(GEO_1)
dim(GEO_1) 
dim(data) 
HCC_sample<-colnames(data) 
GEO_2<-GEO_1[,HCC_sample] 
dim(GEO_2) 
GEO_3<-t(GEO_2) 
#fix(GEO_3)
dim(GEO_3) 
GEO_3.class <- rownames(GEO_3)
GEO_3.pca <- prcomp(GEO_3, scale. = TRUE)                             
abc<-predict(GEO_3.pca)
dim(abc) 
write.table(abc,file="5_GEO-PCA-result_newTab.xls",quote=F,sep="\t")    
library(ggplot2)
cluster_2a=read.table("1_GEO-two-group_cluster(PCA).txt",sep="\t",header=F,check.names=F)     
group=paste0("cluster",as.vector(cluster_2a[,2]))
pcaPredict_2a=predict(GEO_3.pca)
#fix(pcaPredict_2a)
dim(pcaPredict_2a) 
PCA_2a = data.frame(PCA1 = pcaPredict_2a[,1], PCA2 = pcaPredict_2a[,2],group=group)
tiff(file="6_GEO-PCA_two-group.tiff",width =20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data = PCA_2a, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
library(ggplot2)
cluster_3a=read.table("3_GEO-three-group_cluster(PCA).txt",sep="\t",header=F,check.names=F)     
group=paste0("cluster",as.vector(cluster_3a[,2]))
pcaPredict_3a=predict(GEO_3.pca)
#fix(pcaPredict_3a)
dim(pcaPredict_3a) 
PCA_3a = data.frame(PCA1 = pcaPredict_3a[,1], PCA2 = pcaPredict_3a[,2],group=group)
tiff(file="7_GEO-PCA_three-group.tiff",width =20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data = PCA_3a, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
#install.packages("survival")
rm(list=ls())  
library(survival)
cluster_1=read.table("0_GEO-cluster-surv-trait-enrichment.txt",header=T,sep="\t",check.names=F,row.names=1)
##fix(cluster_1)
dim(cluster_1) 
cli_1=read.table("0_GSE14520_OS_valid-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
##fix(cli_1)
dim(cli_1)  
colnames(cli_1)
cli_2<-cli_1[,c("os_time", "os_stat")] 
##fix(cli_2)
dim(cli_2)  
cli_3<-cli_2[row.names(cluster_1),] 
dim(cli_3)
##fix(cli_3) 
aaa <- data.frame(cli_3,cluster_1[match(rownames(cluster_1),rownames(cli_3)),])
dim(aaa)
##fix(aaa)
colnames(aaa)[3]="cluster"
dim(aaa)
##fix(aaa) 
write.table(aaa,file="1_GEO-two-cluster-OS-comb.txt",sep="\t",quote=F)
clusterNum=2    
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }

fit_2 <- survfit(Surv(os_time, os_stat) ~ cluster, data = aaa)
#summary(fit)
tiff(file="2_GEO_cluster-OS-differ(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_2 , 
     lwd=4, 
     col=c("#FF1493","#0066FF"),
     xlab="Time (year)",  
     mark.time=T,
     bg="yellow",
     ylab="Survival rate",
     
     main=paste("GEO-OS (p=", pValue ,")",sep=""))
legend("topright", 
       paste0("cluster",1:clusterNum), 
       lwd=4, 
       bty = "n", 
       col=c("#0066FF","#FF1493"))
dev.off()
library(survival)
library("survminer") 
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit2a <- survfit(Surv(os_time, os_stat) ~ cluster, data = aaa)
tiff(file="2_GEO_cluster-OS-differ(ggsurvplot).tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) #设置图的大小
surPlot=ggsurvplot(fit2a, 
                   data=aaa,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("cluster1", "cluster2"),
                   legend.title="GEO_OS",
                   xlab="os_time(years)",
                   break.os_time.by = 1,
                   risk.table.title="",
                   palette=c("#0066FF","#2E8B57"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
cluster_2=read.table("0_GEO-three-group_cluster(surv-trait-enrichment).txt",header=T,sep="\t",check.names=F, row.names=1)
#fix(cluster_2)
dim(cluster_2) 
cli_1=read.table("0_GSE14520_OS_valid-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1)  
colnames(cli_1)
cli_2<-cli_1[,c("os_time", "os_stat")] 
#fix(cli_2)
dim(cli_2)  
cli_3<-cli_2[row.names(cluster_2),] 
dim(cli_3)
#fix(cli_3) 
bbb <- data.frame(cli_3,cluster_2[match(rownames(cluster_2),rownames(cli_3)),])
dim(bbb)
#fix(bbb)
colnames(bbb)[3]="cluster"
dim(bbb)
#fix(bbb)
write.table(aaa,file="3_GEO-three-group_OS-surv-comb.txt",sep="\t",quote=F)
clusterNum=3    
diff=survdiff(Surv(os_time, os_stat) ~cluster,data = bbb) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }
fit_3 <- survfit(Surv(os_time, os_stat) ~ cluster, data = bbb)
summary(fit_3)
tiff(file="4_GEO_three-group-OS-surv-diff(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_3 , 
     lwd=4, 
     col=c("#FF1493","#0066FF","#006600"),
     xlab="Time (year)",  
     mark.time=T,
     ylab="Survival rate",
     main=paste("GEO-OS (p=", pValue ,")",sep=""))
legend( "bottomleft",    #"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
        paste0("cluster",1:clusterNum), 
        lwd=4, 
        bty = "n", 
        col=c("#0066FF","#FF1493","#006600"))
dev.off()
cluster_1=read.table("0_GEO-cluster-surv-trait-enrichment.txt",header=T,sep="\t",check.names=F,row.names=1)
#fix(cluster_1)
dim(cluster_1)  

cli_1=read.table("0_GSE14520_RFS_valid-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
colnames(cli_1)
cli_2<-cli_1[,c("rfs_time", "rfs_stat")] 
#fix(cli_2)
dim(cli_2) 
cli_3<-cli_2[row.names(cluster_1),] 
dim(cli_3)
#fix(cli_3)
aaa <- data.frame(cli_3,cluster_1[match(rownames(cluster_1),rownames(cli_3)),])
dim(aaa)
#fix(aaa)
colnames(aaa)[3]="cluster"
dim(aaa)
#fix(aaa) 
write.table(aaa,file="5_GEO-tumor-two-subtype_RFS-surv-comb.txt",sep="\t",quote=F)
clusterNum=2    
diff=survdiff(Surv(rfs_time, rfs_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }
fit_2 <- survfit(Surv(rfs_time, rfs_stat) ~ cluster, data = aaa)
#summary(fit)
tiff(file="6_GEO_tumor-two-subtype-RFS-surv-diff(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_2 , 
     lwd=4, 
     col=c("#FF1493","#0066FF"),
     xlab="Time (year)",  
     mark.time=T,
     bg="yellow",
     ylab="Survival rate",
          main=paste("GEO-RFS (p=", pValue ,")",sep=""))
legend("topright", 
       paste0("cluster",1:clusterNum), 
       lwd=4, 
       bty = "n", 
       col=c("#0066FF","#FF1493"))
dev.off()
library(survival)
library("survminer") 
diff=survdiff(Surv(rfs_time, rfs_stat) ~cluster,data = aaa) 
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit2a <- survfit(Surv(rfs_time, rfs_stat) ~ cluster, data = aaa)
tiff(file="6_GEO_tumor-two-subtype-RFS-surv-diff(ggsurvplot).tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
surPlot=ggsurvplot(fit2a, 
                   data=aaa,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("cluster1", "cluster2"),
                   legend.title="GEO_RFS",
                   xlab="rfs_time(years)",
                   break.rfs_time.by = 1,
                   risk.table.title="",
                   palette=c("#0066FF","#2E8B57"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
cluster_2=read.table("0_GEO-three-group_cluster(surv-trait-enrichment).txt",header=T,sep="\t",check.names=F, row.names=1)
#fix(cluster_2)
dim(cluster_2) 
cli_1=read.table("0_GSE14520_RFS_valid-riskscore.txt",header=T,sep="\t",check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1)  
colnames(cli_1)
cli_2<-cli_1[,c("rfs_time", "rfs_stat")] 
#fix(cli_2)
dim(cli_2) 
cli_3<-cli_2[row.names(cluster_2),] 
dim(cli_3)
#fix(cli_3)
bbb <- data.frame(cli_3,cluster_2[match(rownames(cluster_2),rownames(cli_3)),])
dim(bbb)
#fix(bbb)
colnames(bbb)[3]="cluster"
dim(bbb)
##fix(bbb) 
write.table(aaa,file="7_GEO-three-group_RFS-surv-comb.txt",sep="\t",quote=F)
clusterNum=3    
diff=survdiff(Surv(rfs_time, rfs_stat) ~cluster,data = bbb) 
pValue=1-pchisq(diff$chisq,df=clusterNum-1) 
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}  else{ pValue=round(pValue,3)    }
fit_3 <- survfit(Surv(rfs_time, rfs_stat) ~ cluster, data = bbb)
summary(fit_3)
tiff(file="8_GEO_three-group-RFS-surv-diff(plot).tiff",width =10,height = 10,units ="cm",compression="lzw",bg="white",res=300)
plot(fit_3 , 
     lwd=4, 
     col=c("#FF1493","#0066FF","#006600"),
     xlab="Time (year)",  
     mark.time=T,
     ylab="Survival rate",
     main=paste("GEO-RFS (p=", pValue ,")",sep=""))
legend( "bottomleft",    #"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
        paste0("cluster",1:clusterNum), 
        lwd=4, 
        bty = "n", 
        col=c("#0066FF","#FF1493","#006600"))
dev.off()
rm(list=ls())  
GEO_OS_cli<-"0_GEO_OS_model-gene_riskscore_chr.txt"
cluster_1=read.table("0_GEO-cluster-surv-trait-enrichment.txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
#fix(cluster_2)
dim(cluster_2) 
cli_1=read.table(GEO_OS_cli,sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
#fix(cli_2)
dim(cli_2) 
cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
#fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
#fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
#fix(cluster_cli_2) 
write.table(cluster_cli_2, file="1_GEO_OS_cluster（twogroup）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_2"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
GEO_OS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_2","gender", "age","stage")])  
#fix(GEO_OS_clin_rt)
dim(GEO_OS_clin_rt)  
colnames(GEO_OS_clin_rt)
write.table(GEO_OS_clin_rt, file="2_GEO_OS_cluster（twogroup）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
GEO_OS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
#fix(GEO_OS_cluster_rt)
dim(GEO_OS_cluster_rt) 
colnames(GEO_OS_cluster_rt)
write.table(GEO_OS_cluster_rt, file="2_GEO_OS_cluster（twogroup）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
field="cluster_2"  
flag1="1"
flag2="2"
trainFlag=GEO_OS_clin_rt[GEO_OS_clin_rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=GEO_OS_clin_rt[GEO_OS_clin_rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)
#fix(newTable) 
newLabels=c("cluster_2")
for(i in 2:(ncol(GEO_OS_clin_rt)) ){ 
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=fisher.test(tableStat)  
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels 
#fix(newTable)
colnames(GEO_OS_clin_rt)=newLabels
#fix(GEO_OS_clin_rt) 
dim(GEO_OS_clin_rt)
write.table(GEO_OS_clin_rt,file="3_GEO_OS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",quote=F)
#install.packages("pheatmap")
library(pheatmap)
rt=read.table("2_GEO_OS_cluster（twogroup）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("3_GEO_OS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)
Type=Type[order(Type$cluster_2),] 
#fix(Type)
rt_1=rt[,row.names(Type)] 
#fix(rt_1)
GEO_OS_cluster_rt_1=read.table("2_GEO_OS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_OS_cluster_rt_1)
gxj_table=table(GEO_OS_cluster_rt_1[,c("stage","cluster_2")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="4_GEO_OS_clustertwogroup_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("GEO-OS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
GEO_OS_cluster_rt_2=read.table("2_GEO_OS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_OS_cluster_rt_2)
GEO_OS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(GEO_OS_cluster_rt_2)){
  nameStat=colnames(GEO_OS_cluster_rt_2)[i]
  tableStat=table(GEO_OS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  GEO_OS_cluster_rt_2StatOut=rbind(GEO_OS_cluster_rt_2StatOut,tableStatOut)}
gxj_StatOut=GEO_OS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="5_GEO_OS_cluster（twogroup）-num.xls",sep="\t",quote=F,row.names=F)
rm(list=ls())  
subtype
cluster_1=read.table("0_GEO-three-group_cluster(surv-trait-enrichment).txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
#fix(cluster_2)
dim(cluster_2) 
GEO_OS_cli<-"0_GEO_OS_model-gene_riskscore_chr.txt"
cli_1=read.table(GEO_OS_cli,sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
#fix(cli_2)
dim(cli_2) 
cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
#fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
#fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
#fix(cluster_cli_2) 
write.table(cluster_cli_2, file="6_GEO_OS__cluster（three-group）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_3"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
GEO_OS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_3","gender", "age","stage")])  
#fix(GEO_OS_clin_rt)
dim(GEO_OS_clin_rt)  
colnames(GEO_OS_clin_rt)
write.table(GEO_OS_clin_rt, file="7_GEO_OS__cluster（three-group）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
GEO_OS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
#fix(GEO_OS_cluster_rt)
dim(GEO_OS_cluster_rt) 
colnames(GEO_OS_cluster_rt)
write.table(GEO_OS_cluster_rt, file="8_GEO_OS__cluster（three-group）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
#install.packages("pheatmap")
library(pheatmap)
rt=read.table("8_GEO_OS__cluster（three-group）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("7_GEO_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)
Type=Type[order(Type$cluster_3),] 
#fix(Type)
rt_1=rt[,row.names(Type)] 
#fix(rt_1)
GEO_OS_cluster_rt_1=read.table("7_GEO_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_OS_cluster_rt_1)
gxj_table=table(GEO_OS_cluster_rt_1[,c("stage","cluster_3")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="9_GEO_OS_cluster（three-group）_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("GEO-OS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
GEO_OS_cluster_rt_2=read.table("7_GEO_OS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_OS_cluster_rt_2)
GEO_OS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(GEO_OS_cluster_rt_2)){
  nameStat=colnames(GEO_OS_cluster_rt_2)[i]
  tableStat=table(GEO_OS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  GEO_OS_cluster_rt_2StatOut=rbind(GEO_OS_cluster_rt_2StatOut,tableStatOut)}

gxj_StatOut=GEO_OS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="10_GEO_OS__cluster（three-group）-num.xls",sep="\t",quote=F,row.names=F)
GEO_RFS_cli<-"0_GEO_RFS_model-gene_riskscore_chr.txt"
cluster_1=read.table("0_GEO-cluster-surv-trait-enrichment.txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
#fix(cluster_2)
dim(cluster_2) 
cli_1=read.table(GEO_RFS_cli,sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
#fix(cli_2)
dim(cli_2) 
cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
#fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
#fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
#fix(cluster_cli_2) 
write.table(cluster_cli_2, file="1_GEO_RFS_cluster（twogroup）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_2"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "RFS_time"     "RFS_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
GEO_RFS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_2","gender", "age","stage")])  
#fix(GEO_RFS_clin_rt)
dim(GEO_RFS_clin_rt)  
colnames(GEO_RFS_clin_rt)
write.table(GEO_RFS_clin_rt, file="2_GEO_RFS_cluster（twogroup）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
GEO_RFS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
#fix(GEO_RFS_cluster_rt)
dim(GEO_RFS_cluster_rt) 
colnames(GEO_RFS_cluster_rt)
write.table(GEO_RFS_cluster_rt, file="2_GEO_RFS_cluster（twogroup）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
field="cluster_2"  
flag1="1"
flag2="2"
trainFlag=GEO_RFS_clin_rt[GEO_RFS_clin_rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=GEO_RFS_clin_rt[GEO_RFS_clin_rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)
newLabels=c("cluster_2")
for(i in 2:(ncol(GEO_RFS_clin_rt)) ){ 
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=fisher.test(tableStat)  
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels 
colnames(GEO_RFS_clin_rt)=newLabels
dim(GEO_RFS_clin_rt)
write.table(GEO_RFS_clin_rt,file="3_GEO_RFS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",quote=F)
library(pheatmap)

rt=read.table("2_GEO_RFS_cluster（twogroup）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("3_GEO_RFS_cluster（twogroup）_chr-cli-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)
Type=Type[order(Type$cluster_2),] 
#fix(Type)
rt_1=rt[,row.names(Type)] 
#fix(rt_1)
GEO_RFS_cluster_rt_1=read.table("2_GEO_RFS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_RFS_cluster_rt_1)
gxj_table=table(GEO_RFS_cluster_rt_1[,c("stage","cluster_2")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="4_GEO_RFS_clustertwogroup_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("GEO-RFS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
GEO_RFS_cluster_rt_2=read.table("2_GEO_RFS_cluster（twogroup）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_RFS_cluster_rt_2)
GEO_RFS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(GEO_RFS_cluster_rt_2)){
  nameStat=colnames(GEO_RFS_cluster_rt_2)[i]
  tableStat=table(GEO_RFS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  GEO_RFS_cluster_rt_2StatOut=rbind(GEO_RFS_cluster_rt_2StatOut,tableStatOut)}

gxj_StatOut=GEO_RFS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="5_GEO_RFS_cluster（twogroup）-num.xls",sep="\t",quote=F,row.names=F)
rm(list=ls())  
subtype
cluster_1=read.table("0_GEO-three-group_cluster(surv-trait-enrichment).txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(cluster_1)
dim(cluster_1)  
ID<-row.names(cluster_1)  
cluster_2<-data.frame(ID, cluster_1, check.names=FALSE)
#fix(cluster_2)
dim(cluster_2) 
GEO_RFS_cli<-"0_GEO_RFS_model-gene_riskscore_chr.txt"
cli_1=read.table(GEO_RFS_cli,sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
ID<-row.names(cli_1)
cli_2<-data.frame(ID, cli_1, check.names=FALSE)
#fix(cli_2)
dim(cli_2) 
cluster_cli_1<-merge(cluster_2, cli_2,  by.x="ID", by.y="ID")   
dim(cluster_cli_1)
#fix(cluster_cli_1)
row.names(cluster_cli_1)<-cluster_cli_1[,1] 
#fix(cluster_cli_1) 
cluster_cli_2<-cluster_cli_1[,-1] 
#fix(cluster_cli_2) 
write.table(cluster_cli_2, file="6_GEO_RFS__cluster（three-group）_chr_model-gene-exp.txt", quote=FALSE,  sep="\t")   
colnames(cluster_cli_2) 
#[1]      "cluster_3"   "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "RFS_time"     "RFS_stat"    
#[11] "KPNA2"       "FMO3"        "SPP1"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk"   
GEO_RFS_clin_rt<-data.frame(cluster_cli_2[,c( "cluster_3","gender", "age","stage")])  
#fix(GEO_RFS_clin_rt)
dim(GEO_RFS_clin_rt)  
colnames(GEO_RFS_clin_rt)
write.table(GEO_RFS_clin_rt, file="7_GEO_RFS__cluster（three-group）_chr-cli-no-sig.txt", quote=FALSE,  sep="\t") 
GEO_RFS_cluster_rt<-data.frame(cluster_cli_2[,c("KPNA2","FMO3","SPP1","KIF2C","LPCAT1")])  
#fix(GEO_RFS_cluster_rt)
dim(GEO_RFS_cluster_rt) 
colnames(GEO_RFS_cluster_rt)
write.table(GEO_RFS_cluster_rt, file="8_GEO_RFS__cluster（three-group）_model-gene-exp.txt", quote=FALSE,  sep="\t") 
#install.packages("pheatmap")
library(pheatmap)
rt=read.table("8_GEO_RFS__cluster（three-group）_model-gene-exp.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt<-t(rt)
#fix(rt)
dim(rt) 
Type=read.table("7_GEO_RFS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(Type)
dim(Type)
Type=Type[order(Type$cluster_3),] 
#fix(Type)
rt_1=rt[,row.names(Type)] 
#fix(rt_1)
GEO_RFS_cluster_rt_1=read.table("7_GEO_RFS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_RFS_cluster_rt_1)
gxj_table=table(GEO_RFS_cluster_rt_1[,c("stage","cluster_3")]) 
pStat=fisher.test(gxj_table)  
pStat$p.value
tiff(file="9_GEO_RFS_cluster（three-group）_-exp_cli-cor-heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt_1, annotation=Type, 
         main=paste0("GEO-RFS-stage: p=",round (pStat$p.value,digits=7)),
         color = colorRampPalette(c("#6A5ACD", "white", "#FF0066"))(50),
         cluster_cols =F, 
         fontsize=8, 
         fontsize_row=8,
         scale="row", 
         show_colnames=F,
         fontsize_col=3)
dev.off()
GEO_RFS_cluster_rt_2=read.table("7_GEO_RFS__cluster（three-group）_chr-cli-no-sig.txt",sep="\t",header=T,row.names=1,check.names=F)
#fix(GEO_RFS_cluster_rt_2)
GEO_RFS_cluster_rt_2StatOut=data.frame()
for(i in 1:ncol(GEO_RFS_cluster_rt_2)){
  nameStat=colnames(GEO_RFS_cluster_rt_2)[i]
  tableStat=table(GEO_RFS_cluster_rt_2[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  GEO_RFS_cluster_rt_2StatOut=rbind(GEO_RFS_cluster_rt_2StatOut,tableStatOut)}
gxj_StatOut=GEO_RFS_cluster_rt_2StatOut[,-3]
write.table(gxj_StatOut,file="10_GEO_RFS__cluster（three-group）-num.xls",sep="\t",quote=F,row.names=F)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install('org.Hs.eg.db')
#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationDbi")
#biocLite("org.Hs.eg.db")
#install.packages('ggplot2')
#install.packages('pheatmap')
#install.packages('colorspace')
#install.packages('labeling')
rm(list=ls())  
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0_GEO-cluster-surv-trait-enrichment.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1)  
com = subset(group_1, cluster_2 =="1"|cluster_2  =="2",  select=cluster_2 )  
dim(com)   
#fix(com)   
colnames(com)[1]<-c("group")
head(com,5)
gene_1<-read.table("0_GSE14520_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
#fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
#fix(eSet)   。
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
#fix(group_list) 
rownames(group_list) = rownames(com)        
colnames(group_list)[1] = "description"   
dim(group_list) 
#fix(group_list) 
design = model.matrix(~ description + 0, group_list)   
#fix(design)
colnames(design) = levels(as.factor(c("cluster1","cluster2")))   
factor(design) 
head(design)  
cont.matrix = makeContrasts(cluster2-cluster1, levels=design)   
fdr = 0.05
logFC = 1
fit = lmFit(eSet, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="BH", sort.by="logFC", resort.by = "P" ,n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "Pos" 
tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Neg"  
tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "NO"
colnames(tT)[ncol(tT)] = "expression"
dim(tT)
write.table(tT,file="1_GEO_cluster2-cluster1_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
model_gene_1<-read.table("0_lasso_(geneCoef).txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(model_gene_1)  
#fix(model_gene_1)  
model_gene<-row.names(model_gene_1) 
model_gene
diff_gene<-row.names(deg) 
diff_gene[1:10]  
aa_1<-setdiff(diff_gene, model_gene)  
dim(deg)  
deg_1<-deg[aa_1,] 
dim(deg_1)    
write.table(deg_1,file="2_GEO_cluster2-cluster1_diff-gene.xls",sep="\t",quote=F)
write.table(deg_1,file="2_GEO_cluster2-cluster1_diff-gene.txt",sep="\t",quote=F)
group_reg <- factor(deg_1$expression)
labels=levels(group_reg)
gxj<-summary(group_reg)
gxj
#fix(tT)
dim(tT)
aa_2<-setdiff(row.names(tT), model_gene) 
tT_1<-tT[aa_2,]
#fix(tT_1)
dim(tT_1)
temp1 = tT_1[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("Pos","Neg","NO"), order=T)
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+ 
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080" ))+  
  labs(x=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]),y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  #xlim(-5,5)+  
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = "#FFE1FF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "#FFE1FF"))
tiff(file="3_GEO_cluster-volcano_style1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_GEO_cluster-volcano_style2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+  
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080"))
dev.off()
temp2 = data.frame(log10(rowMeans(eSet))) 
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT_1),rownames(temp2)),] 
temp2 = data.frame(temp2,tT_1[,c("logFC","expression")])
colnames(temp2)=c("log10Exp","logFC","expression")
temp2$expression=factor(temp2$expression, levels=c("Pos","Neg","NO"), order=T)
P_MA=ggplot(temp2,aes(x=temp2$log10Exp,y=temp2$logFC))+
  geom_point(alpha=0.8, size=6, 
             aes(color=temp2$expression))+
  scale_color_manual(values =c("Pos" = "#FF4500", "NO" = "gray","Neg" = "#483D8B"))+
  labs(x="log10Exp",y=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+
  # xlim(0.75,1.15)+  
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "black"))
tiff(file="4_GEO_cluster-MA.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)  
print(P_MA)
dev.off()
dim(eSet) 
#fix(eSet)
DEG_list = rownames(deg_1)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),] 
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_GEO_cluster-all-gene-heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset, color = colorRampPalette(c('blue2',  'white',  'red'))(10000),
                   cluster_rows = T, cluster_cols=T,                  
                   scale ="row" , 
                   show_rownames=F,show_colnames=F, 
                   treeheight_row = 0,  
                   clustering_distance_rows = "correlation",  
                   clustering_distance_cols = "euclidean",   
                   border_color = "grey",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()
target_gene<-row.names(deg_1)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_GEO_cluster-target-gene-heatmap_(TOP10-remove-hub).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset_2, color = colorRampPalette(c('#000066',  'white',  '#FF0000'))(60), 
                   cluster_rows = T, #clustering_distance_rows = "correlation", 
                   cluster_cols=T,   clustering_distance_cols = "euclidean",          
                   scale ="row" , 
                   show_rownames=T,  
                   show_colnames= F, 
                   treeheight_row = 0,  
                   border_color = "White",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()

#KPNA2; FMO3; SPP1; KIF2C; LPCAT1
model_gene_1<-read.table("0_lasso_(geneCoef).txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(model_gene_1)  
#fix(model_gene_1)  
deg_eset_2<-eSet[row.names(model_gene_1),]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_GEO_cluster-model-heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset_2, color = colorRampPalette(c('#000066',  'white',  '#FF0000'))(60), 
                   cluster_rows = T, #clustering_distance_rows = "correlation", 
                   cluster_cols=T,   clustering_distance_cols = "euclidean",          
                   scale ="row" , 
                   show_rownames=T,  
                   show_colnames= F, 
                   treeheight_row = 0,  
                   border_color = "White",
                   #display_numbers = TRUE, number_color = "blue",  
                   clustering_method = "average")
dev.off()
rm(list=ls())  

library("colorspace")
library("stringi")
library("rlang")
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
library("pathview")
library("GOplot")
clx_kegg1="9_FC_map"   
dir.create(paste(workdir,clx_kegg1,sep="/"))  
clx_kegg2="10_Pvalue_map"   
dir.create(paste(workdir,clx_kegg2,sep="/"))  
gxj <- read.table("0_GEO_cluster2-cluster1_diff-gene.txt",header=T,comment.char = "",check.names=F,sep="\t", row.names=1 ) 
dim(gxj)
#fix(gxj)
colnames(gxj)
gxj$GeneSymbol<-row.names(gxj)  
transfered_ID <- bitr(gxj$GeneSymbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
transfered_ID<-data.frame(transfered_ID) 
#fix(transfered_ID)
Final_ID_acorss<- merge(gxj, transfered_ID,  by.x="GeneSymbol", by.y = "SYMBOL") 
dim(Final_ID_acorss)  
#fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"ENTREZID"])==F,]   
dim(Final_ID_acorss)  
aaa<-duplicated(Final_ID_acorss$ENTREZID)   
Final_ID_acorss<-Final_ID_acorss[!aaa,]
dim(Final_ID_acorss) 
#fix(Final_ID_acorss)
write.table(Final_ID_acorss, file = "01_GEO_ID_change.txt",sep="\t", row.names =F, quote = F)
gxj <- read.table("01_GEO_ID_change.txt",header=T,comment.char = "",check.names=F,sep="\t")
#fix(gxj)
DEG_list <- gxj$ENTREZID
kk <- enrichGO(gene = DEG_list,
               keyType ="ENTREZID",  
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="01_GEO_all_GO.txt",sep="\t",quote=F,row.names = F)
pdf(file="01_GEO_GO_barplot.pdf",width = 30,height = 15)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="01_GEO_GO_bubble.pdf",width =30,height = 15)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
ego_MF <- enrichGO(gene          = DEG_list ,
                   keyType       ="ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "MF",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,  
                   readable      = TRUE)  
write.table(ego_MF, file = "1_GO_MF_enrichment.txt",sep="\t", row.names =F, quote = F)
#ego_results<-summary(ego)
tiff(file="2_GEO_GO_MF_barplot_top10.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GEO_GO_MF_barplot_top20.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GEO_GO_MF_dotplot_top10.tiff",width = 100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=10)
dev.off()
tiff(file="3_GEO_GO_MF_dotplot_top20.tiff",width =100,height = 50,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=20)
dev.off()
tiff(file="4_GEO_GO_MF_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF)
dev.off()
tiff(file="4_GEO_GO_MF_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GEO_GO_MF_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GEO_GO_MF_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_MF_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GEO_GO_MF_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_MF_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GEO_GO_MF_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_MF,showCategory = 20, foldChange =NULL)
dev.off()
tiff(file="7_GEO_GO_MF_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_MF, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
ego_BP <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "BP",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.table(ego_BP, file = "1_GEO_GO_BP_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GEO_GO_BP_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GEO_GO_BP_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GEO_GO_BP_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=10)
dev.off()
tiff(file="3_GEO_GO_BP_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=20)
dev.off()
tiff(file="4_GEO_GO_BP_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP)
dev.off()
tiff(file="4_GEO_GO_BP_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GEO_GO_BP_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GEO_GO_BP_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_BP_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GEO_GO_BP_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_BP_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GEO_GO_BP_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_BP,showCategory = 20, foldChange =NULL)
dev.off()
tiff(file="7_GEO_GO_BP_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_BP, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
ego_CC <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "CC",  
                   pAdjustMethod = "BH",
                   #universe = names(geneList), 
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE) 
write.table(ego_CC, file = "1_GO_CC_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GEO_GO_CC_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GEO_CC_CC_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GEO_GO_CC_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=10)
dev.off()
tiff(file="3_GEO_GO_CC_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=20)
dev.off()
tiff(file="4_GEO_GO_CC_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC)
dev.off()
tiff(file="4_GEO_GO_CC_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GEO_GO_CC_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GEO_GO_CC_cnetplot_style1-circle-with-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
#cnetplot(ego_CC, showCategory = 5, categorySize="pvalue", foldChange=geneList)
cnetplot(ego_CC, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_CC_cnetplot_style2-circle-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5,  foldChange = NULL, layout = "kk", circular =T, colorEdge = T, node_label = "none")
dev.off()
tiff(file="5_GEO_GO_CC_cnetplot_style3-tree-with-note）.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "all")
dev.off()
tiff(file="5_GEO_GO_CC_cnetplot_style4tres-without-note.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC, showCategory = 5, foldChange = NULL, layout = "kk", circular =F, colorEdge = T, node_label = "none")
dev.off()
tiff(file="6_GEO_GO_CC_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_CC,showCategory = 20, foldChange = NULL)
dev.off()
tiff(file="7_GEO_GO_CC_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_CC, showCategory = 30, color = "p.adjust", layout = "kk")
dev.off()
library("pathview")
kegg <- enrichKEGG(gene         = gxj$ENTREZID,
                   organism     = 'hsa',   
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
write.table(kegg, file = "8_GEO_KEGG_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="9_GEO_KEGG_barplot.tiff",width = 35,height = 22,units ="cm",compression="lzw",bg="white",res=300)
barplot(kegg, drop = T, showCategory = 20,color = "pvalue")
dev.off()

tiff(file="10_GEO_KEGG_dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)  
dotplot(kegg,title="Enrichment KEGG_dot")
dev.off()
setwd(paste(workdir,clx_kegg1,sep="/")) 
map_ids = as.matrix(gxj['logFC'])
for (pathway_id in kegg$ID){
  map = pathview(gene.data  = map_ids,
                 pathway.id = pathway_id,
                 species    = "hsa", kegg.native = TRUE)}
setwd(paste(workdir,clx_kegg2,sep="/")) 
grd<-gxj$P.Value   
names(grd)<-gxj$ENTREZID 
for(i in kegg$ID){
  pv.out <- pathview(gene.data = -log10(grd),  
                     pathway.id = i, 
                     species = "hsa", 
                     out.suffix = "pathview",
                     limit=list(gene=5, cpd=5)  
  )}
rm(list=ls())  

library("colorspace")
library("stringi")
library("rlang")
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOplot")
library("stringr")
gxj_1<- read.table("0_GEO_cluster2-cluster1_diff-gene.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(gxj_1)
fix(gxj_1)
gxj<-data.frame(ID=row.names(gxj_1), logFC=gxj_1[,c(1)]) 
dim(gxj)
#fix(gxj)
transfered_ID <- bitr(gxj$ID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL" ), OrgDb="org.Hs.eg.db")
transfered_ID<-data.frame(transfered_ID) 
##fix(transfered_ID)
Final_ID_acorss<- merge(gxj, transfered_ID,  by.x="ID", by.y = "SYMBOL") 
dim(Final_ID_acorss)  
##fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"ENSEMBL"])==F,]   
aaa<-duplicated(Final_ID_acorss$ENSEMBL)   
Final_ID_acorss<-Final_ID_acorss[!aaa,]
dim(Final_ID_acorss) 
##fix(Final_ID_acorss)
write.table(Final_ID_acorss, file = "01_ID_change.txt",sep="\t", row.names =F, quote = F)
gxj_GO <- read.table("01_ID_change.txt",header=T,comment.char = "",check.names=F,sep="\t")
#fix(gxj_GO)
ego_MF <- enrichGO(gene          = gxj_GO$ENSEMBL,
                   keyType       ="ENSEMBL", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "MF",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,  
                   readable      = TRUE)  
write.table(ego_MF, file = "1_GO_MF_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_BP  <- enrichGO(gene          = gxj_GO$ENSEMBL,
                    keyType       ="ENSEMBL", 
                    OrgDb         = org.Hs.eg.db, 
                    ont           = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,  
                    readable      = TRUE)  
write.table(ego_BP , file = "1_GO_BP_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_CC  <- enrichGO(gene          = gxj_GO$ENSEMBL,
                    keyType       ="ENSEMBL", 
                    OrgDb         = org.Hs.eg.db, 
                    ont           = "CC",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,  
                    readable      = TRUE)  
write.table(ego_CC, file = "1_GO_CC_enrichment.txt",sep="\t", row.names =F, quote = F)
ego_MF<-data.frame(ego_MF)
dim(ego_MF)
#fix(ego_MF)
colnames(ego_MF)
GO_MF<-ego_MF[, c(1,2,8,6)]  
dim(GO_MF)
#fix(GO_MF)
GO_MF$geneID<-str_replace_all (GO_MF$geneID, "/", ",")  
names(GO_MF)<-c("ID", "Term", "Genes", "adj_pval")  
GO_MF$Category<-"MF"  
dim(GO_MF)
#fix(GO_MF)
ego_BP<-data.frame(ego_BP)
dim(ego_BP)
#fix(ego_BP)
colnames(ego_BP)
GO_BP<-ego_BP[, c(1,2,8,6)]  
dim(GO_BP)
#fix(GO_BP)
GO_BP $geneID<-str_replace_all (GO_BP$geneID, "/", ",")  
names(GO_BP)<-c("ID", "Term", "Genes", "adj_pval")  
GO_BP$Category<-"BP"  
dim(GO_BP)
#fix(GO_BP)
ego_CC<-data.frame(ego_CC)
dim(ego_CC)
#fix(ego_CC)
colnames(ego_CC)
GO_CC<-ego_CC[ c(1,2,8,6)]  
dim(GO_CC)
#fix(GO_CC)
GO_CC$geneID<-str_replace_all (GO_CC$geneID, "/", ",")  
names(GO_CC)<-c("ID", "Term", "Genes", "adj_pval")  
GO_CC$Category<-"CC"  
dim(GO_CC)
#fix(GO_CC)
enrichment_MF_BP_CC_1<-rbind(GO_MF,GO_BP, GO_CC)   
dim(enrichment_MF_BP_CC_1)
fix(enrichment_MF_BP_CC_1)
write.table(enrichment_MF_BP_CC_1, file = "2_enrichment_MF_BP_CC(all_inf).xls",sep="\t", row.names =F, quote = F)
enrichment_MF_BP_CC_2<-rbind(GO_MF[1:5,],GO_BP[1:5,], GO_CC[1:5,])   
dim(enrichment_MF_BP_CC_2)
fix(enrichment_MF_BP_CC_2)
write.table(enrichment_MF_BP_CC_2, file = "2_enrichment_MF_BP_CC（top5）.txt",sep="\t", row.names =F, quote = F)
  enrichment_MF_BP_CC_top5<- read.table("2_enrichment_MF_BP_CC（top5）.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(enrichment_MF_BP_CC_top5)
fix(enrichment_MF_BP_CC_top5)
circ<- circle_dat(enrichment_MF_BP_CC_top5, gxj) 
dim(circ)
fix(circ)  
tiff(file="1_GOplot_1.tiff",width =70,height =50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, labels = 2)  
dev.off()
reduced_circ <- reduce_overlap(circ, overlap = 0.75) # Reduce redundant terms with a gene overlap >= 0.75...
tiff(file="1_GOplot_1_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, labels = 3) 
dev.off()
tiff(file="1_GOplot_2.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3) 
dev.off()
tiff(file="1_GOplot_2_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3) 
dev.off()
tiff(file="1_GOplot_3.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)   
dev.off()
tiff(file="1_GOplot_3_small.tiff",width =70,height = 50,units ="cm",compression="lzw",bg="white",res=300)
GOBubble(reduced_circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)   
dev.off()
tiff(file="2_GOplot_BP_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'BP'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)   
dev.off()  
tiff(file="2_GOplot_BP_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'BP'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T) 
dev.off()
tiff(file="2_GOplot_CC_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'CC'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)  
dev.off()  
tiff(file="2_GOplot_CC_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'CC'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T) 
dev.off()
tiff(file="2_GOplot_MF_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'MF'),title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)  
dev.off()  
tiff(file="2_GOplot_MF_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(subset(circ, category == 'MF'),title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)  
dev.off()
tiff(file="2_GOplot_all_1.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('green', 'white','blue'),order.by.zscore=F)
dev.off()
tiff(file="2_GOplot_all_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('red', 'white', 'blue'),order.by.zscore=T)+
  theme(axis.text.x = element_text(size=20))
dev.off()
tiff(file="3_GOplot_.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
IDs <- c("GO:0062023", "GO:0062023", "GO:0048037", "GO:0019318") 
GOCircle(circ, nsub = IDs, zsc.col = c('yellow', 'black', 'cyan'),lfc.col = c('purple', 'orange'))
dev.off()
tiff(file="3_GOplot_2.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
GOCircle(circ, nsub = 6, zsc.col = c('yellow', 'red', 'cyan'),lfc.col = c('blue', 'green'))   
dev.off() 
#fix(gxj)
target_gene<-gxj 
target_process<-unique(enrichment_MF_BP_CC_top5$Term)    
tiff(file="4-chord_all-more.tiff",width =45,height = 50,units ="cm",compression="lzw",bg="white",res=300)
chord <- chord_dat(circ, target_gene, target_process)
GOChord(chord, gene.order = 'logFC', gene.space = 0.15, gene.size = 10)  
dev.off()
tiff(file="4-chord_all-less.tiff",width =45,height =50,units ="cm",compression="lzw",bg="white",res=300)
chord <- chord_dat(circ, target_gene, target_process)
GOChord(chord, limit = c(2, 5), gene.order = 'logFC', gene.space = 0.20, 
        lfc.col=c('red','gray','cyan'),
        gene.size = 6)  
dev.off()    
target_process_2<-c( "endoplasmic reticulum lumen",
                     "collagen-containing extracellular matrix",
                     "integral component of endoplasmic reticulum membrane")  
tiff(file="4-chord_GO_Okokok.tiff",width =45,height =50,units ="cm",compression="lzw",bg="white",res=300)
chord_1 <- chord_dat(circ, target_gene, target_process_2)
GOChord(chord_1, space = 0.02, 
        lfc.col=c('#FF0088', 'white','green'),
        #ribbon.col=brewer.pal(length(target_process_2), "Set2"),
        ribbon.col=colorRampPalette(c("orange", "#3399FF"))(length(target_process_2)),
        gene.order = 'logFC', lfc.min=-3, lfc.max=3,  
        gene.space = 0.20, gene.size = 6,border.size=1,
        process.label=15)
dev.off()
tiff(file="5-heatmap-GOHeat_Count_all.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
abc<-data.frame (dim(chord)) 
chord_no_logFC<-chord[,-abc[2,]] 
GOHeat(chord_no_logFC, nlfc = 0)  #the chord object without logFC column to create the heatmap
dev.off()
tiff(file="5-heatmap-GOHeat_logFC_all.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
GOHeat(chord, nlfc = 1, fill.col = c('red', 'white', 'green'))   #create the heatmap with logFC values and user-defined colour scale
dev.off()
tiff(file="5-heatmap-GOHeat_Count_term.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
abc_1<-data.frame (dim(chord_1))  
chord_no_logFC_1<-chord_1[,-abc_1[2,]] 
GOHeat(chord_no_logFC_1, nlfc = 0)  #the chord object without logFC column to create the heatmap
dev.off()
tiff(file="5-heatmap-GOHeat_logFC_term.tiff",width =120,height = 100,units ="cm",compression="lzw",bg="white",res=300)
GOHeat(chord_1, nlfc = 1, fill.col = c('red', 'white', 'green'))  
dev.off()
tiff(file="6-cluster-GOCluster_logFCcluster_all.tiff",width =120,height = 120,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process, clust.by = 'logFC', term.width = 0.5) 
dev.off()
tiff(file="6-cluster-GOCluster_logFCcluster_all.tiff",width =120,height = 120,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process)  
dev.off()
target_process_2
tiff(file="6-cluster-GOCluster_logFCcluster_.tiff",width =70,height = 80,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process_2, clust.by = 'logFC', term.width = 2)   
dev.off()
tiff(file="6-cluster-GOCluster_termcluster_.tiff",width =70,height = 80,units ="cm",compression="lzw",bg="white",res=300)
GOCluster(circ, target_process_2, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
dev.off()
tiff(file="7-veen-GOVenn.tiff",width =40,height = 30,units ="cm",compression="lzw",bg="white",res=300)
l1 <- subset(circ, term == "collagen-containing extracellular matrix", c(genes,logFC))
l2 <- subset(circ, term == "integral component of endoplasmic reticulum membrane",  c(genes,logFC))
l3 <- subset(circ, term ==  "monooxygenase activity", c(genes,logFC))
GOVenn(l1,l2,l3, label = c("collagen-containing extracellular matrix", 
                           "integral component of endoplasmic reticulum membrane",
                           "monooxygenase activity"))
dev.off()
gxj <- read.table("TCGA-GEO-veen62.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(gxj)
##fix(gxj)
colnames(gxj)[c(1)]=c("GeneSymbol") 
#fix(gxj)
transfered_ID <- bitr(gxj$GeneSymbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
transfered_ID<-data.frame(transfered_ID) 
#fix(transfered_ID)
Final_ID_acorss<- merge(gxj, transfered_ID,  by.x="GeneSymbol", by.y = "SYMBOL") 
dim(Final_ID_acorss)  
#fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"ENTREZID"])==F,]   
dim(Final_ID_acorss)  
aaa<-duplicated(Final_ID_acorss$ENTREZID)   
Final_ID_acorss<-Final_ID_acorss[!aaa,]
dim(Final_ID_acorss) 
#fix(Final_ID_acorss)
write.table(Final_ID_acorss, file = "01_ID_change.txt",sep="\t", row.names =F, quote = F)
gxj <- read.table("01_ID_change.txt",header=T,comment.char = "",check.names=F,sep="\t")
#fix(gxj)
DEG_list <- gxj$ENTREZID
kk <- enrichGO(gene = DEG_list ,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="all_GO.txt",sep="\t",quote=F,row.names = F)
pdf(file="all_barplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="all_bubble.pdf",width = 10,height = 8)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
ego_MF <- enrichGO(gene          = DEG_list ,
                   keyType       ="ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "MF",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,  
                 readable      = TRUE)  
write.table(ego_MF, file = "1_GO_MF_enrichment.txt",sep="\t", row.names =F, quote = F)
#ego_results<-summary(ego)
tiff(file="2_GO_MF_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GO_MF_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_MF, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_MF_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=10)
dev.off()
tiff(file="3_GO_MF_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_MF,showCategory=20)
dev.off()
tiff(file="4_GO_MF_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF)
dev.off()
tiff(file="4_GO_MF_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_MF_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_MF_cnetplot_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF,   circular = T, colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_MF_cnetplot_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF,   circular = T , colorEdge = T, node_label=F)
dev.off()
tiff(file="5_GO_MF_cnetplot_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF,   circular = F , colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_MF_cnetplot_style4.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_MF,   circular = F , colorEdge = T, node_label=F)
dev.off()
#data(geneList)
#head(geneList)
tiff(file="6_GO_MF_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_MF)
dev.off()
tiff(file="7_GO_MF_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_MF)
dev.off()
ego_BP <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "BP",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
write.table(ego_BP, file = "1_GO_BP_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GO_BP_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_GO_BP_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_BP, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_BP_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=10)
dev.off()
tiff(file="3_GO_BP_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_BP,showCategory=20)
dev.off()
tiff(file="4_GO_BP_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP)
dev.off()
tiff(file="4_GO_BP_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_BP_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_BP_cnetplot_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP,   circular = T, colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_BP_cnetplot_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP,   circular = T , colorEdge = T, node_label=F)
dev.off()
tiff(file="5_GO_BP_cnetplot_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP,   circular = F , colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_BP_cnetplot_style4.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_BP,   circular = F , colorEdge = T, node_label=F)
dev.off()
tiff(file="6_GO_BP_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_BP)
dev.off()
tiff(file="7_GO_BP_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_BP)
dev.off()
ego_CC <- enrichGO(gene          = DEG_list ,
                   keyType       = "ENTREZID", 
                   OrgDb         = org.Hs.eg.db, 
                   ont           = "CC",  
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE) 
write.table(ego_CC, file = "1_GO_CC_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="2_GO_CC_barplot_top10.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=10, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="2_CC_CC_barplot_top20.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(ego_CC, drop = F, showCategory=20, x = "GeneRatio",col ="p.adjust", border ="red")
dev.off()
tiff(file="3_GO_CC_dotplot_top10.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=10)
dev.off()
tiff(file="3_GO_CC_dotplot_top20.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(ego_CC,showCategory=20)
dev.off()
tiff(file="4_GO_CC_topgo_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC)
dev.off()
tiff(file="4_GO_CC_topgo_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T) 
dev.off()
tiff(file="4_GO_CC_topgo_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = F,useInfo=F)
dev.off()
tiff(file="5_GO_CC_cnetplot_style1.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC,   circular = T, colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_CC_cnetplot_style2.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC,   circular = T , colorEdge = T, node_label=F)
dev.off()
tiff(file="5_GO_CC_cnetplot_style3.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC,   circular = F , colorEdge = T, node_label=T)
dev.off()
tiff(file="5_GO_CC_cnetplot_style4.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
cnetplot(ego_CC,   circular = F , colorEdge = T, node_label=F)
dev.off()
tiff(file="6_GO_CC_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_CC)
dev.off()
tiff(file="7_GO_CC_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_CC)
dev.off()
library("pathview")
kegg <- enrichKEGG(gene         = gxj$ENTREZID,
                   organism     = 'hsa',   
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
write.table(kegg, file = "8_KEGG_enrichment.txt",sep="\t", row.names =F, quote = F)
tiff(file="9_KEGG_barplot.tiff",width = 35,height = 22,units ="cm",compression="lzw",bg="white",res=300)
barplot(kegg, drop = T, showCategory = 20,color = "pvalue")
dev.off()
tiff(file="10_KEGG_dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)  
dotplot(kegg,title="Enrichment KEGG_dot")
dev.off()








#(16) immune cell infiltration analysis#######################################################################
rm(list=ls())  
options(stringsAsFactors = F)
source("0_CIBERSORT(no change).R")
group_1<-read.table("0_TCGA_training-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1) 
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
#fix(com)  
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_TCGA-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_TCGA-LIHC_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
#fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
#fix(eSet)   
library(limma)
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,] 
dim(exp_data1)    
#fix(exp_data1)  #
write.table(exp_data1, file="2_CIBERSORT_exp_TCGA.txt", quote=FALSE,  sep="\t") 
LM22.file <- "0_LM22.txt"
exp.file <- "2_CIBERSORT_exp_TCGA.txt"
TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN =  T)
score_1<-read.table("000-vip_CIBERSORT-Results.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 343  25
#fix(score_1) 
score_2<-data.frame(ID=row.names(score_1), score_1, check.names=F)
dim(score_2) #[1] 343  26
score_3<-score_2[,c(1:23)] 
dim(score_3) #[1] 343  23
#fix(score_3) 
group_1<-read.table("0_TCGA_training-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  #343  9
#fix(group_1) 
H_L_1 = subset(group_1, risk =="high"|risk =="low",  select=risk)  
dim(H_L_1)   #343  1
#fix(H_L_1)   
H_L_2<-data.frame(ID=row.names(H_L_1), H_L_1, check.names=F)
#fix(H_L_2) 
dim(H_L_2) #[1] 343   2
risk_score_1<-merge(H_L_2,score_3, by.x="ID", by.y="ID")
#fix(risk_score_1)
dim(risk_score_1) #[1] 343  24
write.table(risk_score_1, file="3CIBERSORT-risk_TCGA.txt", quote=FALSE,  sep="\t") 
risk_score_2<-read.table("3CIBERSORT-risk_TCGA.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(risk_score_2)  #[1] 343  24
#fix(risk_score_2) 
row.names(risk_score_2)<-risk_score_2[,1]
#fix(risk_score_2)
risk_score_3<-risk_score_2[,-1] 
#fix(risk_score_3)  
dim(risk_score_3)  #[1] 343  24
risk_score_3$risk<-as.factor(risk_score_3$risk)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
x   #[1] "risk"
colnames(risk_score_3)[1]="Type"  
data_1=melt(risk_score_3,id.vars=c("Type"))
#fix(data_1)
colnames(data_1)=c("Type","immunocells","percentage")
colnames(data_1)
#fix(data_1)
#boxplot
tiff(file="3_CIBERSORT_all_immuno-cells_perc-TCGA_risk_ggboxplot-t.test-no-p.tiff", width =25,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="TCGA-High(171)_low(172)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="4_CIBERSORT_all_immuno-cells_perc-TCGA_risk_ggboxplot-t.testwith-P.tiff", width =100,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="TCGA-High(171)_low(172)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_grey(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
print(p1)
dev.off()
tiff(file="5_CIBERSORT_all_immuno-cells_perc-TCGA_risk_ggplot-t.test-no-p-ok-ok.tiff",width =30,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(immunocells,percentage,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="immunocells", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
colnames(risk_score_3)
#[1] "Type"                         "B cells naive"                "B cells memory"              
#[4] "Plasma cells"                 "T cells CD8"                  "T cells CD4 naive"           
#[7] "T cells CD4 memory resting"   "T cells CD4 memory activated" "T cells follicular helper"   
#[10] "T cells regulatory (Tregs)"   "T cells gamma delta"          "NK cells resting"            
#[13] "NK cells activated"           "Monocytes"                    "Macrophages M0"              
#[16] "Macrophages M1"               "Macrophages M2"               "Dendritic cells resting"     
#[19] "Dendritic cells activated"    "Mast cells resting"           "Mast cells activated"        
#[22] "Eosinophils"                  "Neutrophils"  
risk_score_4<-risk_score_3[,c("Type" ,"Macrophages M0" , "Macrophages M1","Macrophages M2")]
risk_score_4$Type<-as.factor(risk_score_4$Type)
#fix(risk_score_4)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_4)[1]
x   #[1] "risk"
colnames(risk_score_4)[1]="Type"  
data_2=melt(risk_score_4,id.vars=c("Type"))
#fix(data_2)
colnames(data_2)=c("Type","immunocells","percentage")
colnames(data_2)
#fix(data_2)
tiff(file="6_CIBERSORT_immu-perce-TCGA_risk_ggboxplot-t.testwith-P-style1.tiff", width =24,height =12,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_2, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="TCGA-High(171)_low(172)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
library("ggplot2")
tiff(file="7_CIBERSORT_immu-perce-TCGA_risk_ggplot-t.test-no-p-style2-OK.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_2,aes(immunocells,percentage,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_cleveland() + 
  labs(x="immunocells", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(45)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
tiff(file="7_CIBERSORT_immu-perce-TCGA_risk_ggplot-t.testwith-P-style2-OK.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_2,aes(immunocells,percentage,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_cleveland() + 
  labs(x="immunocells", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(45)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
dev.off()
library(readxl) 
library(ggplot2)
library(tidyr)
score_1<-read.table("000-vip_CIBERSORT-Results.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 343  25
#fix(score_1) 
score_4<-score_1[,1:22]
#fix(score_4)
dim(score_4) 
H_1 = subset(group_1, risk =="high",  select=risk)  
dim(H_1)   #171  1
#fix(H_1)  
H_1_CIBERSORT<-score_4[rownames(H_1),] 
dim(H_1_CIBERSORT)    #171  22
#fix(H_1_CIBERSORT) 
colnames(H_1_CIBERSORT) 
case_h<-t(H_1_CIBERSORT)
col=rainbow(nrow(case_h),s=0.7,v=0.7)
tiff(file="8_TCGA-high-risk_%(barplot)-style1-used.tiff",width =50,height = 25,units ="cm",compression="lzw",bg="white",res=300)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(case_h,col=col,yaxt="n",ylab="TCGA-high-risk（171）",xaxt="n",cex.lab=1.8           )
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(case_h),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(case_h[,ncol(case_h)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(case_h),col=col,pch=15,bty="n",cex=1.3)
dev.off()
case_h_a<-data.frame(ID=row.names(case_h), case_h, check.names=F)
library(readxl) 
library(ggplot2)
library(tidyr)
case_h_b<-gather(case_h_a, E1, E2, -ID)
#fix(case_h_b)
dim(case_h_b)
tiff(file="8_TCGA-high-risk_%(ggplot)-style2.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_h_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "TCGA-high-risk（171）", y = "%", fill = "immuno-cells")
dev.off()
library(RColorBrewer)
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01)) 
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
tiff(file="8_TCGA-high-risk_%(ggplot)-style3.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_h_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "TCGA-high-risk（171）", y = "%", fill = "immuno-cells") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
dev.off()
L_1 = subset(group_1, risk =="low",  select=risk)  
dim(L_1)   
#fix(L_1)  
L_1_CIBERSORT<-score_4[rownames(L_1),] 
dim(L_1_CIBERSORT)    #171  22
#fix(L_1_CIBERSORT) 
colnames(L_1_CIBERSORT) 
case_l<-t(L_1_CIBERSORT)
col=rainbow(nrow(case_h),s=0.7,v=0.7)
tiff(file="9_TCGA-low-risk-immuno-style1.tiff",width =50,height = 25,units ="cm",compression="lzw",bg="white",res=300)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(case_l,col=col,yaxt="n",ylab="TCGA-low-risk（172）",xaxt="n",cex.lab=1.8           )
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(case_l),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(case_l[,ncol(case_l)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(case_l),col=col,pch=15,bty="n",cex=1.3)
dev.off()
case_l_a<-data.frame(ID=row.names(case_l), case_l, check.names=F)
library(readxl) 
library(ggplot2)
library(tidyr)
case_l_b<-gather(case_l_a, E1, E2, -ID)
dim(case_l_b)
tiff(file="9_TCGA-low-risk-immuno-style2.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_l_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "TCGA-low-risk（172）", y = "%", fill = "immuno-cells")
dev.off()
library(RColorBrewer)
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01)) 
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
tiff(file="9_TCGA-low-risk-immuno-style3.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_l_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "TCGA-low-risk（172）", y = "%", fill = "immuno-cells") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
dev.off()
rm(list=ls())  
group_1<-read.table("0_TCGA_training-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1) 
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
fix(com)  
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_TCGA-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_TCGA-LIHC_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
fix(eSet)   
library(limma)
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,] 
dim(exp_data1)    
fix(exp_data1)  #
write.table(exp_data1, file="2_ESTIMATE_exp_TCGA.txt", quote=FALSE,  sep="\t") 
library(estimate)
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packTypes("estimate", repos=rforge, dependencies=TRUE)
filterCommonGenes(input.f="2_ESTIMATE_exp_TCGA.txt", 
                  output.f="3_commonGenes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "3_commonGenes.gct",
              output.ds="4_estimateScore.gct", 
              platform="illumina") 
scores=read.table("4_estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
fix(out)
dim(out)
out_1<-out[-1,]
fix(out_1)
dim(out_1) 
colnames(out_1) 
write.table(out_1,file="5_TCGA_ESTIMATE_scores.txt",sep="\t",quote=F)
score_1<-read.table("5_TCGA_ESTIMATE_scores.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 343 3
fix(score_1)   
score_2<-data.frame(ID=row.names(score_1), score_1, check.names=F)
fix(score_2) 
dim(score_2) #[1] 343  4
group_1<-read.table("0_TCGA_training-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  #343  9
fix(group_1) 
H_L_1 = subset(group_1, risk =="high"|risk =="low",  select=risk)  
dim(H_L_1)   #343  1
fix(H_L_1)   
H_L_2<-data.frame(ID=row.names(H_L_1), H_L_1, check.names=F)
fix(H_L_2) 
dim(H_L_2) #[1] 343   2
risk_score_1<-merge(H_L_2,score_2, by.x="ID", by.y="ID")
fix(risk_score_1)
dim(risk_score_1) #[1] 343  5  
colnames(risk_score_1) #[1] "ID" "risk" "StromalScore" "ImmuneScore"  "ESTIMATEScore"
row.names(risk_score_1)<-risk_score_1[,1]
risk_score_2<-risk_score_1[,-1]
fix(risk_score_2) 
write.table(risk_score_2, file="6_ESTIMATE_score_riskgroup_TCGA.txt", quote=FALSE,  sep="\t") 
risk_score_3<-read.table("6_ESTIMATE_score_riskgroup_TCGA.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
fix(risk_score_3)
dim(risk_score_3) 
risk_score_3$risk<-as.factor(risk_score_3$risk)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
x   #[1] "risk"
colnames(risk_score_3)[1]="Type"  
data_1=melt(risk_score_3,id.vars=c("Type"))
fix(data_1)
colnames(data_1)=c("Type","score_type","value")
colnames(data_1)
fix(data_1)
#boxplot
tiff(file="7_ESTIMATE_purity-TCGA_risk_ggboxplot-t.test-no-p.tiff", width =25,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="score_type", y="value", color = "Type", fill="Type",
            ylab="ESTIMATE_purity",  xlab="TCGA-High(171)_low(172)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="8_ESTIMATE_purity-TCGA_risk_ggboxplot-t.testwith-P.tiff", width =30,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="score_type", y="value", color = "Type", fill="Type",
            ylab="ESTIMATE_purity",  xlab="TCGA-High(171)_low(172)-t.test",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_grey(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
print(p1)
dev.off()
tiff(file="9_ESTIMATE_purity-TCGA_risk_ggplot-t.test-no-p-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="score_type", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
tiff(file="10_ESTIMATE_purity-TCGA_risk_ggplot-t.testwith-P-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="score_type", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
dev.off()
tiff(file="11_ESTIMATE_purity-TCGA_risk_ggplot-t.testwith-P-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_cleveland() + 
  labs(x="score_type", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(45)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
library("ggpubr")
fix(risk_score_3)
risk_score_3$Type=factor(risk_score_3$Type)
tiff(file="12_ESTIMATE_StromalScore-TCGA_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="StromalScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="StromalScore",  legend.title="TCGA",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="12_ESTIMATE_StromalScore-TCGA_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="StromalScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="StromalScore",  legend.title="TCGA",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
tiff(file="13_ESTIMATE_ImmuneScore-TCGA_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ImmuneScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="ImmuneScore",  legend.title="TCGA",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="13_ESTIMATE_ImmuneScore-TCGA_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ImmuneScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="ImmuneScore",  legend.title="TCGA",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
tiff(file="14_ESTIMATE_ESTIMATEScore-TCGA_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ESTIMATEScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="ESTIMATEScore",  legend.title="TCGA",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="14_ESTIMATE_ESTIMATEScore-TCGA_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ESTIMATEScore", fill = "Type",
         xlab="TCGA-high-low-risk",  ylab="ESTIMATEScore",  legend.title="TCGA",   ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
rm(list=ls())  
library(GSVA)
library(limma)
library(GSEABase)
group_1<-read.table("0_TCGA_training-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
##fix(group_1) 
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
##fix(com)  
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_TCGA-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_TCGA-LIHC_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
##fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
##fix(eSet)   
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,]  
dim(exp_data1)    
#fix(exp_data1) 
write.table(exp_data1, file="1_TCGA-LIHC_BATCH-3.txt", quote=FALSE,  sep="\t") 
geneSet=getGmt("0-immune cell.gmt", geneIdType=SymbolIdentifier())
ssGSEA_Score=gsva(exp_data1, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
norm_ssGSEA_Score=normalize(ssGSEA_Score)  
norm_ssGSEA_Score=rbind(id=colnames(norm_ssGSEA_Score),norm_ssGSEA_Score)
dim(norm_ssGSEA_Score)    
#fix(norm_ssGSEA_Score) 
norm_ssGSEA_Score<-norm_ssGSEA_Score[-1,]
dim(norm_ssGSEA_Score)    
#fix(norm_ssGSEA_Score)

ssGSEA_1<-apply(norm_ssGSEA_Score,2,as.numeric)
dim(ssGSEA_1)   
#fix(ssGSEA_1)
row.names(ssGSEA_1)<-row.names(norm_ssGSEA_Score)
dim(ssGSEA_1)    
#fix(ssGSEA_1)
ssGSEA_2<-ssGSEA_1[rowMeans(ssGSEA_1)>0,]  
dim(ssGSEA_2)    
#fix(ssGSEA_2) 
ssGSEA_3<-na.omit(ssGSEA_2)  
#fix(ssGSEA_3)
dim(ssGSEA_3) 
write.table(ssGSEA_3,file="2_norm_ssGSEA_Score.txt",sep="\t",quote=F) 

risk=read.table("1_TCGA-LIHC_risk-group.txt",sep="\t",header=T,check.names=F)
dim(risk)    
#fix(risk) 
ID<-row.names(risk)
risk2<-data.frame(ID,risk)
norm_ssGSEA_Score=read.table("2_norm_ssGSEA_Score.txt",sep="\t",header=T,row.names=1,check.names=F)
dim(norm_ssGSEA_Score)    
#fix(norm_ssGSEA_Score) 
rownames(norm_ssGSEA_Score)
score<-t(norm_ssGSEA_Score)
score_1<-data.frame(ID=row.names(score), score, check.names=F)
#fix(score_1)
dim(score_1)
risk_score_1<-merge(risk2,score_1, by.x="ID", by.y="ID")
#fix(risk_score_1) 
dim(risk_score_1)
risk_score_2=risk_score_1[order(risk_score_1$group),]  
#fix(risk_score_2) 
dim(risk_score_2)
write.table(risk_score_2,file="3_risk-immuno-comb.txt",sep="\t",quote=F, row.names=F)
risk_score_3=read.table("3_risk-immuno-comb.txt",sep="\t",header=T,check.names=F, row.names=1)
dim(risk_score_3)    
#fix(risk_score_3) 
colnames(risk_score_3)[1]<-"risk"
#fix(risk_score_3) 
score_2<-risk_score_3[,-1]
#fix(score_2) 
score_3<-t(score_2)
#fix(score_3) 
risk3<-as.data.frame(risk_score_3[,c("risk")])
#fix(risk3)
row.names(risk3)<-row.names(risk_score_3)
colnames(risk3)<-"risk"
#fix(risk3)
library(pheatmap)
tiff(file="4_TCGA-LIHC_risk-immuno_heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_3, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "TCGA-OS-risk-28-immune-cells-heatmap")
dev.off()
risk_score_3$group<-as.factor(risk_score_3$group)
#fix(risk_score_3)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
colnames(risk_score_3)[1]="Type"
data_1=melt(risk_score_3,id.vars=c("Type"))
#fix(data_1)
colnames(data_1)=c("Type","immunocells","level")
colnames(data_1)
#fix(data_1)
#boxplot
tiff(file="5_all-immuno-TCGA-RISK-ggboxplot-t.test-no-p.tiff", width =120,height = 30,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="level", color = "Type", fill="Type",
            ylab="TCGA_immunocells level",  xlab="",
            legend.title=x, 
            # palette = "nejm", 
            palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            # notch = TRUE,
            # theme_gray(), 
            width=0.1)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="5_all-immuno-TCGA-RISK-ggboxplot-t.testwith-P.tiff", width =100,height =25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="level", color = "Type", fill="Type",
            ylab="TCGA_immunocells level",  xlab="",
            legend.title=x, 
            #palette = "nejm", 
            palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            # notch = TRUE,
            # theme_gray(), 
            width=0.1)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
tiff(file="5_all_immuno_TCGA-RISK-ggplot-t.test-no-p-ok.tiff",width =50,height = 25,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(immunocells,level,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x="immunocells", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
colnames(risk_score_3)
#fix(score_3) 
score_4<-score_3[c(2:4,6,13,14,16:20,22:28),]
rownames(score_4)
dim(score_4) #[1]  18 343
risk_score_4<-risk_score_3[,c("Type", rownames(score_4))]
risk_score_4$Type<-as.factor(risk_score_4$Type)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_4)[1]
colnames(risk_score_4)[1]="Type"
data_2=melt(risk_score_4,id.vars=c("Type"))
#fix(data_2)
colnames(data_2)=c("Type","immunocells","level")
colnames(data_2)
#fix(data_2)
tiff(file="6_sig-immuno-TCGA-risk-ggplot-t.test-no-p-ok.tiff",width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_2,aes(immunocells,level,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x="immunocells", y = "TCGA-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T )
dev.off()
library(pheatmap)
tiff(file="7_TCGA-LIHC_risk-sig-stat-immu-heatmap-style1.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_4, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "TCGA-OS-risk-heatmap")
dev.off()
##fix(score_3) 
score_5<-score_3[c(2:3,16,26,28),]
dim(score_5)
library(pheatmap)
tiff(file="7_TCGA-LIHCrisk-sig-positive-immune_heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_5, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "TCGA-OS-risk--heatmap")
dev.off()
rm(list=ls())  
options(stringsAsFactors = F)
source("0_CIBERSORT(no change).R")
group_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1)  
com = subset(group_1, risk =="low"|risk =="high",  select=risk) 
dim(com)   
#fix(com)  
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_GEO-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_GSE14520_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
#fix(eSet)   
library(limma)
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,] 
dim(exp_data1)    
#fix(exp_data1) 
write.table(exp_data1, file="2_CIBERSORT_exp_GEO.txt", quote=FALSE,  sep="\t") 
LM22.file <- "0_LM22.txt"
exp.file <- "2_CIBERSORT_exp_GEO.txt"
TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN =  T)
score_1<-read.table("000-vip_CIBERSORT-Results.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 242  25
#fix(score_1)  
score_2<-data.frame(ID=row.names(score_1), score_1, check.names=F)
#fix(score_2) 
dim(score_2) #[1] 242  26
score_3<-score_2[,c(1:23)] 
dim(score_3) #[1] 242  23
#fix(score_3) 
group_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  #242  9
#fix(group_1)  
H_L_1 = subset(group_1, risk =="high"|risk =="low",  select=risk)  
dim(H_L_1)   #242  1
#fix(H_L_1) 
H_L_2<-data.frame(ID=row.names(H_L_1), H_L_1, check.names=F)
#fix(H_L_2) 
dim(H_L_2) #[1] 242   2
risk_score_1<-merge(H_L_2,score_3, by.x="ID", by.y="ID")
#fix(risk_score_1)
dim(risk_score_1) #[1] 242  24
write.table(risk_score_1, file="3CIBERSORT-risk_GEO.txt", quote=FALSE,  sep="\t") 
risk_score_2<-read.table("3CIBERSORT-risk_GEO.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(risk_score_2)  #[1] 242  24
#fix(risk_score_2) 
row.names(risk_score_2)<-risk_score_2[,1]
#fix(risk_score_2)
risk_score_3<-risk_score_2[,-1] 
#fix(risk_score_3)  
dim(risk_score_3)  #[1] 242  24
risk_score_3$risk<-as.factor(risk_score_3$risk)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
x   #[1] "risk"
colnames(risk_score_3)[1]="Type"  
data_1=melt(risk_score_3,id.vars=c("Type"))
#fix(data_1)
colnames(data_1)=c("Type","immunocells","percentage")
colnames(data_1)
#fix(data_1)
#boxplot
tiff(file="3_CIBERSORT_all_immuno-cells_perc-GEO_risk_ggboxplot-t.test-no-p.tiff", width =25,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="GEO-High(195)_low(47)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="4_CIBERSORT_all_immuno-cells_perc-GEO_risk_ggboxplot-t.testwith-P.tiff", width =100,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="GEO-High(195)_low(47)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_grey(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
print(p1)
dev.off()
tiff(file="5_CIBERSORT_all_immuno-cells_perc-GEO_risk_ggplot-t.test-no-p-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(immunocells,percentage,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="immunocells", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()

colnames(risk_score_3)
#[1] "Type"                         "B cells naive"                "B cells memory"              
#[4] "Plasma cells"                 "T cells CD8"                  "T cells CD4 naive"           
#[7] "T cells CD4 memory resting"   "T cells CD4 memory activated" "T cells follicular helper"   
#[10] "T cells regulatory (Tregs)"   "T cells gamma delta"          "NK cells resting"            
#[13] "NK cells activated"           "Monocytes"                    "Macrophages M0"              
#[16] "Macrophages M1"               "Macrophages M2"               "Dendritic cells resting"     
#[19] "Dendritic cells activated"    "Mast cells resting"           "Mast cells activated"        
#[22] "Eosinophils"                  "Neutrophils"  
risk_score_4<-risk_score_3[,c("Type" ,"Plasma cells" , "T cells CD4 memory resting","Macrophages M0")]
risk_score_4$Type<-as.factor(risk_score_4$Type)
#fix(risk_score_4)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_4)[1]
x   #[1] "risk"
colnames(risk_score_4)[1]="Type"  
data_2=melt(risk_score_4,id.vars=c("Type"))
#fix(data_2)
colnames(data_2)=c("Type","immunocells","percentage")
colnames(data_2)
#fix(data_2)
tiff(file="6_CIBERSORT_immu-perce-GEO_risk_ggboxplot-t.testwith-P-style1.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_2, x="immunocells", y="percentage", color = "Type", fill="Type",
            ylab="CIBERSORT_immuno-cells_perc",  xlab="GEO-High(195)_low(47)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
library("ggplot2")
tiff(file="7_CIBERSORT_immu-perce-GEO_risk_ggplot-t.testwith-P-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_2,aes(immunocells,percentage,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_cleveland() + 
  labs(x="immunocells", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(45)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
library(readxl) 
library(ggplot2)
library(tidyr)
score_1<-read.table("000-vip_CIBERSORT-Results.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 242  25
#fix(score_1)  
score_4<-score_1[,1:22]
#fix(score_4)
dim(score_4) 
H_1 = subset(group_1, risk =="high",  select=risk)  
dim(H_1)   #195  1
#fix(H_1)   
H_1_CIBERSORT<-score_4[rownames(H_1),] 
dim(H_1_CIBERSORT)    #195  22
#fix(H_1_CIBERSORT)  
colnames(H_1_CIBERSORT) 
case_h<-t(H_1_CIBERSORT)
#fix(case_h) 
col=rainbow(nrow(case_h),s=0.7,v=0.7)
tiff(file="8_GEO-high-risk_%(barplot)-style1.tiff",width =50,height = 25,units ="cm",compression="lzw",bg="white",res=300)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(case_h,col=col,yaxt="n",ylab="GEO-high-risk（195）",xaxt="n",cex.lab=1.8           )
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(case_h),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(case_h[,ncol(case_h)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(case_h),col=col,pch=15,bty="n",cex=1.3)
dev.off()
case_h_a<-data.frame(ID=row.names(case_h), case_h, check.names=F)
library(readxl) 
library(ggplot2)
library(tidyr)
case_h_b<-gather(case_h_a, E1, E2, -ID)
dim(case_h_b)
tiff(file="8_GEO-high-risk_%(ggplot)-style2.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_h_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "GEO-high-risk（195）", y = "%", fill = "immuno-cells")
dev.off()
library(RColorBrewer)
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01)) 
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

tiff(file="8_GEO-high-risk_%(ggplot)-style3.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_h_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "GEO-high-risk（195）", y = "%", fill = "immuno-cells") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
dev.off()
L_1 = subset(group_1, risk =="low",  select=risk)  
dim(L_1)   
#fix(L_1)  
L_1_CIBERSORT<-score_4[rownames(L_1),] 
dim(L_1_CIBERSORT)    #195  22
#fix(L_1_CIBERSORT)  
colnames(L_1_CIBERSORT) 
case_l<-t(L_1_CIBERSORT)
#fix(case_l) 
col=rainbow(nrow(case_h),s=0.7,v=0.7)
tiff(file="9_GEO-low-risk-immuno-style1.tiff",width =50,height = 25,units ="cm",compression="lzw",bg="white",res=300)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(case_l,col=col,yaxt="n",ylab="GEO-low-risk（47）",xaxt="n",cex.lab=1.8           )
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(case_l),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(case_l[,ncol(case_l)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(case_l),col=col,pch=15,bty="n",cex=1.3)
dev.off()
case_l_a<-data.frame(ID=row.names(case_l), case_l, check.names=F)
library(readxl) 
library(ggplot2)
library(tidyr)
case_l_b<-gather(case_l_a, E1, E2, -ID)
#fix(case_l_b)
dim(case_l_b)
tiff(file="9_GEO-low-risk-immuno-style2.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_l_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "GEO-low-risk（47）", y = "%", fill = "immuno-cells")
dev.off()
library(RColorBrewer)
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01)) 
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
tiff(file="9_GEO-low-risk-immuno-style3.tiff",width =100,height = 30,units ="cm",compression="lzw",bg="white",res=300)
ggplot(case_l_b) +
  geom_bar(aes(x = E1, y = E2, fill =ID ),
           stat = "identity") +
  labs(x = "GEO-low-risk（47）", y = "%", fill = "immuno-cells") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
dev.off()
group_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1)  
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
#fix(com)  
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_GEO-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_GSE14520_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
#fix(eSet)   
library(limma)
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,] 
dim(exp_data1)    
#fix(exp_data1) 
write.table(exp_data1, file="2_ESTIMATE_exp_GEO.txt", quote=FALSE,  sep="\t") 
library(estimate)
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packTypes("estimate", repos=rforge, dependencies=TRUE)
filterCommonGenes(input.f="2_ESTIMATE_exp_GEO.txt", 
                  output.f="3_commonGenes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "3_commonGenes.gct",
              output.ds="4_estimateScore.gct", 
              platform="illumina") 
scores=read.table("4_estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
#fix(out)
dim(out) #[1] 243   3
out_1<-out[-1,]
#fix(out_1)
dim(out_1)
colnames(out_1) 
write.table(out_1,file="5_GEO_ESTIMATE_scores.txt",sep="\t",quote=F)
score_1<-read.table("5_GEO_ESTIMATE_scores.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(score_1)  #[1] 242 3
#fix(score_1)  
score_2<-data.frame(ID=row.names(score_1), score_1, check.names=F)
#fix(score_2) 
dim(score_2) #[1] 242  4
group_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  #242  9
#fix(group_1)  
H_L_1 = subset(group_1, risk =="high"|risk =="low",  select=risk)  
dim(H_L_1)   #242  1
#fix(H_L_1) 
H_L_2<-data.frame(ID=row.names(H_L_1), H_L_1, check.names=F)
#fix(H_L_2) 
dim(H_L_2) #[1] 242   2
risk_score_1<-merge(H_L_2,score_2, by.x="ID", by.y="ID")
#fix(risk_score_1)
dim(risk_score_1) #[1] 242  5  
colnames(risk_score_1) #[1] "ID" "risk" "StromalScore" "ImmuneScore"  "ESTIMATEScore"
row.names(risk_score_1)<-risk_score_1[,1]
risk_score_2<-risk_score_1[,-1]
#fix(risk_score_2) 
write.table(risk_score_2, file="6_ESTIMATE_score_riskgroup_GEO.txt", quote=FALSE,  sep="\t") 
risk_score_3<-read.table("6_ESTIMATE_score_riskgroup_GEO.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
#fix(risk_score_3)
dim(risk_score_3) 
risk_score_3$risk<-as.factor(risk_score_3$risk)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
x   #[1] "risk"
colnames(risk_score_3)[1]="Type"  
data_1=melt(risk_score_3,id.vars=c("Type"))
#fix(data_1)
colnames(data_1)=c("Type","score_type","value")
colnames(data_1)
#fix(data_1)
#boxplot
tiff(file="7_ESTIMATE_purity-GEO_risk_ggboxplot-t.test-no-p.tiff", width =25,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="score_type", y="value", color = "Type", fill="Type",
            ylab="ESTIMATE_purity",  xlab="GEO-High(171)_low(172)",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_bw(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="8_ESTIMATE_purity-GEO_risk_ggboxplot-t.testwith-P.tiff", width =30,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="score_type", y="value", color = "Type", fill="Type",
            ylab="ESTIMATE_purity",  xlab="GEO-High(171)_low(172)-t.test",
            legend.title=x, 
            palette = "npg", 
            # palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            ggtheme = theme_grey(),
            #size=0.8,
            width=0.25)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
print(p1)
dev.off()
tiff(file="9_ESTIMATE_purity-GEO_risk_ggplot-t.test-no-p-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="score_type", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
tiff(file="10_ESTIMATE_purity-GEO_risk_ggplot-t.testwith-P-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_gray() + 
  labs(x="score_type", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.format", hide.ns =T)
dev.off()
tiff(file="11_ESTIMATE_purity-GEO_risk_ggplot-t.testwith-P-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(score_type,value,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_cleveland() + 
  labs(x="score_type", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(45)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
library("ggpubr")
#fix(risk_score_3)
risk_score_3$Type=factor(risk_score_3$Type)
tiff(file="12_ESTIMATE_StromalScore-GEO_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="StromalScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="StromalScore",  legend.title="GEO",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="12_ESTIMATE_StromalScore-GEO_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="StromalScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="StromalScore",  legend.title="GEO",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
tiff(file="13_ESTIMATE_ImmuneScore-GEO_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ImmuneScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="ImmuneScore",  legend.title="GEO",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="13_ESTIMATE_ImmuneScore-GEO_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ImmuneScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="ImmuneScore",  legend.title="GEO",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
tiff(file="14_ESTIMATE_ESTIMATEScore-GEO_risk_ggplot-t.test(no-p)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ESTIMATEScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="ESTIMATEScore",  legend.title="GEO",     ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white")) 
#  stat_compare_means(method = "wilcox.test", label.y = 11, label.x =2)  + 
#  stat_compare_means(method = "t.test", label.y = 10, label.x = 1)    
dev.off()
tiff(file="14_ESTIMATE_ImmuneScore-GEO_risk_ggplot-t.test(P)-style2.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(risk_score_3, x="Type", y="ESTIMATEScore", fill = "Type",
         xlab="GEO-high-low-risk",  ylab="ESTIMATEScore",  legend.title="GEO",   ggtheme = theme_grey(),
         add = "boxplot", add.params = list(fill="white"))+
  stat_compare_means(method = "wilcox.test", label.y = 1400, label.x =1.8)  + 
  stat_compare_means(method = "t.test", label.y = 1500, label.x = 0.9)    
dev.off()
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GSVA")   
#GSVA: The Gene Set Variation Analysis package for microarray and RNA-seq data
#[1] "Activated B cell"               "Activated CD4 T cell"           "Activated CD8 T cell"           "Central memory CD4 T cell"      "Central memory CD8 T cell"     
#[6] "Effector memeory CD4 T cell"    "Effector memeory CD8 T cell"    "Gamma delta T cell"             "Immature  B cell"               "Memory B cell"                 
#[11] "Regulatory T cell"              "T follicular helper cell"       "Type 1 T helper cell"           "Type 17 T helper cell"          "Type 2 T helper cell"          
#[16] "Activated dendritic cell"       "CD56bright natural killer cell" "CD56dim natural killer cell"    "Eosinophil"                     "Immature dendritic cell"       
#[21] "Macrophage"                     "Mast cell"                      "MDSC"                           "Monocyte"                       "Natural killer cell"           
#[26] "Natural killer T cell"          "Neutrophil"                     "Plasmacytoid dendritic cell"   

anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell',
                'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 
                'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell',
                'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 
               'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
rm(list=ls())  
library(GSVA)
library(limma)
library(GSEABase)
group_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
#fix(group_1)  
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
#fix(com)   
colnames(com)[1]<-c("group")
head(com,5)
write.table(com, file="1_GEO-LIHC_risk-group.txt", quote=FALSE,  sep="\t") 
gene_1<-read.table("0_GSE14520_batch-2.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
#fix(eSet)   
exp_data=as.matrix(eSet)
exp_data1<-avereps(exp_data)
exp_data1<-exp_data[rowMeans(exp_data1)>0,]  
dim(exp_data1)    
#fix(exp_data1) 
write.table(exp_data1, file="1_GEO-LIHC_BATCH-3.txt", quote=FALSE,  sep="\t") 
geneSet=getGmt("0-immune cell.gmt", geneIdType=SymbolIdentifier())
ssGSEA_Score=gsva(exp_data1, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))} 
norm_ssGSEA_Score=normalize(ssGSEA_Score) 
norm_ssGSEA_Score=rbind(id=colnames(norm_ssGSEA_Score),norm_ssGSEA_Score)
dim(norm_ssGSEA_Score)    
#fix(norm_ssGSEA_Score) 
norm_ssGSEA_Score<-norm_ssGSEA_Score[-1,]
dim(norm_ssGSEA_Score)   
ssGSEA_1<-apply(norm_ssGSEA_Score,2,as.numeric)  
dim(ssGSEA_1)    
#fix(ssGSEA_1)
row.names(ssGSEA_1)<-row.names(norm_ssGSEA_Score)
dim(ssGSEA_1)   
#fix(ssGSEA_1)
ssGSEA_2<-ssGSEA_1[rowMeans(ssGSEA_1)>0,] 
dim(ssGSEA_2)    
#fix(ssGSEA_2) 
ssGSEA_3<-na.omit(ssGSEA_2)  
#fix(ssGSEA_3)
dim(ssGSEA_3) 
write.table(ssGSEA_3,file="2_norm_ssGSEA_Score.txt",sep="\t",quote=F) 
risk=read.table("1_GEO-LIHC_risk-group.txt",sep="\t",header=T,check.names=F, row.names=1)
dim(risk)    
#fix(risk)
ID<-row.names(risk)
risk2<-data.frame(ID,risk)
#fix(risk2)
norm_ssGSEA_Score=read.table("2_norm_ssGSEA_Score.txt",sep="\t",header=T,row.names=1,check.names=F)#读取-result文件
dim(norm_ssGSEA_Score)    
#fix(norm_ssGSEA_Score) 
rownames(norm_ssGSEA_Score)
score<-t(norm_ssGSEA_Score)
#fix(score) 
dim(score)
score_1<-data.frame(ID=row.names(score), score, check.names=F)
#fix(score_1) 
dim(score_1)
risk_score_1<-merge(risk2,score_1, by.x="ID", by.y="ID")
#fix(risk_score_1) 
dim(risk_score_1)
risk_score_2=risk_score_1[order(risk_score_1$group),]  
#fix(risk_score_2) 
dim(risk_score_2)
write.table(risk_score_2,file="3_risk-immuno-comb.txt",sep="\t",quote=F, row.names=F) 
risk_score_3=read.table("3_risk-immuno-comb.txt",sep="\t",header=T,check.names=F, row.names=1)
dim(risk_score_3)    
#fix(risk_score_3) 
colnames(risk_score_3)[1]<-"risk"
#fix(risk_score_3) 
score_2<-risk_score_3[,-1]
#fix(score_2) 
score_3<-t(score_2)
#fix(score_3) 
risk3<-as.data.frame(risk_score_3[,c("risk")])
#fix(risk3)
row.names(risk3)<-row.names(risk_score_3)
colnames(risk3)<-"risk"
#fix(risk3)
library(pheatmap)
tiff(file="4_GEO-LIHC_risk-immuno_heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_3, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "GEO-OS-risk-28-immune-cells-heatmap")
dev.off()
risk_score_3$group<-as.factor(risk_score_3$group)
#fix(risk_score_3)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_3)[1]
colnames(risk_score_3)[1]="Type"

data_1=melt(risk_score_3,id.vars=c("Type"))
#fix(data_1)
colnames(data_1)=c("Type","immunocells","level")
colnames(data_1)
#fix(data_1)
#boxplot
tiff(file="5_all-immu-GEO-risk-ggboxplot-t.test-no-p.tiff", width =100,height = 25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="level", color = "Type", fill="Type",
            ylab="GEO_immunocells level",  xlab="",
            legend.title=x, 
            # palette = "nejm", 
            palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            # notch = TRUE,
            # theme_gray(), 
            width=0.1)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
print(p1)
dev.off()
tiff(file="5_all-immu-GEO-risk-ggboxplot-t.testwith-P.tiff", width =100,height =25,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_1, x="immunocells", y="level", color = "Type", fill="Type",
            ylab="GEO_immunocells level",  xlab="",
            legend.title=x, 
            #palette = "nejm", 
            palette = c("#5F9EA0","#FF6633"),
            # error.plot = "pointrange",
            #add = "none",
            # notch = TRUE,
            # theme_gray(), 
            width=0.1)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
tiff(file="5_all-immu-GEO-risk-ggplot-t.test-no-p-ok.tiff",width =18,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_1,aes(immunocells,level,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x="immunocells", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T)
dev.off()
colnames(risk_score_3)
#fix(score_3) 
score_4<-score_3[c(2:4,6,13,14,16:20,22:28),]
rownames(score_4)
dim(score_4) #[1]  18 242
risk_score_4<-risk_score_3[,c("Type", rownames(score_4))]
risk_score_4$Type<-as.factor(risk_score_4$Type)
library(reshape2)
library(ggpubr)
x=colnames(risk_score_4)[1]
colnames(risk_score_4)[1]="Type"
data_2=melt(risk_score_4,id.vars=c("Type"))
#fix(data_2)
colnames(data_2)=c("Type","immunocells","level")
colnames(data_2)
#fix(data_2)
tiff(file="6_sig-immu-GEO-risk-ggplot-t.test-no-p-ok.tiff",width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data_2,aes(immunocells,level,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x="immunocells", y = "GEO-Estiamted Proportion") +
  # theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())+
  rotate_x_text(90)+  
  stat_compare_means(aes(group=Type), method="t.test", label = "p.signif", hide.ns =T )
dev.off()
library(pheatmap)
tiff(file="7_GEO-LIHC_risk-sig-stat-immu-heatmap-style1.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_4, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "GEO-OS-risk-heatmap")
dev.off()
score_5<-score_3[c(2:3,16,26,28),]
dim(score_5)
library(pheatmap)
tiff(file="7_GEO-LIHCrisk-sig-positive-immune_heatmap.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(score_5, annotation=risk3, 
         color =colorRampPalette(color.key)(50),
         cluster_cols =F, 
         cluster_rows=T, 
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=T,
         fontsize_col=3,
         main = "GEO-OS-risk--heatmap")
dev.off()










#(17) clinical information sorting-exp-combi####################################################################
rm(list=ls())  
gxj_1<-read.table("0-TCGA-LIHC-PEACH-2022-ok.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gxj_1)  #[1] 364  46
fix(gxj_1)  
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) 
fix(gxj_2)
colnames(gxj_2)
str(gxj_2)
gxj_2$eastern_cancer_oncology_group  <- as.factor(gxj_2$eastern_cancer_oncology_group)
str(gxj_2)
gxj_3<-summary(gxj_2)
gxj_3
gxj_4<-as.data.frame(gxj_3)
fix(gxj_4) 
write.table(gxj_4,file="1_TCGA-LIHCsurv-count-mean.xls",sep="\t",quote=F,row.names=F)
colnames(gxj_2)
gxj_6<-data.frame(gxj_2[,(2:30)]) 
fix(gxj_6)
dim(gxj_6)
clx=data.frame()
for(i in 1:ncol(gxj_6)){
  nameStat=colnames(gxj_6)[i]
  tableStat=table(gxj_6[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx=rbind(clx,tableStatOut)}
fix(clx)
dim(clx)
write.table(clx,file="2_direct-one-class-LICH-COUNT-OK.xls",sep="\t",quote=F,row.names=F)
clx_1<-read.table("0_TCGA_training-riskscore.txt",sep='\t',quote="", comment.char="",header=T,row.names=1)  
colnames(clx_1) #[1] "os_time"   "os_stat"   "KPNA2"     "FMO3"      "SPP1"      "KIF2C"     "LPCAT1"    "riskScore" "risk" 
dim(clx_1)  
fix(clx_1)  
clx_2<-data.frame(ID=rownames(clx_1), clx_1)
dim(clx_2) 
fix(clx_2)
gxj_1<-read.table("0-TCGA-LIHC-PEACH-2022-ok.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gxj_1)  #[1] 364  50
fix(gxj_1)  
str(gxj_1)
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) #
fix(gxj_2)
all_1<-merge(gxj_2,clx_2,  by.x="ID", by.y = "ID")  
dim(all_1)    #[1] 336  60
fix(all_1)  
write.table(all_1, file="3_TCGA_LIHC_clinial_OS_hubexp_riskscorecomb.txt", quote=FALSE, row.names=F, sep="\t")  
all_1<-read.table("3_TCGA_LIHC_clinial_OS_hubexp_riskscorecomb.txt",sep='\t',quote="", comment.char="",header=T,row.names=1)  
dim(all_1)  
fix(all_1)   
colnames(all_1)
str(all_1)
all_2<-all_1[,c(1:29,55)] 
str(all_2)
all_2a<-all_2[all_2$risk %in% c("high") , ] 
dim(all_2a)
fix(all_2a)
clx_a=data.frame()
for(i in 1:ncol(all_2a)){
  nameStat=colnames(all_2a)[i]
  tableStat=table(all_2a[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx_a=rbind(clx_a,tableStatOut)}
fix(clx_a)
dim(clx_a)# [1] 411   4
colnames(clx_a) #[1] "Covariates" "Type"       "Total"      "tableStat" 
all_2b<-all_2[all_2$risk %in% c("low") , ] 
dim(all_2b)
clx_b=data.frame()
for(i in 1:ncol(all_2b)){
  nameStat=colnames(all_2b)[i]
  tableStat=table(all_2b[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx_b=rbind(clx_b,tableStatOut)}
fix(clx_b)
dim(clx_b)# [1] 411   4
colnames(clx_b) #[1] "Covariates" "Type"       "Total"      "tableStat" 
clx_a$risk<-c("high") 
colnames(clx_a) #[1] "Covariates" "Type"       "Total"      "tableStat"  "group" 
clx_b$risk<-c("low")  
colnames(clx_b) #[1] "Covariates" "Type"       "Total"      "tableStat"  "group" 
clx_all<-rbind(clx_a,clx_b)
dim(clx_all) #[1] 822   5
fix(clx_all)
write.table(clx_all,file="4_TCGA-LIHC-risk-two-class-count-OK.xls",sep="\t",quote=F,row.names=F)
#install.packages("ggpubr")
rm(list=ls())  
gxj_1<-read.table("0_TCGA_LIHC_clinial_OS_hubexp_riskscorecomb.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gxj_1)  #[1] 336  55
fix(gxj_1)  
colnames(gxj_1)
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) 
fix(gxj_2)
#
var1<-"pathologic_stage"  
var2<-"LPCAT1"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
S1<-summary(gxj_3[,"group"])[1]
S2<-summary(gxj_3[,"group"])[2]
S3<-summary(gxj_3[,"group"])[3]
S4<-summary(gxj_3[,"group"])[4]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
tiff(file=paste0("1_", var2,"and",var1,"ggboxplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p2=ggboxplot(gxj_3, x="group", y="index", color="group",
             xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
             ylab=paste0("S1:",S1,"; S2:",S2,"; S3:",S3,"; S4:",S4 ),
             legend.title="group",
             palette ="png", #"png" "jco"， orc("blue","red") 
             add = "jitter") +
  #stat_compare_means( method = "anova", label.y = max(gxj_3$FANCA_exp)+3, label.x =2)  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
p3=p2+rotate_x_text(45)  
print(p3)
dev.off()
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
M0<-summary(gxj_3[,"group"])[1]
M1<-summary(gxj_3[,"group"])[2]
MX<-summary(gxj_3[,"group"])[3]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
tiff(file=paste0("1_", var2,"and",var1,"ggboxplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p2=ggboxplot(gxj_3, x="group", y="index", color="group",
             xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
             ylab=paste0("M0:",M0,"; M1:",M1,"; MX:",MX),
             legend.title="group",
             palette ="png", #"png" "jco"， orc("blue","red") 
             add = "jitter") +
  #stat_compare_means( method = "anova", label.y = max(gxj_3$FANCA_exp)+3, label.x =2)  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
p3=p2+rotate_x_text(45)  
print(p3)
dev.off()
var1<-"N"  
var2<-"LPCAT1"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
N0<-summary(gxj_3[,"group"])[1]
N1<-summary(gxj_3[,"group"])[2]
NX<-summary(gxj_3[,"group"])[3]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
tiff(file=paste0("1_", var2,"and",var1,"ggboxplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p2=ggboxplot(gxj_3, x="group", y="index", color="group",
             xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
             ylab=paste0("N0:",N0,"; N1:",N1,"; NX:",NX),
             legend.title="group",
             palette ="png", #"png" "jco"， orc("blue","red") 
             add = "jitter") +
  #stat_compare_means( method = "anova", label.y = max(gxj_3$FANCA_exp)+3, label.x =2)  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
p3=p2+rotate_x_text(45)  
print(p3)
dev.off()
var1<-"T"  
var2<-"FMO3"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
T1<-summary(gxj_3[,"group"])[1]
T2<-summary(gxj_3[,"group"])[2]
T3<-summary(gxj_3[,"group"])[3]
T4<-summary(gxj_3[,"group"])[4]
TX<-summary(gxj_3[,"group"])[5]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
tiff(file=paste0("1_", var2,"and",var1,"ggboxplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p2=ggboxplot(gxj_3, x="group", y="index", color="group",
             xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
             ylab=paste0("T1:",T1,"; T2:",T2,"; T3:",T3," T4:",T4,"; TX:",TX),
             legend.title="group",
             palette ="png", #"png" "jco"， orc("blue","red") 
             add = "jitter") +
  #stat_compare_means( method = "anova", label.y = max(gxj_3$FANCA_exp)+3, label.x =2)  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
p3=p2+rotate_x_text(45)  
print(p3)
dev.off()
var1<-"histologic_grade"  
var2<-"FMO3"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
G1<-summary(gxj_3[,"group"])[1]
G2<-summary(gxj_3[,"group"])[2]
G3<-summary(gxj_3[,"group"])[3]
G4<-summary(gxj_3[,"group"])[4]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
tiff(file=paste0("2_", var2,"and",var1,"ggviolin.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
ggviolin(gxj_3, x="group", y="index", fill = "group",
         xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
         ylab=paste0("G1:",G1,"; G2:",G2,"; G3:",G3,"; G4:",G4),
         legend.title="group",
         add = "boxplot" , 
         size=0.5,
         add.params = list(fill="#FFFF66")
)+ stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
dev.off()
var1<-"eastern_cancer_oncology_group"  
var2<-"LPCAT1"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
E1<-summary(gxj_3[,"group"])[1]
E2<-summary(gxj_3[,"group"])[2]
E3<-summary(gxj_3[,"group"])[3]
E4<-summary(gxj_3[,"group"])[4]
E5<-summary(gxj_3[,"group"])[4]
group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  
library("ggpubr")
tiff(file=paste0("3_", var2,"and",var1,"ggbarplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p=ggbarplot(gxj_3, x = "group", y = "index",  
            xlab =paste0("E1:",E1,"; E2:",E2,"; E3:",E3," E4:",E4,"; E5:",E5),
            ylab =paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
            fill = "group", 
            color="#666666",size=0.6,  
            error.plot = "upper_errorbar",
            add = c("mean_se"), 
            add.params = list(size=1.2,width=0.5), 
            palette = "png",   
            x.text.angle=45, 
            ggtheme = theme_gray(),
            # size=0.5,
            position = position_dodge())
p=p+ stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
print(p)
dev.off()
var1<-"vascular_tumor_cell_type"  
var2<-c( "FMO3", "KPNA2","SPP1","KIF2C", "LPCAT1" )  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c(var1, var2)])
colnames(gxj_3) #[1]  "vascular_tumor_cell_type" "FMO3"  "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1" 
gxj_3<-na.omit(gxj_3) 
colnames(gxj_3)[1]<-c("group")
colnames(gxj_3) # [1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
summary(gxj_3[,"group"])
Macro<-summary(gxj_3[,"group"])[1]
Micro<-summary(gxj_3[,"group"])[2]
None<-summary(gxj_3[,"group"])[3]
colnames(gxj_3) #[1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
x=colnames(gxj_3)[1] 
colnames(gxj_3)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_3)[2:ncol(gxj_3)]){
  gxj_3a=gxj_3[,c(gene,"Type")]
  colnames(gxj_3a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_3a$Type)))>2){ test=kruskal.test(expression ~ Type, data = gxj_3a)  
  p=test$p.value  }
  else{test=wilcox.test(expression ~ Type, data = gxj_3a)
  p=test$p.value  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)}
colnames(gxj_3)=paste0(colnames(gxj_3),geneSig)
library(reshape2)
data=melt(gxj_3,id.vars=c("Type"))
#fix(data)
colnames(data)=c("Type","Gene","Expression")
#
library(ggplot2)
p2=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = paste0("Macro:",Macro,"; Micro:",Micro,"; None:",None), 
       y = "hub Gene expression_kruskal")+
  geom_boxplot()+ 
  facet_wrap(~Gene,nrow =1)+ theme_gray()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff(file=paste0("4_", var1,"ggplot-three-group.tiff"), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=600)
print(p2)
dev.off()
var1<-"adjacent_hepatic_tissue_inflammation"  
var2<-c( "FMO3", "KPNA2","SPP1","KIF2C", "LPCAT1" )  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c(var1, var2)])
colnames(gxj_3) #[1]  "adjacent_hepatic_tissue_inflammation" "FMO3"  "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1" 
gxj_3<-na.omit(gxj_3) 
colnames(gxj_3)[1]<-c("group")
colnames(gxj_3) # [1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
summary(gxj_3[,"group"])
Mild<-summary(gxj_3[,"group"])[1]
None<-summary(gxj_3[,"group"])[2]
Severe<-summary(gxj_3[,"group"])[3]

colnames(gxj_3) #[1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
x=colnames(gxj_3)[1] 
colnames(gxj_3)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_3)[2:ncol(gxj_3)]){
  gxj_3a=gxj_3[,c(gene,"Type")]
  colnames(gxj_3a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_3a$Type)))>2){ test=kruskal.test(expression ~ Type, data = gxj_3a)  
  p=test$p.value  }
  else{test=wilcox.test(expression ~ Type, data = gxj_3a)
  p=test$p.value  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)}
colnames(gxj_3)=paste0(colnames(gxj_3),geneSig)

library(reshape2)
data=melt(gxj_3,id.vars=c("Type"))
#fix(data)
colnames(data)=c("Type","Gene","Expression")
#
library(ggplot2)
p3=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = paste0("Mild:",Mild,"; None:",None,"; Severe:",Severe), 
       y = "hub Gene expression_kruskal")+
  geom_violin()+ geom_boxplot(width=0.2,position=position_dodge(0.9), fill="white")+ 
  facet_wrap(~Gene,nrow =1)+ theme_gray()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff(file=paste0("5_", var1,"group-hubggplot-3.tiff"), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=600)
print(p3)
dev.off()
var1<-"Tumor_status"  
var2<-c( "FMO3", "KPNA2","SPP1","KIF2C", "LPCAT1" )  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c(var1, var2)])
colnames(gxj_3) #[1]  "Tumor_status" "FMO3"  "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1" 
gxj_3<-na.omit(gxj_3) 
colnames(gxj_3)[1]<-c("group")
colnames(gxj_3) # [1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
summary(gxj_3[,"group"])
a<-summary(gxj_3[,"group"])[1]
b<-summary(gxj_3[,"group"])[2]
colnames(gxj_3) #[1] "group"  "FMO3"   "KPNA2"  "SPP1"   "KIF2C"  "LPCAT1"
x=colnames(gxj_3)[1] 
colnames(gxj_3)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_3)[2:ncol(gxj_3)]){
  gxj_3a=gxj_3[,c(gene,"Type")]
  colnames(gxj_3a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_3a$Type)))>2){ test=kruskal.test(expression ~ Type, data = gxj_3a)  
  p=test$p.value  }
  else{test=wilcox.test(expression ~ Type, data = gxj_3a)
  p=test$p.value  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)}
colnames(gxj_3)=paste0(colnames(gxj_3),geneSig)
library(reshape2)
data=melt(gxj_3,id.vars=c("Type"))
#fix(data)
colnames(data)=c("Type","Gene","Expression")
#
library(ggplot2)
p4=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = paste0("Tumor_free:",a,"; With_tumor:",b), 
       y = "hub Gene expression_wilcox")+
  geom_boxplot() +
  geom_jitter(size = 1.0, shape = 19,  stroke = 0.5,width = 0.15, color = c("#FF66FF"),	show.legend = FALSE )+ #for stripchart
  facet_wrap(~Gene,nrow =1)+ theme_gray()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff(file=paste0("6_", var1,"group-hub-ggplot.tiff"), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=600)
print(p4)
dev.off()
#• geom_boxplot() for box plot
#• geom_violin() for violin plot
#• geom_dotplot() for dot plot
#•geom_jitter() for stripchart
#• geom_line() for line plot
#• geom_bar() for bar plot
var1<-"risk"  
var2<-c( "fetoprotein")  #  "weight" "height" "BMI" "albumin" "creatinine" "fetoprotein"  "platelet" "prothrombin_time"  
gxj_3<-data.frame(gxj_2[,c(var1, var2)])
gxj_3<-na.omit(gxj_3) 
#fix(gxj_3)
dim(gxj_3)
colnames(gxj_3)<-c("group","type")
colnames(gxj_3)
summary(gxj_3[,"group"])
a<-summary(gxj_3[,"group"])[1]
b<-summary(gxj_3[,"group"])[2]
p<-compare_means(type~group, data =gxj_3, method = "t.test")  
p_t<-p[1,5]
p_t
p_w<-compare_means(type~group, data =gxj_3, method = "wilcox.test")  
p_w<-p[1,5]
p_w
library(reshape2)
data=melt(gxj_3,id.vars=c("group"))
#fix(data) 
colnames(data)=c("group","index","value")
#fix(data) 
#boxplot
library(ggpubr)
tiff(file=paste0("7_", var1,"group", var2,"ggboxplot-single-index-two-group.tiff"), width =10,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p1=ggboxplot(data, x="index", y="value", color = "group", 
             ylab=paste0(var2,"--t.test:",p_t,"; wilcox:",p_w),
             xlab=paste0("high_risk:",a,"; low_risk:",b),
             legend.title=colnames(gxj_3)[1], 
             ggtheme =theme_grey(), #选择主题theme_gray(), theme_bw(), theme_minimal(), theme_classic(), theme_void(),
             palette ="npg", # c("blue","red")，"npg" "jama"
             width=0.4, add = "jitter")
p2=p1+rotate_x_text(60)
print(p2)
dev.off()
colnames(gxj_2) #
index<-c("FMO3", "KPNA2","SPP1","KIF2C", "LPCAT1", "riskScore" )
rt<-gxj_2[,index] 
#fix(rt) 
dim(rt) #[1] 336   6
cli<-c("risk","gender", "race","ethnicity","histologic_grade","pathologic_stage","M", "N", "T",             
       "residual_tumor", "vascular_tumor_cell_type", "adjacent_hepatic_tissue_inflammation","Tumor_status",
       "vital_status", "child_pugh_classification_grade","eastern_cancer_oncology_group",      
       "fibrosis_ishak_score","age_2" , 
       "Alcohol_consumption", "Hepatitis_B",          
       "Hepatitis_C", "Smoking","Hemochromatosis",  "Non_Alcoholic_Fatty_Liver_Disease",
       "relative_family_cancer_history", "history_of_neoadjuvant_treatment")
Type<-data.frame(gxj_2[,cli]) 
#fix(Type) 
dim(Type)  #[1] 336  28
colnames(Type)
clx<-"risk" 
Type=Type[order(Type[,clx]),]   
#fix(Type) 
rt=rt[row.names(Type),] 
rt<-as.data.frame(t(rt)) 
library(pheatmap)   
#pdf(file="8_hub-clin-heatmap.pdf", height=12,width=12) 
tiff(file="8_hub-clin-heatmap-style1（main）.tiff", width =24,height = 12,units ="cm",compression="lzw",bg="white",res=600)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("#0033CC", "white", "#A20055"))(10),
         cluster_cols =F,    
         scale="row",   
         show_colnames=F,
         fontsize=5,
         fontsize_row=5,
         fontsize_col=5)
dev.off()
#pdf(file="8_hub-clin-heatmap.pdf", height=30,width=60) 
tiff(file="8_hub-clin-heatmap-style2-note.tiff", width =24,height = 30,units ="cm",compression="lzw",bg="white",res=600)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("#0033CC", "white", "#A20055"))(10),
         cluster_cols =F,    
         scale="row",   
         show_colnames=F,
         fontsize=5,
         fontsize_row=5,
         fontsize_col=5)
dev.off()
rm(list=ls())  
gxj_1<-read.table("0_TCGA_LIHC_clinial_OS_hubexp_riskscorecomb.txt",header=T,sep="\t",check.names=F,row.names=1)     
dim(gxj_1)  # [1] 4259  183
#fix(gxj_1) 
#group="race" 
library(survival)
library(survminer)
var1<-"risk" 
var2<- "drug_treatment" 
rt<-gxj_1[,c( "OS_time","OS_status", var2, var1)]
rt<-na.omit(rt)
#fix(rt)
Type1=paste0(rt[,var1], "/", rt[,var2]) 
rt1=cbind(rt, Type1) 
colnames(rt1)[1:2]=c("futime","fustat")
#fix(rt1)
length=length(levels(factor(Type1)))
diff=survdiff(Surv(futime, fustat) ~ Type1,data = rt1)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type1, data = rt1)
surPlot1=ggsurvplot(fit, 
                    data=rt1,
                    conf.int=F,
                    pval=pValue,
                    pval.size=5,
                    legend.title="group",
                    legend.labs=levels(factor(rt1[,"Type1"])),
                    legend = "top", #c("top", "bottom", "left", "right", "none")
                    xlab=paste0(var2,"_Time(years)"),
                    break.time.by = 1,
                    risk.table=T,
                    risk.table.title="",
                    risk.table.height=.3)
tiff(file=paste0("1-TCGA-LIHC-",var1, "-",var2, "_two-index-curve-ggsurvplot.TIFF"), width =18,height =18,units ="cm",compression="lzw",bg="white",res=300)
print(surPlot1)
dev.off()
rm(list=ls())  
gxj_1<-read.table("0_GSE14520_Tumor-cli--sort-2022-OK.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gxj_1)  #[1] 364  46
fix(gxj_1)  
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) 
fix(gxj_2)
colnames(gxj_2)
#[1] "ID"               "Survival.status"  "OS_time"          "Recurr_status"    "Recurr_time"      "gender"           "age"             
#[8] "HBV_viral_status" "multinodular"     "cirrhosis"        "TNM_stage"        "BCLC_staging"     "CLIP_staging"     "ALT"             
#[15] "Main_tumor_size"  "AFP"  
str(gxj_1)
gxj_3<-summary(gxj_1)
gxj_3
gxj_4<-as.data.frame(gxj_3)
fix(gxj_4) 
write.table(gxj_4,file="1_GSE14520surv-count-mean.xls",sep="\t",quote=F,row.names=F)
colnames(gxj_1)
#[1] "Survival status"  "OS_time"          "Recurr_status"    "Recurr_time"      "gender"           "age"              "HBV_viral_status"
#[8] "multinodular"     "cirrhosis"        "TNM_stage"        "BCLC_staging"     "CLIP_staging"     "ALT"              "Main_tumor_size" 
#[15] "AFP" 
gxj_6<-data.frame(gxj_1[,(5:15)]) 
fix(gxj_6)
dim(gxj_6)
clx=data.frame()
for(i in 1:ncol(gxj_6)){
  nameStat=colnames(gxj_6)[i]
  tableStat=table(gxj_6[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx=rbind(clx,tableStatOut)}

fix(clx)
dim(clx)
write.table(clx,file="2_direst-one-class-GSE14520-OK.xls",sep="\t",quote=F,row.names=F)

clx_1<-read.table("0_GSE14520_OS_valid-riskscore.txt",sep='\t',quote="", comment.char="",header=T,row.names=1)  
colnames(clx_1) #[1] "os_time"   "os_stat"   "KPNA2"     "FMO3"      "SPP1"      "KIF2C"     "LPCAT1"    "riskScore" "risk" 
dim(clx_1)   #[1] 242   9
fix(clx_1)   
clx_2<-data.frame(ID=rownames(clx_1), clx_1)
dim(clx_2) 
fix(clx_2)
----
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) 
fix(gxj_2)

all_1<-merge(gxj_2,clx_2,  by.x="ID", by.y = "ID") 
dim(all_1)    #[1] 242  25
fix(all_1)  
write.table(all_1, file="3_GSE14520_clinial_OS_hubexp_riskscorecomb.txt", quote=FALSE, row.names=F, sep="\t")  
all_1<-read.table("3_GSE14520_clinial_OS_hubexp_riskscorecomb (excel-sort).txt",sep='\t',quote="", comment.char="",header=T,row.names=1)  
dim(all_1)  
fix(all_1)   
colnames(all_1)
#[1] "os_time"          "os_stat"          "Recurr_status"    "Recurr_time"      "gender"           "age"              "HBV_viral_status"
#[8] "multinodular"     "cirrhosis"        "TNM_stage"        "BCLC_staging"     "CLIP_staging"     "ALT"              "Main_tumor_size" 
#[15] "AFP"              "KPNA2"            "FMO3"             "SPP1"             "KIF2C"            "LPCAT1"           "riskScore"       
#[22] "risk" 
str(all_1)
all_2<-all_1[,c(4:15,22)] 
str(all_2)
all_2a<-all_2[all_2$risk %in% c("high") , ] 
dim(all_2a)
fix(all_2a)
clx_a=data.frame()
for(i in 1:ncol(all_2a)){
  nameStat=colnames(all_2a)[i]
  tableStat=table(all_2a[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx_a=rbind(clx_a,tableStatOut)}
fix(clx_a)
dim(clx_a)# [1] 411   4
colnames(clx_a) #[1] "Covariates" "Type"       "Total"      "tableStat" 
all_2b<-all_2[all_2$risk %in% c("low") , ] 
dim(all_2b)
clx_b=data.frame()
for(i in 1:ncol(all_2b)){
  nameStat=colnames(all_2b)[i]
  tableStat=table(all_2b[,i])
  tableStatSum=cbind(Total=sum(tableStat),tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio,2)
  tableStatSum[,2]=paste(tableStatSum[,2]," (",tableStatRatio[,2],")",sep="")
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
  clx_b=rbind(clx_b,tableStatOut)}
fix(clx_b)
dim(clx_b)# [1] 411   4
colnames(clx_b) #[1] "Covariates" "Type"       "Total"      "tableStat" 
clx_a$risk<-c("high") 
colnames(clx_a) #[1] "Covariates" "Type"       "Total"      "tableStat"  "group" 
clx_b$risk<-c("low")  
colnames(clx_b) #[1] "Covariates" "Type"       "Total"      "tableStat"  "group" 
clx_all<-rbind(clx_a,clx_b)
dim(clx_all) #[1] 822   5
fix(clx_all)
write.table(clx_all,file="4_GSE14520-risk-two-class-count-OK.xls",sep="\t",quote=F,row.names=F)
rm(list=ls())  
gxj_1<-read.table("0_GSE14520_clinial_OS_hubexp_riskscorecomb (excel-sort).txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gxj_1)  #[1] 336  55
fix(gxj_1)  
colnames(gxj_1)
#[1] "os_time"          "os_stat"          "Recurr_status"    "Recurr_time"      "gender"           "age"              "HBV_viral_status"
#[8] "multinodular"     "cirrhosis"        "TNM_stage"        "BCLC_staging"     "CLIP_staging"     "ALT"              "Main_tumor_size" 
#[15] "AFP"              "KPNA2"            "FMO3"             "SPP1"             "KIF2C"            "LPCAT1"           "riskScore"       
#[22] "risk" 
gxj_2<-data.frame(ID=rownames(gxj_1), gxj_1)
dim(gxj_2) 
fix(gxj_2)
var1<-"TNM_stage"  
var2<-"FMO3"  #  "FMO3"     "KPNA2"  "SPP1"  "KIF2C"  "LPCAT1" 
#fix(gxj_2)
gxj_3<-data.frame(gxj_2[,c("ID", var1, var2)])
colnames(gxj_3)<-c("ID", "group", "index")
gxj_3<-na.omit(gxj_3) 
library(ggpubr)
p_value<-compare_means(index~group, data =gxj_3, method = "kruskal.test")  
p.adj_1<-p_value[1,3]
p.adj_1
p_value_2<-compare_means(index~group, data =gxj_3, method = "anova")  
p.adj_2<-p_value_2[1,3]
p.adj_2
summary(gxj_3[,"group"])
S1<-summary(gxj_3[,"group"])[1]
S2<-summary(gxj_3[,"group"])[2]
S3<-summary(gxj_3[,"group"])[3]

group=levels(factor(gxj_3$group))
gxj_3$group=factor(gxj_3$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}  

tiff(file=paste0("1_", var2,"and",var1,"ggboxplot.tiff"), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=600)
p2=ggboxplot(gxj_3, x="group", y="index", color="group",
             xlab=paste0(var2,"_kruskal:",p.adj_1, "; anova:", p.adj_2),
             ylab=paste0("S1:",S1,"; S2:",S2,"; S3:",S3),
             legend.title="group",
             palette ="png", #"png" "jco"， orc("blue","red") 
             add = "jitter") +
  #stat_compare_means( method = "anova", label.y = max(gxj_3$FANCA_exp)+3, label.x =2)  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test" )
p3=p2+rotate_x_text(45)  
print(p3)
dev.off()

colnames(gxj_2) #
index<-c("FMO3", "KPNA2","SPP1","KIF2C", "LPCAT1" )
rt<-gxj_2[,index] 
#fix(rt) 
dim(rt) #[1] 336   6
cli<-c( "risk","gender" ,"age" , "HBV_viral_status","multinodular" ,
        "cirrhosis" ,"TNM_stage","BCLC_staging","CLIP_staging","ALT","Main_tumor_size","AFP")
Type<-data.frame(gxj_2[,cli]) 
#fix(Type) 
dim(Type)  #[1] 336  28
colnames(Type)
clx<-"risk" 
Type=Type[order(Type[,clx]),]   
#fix(Type) 
rt=rt[row.names(Type),] 
rt<-as.data.frame(t(rt)) 
library(pheatmap)   
#pdf(file="2_hub-clin-heatmap.pdf", height=12,width=12) 
tiff(file="2_hub-clin-heatmap-style1（main）.tiff", width =12,height = 10,units ="cm",compression="lzw",bg="white",res=600)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("#0033CC", "white", "#A20055"))(10),
         cluster_cols =F,    
         scale="row",   
         show_colnames=F,
         fontsize=5,
         fontsize_row=5,
         fontsize_col=5)
dev.off()
#pdf(file="2_hub-clin-heatmap.pdf", height=30,width=60) 
tiff(file="2_hub-clin-heatmap-style2-note.tiff", width =12,height = 15,units ="cm",compression="lzw",bg="white",res=600)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("#0033CC", "white", "#A20055"))(10),
         cluster_cols =F,    
         scale="row",   
         show_colnames=F,
         fontsize=5,
         fontsize_row=5,
         fontsize_col=5)
dev.off()
rm(list=ls())  
gxj_1<-read.table("0_GSE14520_clinial_OS_hubexp_riskscorecomb (excel-sort).txt",header=T,sep="\t",check.names=F,row.names=1)     
dim(gxj_1)  # [1] 4259  183
library(survival)
library(survminer)
var1<-"risk" 
var2<- "gender" 
rt<-gxj_1[,c( "os_time","os_stat", var2, var1)]
rt<-na.omit(rt)
#fix(rt)
Type1=paste0(rt[,var1], "/", rt[,var2]) 
rt1=cbind(rt, Type1) 
colnames(rt1)[1:2]=c("futime","fustat")
#fix(rt1)
length=length(levels(factor(Type1)))
diff=survdiff(Surv(futime, fustat) ~ Type1,data = rt1)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type1, data = rt1)
surPlot1=ggsurvplot(fit, 
                    data=rt1,
                    conf.int=F,
                    pval=pValue,
                    pval.size=5,
                    legend.title="group",
                    legend.labs=levels(factor(rt1[,"Type1"])),
                    legend = "top", #c("top", "bottom", "left", "right", "none")
                    xlab=paste0(var2,"_Time(years)"),
                    break.time.by = 1,
                    risk.table=T,
                    risk.table.title="",
                    risk.table.height=.3)
tiff(file=paste0("1-GSE14520-",var1, "-",var2, "_two-index-curve-ggsurvplot.TIFF"), width =18,height =18,units ="cm",compression="lzw",bg="white",res=300)
print(surPlot1)
dev.off()








































