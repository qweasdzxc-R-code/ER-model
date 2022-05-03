rm(list=ls()) 
work_dir <- "#/#/#" 
setwd(work_dir)

#(12) clusting analysis#######################################################################################
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





