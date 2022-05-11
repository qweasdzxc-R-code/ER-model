
# ER prognostic model  R transcripts

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


#(1) TCGA exp download#######################################################################
library("org.Hs.eg.db")
library("clusterProfiler")
library(SummarizedExperiment)  
library(TCGAbiolinks)  
work_dir <- "#/#/#" 
setwd(work_dir)
TCGAbiolinks:::getGDCprojects()$project_id   
project <- "TCGA-LIHC"   
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
workflow_type <- "HTSeq - Counts"
legacy <- FALSE
DataDirectory <- paste0(work_dir,"/GDC/",gsub("-","_",project))   
FileNameData <- paste0(DataDirectory, "_","RNAseq_HTSeq",".rda")
query <- GDCquery(project = project,
                  data.category = data_category,
                  data.type = data_type, 
                  workflow.type = workflow_type,
                  legacy = legacy)
samplesDown <- getResults(query,cols=c("cases"))
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,typesample = "NT")
GDCdownload(query = query,
            directory = DataDirectory,files.per.chunk=6, method='client')
data <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  DataDirectory,  
                   save.filename = FileNameData)
data_expr <- assay(data)  
fix(data_expr)
dim(data_expr)  
expr_file <- paste0(DataDirectory, "_","All_HTSeq_FPKM-UQ",".txt")
write.table(data_expr, file = expr_file, sep="\t", row.names =T, quote = F)
gene_info = read.table("./ref/gene_info.txt",header = 1)
gene_selected <- row.names(data_expr) %in% gene_info$Gene_id
gene_expr <- data_expr[gene_selected,]
gene_expr_file <- paste0(DataDirectory, "_","Gene_HTSeq_FPKM-UQ",".txt")
write.table(gene_expr, file = gene_expr_file, sep="\t", row.names =T, quote = F)
lncRNA_info = read.table("./ref/lncRNA_info.txt",header = 1)
lncRNA_selected <- row.names(data_expr) %in% lncRNA_info$Gene_id
lncRNA_expr <- data_expr[lncRNA_selected,]
lnc_expr_file <- paste0(DataDirectory, "_","LncRNA_HTSeq_FPKM-UQ",".txt")
write.table(lncRNA_expr, file = lnc_expr_file, sep="\t", row.names =T, quote = F)




#(2) clini download################################################################################
project <- "TCGA-LIHC"
data_category <- "Clinical"
data_type <- "Clinical Supplement"
file_type = "xml"
legacy <- FALSE
setwd(work_dir)
DataDirectory <- paste0(work_dir,"/GDC/",gsub("-","_",project))

query <- GDCquery(project = project,
                  data.category = data_category,
                  data.type = data_type, 
                  file.type = file_type,
                  legacy = legacy)
samplesDown <- getResults(query,cols=c("cases"))
GDCdownload(query = query,directory = DataDirectory,files.per.chunk=6, method='client')
clinical <- GDCprepare_clinic(query, clinical.info = "patient",directory = DataDirectory)
dim(clinical)
fix(clinical)
clinical_file <- paste0(DataDirectory, "_","clinical",".txt")
write.table(clinical, file = clinical_file, row.names = TRUE,  quote=F,sep="\t")
write.csv2(clinical, file = clinical_file, row.names = F, quote = F)
clinical_gxj <-GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")
fix(clinical_gxj )
dim(clinical_gxj )
write.csv(clinical_2, file = 'TCGA-LIHC_Clinical_OK2.csv', row.names =F)


#(3) Target gene-exp+clin combi################################################################################
gene_1<-read.table("00-TCGA_LIHC_All_HTSeq_FPKM-UQ-ok.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)  
gene_2<-log2(gene_1+1)   
dim(gene_2)  
fix(gene_2) 
gene_2$ENSEMBL<-rownames(gene_2)  
library("org.Hs.eg.db")
library("clusterProfiler")  
transfered_ID <- bitr(gene_2$ENSEMBL, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db") 
transfered_ID<-data.frame(transfered_ID) 
fix(transfered_ID)
Final_ID_acorss<- merge(transfered_ID, gene_2,   by.x="ENSEMBL", by.y = "ENSEMBL") 
dim(Final_ID_acorss)  
fix(Final_ID_acorss)  
Final_ID_acorss<-Final_ID_acorss[is.na(Final_ID_acorss[,"SYMBOL"])==F,]  
dim(Final_ID_acorss)  
library('limma')
exp_1<-avereps(Final_ID_acorss[,-c(1:2)],ID=Final_ID_acorss$SYMBOL)
fix(exp_1) 
dim(exp_1)
write.table(exp_1, "1_TCGA_LIHC_All_HTSeq_FPKM-UQ-ok.txt",  sep="\t",  quote=F, row.names=T)  
dim(Final_ID_acorss)  
aaa<-duplicated(Final_ID_acorss$SYMBOL)  
ID_exp_1<-Final_ID_acorss[!aaa,]
dim(ID_exp_1) 
fix(ID_exp_1)
rownames(ID_exp_1)<-ID_exp_1$SYMBOL
ID_exp<-ID_exp_1[,(1:2)] 
dim(ID_exp) 
fix(ID_exp)
write.table(ID_exp, "2_TCGA_LIHC_ENSEMBL-SYMBOL-2.txt",  sep="\t",  quote=F, row.names=F) 
# 1
all_gene_1<-read.table("1_TCGA_LIHC_All_HTSeq_FPKM-UQ-ok.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(all_gene_1)  
fix(all_gene_1) 
SYMBOL<-rownames(all_gene_1)
all_gene_2<-data.frame(SYMBOL, all_gene_1, check.names=FALSE)  
dim(all_gene_2)  
fix(all_gene_2)   
target_gene<-read.table("00-GSEA _ER（1975）.txt",sep='\t',quote="", comment.char="",stringsAsFactors = FALSE, check.names=FALSE,header=T,row.names=1)  
colnames(target_gene)
dim(target_gene)  
fix(target_gene)   
geneset_1<-merge(target_gene,all_gene_2,  by.x="Symbol", by.y = "SYMBOL") 
dim(geneset_1)    
fix(geneset_1)  
write.table(geneset_1, file="3_TCGA_LIHC_ER_log2（merge）.txt", quote=FALSE, row.names=T, sep="\t")  
mRNA_lncRNA<-read.table("3_TCGA_LIHC_ER_log2（merge）.txt",sep="\t",quote="",header =TRUE, comment.char="", check.names=FALSE)
dim(mRNA_lncRNA)  
fix(mRNA_lncRNA)  
gene.t<-t(mRNA_lncRNA)    
gene.t<-as.data.frame(gene.t)  
dim(gene.t)  
fix(gene.t)  
ID_full<-rownames(gene.t)
ID12<-substr(ID_full,1,12)
gene.m<-data.frame(ID_full, ID12, gene.t) 
dim(gene.m)  
fix(gene.m)  
group_1<-as.numeric(substr(gene.m$ID_full,14,14)) 
group<-ifelse(group_1==0, "Tumor", "Normal") 
group_gene_1<-data.frame(group, gene.m)  
dim(group_gene_1)    
fix(group_gene_1)  
group_gene_2<-group_gene_1[,c(1:3)]
dim(group_gene_2)    
fix(group_gene_2)  
clinical_1<-read.table("00-TCGA-LIHC_clinical-sort.txt",sep="\t",quote="",header =TRUE, comment.char="", check.names=FALSE)
dim(clinical_1)   
fix(clinical_1) 
gene.clinical.cross<-merge(group_gene_2,clinical_1,  by.x="ID12", by.y = "sample")  
dim(gene.clinical.cross)    
fix(gene.clinical.cross)  
write.table(gene.clinical.cross, file="4_LIHC_cli_IDfll_Tumor-Normal-group.xls", quote=FALSE, row.names=F, sep="\t")  



#(4) GSE14520 download-sort######################################################################################
library('Biobase')
library('GEOquery')
GSE = "GSE14520"
GPL = "GPL571"
gset = getGEO(GSE, GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx = grep(GPL, attr(gset, "names")) else idx = 1  
gset = gset[[idx]]
exprSet = exprs(gset)
dim(exprSet)
fix(exprSet)
write.table (exprSet, file= paste ("0_",  GSE, "_exp.txt", sep = ""),  quote=F, sep="\t")  
pData = pData(gset) 
dim(pData)
fix(pData)
colnames(pData)
write.table (pData, file= paste ("0_",  GSE, "_cli.txt", sep = ""), quote=F,sep="\t")  
fData = fData(gset)
dim(fData)
fix(fData)
write.table (fData, file=paste ("0_",  GSE, "_", GPL, "_note.txt", sep = ""),  quote=F,sep="\t")  
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library('Biobase')
library('GEOquery')
expression_1<-read.table("0_GSE14520-GPL571_series_matrix.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(expression_1) 
fix(expression_1)  
expression_2<-read.table("0_GSE14520-GPL3921_series_matrix.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(expression_2) 
fix(expression_2)  
expression <- data.frame(expression_1,expression_2[match(rownames(expression_2),rownames(expression_1)),])
dim(expression)  
fix(expression)
write.table(expression, file="1_GSE14520_cohort_1_2（22268probe，488sam）.txt", quote=FALSE,  sep="\t")   
probe_1<-read.table("0_GSE14520_GPL571_note.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(probe_1) 
fix(probe_1)   
probe<-read.table("0_GSE14520_GPL3921_note.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(probe) 
fix(probe)  
head(expression,1)[1:5]  
dim(expression)
n.sample = ncol(expression)   
cols = rainbow(n.sample)
pdf(file="2_boxplot_nom.pdf", width=24, height=18)
par(cex = 0.7)
if(n.sample>40) par(cex = 0.5)
boxplot(expression,col = cols, main = "expression", las = 2)
dev.off()
tiffFile=paste("2_boxplot_nom",".tiff",sep="")
tiff(file=tiffFile,width = 24,height = 18 ,units ="cm",compression="lzw",bg="white",res=600)
par(cex = 0.7)
if(n.sample>40) par(cex = 0.5)
boxplot(expression,col = cols, main = "expression", las = 2)
dev.off()
fix(probe)
colnames(probe)
probe_1<-probe[,c("Gene title","Gene ID")]  
fix(probe_1)
expression_1<-data.frame(probe_ID=row.names(expression),expression)
fix(expression_1)  
exptemp<- merge(probe_1, expression_1,  by.x="Gene title", by.y = "probe_ID")  
dim(exptemp)  
fix(exptemp)
colnames(exptemp)[c(1:2)]<-c("probeID","GeneID")
dim(exptemp)  #22268   490
fix(exptemp)
exptemp = exptemp[!is.na(exptemp$GeneID),]
#exptemp = exptemp[is.na(exptemp[,"GeneID"])==F,]   
#exptemp = exptemp[exptemp[,"GeneID"]!="",]
dim(exptemp)   #22268   490
exptemp=exptemp[!grepl(pattern="///",exptemp$GeneID),]
dim(exptemp) 
fix(exptemp) 
library('limma')
exp_2<-avereps(exptemp[,-c(1:2)],ID=exptemp$GeneID)
fix(exp_2) 
dim(exp_2)
write.table(exp_2,file="3_GSE14520_IDtran-sort-exp.txt",sep="\t",quote=F)
exp_ID<-read.table("3_GSE14520_IDtran-sort-exp.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T)  
dim(exp_ID) 
fix(exp_ID)  
datExp_1<-data.frame(GeneID=row.names(exp_ID),exp_ID) 
dim(datExp_1)
fix(datExp_1)
datExp=melt(datExp_1, id.vars="GeneID",variable.name="sample",value.name = "expression")  
dim(datExp)
fix(datExp)  
head(datExp,3)
P_density_GeneID = ggplot(datExp,aes(x=expression))+
  geom_density(aes(fill=as.character(datExp$sample),color=as.character(datExp$sample)),alpha = 0,size=1,linetype="solid")+
  labs(x="expression")+
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),
        legend.position = 'right',
            axis.line=element_line(size = 1,color="black"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
pdf(file=paste("4_density_1",".PDF", sep=""), width=96, height=36,onefile = FALSE)
print(P_density_GeneID)  
dev.off()
P_density_GeneID = ggplot(datExp,aes(x=expression))+
  geom_density(aes(fill=as.character(datExp$sample),color=as.character(datExp$sample)),alpha = 0,size=0.05,linetype="solid")
pdf(file=paste("4_density_2",".PDF", sep=""), width=96, height=36,onefile = FALSE)
print(P_density_GeneID)  
dev.off()
gene_1<-read.table("3_GSE14520_IDtran-sort-exp.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1) 
gene_2<-t(gene_1)
dim(gene_2)  
fix(gene_2) 
group_1<-read.table("0_GSE14520_cli-nornal-tumor-group.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1) 
Tumor = subset(group_1, group =="Tumor") 
dim(Tumor)   
fix(Tumor)  
head(Tumor,3)
dim(gene_2) 
gene_3<-gene_2[rownames(Tumor),] 
fix(gene_3) 
dim(gene_3)  
write.table(gene_3, file="5_GSE14520_IDtran-sort-exp_Tumor_sample.txt", quote=FALSE,  sep="\t")  
write.table(Tumor, file="5_GSE14520_Tumor_samplecli.txt", quote=FALSE,  sep="\t")   
write.table(gene_2, file="5_GSE14520_IDtran-sort-exp_all_sample.txt", quote=FALSE,  sep="\t")   




#(5) TCGA-GSE14520 batch match#################################################################################
library(limma)
library(sva) 
rt = read.table("0_TCGA_LIHC_ER_log2（merge）.txt",header=T,sep="\t",check.names=F, row.names=1)
fix(rt) 
dim(rt) 
rt=as.matrix(rt)   
rownames(rt)=rt[,1]  
exp=rt[,2:ncol(rt)]  
dimnames=list(rownames(exp),colnames(exp))  
metab=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
fix(metab) 
dim(metab)  
geo_1 = read.table("0_GSE14520_IDtran-sort-exp_all_sample.txt",header=T,sep="\t",check.names=F)
fix(geo_1) 
dim(geo_1) 
geo_2<-t(geo_1)
fix(geo_2) 
dim(geo_2)  
sameGene=intersect(row.names(metab),row.names(geo_2))  
length(sameGene)  
metabOut=metab[sameGene,] 
dim(metabOut)
geoOut=geo_2[sameGene,] 
dim(geoOut)
fix(geoOut)
fix(metabOut)
all=cbind(metabOut,geoOut)  
fix(all)  
dim(all)
clx<- apply(all, 1, function(row) all(row !=0 ))  
all_1<-all[clx,]
fix(all_1)
dim(all_1)  
all_2<-na.omit(all_1)  
fix(all_2)
dim(all_2) 
batchType=c(rep(1,ncol(metabOut)),rep(2,ncol(geoOut))) 
outTab=ComBat(all_1, batchType,par.prior=TRUE)
fix(outTab)
dim(outTab)  
write.table(outTab,file="1_TCGA-LIHC-GSE14520_batch.txt",sep="\t",quote=F)
metabOut_1=outTab[,colnames(metabOut)] 
metabOut_2=rbind(ID=colnames(metabOut),metabOut_1)  
dim(metabOut_2)
metabOut_3<-metabOut_2[-1,] 
fix(metabOut_3)
dim(metabOut_3)  
write.table(metabOut_3,file="2_TCGA-LIHC_batch.txt",sep="\t",quote=F)
geoOut_1=outTab[,colnames(geoOut)]
geoOut_2=rbind(ID=colnames(geoOut),geoOut_1) 
dim(geoOut_2)
geoOut_3<-geoOut_2[-1,] 
fix(geoOut_3)
dim(geoOut_3) 
write.table(geoOut_3,file="3_GSE14520_batch.txt",sep="\t",quote=F)


#(6) TCGA-Normal-Tumor-diff##################################################################################
rm(list=ls()) 
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0_TCGA-LIHC-Normal-HCCgroup-IDfull.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1) 
com = subset(group_1, group =="Tumor"|group =="Normal",  select=group)  
dim(com)   
fix(com)  
head(com,5)
gene_1<-read.table("1_TCGA-LIHC_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)  
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]  
dim(eSet)   
fix(eSet) 
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
fix(group_list) 
rownames(group_list) = rownames(com)       
colnames(group_list)[1] = "description"  
dim(group_list) 
fix(group_list) 
design = model.matrix(~ description + 0, group_list)  
fix(design)
colnames(design) = levels(as.factor(c("Tumor","Normal")))   
factor(design)  
head(design) 
cont.matrix = makeContrasts(Tumor-Normal, levels=design)   
fdr = 0.05
logFC = 1
fit = lmFit(eSet, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="BH", sort.by="logFC", resort.by = "P" ,n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")

tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "High" 
tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Low"  
tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "NO"
colnames(tT)[ncol(tT)] = "expression"
dim(tT)
fix(tT)  
write.table(tT,file="1_TCGA_Tumor_vs._Normal_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
fix(deg)  
write.table(deg,file="2_TCGA_Tumor_vs._Normal_dif_gene.txt",sep="\t",quote=F)
group_reg <- factor(deg$expression)
labels=levels(group_reg)
gxj<-summary(group_reg) 
gxj

temp1 = tT[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("High","Low","NO"), order=T)
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+  
  scale_color_manual(values =c("High" = "#8B0000", "NO" = "gray","Low" = "#008080" ))+   
  labs(x=paste("log2FC", "High:", gxj[1],"Low:", gxj[2]),y="-log10FDR")+
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
tiff(file="3_volcano-style1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_volcano-style2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "High:", gxj[1],"Low:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+   
  scale_color_manual(values =c("High" = "#8B0000", "NO" = "gray","Low" = "#008080"))
dev.off()

temp2 = data.frame(log10(rowMeans(eSet)))
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT),rownames(temp2)),] 
temp2 = data.frame(temp2,tT[,c("logFC","expression")])
colnames(temp2)=c("log10Exp","logFC","expression")
temp2$expression=factor(temp2$expression, levels=c("High","Low","NO"), order=T)
P_MA=ggplot(temp2,aes(x=temp2$log10Exp,y=temp2$logFC))+
  geom_point(alpha=0.8, size=6,  
             aes(color=temp2$expression))+
  scale_color_manual(values =c("High" = "#FF4500", "NO" = "gray","Low" = "#483D8B"))+
  labs(x="log10Exp",y=paste("log2FC", "High:", gxj[1],"Low:", gxj[2]))+
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
tiff(file="4_MA_plot.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)  #也注意投置片的大小
print(P_MA)
dev.off()

dim(eSet) 
#fix(eSet)
DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),]
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_all_gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
target_gene<-row.names(deg)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target_gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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



#(7) GEO-Normal-Tumor-diff############################################################################################################
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0-GSE14520-N-T-pair-232.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1)  
com = subset(group_1, group =="Tumor"|group =="Normal",  select=group)  
dim(com)   
fix(com)  
head(com,5)
gene_1<-read.table("0_GSE14520_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))   
eSet=gene_1[,gsm] 
dim(eSet)  
fix(eSet)  
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
fix(group_list) 
rownames(group_list) = rownames(com)       
colnames(group_list)[1] = "description"   
dim(group_list) 
fix(group_list) 
design = model.matrix(~ description + 0, group_list)  
fix(design)
colnames(design) = levels(as.factor(c("Tumor","Normal")))   
factor(design) 
head(design)  
cont.matrix = makeContrasts(Tumor-Normal, levels=design)  
fdr = 0.05
logFC = 1
fit = lmFit(eSet, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="BH", sort.by="logFC", resort.by = "P" ,n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "High" 
tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Low" 
tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "NO"
colnames(tT)[ncol(tT)] = "expression"
dim(tT)
fix(tT)  
write.table(tT,file="1_GEO_Tumor_vs._Normal_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
fix(deg)  
write.table(deg,file="2_GEO_Tumor_vs._Normal_diff_gene.txt",sep="\t",quote=F)
group_reg <- factor(deg$expression)
labels=levels(group_reg)
gxj<-summary(group_reg) 

temp1 = tT[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("High","Low","NO"), order=T)
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+  
  scale_color_manual(values =c("High" = "#8B0000", "NO" = "gray","Low" = "#008080" ))+  
  labs(x=paste("log2FC", "High:", gxj[1],"Low:", gxj[2]),y="-log10FDR")+
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
tiff(file="3_volcano-stype1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_volcano-stype2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "High:", gxj[1],"Low:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+  
  scale_color_manual(values =c("High" = "#8B0000", "NO" = "gray","Low" = "#008080"))
dev.off()

temp2 = data.frame(log10(rowMeans(eSet)))
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT),rownames(temp2)),] 
temp2 = data.frame(temp2,tT[,c("logFC","expression")])
colnames(temp2)=c("log10Exp","logFC","expression")
temp2$expression=factor(temp2$expression, levels=c("High","Low","NO"), order=T)
P_MA=ggplot(temp2,aes(x=temp2$log10Exp,y=temp2$logFC))+
  geom_point(alpha=0.8, size=6, 
             aes(color=temp2$expression))+
  scale_color_manual(values =c("High" = "#FF4500", "NO" = "gray","Low" = "#483D8B"))+
  labs(x="log10Exp",y=paste("log2FC", "High:", gxj[1],"Low:", gxj[2]))+
  # xlim(0.75,1.15)+ 。
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
tiff(file="4_MA.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)  #也注意投置片的大小
print(P_MA)
dev.off()

dim(eSet) 
#fix(eSet)
DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),] 
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_all-gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset, color = colorRampPalette(c('blue2',  'white',  'red'))(10000), 
                   cluster_rows = T, cluster_cols=T,                  
                   scale ="row" , 
                   show_rownames=F,show_colnames=F, 
                   treeheight_row = 0,  
                   clustering_distance_rows = "correlation",  
                   clustering_distance_cols = "euclidean",  
                   border_color = "grey",
                   clustering_method = "average")
dev.off()
target_gene<-row.names(deg)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target-gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap::pheatmap(deg_eset_2, color = colorRampPalette(c('#000066',  'white',  '#FF0000'))(60), 
                   cluster_rows = T, #clustering_distance_rows = "correlation", 
                   cluster_cols=T,   clustering_distance_cols = "euclidean",                
                   scale ="row" , 
                   show_rownames=T, 
                   show_colnames= F, 
                   treeheight_row = 0,  
                   border_color = "White",
                  clustering_method = "average")
dev.off()



#(8) TCGA-GEO-uni cox##########################################################################################
rm(list=ls()) 
pFilter=0.05       
surv_1=read.table("0_GSE14520_OS_survival_inf.txt",header=T,sep="\t",check.names=F,row.names=1)   
dim(surv_1)
ID_full<-rownames(surv_1)
surv_2<-data.frame(ID_full, surv_1) 
dim(surv_2)  
exp_1=read.table("0_GSE14520_batch.txt",header=T,sep="\t",check.names=F,row.names=1)     
dim(exp_1)
exp_2<-t(exp_1)
dim(exp_2)  
ID_full_2<-rownames(exp_2)
exp_3<-data.frame(ID_full_2, exp_2)  
dim(exp_3)  
exp_surv_1<-merge(surv_2,exp_3,  by.x="ID_full", by.y = "ID_full_2") 
dim(exp_surv_1)    
write.table(exp_surv_1, file="1_GSE14520_OS_exp_surv.txt", quote=FALSE, row.names=T, sep="\t")  
gxj<-read.table("1_GSE14520_OS_exp_surv.txt",header=T,sep="\t",check.names=F,row.names=1)    
dim(gxj) 
rownames(gxj)<-gxj$ID_full 
gxj_1<-gxj[,-1] 
dim(gxj_1)
outTab=data.frame()
sigGenes=c("month_os","status_os")
for(i in colnames(gxj_1[,3:ncol(gxj_1)])){
  cox <- coxph(Surv(month_os, status_os) ~ gxj_1[,i], data = gxj_1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )  }}
dim(outTab)
table(is.na(outTab)) 
outTab[!complete.cases(outTab),] 
outTab<- na.omit(outTab) 
dim(outTab)  
file_os<-paste("2_GAS14520_OS_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene.txt")
write.table(outTab,file=file_os,sep="\t",row.names=F,quote=F)
uniSigExp=gxj_1[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
dim(uniSigExp)
write.table(uniSigExp,file=paste("3_GAS14520_OS_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene_diff_gene_surv_exp.txt"),sep="\t",row.names=F,quote=F)
rt <- read.table(file_os,header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
tiff(file=paste("4_OS_plot_forest_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-double-colorok.tiff"),width =60,height =120,units ="cm",compression="lzw",bg="white",res=300)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8  
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=5,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=5,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab=paste("Hazard ratio_plot_gene：",dim(outTab)[1]))
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5) 
#arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')  
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.8)  
axis(1)
dev.off()
library(forestplot)  
options(forestplot_new_page = FALSE)  
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")              
rt=read.table(file_os, header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file=paste("4_OS_forestplot_forest_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene.tiff"),width =28,height = 42,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,  
           xlab=paste("Hazard ratio_foresplot_gene：",dim(outTab)[1]))
dev.off()
surv_2=read.table("0_GSE14520_RFS_survival_inf.txt",header=T,sep="\t",check.names=F,row.names=1)  
dim(surv_2)
ID_full<-rownames(surv_2)
surv_3<-data.frame(ID_full, surv_2) 
dim(surv_3)  
exp_1=read.table("0_GSE14520_batch.txt",header=T,sep="\t",check.names=F,row.names=1)    
dim(exp_1)
exp_2<-t(exp_1)
dim(exp_2)  
ID_full_2<-rownames(exp_2)
exp_3<-data.frame(ID_full_2, exp_2)  
exp_surv_RFS<-merge(surv_3,exp_3,  by.x="ID_full", by.y = "ID_full_2") 
dim(exp_surv_RFS)    
write.table(exp_surv_RFS, file="1_GAS14520_RFS_exp_surv.txt", quote=FALSE, row.names=T, sep="\t")   
gxj<-read.table("1_GAS14520_RFS_exp_surv.txt",header=T,sep="\t",check.names=F,row.names=1)    
dim(gxj)
rownames(gxj)<-gxj$ID_full 
gxj_1<-gxj[,-1]  
dim(gxj_1)
outTab=data.frame()
sigGenes=c("month_rfs","status_rfs")
for(i in colnames(gxj_1[,3:ncol(gxj_1)])){
  cox <- coxph(Surv(month_rfs, status_rfs) ~ gxj_1[,i], data = gxj_1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )  }}
dim(outTab)
table(is.na(outTab)) 
outTab[!complete.cases(outTab),] 
outTab<- na.omit(outTab) 
dim(outTab)  
file_rfs<-paste("2_GAS14520_rfs_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene.txt")
write.table(outTab,file=file_rfs,sep="\t",row.names=F,quote=F)
uniSigExp=gxj_1[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
dim(uniSigExp) 
write.table(uniSigExp,file=paste("3_GAS14520_rfs_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene_diff_gene_surv_exp.txt"),sep="\t",row.names=F,quote=F)
rt <- read.table(file_rfs,header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
tiff(file=paste("5_RFS_plot_forest_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-double-colorOK.tiff"),width =60,height =120,units ="cm",compression="lzw",bg="white",res=300)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8  
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=5,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=5,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab=paste("Hazard ratio_plot_gene：",dim(outTab)[1]))
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
#arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')  
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.8)  
axis(1)
dev.off()
————————————————
library(forestplot)  
options(forestplot_new_page = FALSE)  
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")             
rt=read.table(file_rfs, header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )         
tiff(file=paste("5_RFS_forestplot_forest_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene.tiff"),width =28,height = 42,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,  
           xlab=paste("Hazard ratio_foresplot_gene：",dim(outTab)[1]))
dev.off()
rm(list=ls()) 
library(survival)
surv_1=read.table("0-TCGA-LIHC-OS-survival_inf.txt",header=T,sep="\t",check.names=F,row.names=1)    
dim(surv_1)
ID_full<-rownames(surv_1)
surv_2<-data.frame(ID_full, surv_1)  
dim(surv_2)  
exp_1=read.table("0_TCGA-LIHC_batch.txt",header=T,sep="\t",check.names=F,row.names=1)   
dim(exp_1)
exp_2<-t(exp_1)
dim(exp_2)  
ID_full_2<-rownames(exp_2)
exp_3<-data.frame(ID_full_2, exp_2) 
dim(exp_3)  
exp_surv_1<-merge(surv_2,exp_3,  by.x="ID_full", by.y = "ID_full_2") 
dim(exp_surv_1)    
write.table(exp_surv_1, file="1_TCGA-LIHC_exp_surv.txt", quote=FALSE, row.names=T, sep="\t")  
gxj<-read.table("1_TCGA-LIHC_exp_surv.txt",header=T,sep="\t",check.names=F,row.names=1)    
dim(gxj)
rownames(gxj)<-gxj$ID_full  
gxj_1<-gxj[,-1]  
dim(gxj_1)
pFilter=0.05       
outTab=data.frame()
sigGenes=c("time","status")
for(i in colnames(gxj_1[,3:ncol(gxj_1)])){
  cox <- coxph(Surv(time, status) ~ gxj_1[,i], data = gxj_1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )  }}
dim(outTab)
table(is.na(outTab)) 
outTab[!complete.cases(outTab),] 
outTab<- na.omit(outTab) 
dim(outTab)  
write.table(outTab,file=paste("2_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-OK.txt"),sep="\t",row.names=F,quote=F)
uniSigExp=gxj_1[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
#fix(uniSigExp)
dim(uniSigExp)  
write.table(uniSigExp,file=paste("3_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene_diff_gene_surv_exp-OK.txt"),sep="\t",row.names=F,quote=F)
rt <- read.table(paste("2_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-OK.txt"),header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
tiff(file=paste("4_plot_forest_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-OK-double-col-OK.tiff"),width =60,height =120,units ="cm",compression="lzw",bg="white",res=300)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8 
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=5,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=5,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab=paste("Hazard ratio_plot_gene：",dim(outTab)[1]))
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5) 
#arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')  
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.8) 
axis(1)
dev.off()
————————————————
library(forestplot)  
options(forestplot_new_page = FALSE)  
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")              
rt=read.table(paste("2_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-OK.txt"),header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )         
tiff(file=paste("5_forestplot_forest_TCGA-LIHC_UniCox_P_",pFilter,"_",dim(outTab)[1],"prog_gene-OK.tiff"),width =18,height = 42,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4, 
           xlab=paste("Hazard ratio_foresplot_gene：",dim(outTab)[1]))
dev.off()





#(9) Lasso modeling############################################################################################
library("glmnet")
library(survival)
lasso_rt_1=read.table("0_TCGA-LIHC_UniCox_P_ 0.05_prog_gene-exp.txt",header=T,sep="\t",row.names=1)    
dim(lasso_rt_1)
#fix(lasso_rt_1)
gene1<-read.table("0_TCGA_Tumor_vs._Normal_diff_gene.txt",sep='\t',quote="", comment.char="", check.names=FALSE,row.names=1) 
dim(gene1)
#fix(gene1) 
a<-row.names(gene1)
gene2<-read.table("0_TCGA-LIHC_UniCox_P_0.05_496_prog_gene.txt",sep='\t',quote="", head=1, comment.char="", check.names=FALSE,row.names=1)  
dim(gene2)
#fix(gene2)  
b<-row.names(gene2)
#fix(Up_gene3) 
common_genes_2<-Reduce(intersect,list(a, b)) 
length(common_genes_2)  
write.table(common_genes_2, file = "1_TCGA-unicox-TCGAdiff_gene-com-gene.txt",sep="\t", row.names =F, quote = F)
library(venn)
library(ggplot2)
tiff(file="1_venn_TCGA-unicox-TCGAdiff_gene-com-gene.tiff", width =35,height = 35,units ="cm",compression="lzw",bg="white",res=300)
venn(x=list(TCGA_dif=a, TCGA_unicox=b),
     #zcolor = "style",
     zcolor = c("#FFCC33", "#00FFFF"),
     opacity = 0.5, 
     ilcs = 2.6,  
     sncs = 1.85, 
     borders = F, box = F, ggplot = F 
     #main="common_UP_genes", main.cex = 2, main.fontface = 2 
)
dev.off()
#1）GEOdiff_gene
gene3<-read.table("0_GEO_Tumor_vs._Normal_diff_gene-un-used.txt",sep='\t',quote="", comment.char="", check.names=FALSE,row.names=1) 
dim(gene3)
#fix(gene3) 
c<-row.names(gene3)
common_genes_3<-Reduce(intersect,list(a, b)) 
length(common_genes_3) 
write.table(common_genes_3, file = "1_TCGA-unicox-TCGA-GEOdiff_gene-un-used.txt",sep="\t", row.names =F, quote = F)
library(venn)
library(ggplot2)
tiff(file="1_venn_TCGA-unicox-TCGA-GEOdiff_gene.tiff", width =35,height = 35,units ="cm",compression="lzw",bg="white",res=300)
venn(x=list(TCGA_dif=a, TCGA_unicox=b, GEO_dif=c ),
     #zcolor = "style",
     zcolor = c("#FFCC33", "#00FFFF", "#9900FF"), 
     opacity = 0.5,
     ilcs = 2.6,  
     sncs = 1.85, 
     borders = F, box = F, ggplot = F 
     #main="common_UP_genes", main.cex = 2, main.fontface = 2
)
dev.off()
com_1=read.table("000_COM-used.txt",header=F,sep="\t",row.names=1)    
# fix(com_1)
dim(com_1)
common_genes_5<-row.names(com_1)
common_genes_5
lasso_rt<-lasso_rt_1[,c("time", "status",common_genes_5)]  
colnames(lasso_rt)[c(1:2)]=c("os_time", "os_stat")  
#fix(lasso_rt)
dim(lasso_rt)
lasso_rt$os_time=lasso_rt$os_time/365  
#fix(lasso_rt)
dim(lasso_rt)
#(1)
#for (a in names(lasso_rt)[c(3:ncol(lasso_rt))]){ lasso_rt[,a] <- as.factor(lasso_rt[,a])}
x=as.matrix(lasso_rt[,c(3:ncol(lasso_rt))])  
head(x)
y=data.matrix(Surv(lasso_rt$os_time,lasso_rt$os_stat))
head(y)
fit=glmnet(x, y, family = "cox")  
cvfit=cv.glmnet(x, y, family="cox")
#(2) 
tiff(file="2_lasso_cvfit.tiff",width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")  
dev.off()
tiff(file="2_lasso_lambda.tiff",width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300)
plot(fit, xvar = "lambda", label = TRUE)  
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
#fix(geneCoef)
write.table(geneCoef,file="2_lasso_geneCoef.txt",sep="\t",quote=F,row.names=F)
trainFinalGeneExp=lasso_rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("os_time","os_stat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))  
outTab=cbind(lasso_rt[,outCol],riskScore=as.vector(trainScore),risk)
final_outTab<-cbind(id=rownames(outTab),outTab)
##fix(final_outTab)
dim(final_outTab)
write.table(final_outTab,file="3_TCGA_training_risk.txt",sep="\t",quote=F,row.names=F)
surv_1=read.table("0_GSE14520_OS_surv-inf（year）.txt",header=T,sep="\t",check.names=F,row.names=1)   
dim(surv_1)
ID<-rownames(surv_1)
surv_2<-data.frame(ID, surv_1)
dim(surv_2)  
exp_1=read.table("0_GSE14520_batch.txt",header=T,sep="\t",check.names=F,row.names=1)     
dim(exp_1)
exp_2<-t(exp_1)
dim(exp_2) 
ID_2<-rownames(exp_2)
exp_3<-data.frame(ID_2, exp_2) 
dim(exp_3)  
exp_surv_1<-merge(surv_2,exp_3,  by.x="ID", by.y = "ID_2") 
dim(exp_surv_1)    
##fix(exp_surv_1)  
write.table(exp_surv_1, file="4_GSE14520_exp_OSsurv.txt", quote=FALSE, row.names=F, sep="\t")   
surv_3=read.table("0_GSE14520_RFS_surv-inf（year）.txt",header=T,sep="\t",check.names=F,row.names=1)     
dim(surv_3)
ID_3<-rownames(surv_3)
surv_4<-data.frame(ID_3, surv_3) 
dim(surv_4)  
##fix(surv_4)
exp_surv_2<-merge(surv_4,exp_3,  by.x="ID_3", by.y = "ID_2")  
dim(exp_surv_2)    
##fix(exp_surv_2)  
write.table(exp_surv_2, file="5_GSE14520_exp_RFSsurv.txt", quote=FALSE, row.names=F, sep="\t")  
#（1）OS
geo_os_1=read.table("4_GSE14520_exp_OSsurv.txt",header=T,sep="\t",row.names=1)
##fix(geo_os_1)
dim(geo_os_1)  
colnames(geo_os_1)[c(1:2)]=c("os_time", "os_stat")  
dim(geo_os_1)
geo_os<-na.omit(geo_os_1)  
dim(geo_os)
##fix(geo_os)
geo_os[,3:ncol(geo_os)][geo_os[,3:ncol(geo_os)]<0]=0  
#geo_os$os_time=geo_os$os_time/12 
length(lassoGene)
lassoGene[1:6]  
testFinalGeneExp=geo_os[,lassoGene] 
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("os_time","os_stat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low")) 
outTab=cbind(geo_os[,outCol],riskScore=as.vector(testScore),risk)
outTab_1<-cbind(id=rownames(outTab),outTab)
dim(outTab_1)
write.table(outTab_1,file="6_GSE14520_OS_valid--riskscore.txt",sep="\t",quote=F,row.names=F)
#（2）RFS
geo_rfs_1=read.table("5_GSE14520_exp_RFSsurv.txt",header=T,sep="\t",row.names=1)
#fix(geo_rfs_1)
dim(geo_rfs_1)  
colnames(geo_rfs_1)[c(1:2)]=c("rfs_time","rfs_stat" )  
#fix(geo_rfs_1)
dim(geo_rfs_1)
geo_rfs<-na.omit(geo_rfs_1)  
dim(geo_rfs)
#fix(geo_rfs)
geo_rfs[,3:ncol(geo_rfs)][geo_rfs[,3:ncol(geo_rfs)]<0]=0  
#geo_rfs$os_time=geo_rfs$os_time/12   
lassoGene[1:6] 
summary(lassoGene)
testFinalGeneExp=geo_rfs[,lassoGene] 
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("rfs_time","rfs_stat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low")) 
outTab=cbind(geo_rfs[,outCol],riskScore=as.vector(testScore),risk)
outTab_1<-cbind(id=rownames(outTab),outTab)
#fix(outTab_1)
dim(outTab_1)
write.table(outTab_1,file="7_GSE14520_RFS_valid--riskscore.txt",sep="\t",quote=F,row.names=F)
library(survival)
library("survminer")
TCGA_rt=read.table("3_TCGA_training_risk.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(TCGA_rt)
dim(TCGA_rt)
#surv curve 
diff=survdiff(Surv(os_time, os_stat) ~risk,data = TCGA_rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(os_time, os_stat) ~ risk, data = TCGA_rt)
tiff(file="8_TCGA_risk_surv curve.tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
surPlot=ggsurvplot(fit, 
                   data=TCGA_rt,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="TCGA_Risk",
                   xlab="os_time(years)",
                   break.os_time.by = 1,
                   risk.table.title="",
                   palette=c("#FF00FF", "#008080"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
#risk curve
TCGA_rt=TCGA_rt[order(TCGA_rt$riskScore),]   
riskClass=TCGA_rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=TCGA_rt[,"riskScore"]
line[line>10]=10
tiff(file="9_TCGA_risk curve.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="TCGA_Risk score",
     col=c(rep("#008080",lowLength),rep("#FF00FF",highLength)) )
abline(h=median(TCGA_rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("#FF00FF","#008080"),cex=1.2)
dev.off()
#3：surv
color=as.vector(TCGA_rt$os_stat)
color[color==1]="#4169E1"
color[color==0]="#CD5C5C"
tiff(file="10_TCGA_survstatus.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
plot(TCGA_rt$os_time, pch=19,
     xlab="Patients (increasing risk socre)", ylab="TCGA_Survival os_time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#4169E1","#CD5C5C"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()
#4：risk
library(pheatmap)
TCGA_rt1=TCGA_rt[c(3:(ncol(TCGA_rt)-2))]
TCGA_rt1=log2(TCGA_rt1+1)
TCGA_rt1=t(TCGA_rt1)
##fix(TCGA_rt1)
dim(TCGA_rt1) 
annotation=data.frame(type=TCGA_rt[,ncol(TCGA_rt)])
rownames(annotation)=rownames(TCGA_rt)
##fix(annotation)
tiff(file="11_TCGA_heatmap.tiff",width =30,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
pheatmap(TCGA_rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         main="TCGA",
         color = colorRampPalette(c("#32CD32", "white", "#CC0066"))(200) )
dev.off()
#(2)GEO-OS 
GEO_os_rt=read.table("6_GSE14520_OS_valid--riskscore.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(GEO_os_rt)
dim(GEO_os_rt)
#1： surv curve 
diff=survdiff(Surv(os_time, os_stat) ~risk,data = GEO_os_rt) 
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(os_time, os_stat) ~ risk, data = GEO_os_rt)
tiff(file="12_GSE_OSrisk_surv curve.tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
surPlot=ggsurvplot(fit, 
                   data=GEO_os_rt,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="GEO_OS_Risk",
                   xlab="os_time(years)",
                   break.os_time.by = 1,
                   risk.table.title="",
                   palette=c("#FF00FF", "#008080"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
#2： risk curve
GEO_os_rt=GEO_os_rt[order(GEO_os_rt$riskScore),]    #按照riskScore样品排序
riskClass=GEO_os_rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=GEO_os_rt[,"riskScore"]
line[line>10]=10
tiff(file="13_GSE_OSrisk curve.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="GEO_OS_Risk score",
     col=c(rep("#008080",lowLength),rep("#FF00FF",highLength)) )
abline(h=median(GEO_os_rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("#FF00FF","#008080"),cex=1.2)
dev.off()
#3： survstatus
color=as.vector(GEO_os_rt$os_stat)
color[color==1]="#4169E1"
color[color==0]="#CD5C5C"
tiff(file="14_GSE_OSrisksurvstatus.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
plot(GEO_os_rt$os_time, pch=19,
     xlab="Patients (increasing risk socre)", ylab="GEO_OS_Survival os_time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#4169E1","#CD5C5C"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()
#4： riskheatmap
library(pheatmap)
GEO_os_rt1=GEO_os_rt[c(3:(ncol(GEO_os_rt)-2))]
GEO_os_rt1=log2(GEO_os_rt1+1)
GEO_os_rt1=t(GEO_os_rt1)
annotation=data.frame(type=GEO_os_rt[,ncol(GEO_os_rt)])
rownames(annotation)=rownames(GEO_os_rt)
tiff(file="15_GSE_OSriskheatmap.tiff",width =30,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
pheatmap(GEO_os_rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         main="GEO_OS",
         color = colorRampPalette(c("#32CD32", "white", "#CC0066"))(200) )
dev.off()
#(3)GEO-RFS
GEO_rfs_rt=read.table("7_GSE14520_RFS_valid--riskscore.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(GEO_rfs_rt)
dim(GEO_rfs_rt)
#1： surv curve 
diff=survdiff(Surv(rfs_time, rfs_stat) ~risk,data = GEO_rfs_rt) 
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(rfs_time, rfs_stat) ~ risk, data = GEO_rfs_rt)
tiff(file="16_GSE_RFSrisk_surv curve.tiff",width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
surPlot=ggsurvplot(fit, 
                   data=GEO_rfs_rt,
                   conf.int=TRUE,
                   pval=paste0("p=",pValue),
                   pval.size=4,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="GEO_RFS_Risk",
                   xlab="rfs_time(years)",
                   break.rfs_time.by = 1,
                   risk.table.title="",
                   palette=c("#FF00FF", "#008080"),
                   risk.table.height=.25)
print(surPlot)
dev.off()
#2： risk curve
GEO_rfs_rt=GEO_rfs_rt[order(GEO_rfs_rt$riskScore),]    
riskClass=GEO_rfs_rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=GEO_rfs_rt[,"riskScore"]
line[line>10]=10
tiff(file="17_GSE_RFSrisk curve.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="GEO_RFS_Risk score",
     col=c(rep("#008080",lowLength),rep("#FF00FF",highLength)) )
abline(h=median(GEO_rfs_rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("#FF00FF","#008080"),cex=1.2)
dev.off()
#3： survstatus
color=as.vector(GEO_rfs_rt$rfs_stat)
color[color==1]="#4169E1"
color[color==0]="#CD5C5C"
tiff(file="18_GSE_RFSrisksurvstatus.tiff",width =25,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
plot(GEO_rfs_rt$rfs_time, pch=19,
     xlab="Patients (increasing risk socre)", ylab="GEO_RFS_Survival rfs_time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#4169E1","#CD5C5C"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()
#4： riskheatmap
library(pheatmap)
GEO_rfs_rt1=GEO_rfs_rt[c(3:(ncol(GEO_rfs_rt)-2))]
GEO_rfs_rt1=log2(GEO_rfs_rt1+1)
GEO_rfs_rt1=t(GEO_rfs_rt1)
annotation=data.frame(type=GEO_rfs_rt[,ncol(GEO_rfs_rt)])
rownames(annotation)=rownames(GEO_rfs_rt)
tiff(file="19_GSE_RFSriskheatmap.tiff",width =30,height = 15,units ="cm",compression="lzw",bg="white",res=300) 
pheatmap(GEO_rfs_rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         main="GEO_RFS",
         color = colorRampPalette(c("#32CD32", "white", "#CC0066"))(200) )
dev.off()





#(10) uni/multi cox survival -KM##################################################################################################
library(survival)
risk=read.table("0_TCGA_ training-riskscore.txt",header=T,sep="\t",check.names=F,row.names=1)       
dim(risk)
#fix(risk)
cli=read.table("0_TCGA-LIHC-surv-sum.txt",sep="\t",check.names=F,header=T,row.names=1)     
dim(cli)
#fix(cli)
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
TCGA_rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dim(TCGA_rt)
#fix(TCGA_rt)
write.table(TCGA_rt,file="1_TCGA-LIHC_ training-riskscore_surv-sum.txt",sep="\t",quote=F)
#（1）unicox
uniTab=data.frame()
for(i in colnames(TCGA_rt[,3:ncol(TCGA_rt)])){
  cox <- coxph(Surv(futime, fustat) ~ TCGA_rt[,i], data = TCGA_rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])  )}
#fix(uniTab)
write.table(uniTab,file="2_TCGA_uniCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_1 <- read.table("2_TCGA_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_1)
hr <- sprintf("%.3f",plot_rt_1$"HR")
hrLow  <- sprintf("%.3f",plot_rt_1$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_1$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_1$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_1$pvalue))
tiff(file="3_TCGA_plot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="1_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_1)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="uni-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5) 
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
library(forestplot)   
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")              
forestplot_rt_1=read.table("2_TCGA_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_1)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="4_TCGA_forestplot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="uni-cox-Hazard ratio-forestplot")
dev.off()
#（1）multi
multiCox=coxph(Surv(futime, fustat) ~ ., data = TCGA_rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
#fix(multiTab)
write.table(multiTab,file="5_TCGA_multiCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_2 <- read.table("5_TCGA_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_2)
hr <- sprintf("%.3f",plot_rt_2$"HR")
hrLow  <- sprintf("%.3f",plot_rt_2$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_2$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_2$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_2$pvalue))
 tiff(file="6_TCGA_plot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="3_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_2)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="multi-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
#install.packages('forestplot')
library(forestplot)   
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")             
forestplot_rt_2=read.table("5_TCGA_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_2)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="7_TCGA_forestplot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="multi-cox-Hazard ratio-forestplot")
dev.off()
#install.packages('survival')
risk=read.table("0_GSE14520_OS_ valid-riskscore.txt",header=T,sep="\t",check.names=F,row.names=1)        
dim(risk)
#fix(risk)
cli=read.table("0_GSE14520_treat-surv-sum.txt",sep="\t",check.names=F,header=T,row.names=1)     
dim(cli)
#fix(cli)
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
GEO_rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dim(GEO_rt)
#fix(GEO_rt)
write.table(GEO_rt,file="1_GSE14520_ valid-riskscore_OSsurv-sum.txt",sep="\t",quote=F)
#（1）unitreat
uniTab=data.frame()
for(i in colnames(GEO_rt[,3:ncol(GEO_rt)])){
  cox <- coxph(Surv(futime, fustat) ~ GEO_rt[,i], data = GEO_rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])  )}
#fix(uniTab)
write.table(uniTab,file="2_GEO_OS_uniCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_1 <- read.table("2_GEO_OS_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_1)
hr <- sprintf("%.3f",plot_rt_1$"HR")
hrLow  <- sprintf("%.3f",plot_rt_1$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_1$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_1$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_1$pvalue))
tiff(file="3_GEO_OS_plot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="1_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_1)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
 xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="uni-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
 #install.packages('forestplot')
library(forestplot)   
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")              
forestplot_rt_1=read.table("2_GEO_OS_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_1)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="4_GEO_OS_forestplot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="uni-cox-Hazard ratio-forestplot")
dev.off()
#multicox
#（1）multi
multiCox=coxph(Surv(futime, fustat) ~ ., data = GEO_rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
#fix(multiTab)
write.table(multiTab,file="5_GEO_OS_multiCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_2 <- read.table("5_GEO_OS_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_2)
hr <- sprintf("%.3f",plot_rt_2$"HR")
hrLow  <- sprintf("%.3f",plot_rt_2$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_2$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_2$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_2$pvalue))
 tiff(file="6_GEO_OS_plot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="3_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_2)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
 xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="multi-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")             
forestplot_rt_2=read.table("5_GEO_OS_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_2)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="7_GEO_OS_forestplot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="multi-cox-Hazard ratio-forestplot")
dev.off()
#GEO-RFS
risk=read.table("0_GSE14520_RFS_valid-riskscore.txt",header=T,sep="\t",check.names=F,row.names=1)        
dim(risk)
#fix(risk)
cli=read.table("0_GSE14520_treat-surv-sum.txt",sep="\t",check.names=F,header=T,row.names=1)     
dim(cli)
#fix(cli)
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
GEO_rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dim(GEO_rt)
#fix(GEO_rt)
write.table(GEO_rt,file="1_GSE14520_ valid-riskscore_RFSsurv-sum.txt",sep="\t",quote=F)
#（1）uni-COX
uniTab=data.frame()
for(i in colnames(GEO_rt[,3:ncol(GEO_rt)])){
  cox <- coxph(Surv(futime, fustat) ~ GEO_rt[,i], data = GEO_rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])  )}
write.table(uniTab,file="2_GEO_RFS_uniCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_1 <- read.table("2_GEO_RFS_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_1)
hr <- sprintf("%.3f",plot_rt_1$"HR")
hrLow  <- sprintf("%.3f",plot_rt_1$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_1$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_1$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_1$pvalue))
 tiff(file="3_GEO_RFS_plot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="1_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_1)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="uni-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
 #install.packages('forestplot')
library(forestplot)   
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")              
forestplot_rt_1=read.table("2_GEO_RFS_uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_1)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="4_GEO_RFS_forestplot_unicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="uni-cox-Hazard ratio-forestplot")
dev.off()
#（1）multi-cox
multiCox=coxph(Surv(futime, fustat) ~ ., data = GEO_rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
#fix(multiTab)
write.table(multiTab,file="5_GEO_RFS_multiCox.txt",sep="\t",row.names=F,quote=F)
plot_rt_2 <- read.table("5_GEO_RFS_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(plot_rt_2)
hr <- sprintf("%.3f",plot_rt_2$"HR")
hrLow  <- sprintf("%.3f",plot_rt_2$"HR.95L")
hrHigh <- sprintf("%.3f",plot_rt_2$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(plot_rt_2$pvalue<0.001, "<0.001", sprintf("%.3f", plot_rt_2$pvalue))
tiff(file="6_GEO_RFS_plot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
#pdf(file="3_unicox_forest-plot.pdf", width = 6,height = 4.5)
n <- nrow(plot_rt_2)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#forest-plot
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="multi-cox-Hazard ratio-plot")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#336666",lwd=5)  
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#CC3399', '#33CCCC')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=3)
axis(1)
dev.off()
 #install.packages('forestplot')
library(forestplot)   
options(forestplot_new_page = FALSE) 
clrs <- fpColors(box="#33CCCC",line="#990099", summary="royalblue")             
forestplot_rt_2=read.table("5_GEO_RFS_multiCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(forestplot_rt_2)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
tiff(file="7_GEO_RFS_forestplot_multicox_forest-plot.tiff",width =18,height = 13,units ="cm",compression="lzw",bg="white",res=300)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(60, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.25,  
           xlab="multi-cox-Hazard ratio-forestplot")
dev.off()
library(survival)
library("survminer")
library("ggplot2")
TCGA_rt=read.table("3_TCGA_training-riskscore.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(TCGA_rt)
dim(TCGA_rt)
colnames(TCGA_rt)[1:2]<-c("time", "stat")  
#fix(TCGA_rt)
colnames(TCGA_rt)
for(gene in colnames(TCGA_rt)[3:(ncol(TCGA_rt)-2)]) 
{
  group <- ifelse(TCGA_rt[[gene]] > median(TCGA_rt[[gene]]),   "high", "low") 
   diff=survdiff(Surv(time, stat) ~ group, data = TCGA_rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(time, stat) ~ group, data = TCGA_rt)
  tiff(file=paste0("1_TCGA_",gene,".tiff"),width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
  surPlot=ggsurvplot(fit,  data=TCGA_rt,  conf.int=TRUE,  pval=paste0("p=",pValue),  pval.size=4,
                     risk.table=TRUE, legend.labs=c("high expression", "low expression"),
                     legend.title="TCGA_Risk",   xlab=paste0("os_time(years)_",gene),
                     break.os_time.by = 1,   risk.table.title="",
                     palette=c("#CD5C5C", "#2E8B57"),  risk.table.height=.25)
  print(surPlot)
  dev.off()
}
#GEO-OS
GEO_OS_rt=read.table("6_GSE14520_OS_valid-riskscore.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(GEO_OS_rt)
dim(GEO_OS_rt)
colnames(GEO_OS_rt)[1:2]<-c("time", "stat")  
#fix(GEO_OS_rt)
colnames(GEO_OS_rt)
for(gene in colnames(GEO_OS_rt)[3:(ncol(GEO_OS_rt)-2)])
{
  group <- ifelse(GEO_OS_rt[[gene]] > median(GEO_OS_rt[[gene]]),   "high", "low")
  diff=survdiff(Surv(time, stat) ~ group, data = GEO_OS_rt)  
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(time, stat) ~ group, data = GEO_OS_rt)
  tiff(file=paste0("2_GEO_OS_",gene,".tiff"),width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
  surPlot=ggsurvplot(fit,  data=GEO_OS_rt,  conf.int=TRUE,  pval=paste0("p=",pValue),  pval.size=4,
                     risk.table=TRUE, legend.labs=c("high expression", "low expression"),
                     legend.title="GEO_OS_Risk",   xlab=paste0("os_time(years)_",gene),
                     break.os_time.by = 1,   risk.table.title="",
                     palette=c("#CD5C5C", "#2E8B57"),  risk.table.height=.25)
  print(surPlot)
  dev.off()
}
#GEO-RFS
GEO_RFS_rt=read.table("7_GSE14520_RFS_valid-riskscore.txt",sep="\t",header=T,row.names=1,check.names=F) 
#fix(GEO_RFS_rt)
dim(GEO_RFS_rt)
colnames(GEO_RFS_rt)[1:2]<-c("time", "stat")  
#fix(GEO_RFS_rt)
colnames(GEO_RFS_rt)
for(gene in colnames(GEO_RFS_rt)[3:(ncol(GEO_RFS_rt)-2)])
{
  group <- ifelse(GEO_RFS_rt[[gene]] > median(GEO_RFS_rt[[gene]]),   "high", "low")
  diff=survdiff(Surv(time, stat) ~ group, data = GEO_RFS_rt)  
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(time, stat) ~ group, data = GEO_RFS_rt)
  tiff(file=paste0("3_GEO_RFS_",gene,".tiff"),width =15,height = 13,units ="cm",compression="lzw",bg="white",res=300) 
  surPlot=ggsurvplot(fit,  data=GEO_RFS_rt,  conf.int=TRUE,  pval=paste0("p=",pValue),  pval.size=4,
                     risk.table=TRUE, legend.labs=c("high expression", "low expression"),
                     legend.title="GEO_RFS_Risk",   xlab=paste0("os_time(years)_",gene),
                     break.os_time.by = 1,   risk.table.title="",
                     palette=c("#CD5C5C", "#2E8B57"),  risk.table.height=.25)
  print(surPlot)
  dev.off()
}









#(11) ROC######################################################################################################

#install.packages("survivalROC")
library(survivalROC)
bioROC=function(riskFile=null,cliFile=null,outFile=null){
  risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)         
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)           
  sameSample=intersect(row.names(cli),row.names(risk))    
  risk=risk[sameSample,]   
  cli=cli[sameSample,]    
  rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])   
  rocCol=rainbow(ncol(rt)-2)    
  aucText=c()
    tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  #pdf(file=outFile,width=6,height=6)
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")  
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate（1year）", ylab="True positive rate",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))  
  abline(0,1)
  # 
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-1)])){
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")   
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)   
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",cliFile="0_TCGA-LIHC-surv-sum.txt",outFile="1_TCGA_LIHC_ROC（1year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="2_GEO_OS_ROC（1year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="3_GEO_RFS_ROC（1year）.TIFF")

bioROC=function(riskFile=null,cliFile=null,outFile=null){
  risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)         
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)           
  sameSample=intersect(row.names(cli),row.names(risk))    
  risk=risk[sameSample,]   
  cli=cli[sameSample,]    
  rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])   
  rocCol=rainbow(ncol(rt)-2)    
  aucText=c()
    # 
  tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  #pdf(file=outFile,width=6,height=6)
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM") 
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate（3year）", ylab="True positive rate",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))  
  abline(0,1)
  # 
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-1)])){
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =3, method="KM") 
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)   
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",cliFile="0_TCGA-LIHC-surv-sum.txt",outFile="4_TCGA_LIHC_ROC（3year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="5_GEO_OS_ROC（3year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="6_GEO_RFS_ROC（3year）.TIFF")
bioROC=function(riskFile=null,cliFile=null,outFile=null){
  risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)         
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)           
  sameSample=intersect(row.names(cli),row.names(risk))    
  risk=risk[sameSample,]   
  cli=cli[sameSample,]    
  rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])   
  rocCol=rainbow(ncol(rt)-2)    
  aucText=c()
  tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  #pdf(file=outFile,width=6,height=6)
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")  
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate（5year）", ylab="True positive rate",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2, fg="black") 
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))  
  abline(0,1)
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-1)])){
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =5, method="KM")   
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)   
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",cliFile="0_TCGA-LIHC-surv-sum.txt",outFile="7_TCGA_LIHC_ROC（5year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="8_GEO_OS_ROC（5year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",cliFile="0_GSE14520_treat-surv-sum.txt",outFile="9_GEO_RFS_ROC（5year）.TIFF")
#install.packages("survivalROC")
#---------------------1year-km-surv
bioROC=function(riskFile=null,cliFile=null,outFile=null){
  rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        year-km-surv
  dim(rt)
  #Fix(rt)
  colnames(rt)[1:2]<-c("futime", "fustat") year-km-surv
  #fix(rt)
  tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  rocCol=rainbow(ncol(rt)-2)   
  aucText=c()
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")  
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate (1year)", ylab="True positive rate (KM) ",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")")) 
  abline(0,1)
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-2)])){  
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")   
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",outFile="1_TCGA_LIHC_ROC-sig-mul-genes（1year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",outFile="2_GSE14520_OS-sig-mul-genes（1year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",outFile="3_GSE14520_RFS-sig-mul-genes（1year）.TIFF")
#---------------------3year-km-surv
bioROC=function(riskFile=null,cliFile=null,outFile=null){
  rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        year-km-surv
  dim(rt)
  #Fix(rt)
  colnames(rt)[1:2]<-c("futime", "fustat") year-km-surv
  #fix(rt)
  tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  rocCol=rainbow(ncol(rt)-2)   
  aucText=c()
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM")  
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate (3year)", ylab="True positive rate (KM) ",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")")) 
  abline(0,1)
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-2)])){  
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =3, method="KM")   
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",outFile="4_TCGA_LIHC_ROC-sig-mul-genes（3year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",outFile="5_GSE14520_OS-sig-mul-genes（3year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",outFile="6_GSE14520_RFS-sig-mul-genes（3year）.TIFF")
#---------------------5year-km-surv
bioROC=function(riskFile=null,cliFile=null,outFile=null){
  rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        year-km-surv
  dim(rt)
  #Fix(rt)
  colnames(rt)[1:2]<-c("futime", "fustat") year-km-surv
  #fix(rt)
  tiff(file=outFile,width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
  rocCol=rainbow(ncol(rt)-2)   
  aucText=c()
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")  
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate (5year)", ylab="True positive rate (KM) ",
       lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")")) 
  abline(0,1)
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-2)])){  
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =5, method="KM")   
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
  dev.off()
}
bioROC(riskFile="3_TCGA_training-riskscore.txt",outFile="7_TCGA_LIHC_ROC-sig-mul-genes（5year）.TIFF")
bioROC(riskFile="6_GSE14520_OS_valid-riskscore.txt",outFile="8_GSE14520_OS-sig-mul-genes（5year）.TIFF")
bioROC(riskFile="7_GSE14520_RFS_valid-riskscore.txt",outFile="9_GSE14520_RFS-sig-mul-genes（5year）.TIFF")
library(survivalROC)
library(timeROC)
library(survival)
library(ggplot2)
##TCGA-OS
riskFile="0_TCGA_training-riskscore.txt"
rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1) 
#fix(rt)
colnames(rt)[1:2]<-c("futime", "fustat")
colnames(rt)
#-----------1：survivalROC
library(survivalROC)
rocCol=rainbow(ncol(rt)-1)   
aucText=c()
tiff(file="1_TCGA_OS_survivalROC_plot（1-3-5year）.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc<-survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")  #KM法或NNE法，预测1，3，5year
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="TCGA-OS-False positive rate", ylab="True positive rate (KM) ",
     lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score_",1," year","(AUC=",sprintf("%.3f",roc$AUC),")")) 
abline(0,1)
j=1
for(i in c(3,5)){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker =rt$riskScore, predict.time =i, method="KM")   
  j=j+1
  aucText=c(aucText,paste0("risk score_",i," year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
dev.off()
library(timeROC)
library(survival)
time_roc_res <- timeROC(T = rt$futime,delta = rt$fustat,  marker = rt$riskScore,  cause = 1,  weighting="marginal",
                        times = c(1, 3, 5),  ROC = TRUE,  iid = TRUE)  
time_roc_res$AUC   
a<-paste0("TCGA-OS-KM-1 year:", sprintf("%.3f", time_roc_res$AUC[[1]]))  
a
b<-paste0("3 year:", sprintf("%.3f", time_roc_res$AUC[[2]]))  
b
c<-paste0("5 year:", sprintf("%.3f", time_roc_res$AUC[[3]]))  
c
confint(time_roc_res, level = 0.95)$CI_AUC  
#time-dependent ROC
tiff(file="2_TCGA_OS_timeROC_plot（1-3-5year）-1.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(time_roc_res, time=1, col = "#FF00CC", title = FALSE, lwd = 5)  
plot(time_roc_res, time=3, add=TRUE, col="#4169E1", lwd = 5) 
plot(time_roc_res, time=5, add=TRUE, col="#EEB422", lwd = 5)
legend("bottomright",c(a,b,c),  col=c("#FF00CC", "#4169E1", "#EEB422"), lty=5, lwd=4, bty="n") 
dev.off()
tiff(file="2_TCGA_OS_timeROC_ggplot（1-3-5year）-2.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
time_ROC_df <- data.frame(  TP_3year = time_roc_res$TP[, 1],  FP_3year = time_roc_res$FP[, 1],
                            TP_5year = time_roc_res$TP[, 2],  FP_5year = time_roc_res$FP[, 2],
                            TP_10year = time_roc_res$TP[, 3],  FP_10year = time_roc_res$FP[, 3])
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 2, color = "#FF00CC") +   
  geom_line(aes(x = FP_5year, y = TP_5year), size = 2, color = "#4169E1") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 2, color = "#EEB422") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 3) +
  theme_bw() +
  annotate("text",   x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#FF00CC"  ) +
  annotate("text",   x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#4169E1"  ) +
  annotate("text",           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#EEB422"  ) +
  labs(x = "TCGA-OS-False positive rate", y = "True positive rate(KM)") +
  theme(    axis.text = element_text(face = "bold", size = 11, color = "black"),
            axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
            axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()
#GEO-OS
riskFile="0_GSE14520_OS_valid-riskscore.txt"
rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1) 
#fix(rt)
colnames(rt)[1:2]<-c("futime", "fustat")
colnames(rt)
library(survivalROC)
rocCol=rainbow(ncol(rt)-1)  
aucText=c()
tiff(file="3_GEO_OS_survivalROC_plot（1-3-5year）.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc<-survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")  #KM法或NNE法，预测1，3，5year
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="GEO-OS-False positive rate", ylab="True positive rate (KM) ",
     lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score_",1," year","(AUC=",sprintf("%.3f",roc$AUC),")")) 
abline(0,1)
j=1
for(i in c(3,5)){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker =rt$riskScore, predict.time =i, method="KM")   
  j=j+1
  aucText=c(aucText,paste0("risk score_",i," year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
dev.off()

library(timeROC)
library(survival)
time_roc_res <- timeROC(T = rt$futime,delta = rt$fustat,  marker = rt$riskScore,  cause = 1,  weighting="marginal",
                        times = c(1, 3, 5),  ROC = TRUE,  iid = TRUE)  
time_roc_res$AUC   
a<-paste0("GEO-OS-KM-1 year:", sprintf("%.3f", time_roc_res$AUC[[1]]))  
a
b<-paste0("3 year:", sprintf("%.3f", time_roc_res$AUC[[2]]))  
b
c<-paste0("5 year:", sprintf("%.3f", time_roc_res$AUC[[3]]))  
c
confint(time_roc_res, level = 0.95)$CI_AUC  
#time-dependent ROC
tiff(file="4_GEO_OS_timeROC_plot（1-3-5year）-1.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(time_roc_res, time=1, col = "#FF00CC", title = FALSE, lwd = 5)  
plot(time_roc_res, time=3, add=TRUE, col="#4169E1", lwd = 5) 
plot(time_roc_res, time=5, add=TRUE, col="#EEB422", lwd = 5)
legend("bottomright",c(a,b,c),  col=c("#FF00CC", "#4169E1", "#EEB422"), lty=5, lwd=4, bty="n") 
dev.off()
tiff(file="4_GEO_OS_timeROC_ggplot（1-3-5year）-2.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
time_ROC_df <- data.frame(  TP_3year = time_roc_res$TP[, 1],  FP_3year = time_roc_res$FP[, 1],
                            TP_5year = time_roc_res$TP[, 2],  FP_5year = time_roc_res$FP[, 2],
                            TP_10year = time_roc_res$TP[, 3],  FP_10year = time_roc_res$FP[, 3])
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 2, color = "#FF00CC") +   
  geom_line(aes(x = FP_5year, y = TP_5year), size = 2, color = "#4169E1") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 2, color = "#EEB422") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 3) +
  theme_bw() +
  annotate("text",   x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#FF00CC"  ) +
  annotate("text",   x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#4169E1"  ) +
  annotate("text",           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#EEB422"  ) +
  labs(x = "GEO-OS-False positive rate", y = "True positive rate(KM)") +
  theme(    axis.text = element_text(face = "bold", size = 11, color = "black"),
            axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
            axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()
#GEO-RFS
riskFile="0_GSE14520_RFS_valid-riskscore.txt"
rt=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1) 
#fix(rt)
colnames(rt)[1:2]<-c("futime", "fustat")
colnames(rt)

library(survivalROC)
rocCol=rainbow(ncol(rt)-1)   
aucText=c()
tiff(file="5_GEO_RFS_survivalROC_plot（1-3-5year）.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc<-survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")  #KM法或NNE法，预测1，3，5year
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="GEO-RFS-False positive rate", ylab="True positive rate (KM) ",
     lwd = 5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score_",1," year","(AUC=",sprintf("%.3f",roc$AUC),")")) 
abline(0,1)

j=1
for(i in c(3,5)){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker =rt$riskScore, predict.time =i, method="KM")   
  j=j+1
  aucText=c(aucText,paste0("risk score_",i," year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 5)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)  
dev.off()

library(timeROC)
library(survival)
time_roc_res <- timeROC(T = rt$futime,delta = rt$fustat,  marker = rt$riskScore,  cause = 1,  weighting="marginal",
                        times = c(1, 3, 5),  ROC = TRUE,  iid = TRUE)  
time_roc_res$AUC   
a<-paste0("GEO-RFS-KM-1 year:", sprintf("%.3f", time_roc_res$AUC[[1]]))  
a
b<-paste0("3 year:", sprintf("%.3f", time_roc_res$AUC[[2]]))  
b
c<-paste0("5 year:", sprintf("%.3f", time_roc_res$AUC[[3]]))  
c
confint(time_roc_res, level = 0.95)$CI_AUC  
#time-dependent ROC
tiff(file="6_GEO_RFS_timeROC_plot（1-3-5year）-1.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
plot(time_roc_res, time=1, col = "#FF00CC", title = FALSE, lwd = 5)  
plot(time_roc_res, time=3, add=TRUE, col="#4169E1", lwd = 5) 
plot(time_roc_res, time=5, add=TRUE, col="#EEB422", lwd = 5)
legend("bottomright",c(a,b,c),  col=c("#FF00CC", "#4169E1", "#EEB422"), lty=5, lwd=4, bty="n") 
dev.off()
tiff(file="6_GEO_RFS_timeROC_ggplot（1-3-5year）-2.TIFF",width =16,height = 15,units ="cm",compression="lzw",bg="white",res=300)
time_ROC_df <- data.frame(  TP_3year = time_roc_res$TP[, 1],  FP_3year = time_roc_res$FP[, 1],
                            TP_5year = time_roc_res$TP[, 2],  FP_5year = time_roc_res$FP[, 2],
                            TP_10year = time_roc_res$TP[, 3],  FP_10year = time_roc_res$FP[, 3])
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 2, color = "#FF00CC") +   
  geom_line(aes(x = FP_5year, y = TP_5year), size = 2, color = "#4169E1") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 2, color = "#EEB422") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 3) +
  theme_bw() +
  annotate("text",   x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#FF00CC"  ) +
  annotate("text",   x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#4169E1"  ) +
  annotate("text",           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#EEB422"  ) +
  labs(x = "GEO-RFS-False positive rate", y = "True positive rate(KM)") +
  theme(    axis.text = element_text(face = "bold", size = 11, color = "black"),
            axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
            axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()















#(12) Nomogram################################################################################################

library(survival)  #coxph()
library(rms)  #cph() lrm()
########1TCGA 
rt_model=read.table("0_TCGA-LIHC_training-riskscore-sur-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
dim(rt_model)
rt_model$gender<-factor(rt_model$gender,labels=c('male','female'))  
rt_model$age<-factor(rt_model$age,labels=c('<40','40~60','>=60')) 
rt_model$stage<-factor(rt_model$stage,labels=c('S1','S2','S3')) 
str(rt_model)
write.table(rt_model, file = "1_TCGA-LIHC__training-riskscore-sur-chr.txt",sep="\t", row.names =T, quote = F)
colnames(rt_model)
dd <- datadist(rt_model)
options(datadist="dd")
units(rt_model$futime)<-"year"
f <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_model) 
surv <- Survival(f)
gxj1<-function(x)surv(1, lp=x) 
gxj3<-function(x)surv(3, lp=x) 
gxj5<-function(x)surv(5, lp=x) 
nom <- nomogram(f, fun=list(gxj1, gxj3, gxj5),  
                lp=F, funlabel=c("1-year OS", "3-year OS", "5-year OS"), 
                maxscale=100,  fun.at=c(0.95, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30,0.20,0.10))  
outFile="2_TCGA-LIHC-OS-line-nomogram-style1-.TIFF"
tiff(file=outFile,width =30,height = 30,units ="cm",compression="lzw",bg="white",res=300)
plot(nom)
dev.off()
outFile="2_TCGA-LIHC-OS-line-nomogram-style2-.TIFF"
tiff(file=outFile,width =30,height = 30,units ="cm",compression="lzw",bg="white",res=300)
plot((nom), col.grid=c("#F4A460","#87CEFA"))
dev.off()
#
library(regplot)
f <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_model) 
surv <- Survival(f)
regplot(f,
        title="Nomograme", 
        clickable=F, 
        observation=rt_model[6,], 
        failtime=c(1,3,5), 
        points=TRUE, 
        interval="confidence",
        rank="sd", 
        showP=T, 
        prfail = T, 
        adds=T,
        droplines=T, 
        #  other=(list(bvcol="#FF00CC",sq="#4169E1",obscol="#EEB422"))
)
library(survival)
f1 <- coxph(Surv(futime, fustat) ~ age+gender+stage+riskScore,  data=rt_model) 
sum_f1<-summary(f1)
c_index<-sum_f1$concordance
options(digits = 3)
c_index
library(rms)
f1 <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore,  data=rt_model,x=T, y=T, surv=T) 
v1<-validate(f1, method="boot", B=1000, dxy=T)
Dxy<-v1[rownames(v1)=="Dxy", colnames(v1)=="index.corrected"]
orig_Dxy<-v1[rownames(v1)=="Dxy", colnames(v1)=="index.orig"]
bias_corrected_c_index<-abs(Dxy)/2+0.5
orig_c_index<-abs(orig_Dxy)/2+0.5
orig_c_index 
bias_corrected_c_index 
dd <- datadist(rt_model)
options(datadist="dd")
units(rt_model$futime)<-"year"
fcox_model <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_model) 
f_TCGA_OS <- cph(Surv(futime, fustat) ~ predict(fcox_model, newdata=rt_model), x=T, y=T, surv=T, data=rt_model) 
validate(f_TCGA_OS, method="boot", B=1000, dxy=T)
c <-rcorrcens(Surv(futime, fustat) ~predict(fcox_model, newdata=rt_model), data = rt_model)
c
c_index<-1-c[1]  
c_index  
options(digits = 3)
SE<- c[4]/2  
c_95_low = c_index - SE*1.96
c_95_low
options(digits = 3)
c_95_high = c_index + SE*1.96
c_95_high
TCGA_OS_C_index<-paste(round(c_index, digits=3),"(", round(c_95_low, digits=3),",", round(c_95_high, digits=3), ")")
TCGA_OS_C_index
GEO_OS_rt=read.table("0_GSE14520_valid-riskscore_OSsurv-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
dim(GEO_OS_rt)
GEO_OS_rt$gender<-factor(GEO_OS_rt$gender,labels=c('male','female'))  
GEO_OS_rt$age<-factor(GEO_OS_rt$age,labels=c('<40','40~60','>=60')) 
GEO_OS_rt$stage<-factor(GEO_OS_rt$stage,labels=c('S1','S2','S3')) 
write.table(GEO_OS_rt, file = "4_GEO-OS_valid-riskscore_suv-chr.txt",sep="\t", row.names =T, quote = F)
colnames(GEO_OS_rt) 
dd <- datadist(GEO_OS_rt)
options(datadist="dd")
units(GEO_OS_rt$futime)<-"year"
f_GEO_OS <- cph(Surv(futime, fustat) ~ predict(fcox_model, newdata=GEO_OS_rt), x=T, y=T, surv=T, data=GEO_OS_rt) 
validate(f_GEO_OS, method="boot", B=1000, dxy=T)
c <-rcorrcens(Surv(futime, fustat) ~predict(fcox_model, newdata=GEO_OS_rt), data = GEO_OS_rt)
c
c_index<-1-c[1]  
c_index  
options(digits = 3)
SE<- c[4]/2  
c_95_low = c_index - SE*1.96
options(digits = 3)
c_95_high = c_index + SE*1.96
GEO_OS_C_index<-paste(round(c_index, digits=3),"(", round(c_95_low, digits=3),",", round(c_95_high, digits=3), ")")
GEO_OS_C_index
GEO_RFS_rt=read.table("0_GSE14520_valid-riskscore_RFSsurv-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
#fix(GEO_RFS_rt)
dim(GEO_RFS_rt)
GEO_RFS_rt$gender<-factor(GEO_RFS_rt$gender,labels=c('male','female'))  
GEO_RFS_rt$age<-factor(GEO_RFS_rt$age,labels=c('<40','40~60','>=60')) 
GEO_RFS_rt$stage<-factor(GEO_RFS_rt$stage,labels=c('S1','S2','S3')) 
write.table(GEO_RFS_rt, file = "4_GEO-RFS_valid-riskscore_suv-chr.txt",sep="\t", row.names =T, quote = F)
colnames(GEO_RFS_rt)
dd <- datadist(GEO_RFS_rt)
options(datadist="dd")
units(GEO_RFS_rt$futime)<-"year"
f_GEO_RFS <- cph(Surv(futime, fustat) ~ predict(fcox_model, newdata=GEO_RFS_rt), x=T, y=T, surv=T, data=GEO_RFS_rt) 
validate(f_GEO_RFS, method="boot", B=1000, dxy=T)
c <-rcorrcens(Surv(futime, fustat) ~predict(fcox_model, newdata=GEO_RFS_rt), data = GEO_RFS_rt)
c_index<-1-c[1]  
c_index  
options(digits = 3)
SE<- c[4]/2  
c_95_low = c_index - SE*1.96
options(digits = 3)
c_95_high = c_index + SE*1.96
GEO_RFS_C_index<-paste(round(c_index, digits=3),"(", round(c_95_low, digits=3),",", round(c_95_high, digits=3), ")")
GEO_RFS_C_index
rt_TCGA_OS=read.table("0_TCGA-LIHC_training-riskscore-sur-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
#fix(rt_TCGA_OS)
dd <- datadist(rt_TCGA_OS)
options(datadist="dd")
units(rt_TCGA_OS$futime)<-"year"  
f_model_1 <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_TCGA_OS, time.inc=1) 
dim(rt_TCGA_OS) 
m=round(dim(rt_TCGA_OS)[1]/4, digits = 0)  
m
cal_1 <- calibrate(f_model_1, cmethod="KM", method="boot", u=1, m=80, B=1000)
f_model_3 <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_TCGA_OS, time.inc=3) 
cal_3 <- calibrate(f_model_3, cmethod="KM", method="boot", u=3, m=80, B=1000)
f_model_5 <- cph(Surv(futime, fustat) ~ age+gender+stage+riskScore, x=T, y=T, surv=T, data=rt_TCGA_OS, time.inc=5) 
cal_5 <- calibrate(f_model_5, cmethod="KM", method="boot", u=5, m=80, B=1000)
tiff(file="5_TCGA-OS-training-Nomogram-calibrate_curve（1year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_1,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#191970", 
     col="#191970" ,
     par.corrected=list(pch=8,col="#191970", lty=4, lwd=4),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("TCGA-OS-1-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
# 3year
tiff(file="5_TCGA-OS-test-Nomogram-calibrate_curve（3year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_3,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#CC0000", 
     col="#CC0000" ,
     par.corrected=list(pch=9,col="#CC0000", lty=4, lwd=4),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("TCGA-OS-3-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
# 5year
tiff(file="5_TCGA-OS-test-Nomogram-calibrate_curve（5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_5, 
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col= "#006400", 
     col= "#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=4),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("TCGA-OS-5-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
# 1,3, 5year-all
tiff(file="5_TCGA-OS-test-Nomogram-calibrate_curve_（1-3-5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#191970", 
     col= "#191970" ,
     par.corrected=list(pch=8,col="#191970", lty=4, lwd=9))  
plot(cal_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#CC0000", 
     col= "#CC0000" ,
     par.corrected=list(pch=9,col="#CC0000", lty=4, lwd=9),
     add=T) 
plot(cal_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#006400", 
     col= "#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("1 year","3 year","5 year"),  col=c("#191970", "#CC0000", "#006400"), lty=5, lwd=4, bty="n") 
title(main = list(paste0("TCGA_OS_C_index: ",TCGA_OS_C_index), cex = 0.8,
                  col = "#CC0000", font = 1))
dev.off()

#（2）GEO-OS
rt_GEO_OS=read.table("0_GSE14520_valid-riskscore_OSsurv-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
#fix(rt_GEO_OS)
library(rms)
colnames(rt_GEO_OS)
dd <- datadist(rt_GEO_OS)
options(datadist="dd")
units(rt_GEO_OS$futime)<-"year"


f_GEO_OS_1<- cph(Surv(futime, fustat) ~ predict(f_model_1, newdata=rt_GEO_OS), x=T, y=T, surv=T, data=rt_GEO_OS, time.inc=1) 
dim(rt_GEO_OS) 
m=round(dim(rt_GEO_OS)[1]/4, digits = 0) 
m
cal_GEO_OS_1 <- calibrate(f_GEO_OS_1 , cmethod="KM", method="boot", u=1, m=56, B=1000)

f_GEO_OS_3<- cph(Surv(futime, fustat) ~ predict(f_model_3, newdata=rt_GEO_OS), x=T, y=T, surv=T, data=rt_GEO_OS, time.inc=3) 
cal_GEO_OS_3 <- calibrate(f_GEO_OS_3 , cmethod="KM", method="boot", u=1, m=56, B=1000)

f_GEO_OS_5<- cph(Surv(futime, fustat) ~ predict(f_model_5, newdata=rt_GEO_OS), x=T, y=T, surv=T, data=rt_GEO_OS, time.inc=5) 
cal_GEO_OS_5 <- calibrate(f_GEO_OS_5 , cmethod="KM", method="boot", u=1, m=56, B=1000)
tiff(file="6_GEO-OS-test-Nomogram-calibrate_curve（1year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_OS_1,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#191970", 
     col="#191970" ,
     par.corrected=list(pch=11,col="#191970", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-OS-1-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
tiff(file="6_GEO-OS-test-Nomogram-calibrate_curve（3year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_OS_3,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#CC0000", 
     col="#CC0000" ,
     par.corrected=list(pch=11,col="#CC0000", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-OS-3-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
tiff(file="6_GEO-OS-test-Nomogram-calibrate_curve（5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_OS_5,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#006400", 
     col="#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-OS-5-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
tiff(file="6_GEO-OS-valid-Nomogram-calibrate_curve_（1-3-5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_OS_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#191970", 
     col= "#191970" ,
     xlab="predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     par.corrected=list(pch=8,col="#191970", lty=4, lwd=9))  
plot(cal_GEO_OS_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#CC0000", 
     col= "#CC0000" ,
     par.corrected=list(pch=9,col="#CC0000", lty=4, lwd=9),
     add=T) 
plot(cal_GEO_OS_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#006400", 
     col= "#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("1 year","3 year","5 year"),  col=c("#191970", "#CC0000", "#006400"), lty=5, lwd=4, bty="n") 
title(main = list(paste0("GEO_OS_C_index: ",GEO_OS_C_index), cex = 0.8,
                  col = "#CC0000", font = 1))
dev.off()
#（3）GEO-RFS
rt_GEO_RFS=read.table("0_GSE14520_valid-riskscore_RFSsurv-sum.txt",header=T,sep="\t",check.names=F,row.names=1)        
#fix(rt_GEO_RFS)
dd <- datadist(rt_GEO_RFS)
options(datadist="dd")
units(rt_GEO_RFS$futime)<-"year"
f_GEO_RFS_1<- cph(Surv(futime, fustat) ~ predict(f_model_1, newdata=rt_GEO_RFS), x=T, y=T, surv=T, data=rt_GEO_RFS, time.inc=1) 
dim(rt_GEO_RFS) 
m=round(dim(rt_GEO_RFS)[1]/4, digits = 0) 
m
cal_GEO_RFS_1 <- calibrate(f_GEO_RFS_1 , cmethod="KM", method="boot", u=1, m=56, B=1000)
f_GEO_RFS_3<- cph(Surv(futime, fustat) ~ predict(f_model_1, newdata=rt_GEO_RFS), x=T, y=T, surv=T, data=rt_GEO_RFS, time.inc=3) 
cal_GEO_RFS_3 <- calibrate(f_GEO_RFS_3 , cmethod="KM", method="boot", u=1, m=56, B=1000)
f_GEO_RFS_5<- cph(Surv(futime, fustat) ~ predict(f_model_1, newdata=rt_GEO_RFS), x=T, y=T, surv=T, data=rt_GEO_RFS, time.inc=5) 
cal_GEO_RFS_5 <- calibrate(f_GEO_RFS_5 , cmethod="KM", method="boot", u=1, m=56, B=1000)

tiff(file="7_GEO-RFS-test-Nomogram-calibrate_curve（1year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_RFS_1,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#191970", 
     col="#191970" ,
     par.corrected=list(pch=11,col="#191970", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-RFS-1-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()


tiff(file="7_GEO-RFS-test-Nomogram-calibrate_curve（3year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_RFS_3,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#CC0000", 
     col="#CC0000" ,
     par.corrected=list(pch=11,col="#CC0000", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-RFS-3-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()


tiff(file="7_GEO-RFS-test-Nomogram-calibrate_curve（5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_RFS_5,
     xlim=c(0,1.2),ylim=c(0,1.2), 
     xlab="predicted probability of nonadherence",
     ylab="3 year_Actual diagnosed nonadherence (proportion)",
     #legend=T, 
     conf.int=T, 
     errbar.col="#006400", 
     col="#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=9),  
     subtitles=T  ) 
abline(0,1,lty=2,lwd=1,col="#666666") 
title(main = list("GEO-RFS-5-year", cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()


tiff(file="7_GEO-RFS-valid-Nomogram-calibrate_curve_（1-3-5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_GEO_RFS_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#191970", 
     col= "#191970" ,
     xlab="predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     par.corrected=list(pch=8,col="#191970", lty=4, lwd=9))  
plot(cal_GEO_RFS_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#CC0000", 
     col= "#CC0000" ,
     par.corrected=list(pch=9,col="#CC0000", lty=4, lwd=9),
     add=T) 
plot(cal_GEO_RFS_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#006400", 
     col= "#006400" ,
     par.corrected=list(pch=11,col="#006400", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("1 year","3 year","5 year"),  col=c("#191970", "#CC0000", "#006400"), lty=5, lwd=4, bty="n") 
title(main = list(paste0("GEO_RFS_C_index: ",GEO_RFS_C_index), cex = 0.8, 
                  col = "#CC0000", font = 1))
dev.off()
tiff(file="8_TCGA-OS-GEO-OS-GEO-RFS-comb_curve（1year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#4B0082", 
     col= "#4B0082" ,
     xlab="predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     par.corrected=list(pch=8,col="#4B0082", lty=4, lwd=9))  
plot(cal_GEO_OS_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#FA8072", 
     col= "#FA8072" ,
     par.corrected=list(pch=9,col="#FA8072", lty=4, lwd=9),
     add=T) 
plot(cal_GEO_RFS_1,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#20B2AA", 
     col= "#20B2AA" ,
     par.corrected=list(pch=11,col="#20B2AA", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("TCGA-OS","GEO-OS","GEO-RFS"),  col=c("#4B0082", "#FA8072", "#20B2AA"), lty=5, lwd=4, bty="n") 
title(main = "1-year", cex = 0.8, 
      col = "#CC0000", font = 1)
dev.off()
tiff(file="8_TCGA-OS-GEO-OS-GEO-RFS-comb_curve（3year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#4B0082", 
     col= "#4B0082" ,
     xlab="predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     par.corrected=list(pch=8,col="#4B0082", lty=4, lwd=9))  
plot(cal_GEO_OS_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#FA8072", 
     col= "#FA8072" ,
     par.corrected=list(pch=9,col="#FA8072", lty=4, lwd=9),
     add=T) 
plot(cal_GEO_RFS_3,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#20B2AA", 
     col= "#20B2AA" ,
     par.corrected=list(pch=11,col="#20B2AA", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("TCGA-OS","GEO-OS","GEO-RFS"),  col=c("#4B0082", "#FA8072", "#20B2AA"), lty=5, lwd=4, bty="n") 
title(main = "3-year", cex = 0.8, 
      col = "#CC0000", font = 1)
dev.off()
tiff(file="8_TCGA-OS-GEO-OS-GEO-RFS-comb_curve（5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)                     
plot(cal_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#4B0082", 
     col= "#4B0082" ,
     xlab="predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     par.corrected=list(pch=8,col="#4B0082", lty=4, lwd=9))  
plot(cal_GEO_OS_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#FA8072", 
     col= "#FA8072" ,
     par.corrected=list(pch=9,col="#FA8072", lty=4, lwd=9),
     add=T) 
plot(cal_GEO_RFS_5,xlim=c(0,1.2),ylim=c(0,1.2),
     errbar.col= "#20B2AA", 
     col= "#20B2AA" ,
     par.corrected=list(pch=11,col="#20B2AA", lty=4, lwd=9),
     add=T) 
abline(0,1,lty=2,lwd=1,col="#666666") 
legend("bottomright",c("TCGA-OS","GEO-OS","GEO-RFS"),  col=c("#4B0082", "#FA8072", "#20B2AA"), lty=5, lwd=4, bty="n") 
title(main = "5-year", cex = 0.8, 
      col = "#CC0000", font = 1)
dev.off()
########IDI-NRI
library('survIDINRI')
data_TCGA_OS<-read.table("0_TCGA-LIHC_training-riskscore-sur-sum.txt",sep="\t",check.names=F,header=T,row.names=1) 
#fix(data_TCGA_OS)
outcome=data_TCGA_OS[,c(1,2)] 
covs1<-as.matrix(data_TCGA_OS[,c(-1,-2)]) 
covs0<-as.matrix(data_TCGA_OS[,c(-1,-2, -6)]) 
#1year
t=1
x_1<-IDI.INF(outcome, covs0, covs1, t, npert=1000)  
y_1<-as.data.frame(IDI.INF.OUT(x_1) )
tiff(file="9_NRI_IDI_survIDINRI_（1year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
IDI.INF.GRAPH(x_1, main=paste0("1 year---IDI:",round(y_1[1,1], digits=3)," (p=", round(y_1[1,4], digits=3),")",
 "----NRI:",round(y_1[2,1], digits=3)," (p=", round(y_1[2,4], digits=3),")")  )
dev.off()
#fix(y_1)
#3year
t=3
x_3<-IDI.INF(outcome, covs0, covs1, t, npert=1000)
y_3<-as.data.frame(IDI.INF.OUT(x_3) )
tiff(file="9_NRI_IDI_survIDINRI_（3year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
IDI.INF.GRAPH(x_3, main=paste0("3 year---IDI:",round(y_3[1,1], digits=3)," (p=", round(y_3[1,4], digits=3),")",
                               "----NRI:",round(y_3[2,1], digits=3)," (p=", round(y_3[2,4], digits=3),")")  )
dev.off()
#fix(y_3)
#5year
t=5
x_5<-IDI.INF(outcome, covs0, covs1, t, npert=1000)
y_5<-as.data.frame(IDI.INF.OUT(x_5)) 
tiff(file="9_NRI_IDI_survIDINRI_（5year）.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
IDI.INF.GRAPH(x_5, main=paste0("5 year---IDI:",round(y_5[1,1], digits=3)," (p=", round(y_5[1,4], digits=3),")",
                               "----NRI:",round(y_5[2,1], digits=3)," (p=", round(y_5[2,4], digits=3),")")  )
dev.off()
y_1<-data.frame(y_1, year="1 year")
rownames(y_1)[1:3]<-c("IDI_1","NRI_1","median-diff_1")
y_3<-data.frame(y_3, year="3 year")
rownames(y_3)[1:3]<-c("IDI_3","NRI_3","median-diff_3")
y_5<-data.frame(y_5, year="5 year")
rownames(y_5)[1:3]<-c("IDI_5","NRI_5","median-diff_5")
NRI_IDI<-rbind(y_1, y_3, y_5 )
fix(NRI_IDI)
write.table(NRI_IDI, file = "9_NRI_IDI_survIDINRI_（1_3_5year）.xls",sep="\t", row.names =T, quote = F)

###"nricens"
rm(list=ls()) 
library('nricens')
data_TCGA_OS<-read.table("0_TCGA-LIHC_training-riskscore-sur-sum.txt",sep="\t",check.names=F,header=T,row.names=1) 
fix(data_TCGA_OS)  
time  = data_TCGA_OS$futime
event = data_TCGA_OS$fustat
z.std = as.matrix(subset(data_TCGA_OS, select = c(age, gender, stage)))
z.new = as.matrix(subset(data_TCGA_OS, select = c(age, gender, stage, riskScore)))
## coxph fit
mstd = coxph(Surv(time,event) ~ ., data.frame(time,event,z.std), x=TRUE) 
mnew = coxph(Surv(time,event) ~ ., data.frame(time,event,z.new), x=TRUE) 
tiff(file="10_TCGA-OS-NRI_nricens_category_(0.3-0.6)(1year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=1)
p.new = get.risk.coxph(mnew, t0=1)
NRI_1<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =1, updown = 'category', cut = c(0.3, 0.6),  niter = 1000, point.method = 'km')  
title(main=paste0("1 year_NRI: ", round (NRI_1$nri[1,1], digits=3),"--95%: ", round (NRI_1$nri[1,2], digits=3),"-",round (NRI_1$nri[1,3], digits=3),
                  "---niter=1000; 0.3-0.6"),
      cex = 0.8,  col="#CC0000", font = 1)   
dev.off()

tiff(file="10_TCGA-OS-NRI_nricens_category_(0.3-0.6)(3year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=3)
p.new = get.risk.coxph(mnew, t0=3)
NRI_3<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =3, updown = 'category', cut = c(0.3, 0.6),  niter = 1000, point.method = 'km')  
title(main=paste0("3 year_NRI: ", round (NRI_3$nri[1,1], digits=3),"--95%: ", round (NRI_3$nri[1,2], digits=3),"-",round (NRI_3$nri[1,3], digits=3),
                  "---niter=1000; 0.3-0.6"),
      cex = 0.8,  col="#CC0000", font = 1)   
dev.off()
tiff(file="10_TCGA-OS-NRI_nricens_category_(0.3-0.6)(5year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=5)
p.new = get.risk.coxph(mnew, t0=5)
NRI_5<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =5, updown = 'category', cut = c(0.3, 0.6),  niter = 1000, point.method = 'km')  
title(main=paste0("5 year_NRI: ", round (NRI_5$nri[1,1], digits=3),"--95%: ", round (NRI_5$nri[1,2], digits=3),"-",round (NRI_5$nri[1,3], digits=3),
                  "---niter=1000; 0.3-0.6"),
      cex = 0.8,  col="#CC0000", font = 1)   
dev.off()
tiff(file="10_TCGA-OS-NRI_nricens_diff_(0.05)(1year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=1)
p.new = get.risk.coxph(mnew, t0=1)
NRI_1a<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =1, updown = 'diff', cut = 0.05,  niter = 1000, point.method = 'km')  
title(main=paste0("1 year_NRI: ", round (NRI_1a$nri[1,1], digits=3),"--95%: ", round (NRI_1a$nri[1,2], digits=3),"-",round (NRI_1a$nri[1,3], digits=3),
                  "---niter=1000; diff:0.05"),
      cex = 0.5,  col="#CC0000", font = 1)   
dev.off()
#3yearNRI
tiff(file="10_TCGA-OS-NRI_nricens_diff_(0.05)(3year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=1)
p.new = get.risk.coxph(mnew, t0=1)
NRI_3a<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =1, updown = 'diff', cut = 0.05,  niter = 1000, point.method = 'km')  
title(main=paste0("3 year_NRI: ", round (NRI_3a$nri[1,1], digits=3),"--95%: ", round (NRI_3a$nri[1,2], digits=3),"-",round (NRI_3a$nri[1,3], digits=3),
                  "---niter=1000; diff:0.05"),
      cex = 0.5,  col="#CC0000", font = 1)   
dev.off()
#5yearNRI
tiff(file="10_TCGA-OS-NRI_nricens_diff_(0.05)(5year).TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
p.std = get.risk.coxph(mstd, t0=1)
p.new = get.risk.coxph(mnew, t0=1)
NRI_5a<-nricens(mdl.std = mstd, mdl.new = mnew, t0 =1, updown = 'diff', cut = 0.05,  niter = 1000, point.method = 'km')  
title(main=paste0("5 year_NRI: ", round (NRI_5a$nri[1,1], digits=3),"--95%: ", round (NRI_5a$nri[1,2], digits=3),"-",round (NRI_5a$nri[1,3], digits=3),
                  "---niter=1000; diff:0.05"),
      cex = 0.5,  col="#CC0000", font = 1)   
dev.off()
##
library(ggDCA)
library(rms)
library(foreign)
rm(list=ls()) 
data_TCGA_OS<-read.table("0_TCGA-LIHC_training-riskscore-sur-sum.txt",sep="\t",check.names=F,header=T,row.names=1) 
#fix(data_TCGA_OS)
dim(data_TCGA_OS)
#data_TCGA_OS <- na.omit(data_TCGA_OS)
dd <- datadist(data_TCGA_OS)
options(datadist="dd")
units(data_TCGA_OS$futime)<-"year"
f_DCA_1 <- cph(Surv(futime, fustat) ~ age+gender,   x=T, y=T, surv=T, data=data_TCGA_OS) 
f_DCA_2 <- cph(Surv(futime, fustat) ~ age+gender+stage, x=T, y=T, surv=T,  data_TCGA_OS) 
f_DCA_3 <- cph(Surv(futime, fustat) ~age+gender+stage+riskScore, x=T, y=T, surv=T, data= data_TCGA_OS) 

tiff(file="11_TCGA-OS-DCA_1year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_tcga_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3,  times=1, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_tcga_all, linetype=2)
dev.off()
tiff(file="11_TCGA-OS-DCA_3year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_tcga_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3,  times=3, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_tcga_all, linetype=2)
dev.off()
tiff(file="11_TCGA-OS-DCA_5year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_tcga_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3,  times=5, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_tcga_all, linetype=2)
dev.off()
tiff(file="11_TCGA-OS-DCA_1-3-5year.TIFF",width =35,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_tcga_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3,  times=c(1,3,5), model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_tcga_all)
dev.off()
data_GEO_OS<-read.table("0_GSE14520_valid-riskscore_OSsurv-sum.txt",sep="\t",check.names=F,header=T,row.names=1) 
#fix(data_TCGA_OS)
#dim(data_GEO_OS)
dd <- datadist(data_GEO_OS)
options(datadist="dd")
units(data_GEO_OS$futime)<-"year"
tiff(file="12_GEO-OS-DCA_1year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_OS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_OS,  times=1, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_OS_all, linetype=2)
dev.off()
tiff(file="12_GEO-OS-DCA_3year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_OS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_OS,  times=3, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_OS_all, linetype=2)
dev.off()
tiff(file="12_GEO-OS-DCA_5year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_OS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_OS,  times=5, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_OS_all, linetype=2)
dev.off()
tiff(file="12_GEO-OS-DCA_1-3-5year.TIFF",width =35,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_OS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_OS,  times=c(1,3,5), model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_OS_all)
dev.off()
data_GEO_RFS<-read.table("0_GSE14520_valid-riskscore_RFSsurv-sum.txt",sep="\t",check.names=F,header=T,row.names=1) 
dd <- datadist(data_GEO_RFS)
options(datadist="dd")
units(data_GEO_RFS$futime)<-"year"
tiff(file="13_GEO-RFS-DCA_1year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_RFS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_RFS,  times=1, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_RFS_all, linetype=2)
dev.off()
tiff(file="13_GEO-RFS-DCA_3year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_RFS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_RFS,  times=3, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_RFS_all, linetype=2)
dev.off()
tiff(file="13_GEO-RFS-DCA_5year.TIFF",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_RFS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_RFS,  times=5, model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_RFS_all, linetype=2)
dev.off()
tiff(file="13_GEO-RFS-DCA_1-3-5year.TIFF",width =35,height = 15,units ="cm",compression="lzw",bg="white",res=300)  
d_GEO_RFS_all <- dca(f_DCA_1,f_DCA_2,f_DCA_3, new.data=data_GEO_RFS,  times=c(1,3,5), model.names =c('age+gender','age+gender+stage','age+gender+stage+riskScore')) 
ggplot(d_GEO_RFS_all)
dev.off()





#(13) risk-clinical-differ####################################################################################
library("limma") 
group_1<-read.table("0_TCGA-LIHC-Normal-HCCgroup-IDfull.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1)  
com = subset(group_1, group =="Tumor"|group =="Normal",  select=group, Decrease=T)  
dim(com)   
fix(com)   
library(dplyr)
gxj<-arrange(com, group) 
com<-as.data.frame(gxj)
fix(com)
clx_1<- factor(com$group)
clx_2<-summary(clx_1) 
clx_2
conNum=clx_2 [1]                                                        
conNum
treatNum=clx_2 [2]                                                
treatNum
gene_1<-read.table("0_TCGA-LIHC_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
fix(eSet)   
write.table(eSet,file="1_TCGA-LIHC_batch-order.txt",sep="\t",quote=F)
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum)) 
inputFile="1_TCGA-LIHC_batch-order.txt"     
rt_1=read.table(inputFile,sep="\t",header=T,check.names=F)
fix(rt_1)
dim(rt_1) 
target_1=read.table( "0_lasso_(geneCoef).txt",sep="\t",header=T,check.names=F,row.names = 1)
fix(target_1)
target<-row.names(target_1)
target  
rt<-rt_1[target,]
fix(rt) 
write.table(rt,file="1_TCGA-LIHC_batch(ord_targ).txt",sep="\t",quote=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)  
newGeneLists=c()  
fix(data)
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)  
  pvalue=wilcoxTest$p.value 
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    if(pvalue<0.001){
    newGeneLists=c(newGeneLists,paste0(i,"***"))
  }else if(pvalue<0.01){
    newGeneLists=c(newGeneLists,paste0(i,"**"))
  }else if(pvalue<0.05){
    newGeneLists=c(newGeneLists,paste0(i,"*"))
  }else{
    newGeneLists=c(newGeneLists,i)
  }
}
fix(outTab)
dim(outTab)
write.table(outTab,file="2_all-gene-dif.xls",sep="\t",row.names=F,quote=F)
heatmap=cbind(ID=newGeneLists,data)
fix(heatmap)
dim(heatmap)
write.table(heatmap,file="3_heatmap_geneSigExp.txt",sep="\t",row.names=F,quote=F)
#install.packages("pheatmap")
library(pheatmap)
gxj_exp_1=read.table("3_heatmap_geneSigExp.txt",sep="\t",header=T,row.names=1,check.names=F)
fix(gxj_exp_1)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))  
names(Type)=colnames(gxj_exp_1)
Type=as.data.frame(Type)
tiff(file="4_heatmap-statistical.tiff",width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(gxj_exp_1, 
         annotation=Type, 
         color = colorRampPalette(c("#00FF7F", "white", "#800080"))(200),
         cluster_cols =F,
         show_colnames = F,
         scale="row",  
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)
dev.off()   
library(vioplot)                                                    
normal=conNum                                                           
tumor=treatNum   
gene_number<-ncol(gxj_diff_2)

gxj_diff_1=read.table("1_TCGA-LIHC_batch(ord_targ).txt",sep="\t",header=T,row.names=1,check.names=F)
gxj_diff_2=t(gxj_diff_1) 
#fix(gxj_diff_2)
 tiff(file="5_vioplot.tiff",width =20,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(las=2,mar=c(4,5,3,3))  
x=c(1:ncol(gxj_diff_2))
y=c(1:ncol(gxj_diff_2))
plot(x,y,
     xlim=c(-1,15),ylim=c(min(gxj_diff_2),max(gxj_diff_2)*1.1),  
         main="",xlab=paste("Tumor:",treatNum, "Normal:",conNum), ylab="TCGA-LIHC-Gene expression",
     pch=21,
     cex.lab=1,
     col="white",
     xaxt="n")
for(i in 1:ncol(gxj_diff_2)){
  normalData=gxj_diff_2[1:normal,i] 
  tumorData=gxj_diff_2[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = '#20B2AA')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = '#FF69B4') 
  wilcoxTest=wilcox.test(normalData,tumorData)  
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData)) 
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))  
  text(x=3*(i-1)+0.5, y=mx*1.05, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)  
}
legend("topleft", legend=c("Normal", "Tumor"), fill=c("#20B2AA", "#FF69B4"), cex = 1) 
#axis(side=1, at = c(0,1,3,4,6,7,9,10), labels = F, tick = TRUE)
text(seq(0,12,3), rep(par("usr")[3]-0.9,gene_number-1), xpd = NA,labels=colnames(gxj_diff_2),cex = 1,srt = 45,pos=1)  
dev.off()
#install.packages("ggcorrplot")
library(corrplot)
library(ggcorrplot)
library(ggthemes)
library(psych)
gxj_diff_1=read.table("1_TCGA-LIHC_batch(ord_targ).txt",sep="\t",header=T,row.names=1,check.names=F)
gxj_diff_2=t(gxj_diff_1) 
env.cor <- round(cor(gxj_diff_2, method = "pearson"),3)  
env.cor
library(psych)
cor <- corr.test(gxj_diff_2,gxj_diff_2,method = "pearson",adjust = "BH",ci = F)
cmt<-cor$r
library(ggcorrplot) 
env.p <-cor_pmat(gxj_diff_2,method = "pearson" )
env.p 
cor <- corr.test(gxj_diff_2,gxj_diff_2,method = "spearman",adjust = "BH",ci = F) #
pmt<-cor$p.adj
pmt
res1 <- cor.mtest(gxj_diff_2, conf.level = 0.95)
env.uppCI = round(cor.mtest(gxj_diff_2, conf.level = 0.95)$uppCI,3)
env.uppCI 
env.lowCI = round(cor.mtest(gxj_diff_2, conf.level = 0.95)$lowCI,3)
env.lowCI 
corrplot(corr=cor(gxj_diff_2),
         method = "pie",  # "circle" , "square", "ellipse", "number", "pie", "shade" ， "color". 
         order = "hclust",
         tl.col="black",  
         addCoef.col = "#2200AA",  
         p.mat = res1$p, # p.mat = env.p, p.mat = res1$p
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1.5,
         type = "lower",
         col=colorRampPalette(c("#FFBB00", "white", "#FF0088"))(50)   )
dev.off()

tiff(file="6_genes-cor-P.tiff",width =20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
corrplot(corr=cor(gxj_diff_2),
         method = "pie",  # "circle" , "square", "ellipse", "number", "pie", "shade" ， "color". 
         order = "hclust",
         tl.col="black",  
         addCoef.col = "#2200AA",  
         p.mat = env.p, 
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1.5,
         type = "lower",
         col=colorRampPalette(c("#FFBB00", "white", "#FF0088"))(50)   )
dev.off()
library('pheatmap')
library('ggplot2')
# KPNA2	FMO3	SPP1	KIF2C	LPCAT1	
title_lable<-"GSE14520"
target_gene<-"FMO3"
group_1<-read.table("0-GSE14520-Pair-232.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
com = subset(group_1, group =="Tumor"|group =="Normal") 
dim(com)   
#fix(com)   
head(com,5)
gene_1<-read.table("0_GSE14520_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
gene_2<-t(gene_1) 
dim(gene_2)  
#fix(gene_2) 
gene_3<-gene_2[rownames(com),]  
#fix(gene_3) 
dim(gene_3) 
same<-match(rownames(gene_3),rownames(com))
gene_group <- data.frame(com,gene_3[same,]) 
dim(gene_group)
#fix(gene_group)
write.table(gene_group, file="1_GSE14520_all_expression_group.txt", quote=FALSE,  sep="\t")  
target_gene_group_1<-gene_group[,c("group","pair",target_gene)]
dim(target_gene_group_1)  
#fix(target_gene_group_1) 
write.table(target_gene_group_1, file=paste("1_",title_lable,"_", target_gene,"_group.txt",sep="" ), quote=FALSE,  sep="\t")  
gxj<-read.table(file=paste("1_",title_lable,"_", target_gene,"_group.txt",sep="" ),sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
#fix(gxj)  
dim(gxj)  
head(gxj)[1:3, 1:3]  
names(gxj)=c("group","pair","target_gene") 
head(gxj)[1:3, 1:3] 
library(ggpubr)  
compare_means(target_gene~ group, data =gxj, method = "wilcox.test")  
#2"ggboxplot"
tiff(file=paste("1_",title_lable,"_",target_gene,"wilcoxon.tiff",sep="" ), width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj, x = "group", y = "target_gene", 
          palette = "UChicago",  
          #palette = "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          color = "group", add = "jitter", 
          notch = T, font.x=10,font.xtickslab = 20,font.ytickslab = 20,
          width = 0.1, size = 1.8, 
          xlab = "Group-Wilcoxon",   
          ylab = paste("expression of ", target_gene, sep="") 
) + stat_compare_means(method = "wilcox.test", 
                       label =  "p.format", #label可以选.signif, ..p.signif.., p.format, ..p.format.., p, ..p..
                       #label.x = 0.2,
                       #label.y = 9,   
                       hide.ns = TRUE)
dev.off()
#3------"ggpaired"
tiff(file=paste("2_",title_lable,"_",target_gene, "Wilcoxon.tiff",sep="" ), width =10,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggpaired(gxj, x = "group", y = "target_gene", color = "purple", 
         fill = "group", line.color = "gray", line.size = 0.1,  
         palette = "AAAS",  #"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
         #ggtheme = theme_pubr(), shape='dose', 
         font.xtickslab = 20,font.ytickslab = 20,
         width = 0.3, size = 1.8, 
         xlab = "Group-Wilcoxon",  
         ylab = paste("expression of ", target_gene, sep="") 
)+ stat_compare_means(paired = TRUE,
                      method = "wilcox.test", 
                      label =  "p.format",
                      hide.ns = TRUE)
dev.off()
#4"ggpubr" 
library(ggpubr)
gxj$group<-factor(gxj$group)
##fix(gxj)
tiff(file=paste("3_",title_lable,"_",target_gene,".tiff",sep="" ), width =30,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj, x = "pair", y = "target_gene", 
          fill = "group",  
          color="group",
          palette = "jco",    
          #sort.val = 'desc', 
          position = position_dodge(nrow(gxj)), 
          sort.by.groups = T, 
          #ggtheme = theme_replace(), 
          # lab.col = "black",lab.pos = "out",
          x.text.angle=90)
dev.off()
#5----"ggdotchart"-OK
tiff(file=paste("4_",title_lable,"_",target_gene,"ggdotchart_style1-okokok.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggdotchart(gxj, x = "pair", y = "target_gene", 
           ylab=paste0( target_gene, " expression"),
           #fill = 'group',  
           color='group',
           position = position_dodge(0.3),
           palette = 'nejm' ,  
           #palette = "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
           add = "segments" ,
           dot.size =3, 
           #rotate = T,     
           group = "group", 
           ggtheme = theme_replace() ,       
           add.params = list(color = "group", size = 0.05)  
) 
dev.off()
tiff(file=paste("4_",title_lable,"_",target_gene,"ggdotchart_style2.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggdotchart(gxj, x = "pair", y = "target_gene",
           #fill = 'group',  
           color='group',
           palette = 'jco' , 
           add = "segments" ,
           dot.size = 5, 
           #rotate = T,     
           group = "group", 
           ggtheme =        
             #add.params = list(color = "group", size = 1)  
)  
dev.off()
tiff(file=paste("4_",title_lable,"_",target_gene,"ggdotchart_style3.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggdotchart(gxj, x = "pair", y ="target_gene",
           #fill = 'group',  
           color='group',
           palette = 'jco' , 
           add = "segments" ,
           dot.size = 5, 
           rotate = T,     
           group = "group", 
           ggtheme = theme_pubr()       
           #add.params = list(color = "group", size = 1)  
)  + theme_cleveland()  
dev.off()
#6----vioplot-"ggviolin"   
library(ggpubr)
tiff(file=paste("5_",title_lable,"_",target_gene,"vioplot_style1.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj,x='group',  y="target_gene", 
         fill = 'group', 
         palette = 'jco',    
         #palette = c('#00AFBB','#E7B800','#FC4E07'), 
         add = 'boxplot', 
         error.plot = "crossbar",
         add.params = list(fill='green')) + 
  theme_cleveland()  
dev.off()
tiff(file=paste("5_",title_lable,"_",target_gene,"vioplot_style2.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj,x='group',  y="target_gene", 
         fill = 'group', 
         palette = 'jco',    
         #palette = c('#00AFBB','#E7B800','#FC4E07'), 
         add = 'dotplot', 
         add.params = list(fill='white')) + 
  theme_cleveland()  
dev.off()
tiff(file=paste("5_",title_lable,"_",target_gene,"vioplot_style3.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj,x='group',  y="target_gene", 
         fill = 'group', 
         palette = 'jco',    
         #palette = c('#00AFBB','#E7B800','#FC4E07'), 
         add = c("jitter", "mean_sd"), 
         add.params = list(fill='white')) + 
  theme_cleveland()  
dev.off()

#7----"ggplot"   
#fix(gxj) 
library(ggplot2)
tiff(file=paste("6_",title_lable,"_",target_gene,"ggplot_style1.tiff",sep="" ), width =20,height = 28,units ="cm",compression="lzw",bg="white",res=300)
ggplot(gxj, aes(x =group, y =target_gene, color=group, fill= group)) +
  #geom_point()
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.05)+
  stat_summary(fun = "mean", geom = "crossbar",
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5, col="purple") +
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.5, col="purple")
dev.off()
tiff(file=paste("6_",title_lable,"_",target_gene,"ggplot_style2.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(gxj, aes(x =group, y =target_gene, color=group, fill= group)) +
  geom_point()+
  stat_summary(fun = "mean", geom = "crossbar",
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5, col="purple") +
  stat_summary(fun.data = "mean_se", geom="errorbar", width =0.5, col="purple")
dev.off()
tiff(file=paste("6_",title_lable,"_",target_gene,"ggplot_style3.tiff",sep="" ), width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggplot(gxj, aes(x =group, y =target_gene, color=group, fill= group)) +
  geom_jitter(width = 0.05, na.rm = TRUE) + 
  stat_summary(fun = "mean", geom = "crossbar",
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.4, col="purple") + 
  stat_summary(fun.data = "mean_se", geom="errorbar", width =0.2, col="purple")  
dev.off()
library(ggpubr)  
p_value<-compare_means(target_gene~ group, data =gxj, method = "wilcox.test")  
#fix(p_value)
p_value[1,1]<-target_gene
#fix(p_value)
p_value
write.csv(p_value, file=paste("7_",title_lable,"_", target_gene,"_statistical.csv"),  quote=FALSE)  
rm(list=ls()) 
exp_1=read.table("0_TCGA_taining-riskscore.txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(exp_1)
dim(exp_1)  
cli_1=read.table("0_TCGA-LIHC_taining-riskscore_surv_chr.txt",sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
com<-exp_1[rownames(cli_1),]  
#fix(com) 
dim(com) 
same<-match(rownames(com),rownames(cli_1))  
exp_cli_1 <- data.frame(cli_1,com[same,]) 
dim(exp_cli_1)
#fix(exp_cli_1) 
write.table(exp_cli_1, file="1_TCGA_hub_gene_riskscore_chr.txt", quote=FALSE,  sep="\t")   
colnames(exp_cli_1) # 
#[1] "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"     "KPNA2"       "FMO3"       
#[11] "KIF2C"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk" 
#stage
library(ggplot2)
library(reshape2)
gxj_1<-exp_cli_1[, c("stage","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_1)   #    
x=colnames(gxj_1)[1]
colnames(gxj_1)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_1)[2:ncol(gxj_1)]){
  gxj_1a=gxj_1[,c(gene,"Type")]
  colnames(gxj_1a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_1a$Type)))>2){
    test=kruskal.test(expression ~ Type, data = gxj_1a)
    p=test$p.value
  }else{
    test=wilcox.test(expression ~ Type, data = gxj_1a)
    p=test$p.value
  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)
}
colnames(gxj_1)=paste0(colnames(gxj_1),geneSig)
data=melt(gxj_1,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")

p1=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = "TCGA_stage--kruskal.test", y = "Gene expression")+
  geom_boxplot()+ facet_wrap(~Gene,nrow =1)+ theme_bw()+  #theme_classic() ； theme_bw()
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tiff(file=paste("2_all_gene-TCGA_Stage_box_panel_kruskal.tiff",sep="" ), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(p1)
dev.off()
library("ggpubr")
gxj_2<-exp_cli_1[, c("stage","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
###-----（1）KPNA2
fix(gxj_2) 
group=levels(factor(gxj_2$stage))
gxj_2$stage=factor(gxj_2$stage, levels=group)
comp=combn(group,2)
two_comparisons=list()
for(i in 1:ncol(comp)){two_comparisons[[i]]<-comp[,i]}
tiff(file=paste("3_KPNA2-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KPNA2", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KPNA2",  legend.title="TCGA",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_KPNA2-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KPNA2", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KPNA2",  legend.title="TCGA",    add = "jitter")  
dev.off()
###-----（2）FMO3
tiff(file=paste("3_FMO3-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="FMO3", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="FMO3",  legend.title="TCGA",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_FMO3-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="FMO3", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="FMO3",  legend.title="TCGA",    add = "jitter")  
dev.off()
###-----（3）SPP1
tiff(file=paste("3_SPP1-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="SPP1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="SPP1",  legend.title="TCGA",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_SPP1-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="SPP1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="SPP1",  legend.title="TCGA",    add = "jitter")  
dev.off()
###-----（4）KIF2C
tiff(file=paste("3_KIF2C-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KIF2C", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KIF2C",  legend.title="TCGA",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_KIF2C-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KIF2C", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KIF2C",  legend.title="TCGA",    add = "jitter")  
dev.off()
###-----（5）LPCAT1
tiff(file=paste("3_PCAT1-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="LPCAT1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="LPCAT1",  legend.title="TCGA",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_PCAT1-TCGA_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="LPCAT1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="LPCAT1",  legend.title="TCGA",    add = "jitter")  
dev.off()
###age
library(ggplot2)
library(reshape2)
gxj_3<-exp_cli_1[, c("age","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_3)       
x=colnames(gxj_3)[1]
colnames(gxj_3)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_3)[2:ncol(gxj_3)]){
  gxj_3a=gxj_3[,c(gene,"Type")]
  colnames(gxj_3a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_3a$Type)))>2){
    test=kruskal.test(expression ~ Type, data = gxj_3a)
    p=test$p.value
  }else{
    test=wilcox.test(expression ~ Type, data = gxj_3a)
    p=test$p.value
  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)
}
colnames(gxj_3)=paste0(colnames(gxj_3),geneSig)
data=melt(gxj_3,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")
p2=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = "TCGA_age", y = "Gene expression")+
  geom_violin()+ geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  facet_wrap(~Gene,nrow =1)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff(file=paste("4_all_gene-TCGA_age_vioplot.tiff",sep="" ), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(p2)
dev.off()
###age
library("ggpubr")
gxj_4<-exp_cli_1[, c("age","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
###-----（1）KPNA2
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="KPNA2", fill = "age",
                  xlab="age",  ylab="KPNA2",  legend.title="TCGA",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_KPNA2-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_KPNA2-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="KPNA2",  fill = "age",
         xlab="age",  ylab="KPNA2",  legend.title="TCGA",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（2）FMO3
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="FMO3", fill = "age",
                  xlab="age",  ylab="FMO3",  legend.title="TCGA",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_FMO3-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_FMO3-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="FMO3",  fill = "age",
         xlab="age",  ylab="FMO3",  legend.title="TCGA",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（3）SPP1
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="SPP1", fill = "age",
                  xlab="age",  ylab="SPP1",  legend.title="TCGA",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

tiff(file=paste("5_SPP1-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_SPP1-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="SPP1",  fill = "age",
         xlab="age",  ylab="SPP1",  legend.title="TCGA",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（4）KIF2C
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="KIF2C", fill = "age",
                  xlab="age",  ylab="KIF2C",  legend.title="TCGA",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_KIF2C-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_KIF2C-TCGA_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="KIF2C",  fill = "age",
         xlab="age",  ylab="KIF2C",  legend.title="TCGA",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（5）LPCAT1
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="LPCAT1", fill = "age",
                  xlab="age",  ylab="LPCAT1",  legend.title="TCGA",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file="5_LPCAT1-TCGA_age_.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file="5_LPCAT1-TCGA_age_.tiff" , width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="LPCAT1",  fill = "age",
         xlab="age",  ylab="LPCAT1",  legend.title="TCGA",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
library(ggplot2)
library(reshape2)
library(plyr)
gxj_5<-exp_cli_1[, c("gender","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_5)   #    

gxj_5$gender<-as.factor(gxj_5$gender)
head(gxj_5)
###all_gene
library(reshape2)
library(ggpubr)
x=colnames(gxj_5)[1]
colnames(gxj_5)[1]="Type"
data_5=melt(gxj_5,id.vars=c("Type"))
colnames(data_5)=c("Type","Gene","Expression")
colnames(data_5)
boxplot
tiff(file="6_all_gene-TCGA_gender_ggboxplot.tiff", width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_5, x="Gene", y="Expression", color = "Type", fill="Type",
            ylab="TCGA_Gene expression",  xlab="",
            legend.title=x, 
            palette = "nejm", 
            #palette = c("blue","red"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            # theme_gray(), 
            width=0.5)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
#"KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" 
gxj_5<-exp_cli_1[, c("gender","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_5)   #    
gxj_5$gender<-as.factor(gxj_5$gender)
###-----（1）KPNA2
target_gene<-"KPNA2"
####ggpubr
library(ggpubr)  
p_value<-compare_means(KPNA2~gender, data =gxj_5, method = "t.test")  
p_value
#p_value[1,1]<-"KPNA2"
#fix(p_value)
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
p1
tiff(file="7_KPNA2-TCGA_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "KPNA2", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("TCGA_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2

mean<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","KPNA2","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_KPNA2-TCGA_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=KPNA2, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=KPNA2, ymax=KPNA2+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("TCGA_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（2）FMO3
target_gene<-"FMO3"
####ggpubr
library(ggpubr)  
p_value<-compare_means(FMO3~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_FMO3-TCGA_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "FMO3", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("TCGA_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","FMO3","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_FMO3-TCGA_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=FMO3, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=FMO3, ymax=FMO3+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("TCGA_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（3）SPP1
target_gene<-"SPP1"
####ggpubr
library(ggpubr)  
p_value<-compare_means(SPP1~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_SPP1-TCGA_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "SPP1", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("TCGA_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"), 
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2--
mean<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","SPP1","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_SPP1-TCGA_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=SPP1, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=SPP1, ymax=SPP1+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("TCGA_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（4）KIF2C
target_gene<-"KIF2C"
####ggpubr
library(ggpubr)  
p_value<-compare_means(KIF2C~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_KIF2C-TCGA_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "KIF2C", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("TCGA_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
####ggplot2
mean<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","KIF2C","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_KIF2C-TCGA_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=KIF2C, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=KIF2C, ymax=KIF2C+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("TCGA_t_test_p.adj=",p1,seq=" "))
dev.off()

###-----（5）LPCAT1
target_gene<-"LPCAT1"
####ggpubr
library(ggpubr)  
p_value<-compare_means(LPCAT1~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_LPCAT1-TCGA_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "LPCAT1", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("TCGA_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","LPCAT1","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_LPCAT1-TCGA_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=LPCAT1, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=LPCAT1, ymax=LPCAT1+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("TCGA_t_test_p.adj=",p1,seq=" "))
dev.off()
rm(list=ls()) 
exp_1=read.table("0_GSE14520_OS_valid-group-riskscore.txt",sep="\t",header=T,check.names=F, row.names=1) 
#fix(exp_1)
dim(exp_1)  
cli_1=read.table("0_GEO-OS_valid-group-riskscore_surv_chr.txt",sep="\t",header=T,check.names=F, row.names=1)
#fix(cli_1)
dim(cli_1) 
com<-exp_1[rownames(cli_1),]  
#fix(com) 
dim(com) 
same<-match(rownames(com),rownames(cli_1))  
exp_cli_1 <- data.frame(cli_1,com[same,]) 
dim(exp_cli_1)
#fix(exp_cli_1) 
write.table(exp_cli_1, file="1_GEO_hub_gene_riskscore_chr.txt", quote=FALSE,  sep="\t")   
colnames(exp_cli_1) # 
#[1] "futime"      "fustat"      "gender"      "age"         "stage"       "riskScore"   "os_time"     "os_stat"     "KPNA2"       "FMO3"       
#[11] "KIF2C"        "KIF2C"       "LPCAT1"      "riskScore.1" "risk" 
#stage  
library(ggplot2)
library(reshape2)
gxj_1<-exp_cli_1[, c("stage","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_1)   #    
x=colnames(gxj_1)[1]
colnames(gxj_1)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_1)[2:ncol(gxj_1)]){
  gxj_1a=gxj_1[,c(gene,"Type")]
  colnames(gxj_1a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_1a$Type)))>2){
    test=kruskal.test(expression ~ Type, data = gxj_1a)
    p=test$p.value
  }else{
    test=wilcox.test(expression ~ Type, data = gxj_1a)
    p=test$p.value
  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)
}
colnames(gxj_1)=paste0(colnames(gxj_1),geneSig)
data=melt(gxj_1,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")
p1=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = "GEO_stage--kruskal.test", y = "Gene expression")+
  geom_boxplot()+ facet_wrap(~Gene,nrow =1)+ theme_bw()+  #theme_classic() ； theme_bw()
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff(file=paste("2_all_gene-GEO_Stage_box_panel_kruskal.tiff",sep="" ), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(p1)
dev.off()
###stage
library("ggpubr")
gxj_2<-exp_cli_1[, c("stage","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
###-----（1）KPNA2
fix(gxj_2) 
group=levels(factor(gxj_2$stage))
gxj_2$stage=factor(gxj_2$stage, levels=group)
comp=combn(group,2)
two_comparisons=list()
for(i in 1:ncol(comp)){two_comparisons[[i]]<-comp[,i]}
tiff(file=paste("3_KPNA2-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KPNA2", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KPNA2",  legend.title="GEO",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_KPNA2-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KPNA2", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KPNA2",  legend.title="GEO",    add = "jitter")  
dev.off()
###-----（2）FMO3
tiff(file=paste("3_FMO3-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="FMO3", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="FMO3",  legend.title="GEO",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_FMO3-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="FMO3", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="FMO3",  legend.title="GEO",    add = "jitter")  
dev.off()
###-----（3）SPP1
tiff(file=paste("3_SPP1-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="SPP1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="SPP1",  legend.title="GEO",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_SPP1-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="SPP1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="SPP1",  legend.title="GEO",    add = "jitter")  
dev.off()
###-----（4）KIF2C

tiff(file=paste("3_KIF2C-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KIF2C", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KIF2C",  legend.title="GEO",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_KIF2C-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="KIF2C", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="KIF2C",  legend.title="GEO",    add = "jitter")  
dev.off()
###-----（5）LPCAT1

tiff(file=paste("3_PCAT1-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="LPCAT1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="LPCAT1",  legend.title="GEO",    add = "jitter")+  
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons=two_comparisons )
dev.off()
tiff(file=paste("3_PCAT1-GEO_Stage_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggboxplot(gxj_2, x="stage", y="LPCAT1", color="stage",palette = "npg", wide=2.5,
          xlab="stage",  ylab="LPCAT1",  legend.title="GEO",    add = "jitter")  
dev.off()
##################age--vioplot
library(ggplot2)
library(reshape2)
gxj_3<-exp_cli_1[, c("age","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_3)   #    
x=colnames(gxj_3)[1]
colnames(gxj_3)[1]="Type"
geneSig=c("")
for(gene in colnames(gxj_3)[2:ncol(gxj_3)]){
  gxj_3a=gxj_3[,c(gene,"Type")]
  colnames(gxj_3a)=c("expression","Type")
  p=1
  if(length(levels(factor(gxj_3a$Type)))>2){
    test=kruskal.test(expression ~ Type, data = gxj_3a)
    p=test$p.value
  }else{
    test=wilcox.test(expression ~ Type, data = gxj_3a)
    p=test$p.value
  }
  Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
  geneSig=c(geneSig,Sig)
}
colnames(gxj_3)=paste0(colnames(gxj_3),geneSig)
data=melt(gxj_3,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")
p2=ggplot(data,aes(x=Type,y=Expression,fill=Type))+
  guides(fill=guide_legend(title=x))+
  labs(x = "GEO_age", y = "Gene expression")+
  geom_violin()+ geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  facet_wrap(~Gene,nrow =1)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tiff(file=paste("4_all_gene-GEO_age_vioplot.tiff",sep="" ), width =18,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(p2)
dev.off()
################age_vioplot
library("ggpubr")
gxj_4<-exp_cli_1[, c("age","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
###-----（1）KPNA2
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="KPNA2", fill = "age",
                  xlab="age",  ylab="KPNA2",  legend.title="GEO",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_KPNA2-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_KPNA2-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="KPNA2",  fill = "age",
         xlab="age",  ylab="KPNA2",  legend.title="GEO",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（2）FMO3
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="FMO3", fill = "age",
                  xlab="age",  ylab="FMO3",  legend.title="GEO",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_FMO3-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_FMO3-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="FMO3",  fill = "age",
         xlab="age",  ylab="FMO3",  legend.title="GEO",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（3）SPP1
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="SPP1", fill = "age",
                  xlab="age",  ylab="SPP1",  legend.title="GEO",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_SPP1-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_SPP1-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="SPP1",  fill = "age",
         xlab="age",  ylab="SPP1",  legend.title="GEO",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（4）KIF2C
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="KIF2C", fill = "age",
                  xlab="age",  ylab="KIF2C",  legend.title="GEO",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
tiff(file=paste("5_KIF2C-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file=paste("5_KIF2C-GEO_age_.tiff",sep="" ), width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="KIF2C",  fill = "age",
         xlab="age",  ylab="KIF2C",  legend.title="GEO",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
###-----（5）LPCAT1
group=levels(factor(gxj_4$age))
gxj_4$age=factor(gxj_4$age, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggviolin=ggviolin(gxj_4, x="age", y="LPCAT1", fill = "age",
                  xlab="age",  ylab="LPCAT1",  legend.title="GEO",    
                  add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "kruskal.test", label.y = 11, label.x =2)  + 
  stat_compare_means(method = "anova", label.y = 10, label.x = 1) +   
  stat_compare_means(comparisons = my_comparisons)  
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

tiff(file="5_LPCAT1-GEO_age_.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
print(ggviolin)
dev.off()
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
tiff(file="5_LPCAT1-GEO_age_.tiff" , width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggviolin(gxj_4, x="age", y="LPCAT1",  fill = "age",
         xlab="age",  ylab="LPCAT1",  legend.title="GEO",    
         add = "boxplot", add.params = list(fill="white")) 
dev.off()
library(ggplot2)
library(reshape2)
library(plyr)
gxj_5<-exp_cli_1[, c("gender","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_5)   #    
gxj_5$gender<-as.factor(gxj_5$gender)
head(gxj_5)
library(reshape2)
library(ggpubr)
x=colnames(gxj_5)[1]
colnames(gxj_5)[1]="Type"
data_5=melt(gxj_5,id.vars=c("Type"))
colnames(data_5)=c("Type","Gene","Expression")
colnames(data_5)
boxplot
tiff(file="6_all_gene-GEO_gender_ggboxplot.tiff", width =20,height = 18,units ="cm",compression="lzw",bg="white",res=300)
p=ggboxplot(data_5, x="Gene", y="Expression", color = "Type", fill="Type",
            ylab="GEO_Gene expression",  xlab="",
            legend.title=x, 
            palette = "nejm", 
            #palette = c("blue","red"),
            # error.plot = "pointrange",
            #add = "none",
            notch = TRUE,
            # theme_gray(), 
            width=0.5)
p=p+rotate_x_text(45)  
p1=p+stat_compare_means(aes(group=Type), method="t.test", label = "p.format")
print(p1)
dev.off()
#"KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" 
gxj_5<-exp_cli_1[, c("gender","KPNA2" ,"FMO3","SPP1", "KIF2C","LPCAT1" )]  
dim(gxj_5)   #    

gxj_5$gender<-as.factor(gxj_5$gender)
###-----（1）KPNA2
target_gene<-"KPNA2"
####ggpubr
library(ggpubr)  
p_value<-compare_means(KPNA2~gender, data =gxj_5, method = "t.test")  
p_value
#p_value[1,1]<-"KPNA2"
#fix(p_value)
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
p1
tiff(file="7_KPNA2-GEO_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "KPNA2", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("GEO_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$KPNA2,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","KPNA2","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_KPNA2-GEO_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=KPNA2, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=KPNA2, ymax=KPNA2+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("GEO_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（2）FMO3
target_gene<-"FMO3"
####ggpubr
library(ggpubr)  
p_value<-compare_means(FMO3~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_FMO3-GEO_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "FMO3", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("GEO_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$FMO3,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","FMO3","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_FMO3-GEO_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=FMO3, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=FMO3, ymax=FMO3+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("GEO_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（3）SPP1
target_gene<-"SPP1"
####ggpubr
library(ggpubr)  
p_value<-compare_means(SPP1~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_SPP1-GEO_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "SPP1", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("GEO_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$SPP1,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","SPP1","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_SPP1-GEO_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=SPP1, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=SPP1, ymax=SPP1+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("GEO_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（4）KIF2C
target_gene<-"KIF2C"
####ggpubr
library(ggpubr)  
p_value<-compare_means(KIF2C~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_KIF2C-GEO_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "KIF2C", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("GEO_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$KIF2C,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","KIF2C","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_KIF2C-GEO_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=KIF2C, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=KIF2C, ymax=KIF2C+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("GEO_t_test_p.adj=",p1,seq=" "))
dev.off()
###-----（5）LPCAT1
target_gene<-"LPCAT1"
####ggpubr
library(ggpubr)  
p_value<-compare_means(LPCAT1~gender, data =gxj_5, method = "t.test")  
p.adj<-p_value[1,5]
p1<-round (p.adj,digits=3)
tiff(file="7_LPCAT1-GEO_gender_ggbarplot.tiff", width =12,height = 18,units ="cm",compression="lzw",bg="white",res=300)
ggbarplot(gxj_5, x = "gender", y = "LPCAT1", 
          ylab =paste(target_gene,"mean+sd"), xlab =paste0("GEO_t_test_p.adj=",p1,seq=" "),
          fill = "gender",  
          color="#666666",size=1.2,  
          error.plot = "upper_errorbar",
          add = c("mean_sd"),
          add.params = list(size=1.2,width=0.5), 
          palette = "jco",    # "lancet"  "npg"  "AAAS"  "JCO"   "UCSCGB"   "D3"   "UChicago" "The Simpsons"   "Futurama"  "Rick and Morty GSEA""
          x.text.angle=45,
          ggtheme = theme_gray(),
          position = position_dodge())
dev.off()
#ggplot2
mean<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=mean) 
sd<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=sd)
N<-aggregate(gxj_5$LPCAT1,by=list(gxj_5$gender),FUN=length)
gxj_data<-data.frame(mean,sd=sd$x,N=N$x) 
#fix(gxj_data)
colnames(gxj_data)=c("gender","LPCAT1","sd","N")
gxj_data$se <- gxj_data$sd / sqrt(gxj_data$N)
#fix(gxj_data)
gxj_data
dat<-as.data.frame(gxj_data)
dat$gender =as.factor(dat$gender)
library(ggplot2)
tiff(file="7_LPCAT1-GEO_gender_ggplot.tiff", width =12,height = 12,units ="cm",compression="lzw",bg="white",res=300)
ggplot(data=dat, aes(x=gender, y=LPCAT1, fill=gender)) +
  geom_bar(stat="identity", color="#666666", size=1.3,position=position_dodge()) +
  geom_errorbar(aes(ymin=LPCAT1, ymax=LPCAT1+se), color="#666666", size=1.5, width=0.7)+
  xlab("Gender") + ylab(paste(target_gene,"mean--se")) +
  ggtitle(paste0("GEO_t_test_p.adj=",p1,seq=" "))
dev.off()
library(ggplot2)
df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
#fix(df)
df2 <- data_summary(df, varname="len",  groupnames=c("supp", "dose"))
# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)
head(df2)
library(ggplot2)
# Default bar plot
p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="#666666", size=1.5, 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=len, ymax=len+sd), color="#666666", size=1.5, width=0.7,
                position=position_dodge(.9)) 
print(p)
# Finished bar plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  theme_classic() +
  scale_fill_manual(values=c('#FF0066','#009999'))

mean<-aggregate(df$len,by=list(df$supp,df$dose),FUN=mean) 
sd<-aggregate(df$len,by=list(df$supp,df$dose),FUN=sd)
N<-aggregate(df$len,by=list(df$supp,df$dose),FUN=length)
data<-data.frame(mean,sd=sd$x,N=N$x) 
colnames(data)=c("supp","dose","len","sd","N")
data$se <- data$sd / sqrt(data$N)
data
ggplot(data, aes(x=dose, y=len, fill=supp))+
  geom_bar(position=position_dodge(), color="black",stat="identity",width=.6)+
  geom_errorbar(aes(ymin=len-sd,ymax=len+sd),width=.2, position=position_dodge(.6))
ggplot(data, aes(x=dose, y=len, fill=supp))+
  geom_bar(position=position_dodge(), color="black",stat="identity",width=.6)+
  geom_errorbar(aes(ymin=len, ymax=len+sd),width=.2,position=position_dodge(.6))












#(14) GO-KEGG-GSEA##########################################################################
rm(list=ls()) 
library(limma) 
library(sva) 
rt = read.table("0_TCGA_LIHC_All_HTSeq_FPKM-UQ-ok.txt",header=T,sep="\t",check.names=F, row.names=1)
fix(rt) 
dim(rt) 
geo_1 = read.table("0_GSE14520_ID-ALL-Sample.txt",header=T,sep="\t",check.names=F)
fix(geo_1) 
dim(geo_1)  
geo_2<-t(geo_1)
fix(geo_2) 
dim(geo_2)  
sameGene=intersect(row.names(rt),row.names(geo_2)) 
length(sameGene) 
rtOut=rt[sameGene,] 
dim(rtOut)
geoOut=geo_2[sameGene,] 
dim(geoOut)
fix(geoOut)
all=cbind(rtOut,geoOut)  
fix(all) 
dim(all)
clx<- apply(all, 1, function(row) all(row !=0 ))  
all_1<-all[clx,]
fix(all_1)
dim(all_1)  
all_2<-na.omit(all_1)  
fix(all_2)
dim(all_2)  
all_2<-as.matrix.data.frame(all_2) 
str(all_2)
library(sva)
batchType=c(rep(1,ncol(rtOut)),rep(2,ncol(geoOut))) 
batchType
length(all_2) #912
outTab=ComBat(all_2, batchType,par.prior=TRUE) 
fix(outTab)
dim(outTab)  
write.table(outTab,file="1_TCGA-LIHC-GSE14520_batch.txt",sep="\t",quote=F)
rtOut_1=outTab[,colnames(rtOut)] 
rtOut_2=rbind(ID=colnames(rtOut),rtOut_1)  
#fix(rtOut_2)
dim(rtOut_2)
rtOut_3<-rtOut_2[-1,] 
fix(rtOut_3)
dim(rtOut_3)  
write.table(rtOut_3,file="2_TCGA-LIHC_batch.txt",sep="\t",quote=F)
geoOut_1=outTab[,colnames(geoOut)]
geoOut_2=rbind(ID=colnames(geoOut),geoOut_1)  
#fix(geoOut_2)
dim(geoOut_2)
geoOut_3<-geoOut_2[-1,] 
fix(geoOut_3)
dim(geoOut_3) 
write.table(geoOut_3,file="3_GSE14520_batch.txt",sep="\t",quote=F)
rm(list=ls()) 
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0_TCGA_taining-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1)  
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
fix(com)   
colnames(com)[1]<-c("group")
head(com,5)
gene_1<-read.table("0_TCGA-LIHC_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
fix(eSet)  
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
fix(group_list) 
rownames(group_list) = rownames(com)      
colnames(group_list)[1] = "description"  
dim(group_list) 
fix(group_list) 
design = model.matrix(~ description + 0, group_list)  
fix(design)
colnames(design) = levels(as.factor(c("low","high")))   
factor(design)  
head(design) 
cont.matrix = makeContrasts(high-low, levels=design)   
fdr = 0.05
logFC = 1
fit = lmFit(eSet, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="BH", sort.by="logFC", resort.by = "P" ,n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "Pos"  # high-low,logFC>1,
tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Neg"  # high-low,logFC< 1, 
tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "NO"
colnames(tT)[ncol(tT)] = "expression"
dim(tT)
fix(tT)  
write.table(tT,file="1_TCGA_high-low_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
fix(deg) 
write.table(deg,file="2_TCGA_high-low_diff_gene.xls",sep="\t",quote=F)
group_reg <- factor(deg$expression)
labels=levels(group_reg)
gxj<-summary(group_reg) 
gxj
#1：colcano
temp1 = tT[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("Pos","Neg","NO"), order=T)
#
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+  
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080" ))+   
  labs(x=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]),y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
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
tiff(file="3_colcano_style1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_colcano_style2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+   
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080"))
dev.off()
temp2 = data.frame(log10(rowMeans(eSet)))
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT),rownames(temp2)),] 
temp2 = data.frame(temp2,tT[,c("logFC","expression")])
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
tiff(file="4_MA.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300) 
print(P_MA)
dev.off()
#3：heatmap(all_gene)
dim(eSet)  #
#fix(eSet)
DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),] 
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_all_gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
target_gene<-row.names(deg)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target_gene_heatmap_(TOP10).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
gene_hotmap<-read.table("000-target_gene-remove_hub_gene.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_hotmap)  
fix(gene_hotmap)  
target_gene_2<-rownames(gene_hotmap)
deg_eset_2<-eSet[target_gene_2,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target_gene_heatmap_(remove_hub_gene).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
library('limma')  
library('ggplot2')
library('reshape2')
library('pheatmap')
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggrepel")
library(dplyr)
group_1<-read.table("0_GSE14520_OS_valid-group-riskscore.txt",  sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(group_1)  
fix(group_1)  
com = subset(group_1, risk =="low"|risk =="high",  select=risk)  
dim(com)   
fix(com)   
colnames(com)[1]<-c("group")
head(com,5)
gene_1<-read.table("0_GSE14520_batch.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_1)  
fix(gene_1)   
gsm = match(rownames(com),colnames(gene_1))    
eSet=gene_1[,gsm]   
dim(eSet)   
fix(eSet)  
group_list = data.frame(as.character(com$group)) 
dim(group_list) 
fix(group_list) 
rownames(group_list) = rownames(com)      
colnames(group_list)[1] = "description"  
dim(group_list) 
fix(group_list) 
design = model.matrix(~ description + 0, group_list)  
fix(design)
colnames(design) = levels(as.factor(c("low","high")))   
factor(design)  
head(design) 
cont.matrix = makeContrasts(high-low, levels=design)   
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
fix(tT)  
write.table(tT,file="1_GEO_high-low_all_gene.txt",sep="\t",quote=F)
deg = tT[fdr > tT[,"FDR"]  &  abs(tT[,"logFC"]) >= 1,]
dim(deg)
fix(deg) 
write.table(deg,file="2_GEO_high-low_diff_gene.xls",sep="\t",quote=F)
group_reg <- factor(deg$expression)
labels=levels(group_reg)
gxj<-summary(group_reg) 
gxj
#1：colcano
temp1 = tT[,c("FDR","logFC","expression")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","expression")
temp1$expression=factor(temp1$expression, levels=c("Pos","Neg","NO"), order=T)
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(alpha=0.9, size=5, aes(color=temp1$expression))+  
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080" ))+   
  labs(x=paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]),y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
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
tiff(file="3_colcano_style1.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
print(P_volcano)
dev.off()
tiff(file="3_colcano_style2.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
ggplot(temp1, aes(x=temp1$logFC,y=temp1[,"-log10FDR"],colour=temp1$expression)) +
  xlab(paste("log2FC", "Pos:", gxj[1],"Neg:", gxj[2]))+ 
  ylab('-log10p-value') +
  geom_point(alpha=0.8, size=6)+   
  scale_color_manual(values =c("Pos" = "#8B0000", "NO" = "gray","Neg" = "#008080"))
dev.off()
temp2 = data.frame(log10(rowMeans(eSet)))
rownames(temp2) = rownames(eSet)
temp2 = temp2[match(rownames(tT),rownames(temp2)),] 
temp2 = data.frame(temp2,tT[,c("logFC","expression")])
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
tiff(file="4_MA.tiff", width =25,height = 20,units ="cm",compression="lzw",bg="white",res=300) 
print(P_MA)
dev.off()
#3：heatmap(all_gene)
dim(eSet)  #
#fix(eSet)
DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
head(DEG_list)
deg_eset = eSet[match(DEG_list,rownames(eSet)),] 
dim(deg_eset)
#fix(deg_eset)
tiff(file="5_all_gene_heatmap.tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
#heatmap cor
target_gene<-row.names(deg)[c(1:10)]
target_gene
deg_eset_2<-eSet[target_gene,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target_gene_heatmap_(TOP10).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
gene_hotmap<-read.table("000-target_gene-remove_hub_gene.txt",sep='\t',quote="", comment.char="", check.names=FALSE,header=T,row.names=1)  
dim(gene_hotmap)  
fix(gene_hotmap)  
target_gene_2<-rownames(gene_hotmap)
deg_eset_2<-eSet[target_gene_2,]
dim(deg_eset_2)
#fix(deg_eset_2)
tiff(file="6_target_gene_heatmap_(remove_hub_gene).tiff", width =35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
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
#install.packages("VennDiagram")
rm(list=ls()) 
library(VennDiagram)  
outFile="2-R—intersectGenes.txt"       
files=dir()                      
files=grep("txt$",files,value=T)  
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)        
  geneNames=as.vector(rt[,1])              
  geneNames=gsub("^ | $","",geneNames)     
  uniqGene=unique(geneNames)               
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}
tiff(file="1-R-TCGA-GEO-venn.tiff", width =12,height =12,units ="cm",compression="lzw",bg="white",res=300)
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
grid.draw(venn.plot)
dev.off()
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
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
library("enrichplot")
#GO
gxj <- read.table("TCGA-GEO-130.txt",header=T,comment.char = "",check.names=F,sep="\t")
dim(gxj)
fix(gxj)
colnames(gxj)[c(1)]=c("GeneSymbol")  
fix(gxj)
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
#GO
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
#molecular function (MF)
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
plotGOgraph(ego_MF,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T)  #useInfo就看是否要加GO:
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
tiff(file="6_GO_MF_heatplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
heatplot(ego_MF)
dev.off()
tiff(file="7_GO_MF_emapplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=300)
emapplot(ego_MF)
dev.off()
#biological process (BP)
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
plotGOgraph(ego_BP,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T)  #useInfo就看是否要加GO:
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
#cellular component (CC) 
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
plotGOgraph(ego_CC,firstSigNodes =5, sigForAll = TRUE,useFullNames = T,useInfo=T)  #useInfo就看是否要加GO:
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
#KEGG pathway
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
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
files=grep(".xls",dir(),value=T)                                       
data = lapply(files,read.delim)                                    
names(data) = files
dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                            
gseaCol=c("red","blue","green","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
  theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
  theme(legend.background = element_blank()) + theme(legend.key = element_blank())
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "High expression<----------->Low expression", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "red"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)
gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()
tiff(file="1_TCGA_multipleGSEA.tiff",width =30,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
files=grep(".xls",dir(),value=T)                                        
data = lapply(files,read.delim)                                       
names(data) = files
dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                            
gseaCol=c("red","blue","green","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
  theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
  theme(legend.background = element_blank()) + theme(legend.key = element_blank())
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "High expression<----------->Low expression", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "red"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)
gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()
tiff(file="1_GEO_multipleGSEA.tiff",width =30,height = 15,units ="cm",compression="lzw",bg="white",res=300)
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()

