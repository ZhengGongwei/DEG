# M Acute liver failure, Model group
# CK Acute liver failure, Control group
# GSE217659_RAW

# https://mp.weixin.qq.com/s/lJatjIvzqySEBgsuhh9s_w

setwd("C:/Users/Lenovo/OneDrive/liver/GSE217659_RAW")
#Dclk1

# dat=lapply(files, function(i){
#   read.table(i,sep="\t")
#   })

# extracted_dat <- lapply(dat, function(df) {
#   df <- df[, (ncol(df)-1):ncol(df)]  # 提取最后两列
#   colnames(df) <- c("gene", "count")  # 重命名列名
#   return(df)
# })

# # 使用merge函数将提取后的数据框合并为一个数据框
# merged_dat <- Reduce(function(x, y) merge(x, y, by = c("gene", "count"), all = TRUE), extracted_dat)

files <- list.files(pattern = "txt.gz$")
count_matrix <- NULL
for(i in files){
  temp <- read.table(i, sep = "\t")
  samplename <- strsplit(gsub(".count.txt.gz", "", i), split = "_")[[1]][2]
  colnames(temp) <- c("ENSG", "GeneSymbol", samplename)
  temp <- temp[, c("GeneSymbol", samplename)]
  temp <- temp[!duplicated(temp$GeneSymbol),]

  if(is.null(count_matrix)){
    count_matrix <- temp
  } else {
    count_matrix <- merge(count_matrix, temp, by = "GeneSymbol",all = FALSE)
    count_matrix <- count_matrix[!duplicated(count_matrix$GeneSymbol),]
  }
}
rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]

library("DESeq2")

group <- factor(gsub("(M|CK).*", "\\1", colnames(count_matrix)))
Data <- data.frame(row.names = colnames(count_matrix),group = group)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = Data,
                              design = ~ group)
# 标准化数据
dds2 <- DESeq(dds)
# 从DESeq矩阵中提取差异表达表
tmp <- results(dds2,contrast = c("group","M","CK"))
# > tmp
# log2 fold change (MLE): group M vs CK 
# Wald test p-value: group M vs CK 
# 将差异表达表转化数据框
DEG_DESeq2 <- as.data.frame(tmp)
# 设定阈值
FC <- 2
padj <- 0.05
# 为差异表达表添加"Significant"一列，并暂时全部赋予"normal"，用于后续标差异基因
DEG_DESeq2$Significant <- "normal"
# 筛选上调和下调的的基因，up/down返回值为上调/下调基因的行号
up <- intersect(which(DEG_DESeq2$log2FoldChange > log2(FC) ),
                which(DEG_DESeq2$padj < padj))
                  
down <- intersect(which(DEG_DESeq2$log2FoldChange < (-log2(FC))),
                  which(DEG_DESeq2$padj < padj))
# 根据上下调基因情况将"Significant"列中的"normal"替换为"up"或"down"
DEG_DESeq2$Significant[up] <- "up"
DEG_DESeq2$Significant[down] <- "down"
# 统计频数
table(DEG_DESeq2$Significant)
# 保存R格式
save(DEG_DESeq2, file = "DEG_DESeq2_Modelvscontrol.Rdata")
# 保存xls格式
write.table(DEG_DESeq2,"DEG_DESeq2_Modelvscontrol2.xls",
          sep="\t",quote = F)



# 第二个数据的 Rat数据 分为Bashline和APAP两组 83825_at
# GPL=getGEO(filename = 'GSE205203_family.soft.gz')
# # 提取信息（可以通用）
# gpl=GPL@gpls[[1]]@dataTable@table
# # 83825_at
# # 83825
# # doublecortin-like kinase 1
library(limma)
library(edgeR)
setwd("C:/Users/Lenovo/OneDrive/liver/GSE205203_RAW")
library(GEOquery)
gset <- getGEO(filename='GSE205203-GPL22499_series_matrix.txt.gz',getGPL = F) # Rat Genome
Rat_matrix <- gset@assayData[["exprs"]]
colnames(Rat_matrix) <- gset@phenoData@data[["title"]]

# QC
# par(cex = 0.7)
# n.sample=ncol(Rat_matrix)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# boxplot(Rat_matrix, col = cols,main="expression value",las=2)

samplelist <- c("Rat_Baseline_0h_1","Rat_Baseline_0h_2","Rat_Baseline_0h_3","Rat_Baseline_0h_4","Rat_Baseline_0h_5",
	"Rat_APAP_03h_1","Rat_APAP_03h_2","Rat_APAP_03h_3","Rat_APAP_03h_4","Rat_APAP_03h_5","Rat_APAP_06h_1","Rat_APAP_06h_2",
	"Rat_APAP_06h_3","Rat_APAP_06h_4","Rat_APAP_06h_5","Rat_APAP_09h_1","Rat_APAP_09h_2","Rat_APAP_09h_3","Rat_APAP_09h_4",
	"Rat_APAP_09h_5","Rat_APAP_24h_1","Rat_APAP_24h_2","Rat_APAP_24h_3","Rat_APAP_24h_4","Rat_APAP_24h_5")
Rat_matrix <- Rat_matrix[,samplelist]
group_list <- factor(gsub(".*_(Baseline|APAP)_.*", "\\1", samplelist))

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(Rat_matrix)
design

contrast.matrix<-makeContrasts(APAP-Baseline,levels = design)
contrast.matrix ##这个矩阵声明，我们要把APAP组Baseline跟进行差异分析比较

##step1
fit <- lmFit(Rat_matrix,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
# 保存R格式
save(nrDEG, file = "DEG_limma_APAPvsBaseline_Rat.Rdata")
# 保存xls格式
write.table(nrDEG,"DEG_limma_APAPvsBaseline_Rat.xls",
          sep="\t",quote = F)



# 第二个数据的 Mouse数据 分为Bashline和APAP两组 13175_at

# # GPL29970
# gpl2=GPL@gpls[[2]]@dataTable@table
# 13175_at
# Dclk1
# 13175
# doublecortin-like kinase 1
# 我只要ID和symbol
ids=gpl[,c(1,12)]
# 写出文件
write.table(ids,file = "ids.txt",siep = "\t",row.names=F,col.names = T)

library(limma)
library(edgeR)
setwd("C:/Users/Lenovo/OneDrive/liver/GSE205203_RAW")
library(GEOquery)
gset <- getGEO(filename='GSE205203-GPL29970_series_matrix.txt.gz',getGPL = F) # Mouse Genome
Mouse_matrix <- gset@assayData[["exprs"]]
colnames(Mouse_matrix) <- gset@phenoData@data[["title"]]
# QC
# par(cex = 0.7)
# n.sample=ncol(Mouse_matrix)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# boxplot(Mouse_matrix, col = cols,main="expression value",las=2)

samplelist <- c("Mouse_Baseline_0h_1","Mouse_Baseline_0h_2","Mouse_Baseline_0h_3","Mouse_Baseline_0h_4","Mouse_Baseline_0h_5",
	"Mouse_APAP_03h_1","Mouse_APAP_03h_2","Mouse_APAP_03h_3","Mouse_APAP_03h_4","Mouse_APAP_03h_5","Mouse_APAP_06h_1","Mouse_APAP_06h_2",
	"Mouse_APAP_06h_3","Mouse_APAP_06h_4","Mouse_APAP_06h_5","Mouse_APAP_09h_1","Mouse_APAP_09h_2","Mouse_APAP_09h_3",
	"Mouse_APAP_09h_4","Mouse_APAP_09h_5","Mouse_APAP_24h_1","Mouse_APAP_24h_2","Mouse_APAP_24h_3","Mouse_APAP_24h_4","Mouse_APAP_24h_5")

Mouse_matrix <- Mouse_matrix[,samplelist]
group_list <- factor(gsub(".*_(Baseline|APAP)_.*", "\\1", samplelist))

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(Mouse_matrix)
design

contrast.matrix<-makeContrasts(APAP-Baseline,levels = design)
contrast.matrix ##这个矩阵声明，我们要把APAP组Baseline跟进行差异分析比较

##step1
fit <- lmFit(Mouse_matrix,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
# 保存R格式
save(nrDEG, file = "DEG_limma_APAPvsBaseline_Mouse.Rdata")
# 保存xls格式
write.table(nrDEG,"DEG_limma_APAPvsBaseline_Mouse.xls",sep="\t",quote = F)







# RNA-seq数据
library(EnhancedVolcano)
setwd("C:/Users/Lenovo/OneDrive/liver/GSE217659_RAW")
load("DEG_DESeq2_Modelvscontrol.Rdata")

p1 <- EnhancedVolcano(DEG_DESeq2,
	    lab = rownames(DEG_DESeq2),
	    x = 'log2FoldChange',
	    y = 'pvalue',
	    selectLab = c('Dclk1'),
	    title = 'Model group versus Control group',
	    cutoffLineType = 'twodash',
	    xlab = bquote(~Log[2]~ 'fold change'),
	    pCutoff = 10e-12,
	    FCcutoff = 1.5,
	    pointSize = 2.0,
	    labSize = 6.0,
	    colAlpha = 1,
	    legendLabSize = 12,
	    legendIconSize = 4.0,
	    drawConnectors = TRUE,
	    widthConnectors = 0.75)

ggsave(p1,device = "pdf",file="MvsCk.pdf",width = 7,height = 7)

setwd("C:/Users/Lenovo/OneDrive/liver/GSE205203_RAW")
load("DEG_limma_APAPvsBaseline_Rat.Rdata")

p2 <- EnhancedVolcano(nrDEG,
	    lab = rownames(nrDEG),
	    x = 'logFC',
	    y = 'P.Value',
	    selectLab = c('83825_at'),
	    title = 'Rat APAP group versus Baseline group',
	    cutoffLineType = 'twodash',
	    xlab = bquote(~Log[2]~ 'fold change'),
	    pCutoff = 0.008,
	    FCcutoff = 1,
	    pointSize = 2.0,
	    labSize = 6.0,
	    colAlpha = 1,
	    legendLabSize = 12,
	    legendIconSize = 4.0,
	    drawConnectors = TRUE,
	    widthConnectors = 0.75)

ggsave(p2,device = "pdf",file="Rat_APAP_vs_Baseline.pdf",width = 7,height = 7)

load("DEG_limma_APAPvsBaseline_Mouse.Rdata")
p3 <- EnhancedVolcano(nrDEG,
	    lab = rownames(nrDEG),
	    x = 'logFC',
	    y = 'P.Value',
	    selectLab = c('13175_at'),
	    title = 'Mouse APAP group versus Baseline group',
	    cutoffLineType = 'twodash',
	    xlab = bquote(~Log[2]~ 'fold change'),
	    pCutoff = 0.008,
	    FCcutoff = 1,
	    pointSize = 2.0,
	    labSize = 6.0,
	    colAlpha = 1,
	    legendLabSize = 12,
	    legendIconSize = 4.0,
	    drawConnectors = TRUE,
	    widthConnectors = 0.75)

ggsave(p3,device = "pdf",file="Mouse_APAP_vs_Baseline.pdf",width = 7,height = 7)