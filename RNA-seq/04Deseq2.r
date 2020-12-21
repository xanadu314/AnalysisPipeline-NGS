## import expr
library(reshape2)#载入reshape2这个包，用于后面长宽格式的转化
library(readr)
tmp <- read_table2("tmp.txt", col_names = FALSE) #这里要注意选择分隔符，header的格式，多尝试几次，确定文件导入正确，可以手动导入。
expr <- dcast(tmp,formula = X2~X1)#这里的X2和X1的关系要根据文件格式来确定
head(expr)

## DEG generate
library(DESeq2)
###取出样本count值
count_data<-expr[,2:5]
###把第一列设为行名
row.names(count_data)<-expr[,1]
###构建DESeq2表达矩阵
condition<-factor(c("trt","trt","untrt","untrt"),levels=c("trt","untrt"))#处理组放前面
col_data<-data.frame(row.names=colnames(count_data),condition)
dds<-DESeqDataSetFromMatrix(countData=count_data,colData=col_data,design=~condition)
dds_filter<-dds[rowSums(counts(dds))>1,]#样本基因表达量之和大于1
dds_out<-DESeq(dds_filter)
res<-results(dds_out)

## DES analysis
summary(res)
### ID转换
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol<-mapIds(org.Hs.eg.db,data,keytype="ENSEMBL",column="SYMBOL")
res<-res[orser(res$padj),]
diff_gene<-subset(res,padj <0.05 & (log2FoldChange >1 | log2FoldChange <-1))
write.csv(diff_gene,file="DEG_trt_vs_untrt.csv")

