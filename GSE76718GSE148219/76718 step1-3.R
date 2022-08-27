#step1
# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)
exp <- read.delim("76718featureCount/counts_GSE76718.txt", comment.char="#", header=T, row.names=1 )

name_tmp=c()
for (i in 1:ncol(exp))
{
  name_tmp[i]=unlist(strsplit(colnames(exp)[i],"[.]"))[9]
  }
     
colnames(exp) <- name_tmp  

# 去除前的基因表达矩阵情况
dim(exp)
boxplot(exp)
# 获取分组信息
library("GEOquery")
gse=getGEO("GSE76718")
a=gse[[1]]
pd=pData(a)
library(stringr)
#group_list----
group_list=ifelse(str_detect(pd$title,"Bicuspid"),"Bicuspid",
                           ifelse(str_detect(pd$title,"non-calcified"),"control","Tricuspid"))
group_list = factor(group_list,
                             levels = c("control","Tricuspid","Bicuspid"),ordered = T)                  

# 过滤在至少在75%的样本中都有表达的基因
tmp=as.matrix(exp)
keep_exp <- rowSums(tmp>0) >= floor(0.75*ncol(tmp))
table(keep_exp)
exp <- exp[keep_exp,]
exp[1:4,1:4]

# 加载edgeR包计算counts per millio(cpm) 表达矩阵,并对结果取log2值
rawdata76718=exp  #  保存非标化供后续分析
library(edgeR)
exp <- log2(cpm(exp)+1)
exp[1:6,1:6]
boxplot(exp)
exp76718 <- normalizeBetweenArrays(exp)
exp[1:6,1:6]
boxplot(exp76718)

# 保存表达矩阵和分组结果
save(group_list,exp76718,rawdata76718,pd,file = "Step01-exp76718.Rdata")

