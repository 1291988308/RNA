#step1
# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)
library(GEOquery)

exp <- read.delim("148219featureCount/all.id.txt", comment.char="#", header=T, row.names=1 )
exp <- exp[,-c(1:5)]    ##多了前面五个，不要
name_tmp=c()
for (i in 1:ncol(exp))
{
  name_tmp[i]=unlist(strsplit(colnames(exp)[i],"[.]"))[9]
}

colnames(exp) <- name_tmp  

# 去除前的基因表达矩阵情况
dim(exp)
exp=exp[,c(1:18)]    #取Adult标本

# 获取分组信息
gse=getGEO("GSE148219")
a=gse[[1]]
pd=pData(a)
library(stringr)

###########################################################################
#group_list----          
group_list=ifelse(str_detect(pd$`disease:ch1`,"CTR"),"control",
                      ifelse(str_detect(pd$`disease:ch1`,"BicC"),"Bicuspid","Tricuspid"))
group_list = factor(group_list,
                        levels = c("control","Tricuspid","Bicuspid"),ordered = T)                  

# 过滤在至少在75%的样本中都有表达的基因
tmp=as.matrix(exp)
keep_exp <- rowSums(tmp>0) >= floor(0.75*ncol(tmp))
table(keep_exp)
exp <- exp[keep_exp,]
exp[1:4,1:4]

# 加载edgeR包计算counts per millio(cpm) 表达矩阵,并对结果取log2值
rawdata148219=exp  #  保存非标化供后续分析
library(edgeR)
exp <- log2(cpm(exp)+1)
exp[1:6,1:6]
boxplot(exp)
exp148219 <- normalizeBetweenArrays(exp)
exp[1:6,1:6]
boxplot(exp148219)
# 保存表达矩阵和分组结果
save(group_list,exp148219,rawdata148219,pd,file = "Step01-exp148219.Rdata")






