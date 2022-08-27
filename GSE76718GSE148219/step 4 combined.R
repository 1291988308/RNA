rm(list = ls())
options(stringsAsFactors = F)
if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray")
if(!require(AnnoProbe))devtools::install_github("jmzeng1314/AnnoProbe")
library(GEOquery)
library(stringr)
library(AnnoProbe)
library(ggplot2)
library(tinyarray)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)



######################################################################
#####################################################################
###################################################################
###########数据集准备完毕
rm(list = ls())
library("sva")
library("dplyr")
library("limma")
load(file = "Step01-exp76718.RData")
load(file = "Step01-exp148219.RData")


boxplot(exp76718)
exp76718 <- as.data.frame(exp76718)
exp76718$ENTREZID <- rownames(exp76718)      ##增加ENTREZID一列

boxplot(exp148219)
exp148219 <- as.data.frame(exp148219)
exp148219$ENTREZID <- rownames(exp148219)      ##增加ENTREZID一列



merge_eset=inner_join(exp148219,exp76718,by="ENTREZID")
dim(merge_eset)

table(duplicated(merge_eset$ENTREZID))    ##设置各个数据集的行名，去除重复的基因，如果不设置的话，后面outTab变到deg的时候行名就不一致了
merge_eset <- merge_eset[!duplicated(merge_eset$ENTREZID),]
table(duplicated(merge_eset$ENTREZID)) 
rownames(merge_eset) <- merge_eset$ENTREZID


merge_ids <- merge_eset[19]  ###留取probe ID，symbol，entrezID
##merge_ids <- ids77287[match(merge_ids$symbol,ids77287$symbol),]    ###留取probe ID，symbol，entrezID


merge_eset <- merge_eset[,-19]    ###去掉entrezid
dim(merge_eset)
boxplot(merge_eset)

##整理成sva包所需要的数据类型，首先变成矩阵格式，接着将每个维度的名字变成列表，最后整理成矩阵格式的data，如下：
exp <- as.matrix(merge_eset)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)

##设置需要考虑的批次效应有哪些，一是不同的GSE号，二是不同的normal和钙化样本数量。
batchType=c(rep(1,18),rep(2,27))
modType=c(rep("Bicuspid",1),rep("control",3),rep("Tricuspid",3),rep("Bicuspid",1),rep("Tricuspid",4),rep("Bicuspid",1),rep("control",3),rep("Bicuspid",2),rep("Bicuspid",10),rep("Tricuspid",9),rep("control",8))
mod=model.matrix(~as.factor(modType))

###可以看到先是batchType，在是modType，最后采用model.matrix生成一个mod对象即可。
##一行命令搞定批次去除，主要是Combat命令的使用，将我们的需要考虑的batchType，mod，作为参数填充，即可
outTab=data.frame(ComBat(data,batchType,mod,par.prior = TRUE))
dim(outTab)
outTab[1:4,1:4]
boxplot(outTab)
#####
####??去完批次效应后发现组内仍有不齐，是否需要再次组内normalize？
if(T){outTab =normalizeBetweenArrays(outTab)
outTab[1:4,1:4]
boxplot(outTab)}


##############################
##############################
group_list=c(rep("Bicuspid",1),rep("control",3),rep("Tricuspid",3),rep("Bicuspid",1),rep("Tricuspid",4),rep("Bicuspid",1),rep("control",3),rep("Bicuspid",2),rep("Bicuspid",10),rep("Tricuspid",9),rep("control",8))
group_list
group_list=factor(group_list,levels = c("control","Tricuspid","Bicuspid"),ordered = T)

library(limma)
design=model.matrix(~0+group_list)
colnames(design)=levels(group_list)
px = levels(group_list)
contrast.matrix <- makeContrasts(paste0(px[3],"-",px[1]), 
                                 paste0(px[2],"-",px[1]), 
                                 paste0(px[3],"-",px[2]), 
                                 levels=design)


fit <- lmFit(outTab, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg = lapply(1:ncol(design), function(i){
  topTable(fit2,coef = i,number = Inf)
})
names(deg) = colnames(contrast.matrix)

################################
###################################
for (i in 1:length(deg)) {
  logFC_t=log2(2)
  P.Value_t = 0.01
  k1 = (deg[[i]]$P.Value < P.Value_t)&(deg[[i]]$logFC < -logFC_t)
  k2= (deg[[i]]$P.Value < P.Value_t)&(deg[[i]]$logFC > logFC_t)
  change = ifelse(k1,"down",ifelse(k2,"up","stable"))
  deg[[i]] <- mutate(deg[[i]],change)
}
deg_BV_con<- deg$`Bicuspid-control`
deg_TV_con <- deg$`Tricuspid-control`
deg_BV_TV <- deg$`Bicuspid-Tricuspid`


##########################
######################给DEG加symbol，EntrezID
deg_BV_con$ENTREZID <- rownames(deg_BV_con)
deg_TV_con$ENTREZID <- rownames(deg_TV_con)
deg_BV_TV$ENTREZID <- rownames(deg_BV_TV)

table(deg$`Bicuspid-control`$change)
table(deg$`Tricuspid-control`$change)
table(deg$`Bicuspid-Tricuspid`$change)

# 添加一列gene symbol
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(clusterProfiler)

id2symbol <- bitr(rownames(deg_BV_con), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
head(id2symbol)
symbol <- rep("NA",time=nrow(deg_BV_con))
symbol[match(id2symbol[,1],rownames(deg_BV_con))] <- id2symbol[,2]
deg_BV_con <- cbind(rownames(deg_BV_con),symbol,deg_BV_con)
colnames(deg_BV_con)[1] <- "GeneID"
head(deg_BV_con)

id2symbol <- bitr(rownames(deg_TV_con), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
head(id2symbol)
symbol <- rep("NA",time=nrow(deg_TV_con))
symbol[match(id2symbol[,1],rownames(deg_TV_con))] <- id2symbol[,2]
deg_TV_con <- cbind(rownames(deg_TV_con),symbol,deg_TV_con)
colnames(deg_TV_con)[1] <- "GeneID"
head(deg_TV_con)

id2symbol <- bitr(rownames(deg_BV_TV), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
head(id2symbol)
symbol <- rep("NA",time=nrow(deg_BV_TV))
symbol[match(id2symbol[,1],rownames(deg_BV_TV))] <- id2symbol[,2]
deg_BV_TV <- cbind(rownames(deg_BV_TV),symbol,deg_BV_TV)
colnames(deg_BV_TV)[1] <- "GeneID"
head(deg_BV_TV)


# 保存
write.table(deg_BV_con,"deg_BV_con_all.xls",row.names = F,sep="\t",quote = F)
write.table(deg_BV_con,"deg_TV_con_all.xls",row.names = F,sep="\t",quote = F)
write.table(deg_BV_con,"deg_BV_TV_all.xls",row.names = F,sep="\t",quote = F)



####################################################
#####################################################
#1.volcano---
library(dplyr)
library(ggplot2)
dat  = deg_BV_TV

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 20))
p

######################################################################
######################################################################