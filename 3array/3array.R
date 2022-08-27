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


gse = "GSE77287"
eSet <- getGEO(gse, destdir = '.', getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
#exp = log2(exp+1)
boxplot(exp)
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

eSet[[1]]@annotation
#http://www.bio-info-trainee.com/1399.html
#hgu133plus2
if(F){
  if(!require(hgu133a.db))BiocManager::install("hgu133a.db")
  library(hgu133a.db)
  ls("package:hgu133a.db")
  ids <- toTable(hgu133aSYMBOL)
  head(ids)
}else if(F){
  getGEO(gpl)
  ids = data.table::fread(paste0(gpl,".soft"),header = T,skip = "ID",data.table = F)
  ids = ids[c("ID","Gene Symbol"),]
  colnames(ids) = c("probe_id")
}else if(T){
  ids = idmap(gpl,type = "bioc")
}

#http://www.bio-info-trainee.com/1399.html
#hgu133plus2   ###这个GPL570  如果是GPL16686需要用GPL16686_with_symbol.txt
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
#ls("package:hgu133plus2.db")
#ids <- toTable(hgu133plus2SYMBOL)
GPL16686<-read.table("GPL16686_with_symbol.txt",sep="\t",header = T)
ids<-data.frame(GPL16686$ID,GPL16686$symbol)
head(ids)

#1.加probe_id列
exp <- as.data.frame(exp)
exp <- mutate(exp,probe_id=rownames(exp))
exp <- merge(exp,ids,by.x="probe_id",by.y="GPL16686.ID")
##exp <- inner_join(exp,ids,by="probe_id")   ##其他数据集用此条命令
exp <- exp[!duplicated(exp$GPL16686.symbol),]


##4.加ENTREZID列，后面富集分析要用
s2e <- bitr(unique(exp$GPL16686.symbol), fromType = "SYMBOL",toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
exp <- inner_join(exp,s2e,by=c("GPL16686.symbol"="SYMBOL"))
colnames(exp)[8] <- c("symbol")


exp77287 <- exp
save(exp77287,file = "exp_77287_symbol.RData")
######################################################################
#####################################################################
###################################################################
###########数据集准备完毕
rm(list = ls())
library("sva")
library("dplyr")
library("limma")
load(file = "exp_51472_symbol.RData")
load(file = "exp_12644_symbol.RData")
load(file = "exp_77287_symbol.RData")

if(F){######设置各个数据集的行名，去除重复的基因，如果不设置的话，后面outTab变到deg的时候行名就不一致了
      ######放在后面做了
table(duplicated(exp77287$symbol))
exp77287 <- exp77287[!duplicated(exp77287$symbol),]
rownames(exp77287) <- exp77287$symbol

table(duplicated(exp51472$symbol))
exp51472 <- exp51472[!duplicated(exp51472$symbol),]
rownames(exp51472) <- exp51472$symbol

table(duplicated(exp12644$symbol))
exp12644 <- exp12644[!duplicated(exp12644$symbol),]
rownames(exp12644) <- exp12644$symbol
}


exp51472[1:4,1:4]
ids51472 <- exp51472[,c(16:18)]
exp51472=log2(exp51472[,c(1:15)]+1)
boxplot(exp51472)
exp51472$symbol <- ids51472$symbol      ##增加symbol一列



exp12644[1:4,1:4]
boxplot(exp12644[,c(1:20)])
exp12644 <- normalizeBetweenArrays(exp12644[,c(1:20)])
boxplot(exp12644)
exp12644 <- as.data.frame(exp12644)
exp12644$symbol <- ids51472$symbol      ##增加symbol一列,同一个平台公用一个ids


exp77287[1:4,1:4]
ids77287 <- exp77287[,c(1,8,9)]
boxplot(exp77287[,c(2:7)])
exp77287 <- normalizeBetweenArrays(exp77287[,c(2:7)])
boxplot(exp77287)
exp77287 <- as.data.frame(exp77287)
exp77287$symbol <- ids77287$symbol       ##增加symbol一列


merge_eset1=inner_join(exp51472,exp12644,by="symbol")
merge_eset=inner_join(merge_eset1,exp77287,by="symbol")
dim(merge_eset)
table(duplicated(merge_eset$symbol))    ##设置各个数据集的行名，去除重复的基因，如果不设置的话，后面outTab变到deg的时候行名就不一致了
merge_eset <- merge_eset[!duplicated(merge_eset$symbol),]
table(duplicated(merge_eset$symbol)) 
rownames(merge_eset) <- merge_eset$symbol


merge_ids <- merge_eset[16]  ###留取probe ID，symbol，entrezID
merge_ids <- ids77287[match(merge_ids$symbol,ids77287$symbol),]    ###留取probe ID，symbol，entrezID


merge_eset <- merge_eset[,-c(6:10,16)]    ###去掉symbol和GSE51472中的硬化组
dim(merge_eset)
boxplot(merge_eset)

##整理成sva包所需要的数据类型，首先变成矩阵格式，接着将每个维度的名字变成列表，最后整理成矩阵格式的data，如下：
exp <- as.matrix(merge_eset)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)

##设置需要考虑的批次效应有哪些，一是不同的GSE号，二是不同的normal和钙化样本数量。
batchType=c(rep(1,10),rep(2,20),rep(3,6))
modType=c(rep("normal",5),rep("calcification",5),rep("normal",5),rep("calcification",5),rep("normal",5),rep("calcification",5),rep("calcification",3),rep("normal",3))
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
write.table(outTab,file = "3array_outTab.txt", sep ="", row.names =TRUE, col.names =TRUE, quote =TRUE)

##############################
##############################
group_list=c(rep("normal",5),rep("calcification",5),rep("normal",5),rep("calcification",5),rep("normal",5),rep("calcification",5),rep("calcification",3),rep("normal",3))
group_list
group_list=factor(group_list,levels = c("normal","calcification"),ordered = T)

#####################################
#########做PCA图
dat_pca=as.data.frame(t(outTab))
library(FactoMineR)
library(factoextra) 
# pca
dat.pca <- PCA(dat_pca, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,  labelsize = 5,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = F, # Concentration ellipses
                         legend.title = "Groups",
                         pointsize = 2,
                        repel = TRUE,
                        )+
  labs(x="Dim1(45.3%)",y="Dim2(8.6%)",
       title = "Individuals - PCA",
       subtitle = "normal vs calcification",
     # caption = "Data: NMMAPS",
       tag = "Fig. 1")+
  theme(axis.title = element_text(size = 15, color = "firebrick", face = "italic"),
   legend.text = element_text(size = 15, color = "firebrick", face = "italic"),
   legend.title = element_text(size = 15, color = "firebrick", face = "italic"),
   title = element_text(size = 15, color = "firebrick", face = "italic"))+
scale_color_discrete("Groups") +
  guides(color = guide_legend(override.aes = list(size = 5)))
pca_plot

################################################################################
###PCA修图
pca_plot <- fviz_pca_ind(dat.pca,  labelsize = 5,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = F, # Concentration ellipses
                         legend.title = "Groups",
                         pointsize = 5,
                         repel = TRUE,
)+
  
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 20))




save(pca_plot,file = "pca_plot3array.Rdata")


library(limma)
design=model.matrix(~group_list)
fit=lmFit(outTab,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

################################
###################################


  logFC_t=log2(1.5)
  P.Value_t = 0.01
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 =  (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)



##########################
######################给DEG加symbol，EntrezID
deg$symbol <- rownames(deg)
a <- merge_ids[match(deg$symbol,merge_ids$symbol),]    ########匹配顺序#####很关键
deg$ENTREZID <- a[,3]

table(deg$change)
write.csv(deg,file = "正确deg（log2(1.5)）.csv")
####################################################
#####################################################
#1.volcano---
library(dplyr)
library(ggplot2)
dat  = deg

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

#2 heatmap
cg = deg$symbol[deg$change !="stable"]
length(cg)
n=outTab[cg,]    ####注意是outTab
dim(n)
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   color = colorRampPalette(c("blue","white","red"))(30),
                                   scale = "row",
                                   cluster_cols = T, 
                                   annotation_col=annotation_col,
                                   fontsize=20)) 
heatmap_plot
ggsave(heatmap_plot,filename = "heatmap3array.png")
##################################################################################
##############################################################################
###############################################################################
##GO与kegg富集分析
save(deg,outTab,group_list,logFC_t,P.Value_t,file = "GO-KEGG analysis.Rdata")

rm(list = ls())
load("GO-KEGG analysis.Rdata")


library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")
#source表示运行整个kegg_plot_function.R脚本，里面是一个function
#以up_kegg和down_kegg为输入数据做图

#1.GO database analysis ----

#(1)输入数据
if(F){gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
}

if(T){
  library(readxl)
X70_gene_DEG <- read_excel("70 gene DEG.xlsx")
View(X70_gene_DEG)
X70_gene_DEG <- X70_gene_DEG[,-1]
rownames(X70_gene_DEG) <- X70_gene_DEG$symbol
X70_gene_DEG <- as.data.frame(X70_gene_DEG)
gene_up = X70_gene_DEG[X70_gene_DEG$change == 'up','ENTREZID'] 
gene_down = X70_gene_DEG[X70_gene_DEG$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = X70_gene_DEG[,'ENTREZID']
}


#(2)GO分析，分三部分
#以下步骤耗时很长，实际运行时注意把if后面的括号里F改成T
if(T){
  #细胞组分
  ego_CC <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #生物过程
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #分子功能：
  ego_MF <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  save(ego_CC,ego_BP,ego_MF,file = "ego_3array.Rdata")
}

write.csv(ego_BP,file = "go-BP.csv")
write.csv(ego_CC,file = "go-cc.csv")
write.csv(ego_MF,file = "go-MF.csv")

#(3)可视化
#条带图
barplot(ego_CC,showCategory=10,font.size = 50)
barplot(ego_BP,showCategory=10,font.size = 10)
barplot(ego_MF,showCategory=10)

#气泡图
dotplot(ego_CC)
dotplot(ego_BP)
dotplot(ego_MF)


## 2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
if(F){gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
}
#（2）对上调/下调/所有差异基因进行富集分析
#注意这里又有个F
if(T){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  save(kk.diff,kk.down,kk.up,file = "GSE3arraykegg.Rdata")
}


#(3)从富集结果中提取出结果数据框
kegg_diff_dt <- kk.diff@result
write.csv(kegg_diff_dt,file = "kegg_diff.csv")
#(4)按照pvalue筛选通路
#在enrichkegg时没有设置pvaluecutoff，在此处筛选
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)
#(5)可视化
g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg
kegg <- kegg_plot(up_kegg,down_kegg)
kegg
#g_kegg +scale_y_continuous(labels = c(20,15,10,5,0,5))
ggsave(g_kegg,filename = 'kegg_up_down.png')

