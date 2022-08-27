### step1_get_data
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
library(AnnoProbe)
gse = "GSE102249"
eSet <- getGEO(gse, 
               destdir = '.', 
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
boxplot(exp)
# exp = log2(exp+1)   #注意确认
#library(limma) 
#exp=normalizeBetweenArrays(exp)
#boxplot(exp)

#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation

p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

save(gse,exp,pd,gpl,file = "step1_output.Rdata")
write.table(pd,file = paste0(gse,'_pdata.txt'))

####  step2_group_id
rm(list = ls())
load("step1_output.Rdata")
#2.ids 芯片注释----
gpl

#http://www.bio-info-trainee.com/1399.html
#hgu133plus2
#if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
#library(hgu133plus2.db)
#illuminaHumanv4
if(!require(illuminaHumanv4.db))BiocManager::install("illuminaHumanv4")
library(illuminaHumanv4.db)
ls("package:illuminaHumanv4.db")
ids <- toTable(illuminaHumanv4SYMBOL)
head(ids)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
exp <- as.data.frame(exp)
exp <- mutate(exp,probe_id=rownames(exp))
head(exp)
#2.加symbol列，火山图要用
exp <- inner_join(exp,ids,by="probe_id")
head(exp)
#按照symbol列去重复
exp <- exp[!duplicated(exp$symbol),]
save(exp,ids,file = "step2_output.Rdata")

rm(list = ls())
load("step1_output.Rdata")
load("step2_output.Rdata")

a <- exp[exp$symbol=="XIST",]
b <- exp[exp$symbol=="RPS4Y1",]
a
b
write.csv(a,file = "XIST.csv")
write.csv(b,file = "RPS4Y1.csv")



#male <- c(2,	3,	5,	7,	9,	10,	11,	12,	13,	15,	16,	17,	19,	20,	22,	26,	27,	29,	30,	31,	33,	34,	36,	37,	38,	40,	47,	48,	49,	54,	55,	56,	57,	62,	64,	65,	66,	67,	68,	71,	74,	75,	76,	78,	79,	80,	83,	84,	85,	88,	89,	91,	92,	94,	96,	97,	100,	101,	103,	104,	105,	106,	108,	110,	114,	115,	116,	117,	118,	119,	120,	121,	122,	123,	124,	125,	126,	127,	129,	130,	131,	132,	133,	135,	136,	140,	142,	143,	144,	145,	146,	148,	150,	151,	152,	153,	154,	155,	156,	158,	159,	160,	161,	162,	163,	164,	165,	166,	167,	168,	169,	170,	172,	174,	175,	177,	179,	180,	181,	182)
#female <- c(1,	4,	6,	8,	14,	18,	21,	23,	24,	25,	28,	32,	35,	39,	43,	44,	45,	46,	50,	51,	52,	53,	58,	59,	60,	61,	63,	69,	70,	72,	73,	77,	81,	82,	86,	87,	90,	93,	95,	98,	99,	102,	107,	109,	111,	112,	113,	128,	134,	137,	138,	139,	141,	147,	149,	157,	171,	173,	176,	178,	183,	184,	185,	186,	187,	188,	189,	190,	191,	192,	193,	194,	195,	196,	197,	198,	199,	200,	201,	202,	203,	204,	205,	206,	207,	208,	209,	210,	211,	212,	213,	214,	215,	216,	217,	218,	219,	220,	221,	222,	223,	224,	225,	226,	227,	228,	229,	230,	231,	232,	233,	234,	235,	236,	237,	238,	239,	240)              


group_list=c("female",	"male",	"male",	"female",	"male",	"female",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"male",	"male",	"female",	"male",	"male",	"female",	"male",	"female",	"female",	"female",	"male",	"male",	"female",	"male",	"male",	"male",	"female",	"male",	"male",	"female",	"male",	"male",	"male",	"female",	"male",	"female",	"female",	"female",	"female",	"female",	"female",	"male",	"male",	"male",	"female",	"female",	"female",	"female",	"male",	"male",	"male",	"male",	"female",	"female",	"female",	"female",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"female",	"female",	"male",	"female",	"female",	"male",	"male",	"male",	"female",	"male",	"male",	"male",	"female",	"female",	"male",	"male",	"male",	"female",	"female",	"male",	"male",	"female",	"male",	"male",	"female",	"male",	"female",	"male",	"male",	"female",	"female",	"male",	"male",	"female",	"male",	"male",	"male",	"male",	"female",	"male",	"female",	"male",	"female",	"female",	"female",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"male",	"female",	"female",	"female",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"male",	"female",	"male",	"female",	"male",	"male",	"female",	"male",	"female",	"male",	"male",	"male",	"male",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female",	"female")
group_list=factor(group_list,levels = c("male","female"),ordered = T)
table(group_list)

symbol <- exp$symbol
probe_id <- exp$probe_id
save(group_list,symbol,probe_id,file = "step3_output.Rdata")


### step3_PCA
rm(list = ls())
load("step1_output.Rdata")
load("step3_output.Rdata")
dat=as.data.frame(t(exp))
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra) 
# pca的统一操作走起
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
print(pca_plot)
ggsave(plot = pca_plot,filename = paste0(gse,"PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")

#热图 
cg=names(tail(sort(apply(exp,1,sd)),1000))
n=exp[cg,]

#绘制热图
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(pheatmap)
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row")

dev.off()

#####step4-deg
rm(list = ls())
load("step1_output.Rdata")
load("step3_output.Rdata")


#差异分析，用limma包来做
#需要表达矩阵和group_list，不需要改
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)


#2.ids 芯片注释----
gpl

#http://www.bio-info-trainee.com/1399.html
#hgu133plus2
#if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
#library(hgu133plus2.db)
#illuminaHumanv4
if(!require(illuminaHumanv4.db))BiocManager::install("illuminaHumanv4")
library(illuminaHumanv4.db)
ls("package:illuminaHumanv4.db")
ids <- toTable(illuminaHumanv4SYMBOL)
head(ids)


#2.加symbol列，火山图要用
deg <- inner_join(deg,ids,by="probe_id")
head(deg)
#按照symbol列去重复
deg <- deg[!duplicated(deg$symbol),]


#3.加change列,标记上下调基因
#logFC_t=log2(1.5)
logFC_t=0.2
P.Value_t = 0.05
k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
deg <- deg[!duplicated(deg$symbol),]
save(group_list,deg,logFC_t,P.Value_t,file = "step4_output.Rdata")
table(deg$change)
write.csv(deg,file = "deg_102249.csv")
###############
rm(list = ls()) 
load(file = "step1_output.Rdata")
load(file = "step4_output.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(adj.P.Val))) +
     xlim(-1,1)+
    ylim(0,15
         )+
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


