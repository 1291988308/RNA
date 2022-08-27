#####ceRNA

### step1_get_data
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
library(AnnoProbe)
library(limma)
gse = "GSE87885"
eSet <- getGEO(gse, 
                 destdir = '.', 
                 getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
boxplot(exp)
exp = log2(exp+1)
boxplot(exp)
exp <- normalizeBetweenArrays(exp)
boxplot(exp)
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation
p = identical(rownames(pd),colnames(exp));p
if(!p) {exp = exp[,match(rownames(pd),colnames(exp))]}

save(gse,exp,pd,gpl,file = "step1_output.Rdata")

####  step2_group_id

rm(list = ls())
load("step1_output.Rdata")
group_list = ifelse(str_detect(pd$title,"ZC"),"control","calcified")
group_list=factor(group_list,levels = c("control","calcified"),ordered = T)


#2.ids 芯片注释----
gpl

  ids <- read.csv2("GSE87885_family.soft.gz",
                 encoding="UTF-8",skip = 60,header=T, quote="",
                 fill = T,sep = "\t")
  ids=ids[,-3]
save(group_list,ids,file = "step2_output.Rdata")

### step3_PCA
rm(list = ls())
load("step1_output.Rdata")
load("step2_output.Rdata")
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
################################################PCA修图
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
print(pca_plot)

ggsave(plot = pca_plot,filename = paste0(gse,"PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")
##########################################################################
###step4_deg
rm(list = ls())
load("step1_output.Rdata")
load("step2_output.Rdata")

library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
head(deg)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
#tibble::rownames_to_column()
head(deg)

colnames(ids)=c("probe_id",'ID')
#merge
deg <- inner_join(deg,ids,by="probe_id")
deg <- deg[!duplicated(deg$ID),]
head(deg)
deg=deg[( grepl('hsa-miR',deg$ID)|grepl("hsa-let",deg$ID)),]   #注意名称

#3.加change列：上调或下调，火山图要用
#logFC_t=mean(deg$logFC)+2*sd(deg$logFC)
logFC_t=log2(2)
change=ifelse(deg$P.Value>0.01,'stable', 
              ifelse( deg$logFC >logFC_t,'up', 
                      ifelse( deg$logFC < -logFC_t,'down','stable') )
)
deg <- mutate(deg,change)
head(deg)
table(deg$change)
write.csv(deg,file = "GSE87885DEG.csv")
save(logFC_t,deg,file = "step4_output.Rdata")


###############
rm(list = ls()) 
load(file = "step1_output.Rdata")
load(file = "step4_output.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg
P.Value_t=0.01
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  xlim(-4,4) + ylim(0,4)+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 20))
p

#####################
if(F){
  #自选基因
  for_label <- dat%>% 
    filter(symbol %in% c("TRPM3","SFRP1")) 
}
if(F){
  #p值最小的10个
  for_label <- dat %>% head(10)
}
if(T) {
  #p值最小的前3下调和前3上调
  x1 = dat %>% 
    filter(change == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = ID),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = volcano_plot,filename = paste0(gse,"volcano.png"))

#2.差异基因热图----

load(file = 'step2_output.Rdata')
if(T){
  #全部差异基因
  cg = deg$probe_id[deg$change !="stable"]
  length(cg)
}else{
  #取前30上调和前30下调
  x=deg$logFC[deg$change !="stable"] 
  names(x)=deg$probe_id[deg$change !="stable"] 
  cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))
  length(cg)
}
n=exp[cg,]
dim(n)

#作热图
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   color = colorRampPalette(c("blue","white","red"))(30),
                                   scale = "row",
                                   #cluster_cols = F, 
                                   annotation_col=annotation_col,
                                   fontsize=10)) 
heatmap_plot
ggsave(heatmap_plot,filename = paste0(gse,"heatmap.png"))
load("pca_plot.Rdata")
library(patchwork)
(pca_plot + volcano_plot +heatmap_plot)+ plot_annotation(tag_levels = "A")

