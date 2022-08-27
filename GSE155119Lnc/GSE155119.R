

### step1_get_data
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
library(AnnoProbe)
gse = "GSE155119"
eSet <- getGEO(gse, 
                 destdir = '.', 
                 getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
boxplot(exp)
exp = log2(exp+1)
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
group_list = ifelse(str_detect(pd$title,"Control"),"control","calcified")
group_list=factor(group_list,levels = c("control","calcified"),ordered = T)


#2.ids 芯片注释----
gpl

  ids <- read.csv2("GSE155119_family.soft.gz",
                 encoding="UTF-8",skip = 60,header=T, quote="",
                 fill = T,sep = "\t")
  ids=ids[,-3]

   
   b<- read.table("GPL22120-25936.txt",header=T, 
                 fill = T,sep = "\t")
   
   
   lncr<-b[-c(1,2),]
   table(lncr$TYPE)
   lncr<-lncr[lncr$TYPE=="lncRNA",]
   lncr_s<-lncr[,c(1,4,32)]
   
lncr_s<-lncr_s[,-3]
colnames(lncr_s)<-c("SPOT_ID","ACCESSION")
library("dplyr")
ids2<-inner_join(ids,lncr_s,by="SPOT_ID")   
 
# tmp<-ids2$SPOT_ID
# 
# for (i in 1:nrow(ids))
# {
#   if(grepl("LNCV6",ids$SPOT_ID[i]))
#   { 
#     index<- which(tmp==ids$SPOT_ID[i])
#     if(!(length(index)==0)) {ids$SPOT_ID[i]<-ids2$ACCESSION[index]}
#   }
# }
#   
 save(group_list,ids,ids2,file = "step2_output.Rdata")





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

############################################################PCA修图
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


ggsave(plot = pca_plot,filename = paste0(gse,"PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")

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
deg=deg[( grepl('hsa_circ',deg$ID)|grepl("LNCV6",deg$ID)),]
deg_circ=deg[( grepl('hsa_circ',deg$ID)),]
deg_lncr=deg[(grepl("LNCV6",deg$ID)),]

colnames(ids2)<-c("name1","ID","ACCESSION")
ids2<-ids2[,-1]
deg_lncr=inner_join(deg_lncr,ids2,by="ID")



#3.加change列：上调或下调，火山图要用
#logFC_t=mean(deg$logFC)+2*sd(deg$logFC)
logFC_t=log2(1.5)

change_lncr=ifelse(deg_lncr$P.Value>0.05,'stable', 
              ifelse( deg_lncr$logFC >logFC_t,'up', 
                      ifelse( deg_lncr$logFC < -logFC_t,'down','stable') )
)


deg_lncr <- mutate(deg_lncr,change_lncr)



head(deg_lncr)
table(deg_lncr$change)

save(logFC_t,deg_circ,deg_lncr,change_lncr,file = "step4_output.Rdata")
write.csv(deg_lncr,file = "deg_lncr.csv")
#write.csv(deg_circ,file = "deg_circ.csv")
###############
#####
rm(list = ls()) 
load(file = "step1_output.Rdata")
load(file = "step4_output.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)

if(F){
dat  = deg_circ
P.Value_t=0.05
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=5, 
             aes(color=change_circ)) +
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
}

dat1  = deg_lncr
P.Value_t=0.05
p1 <- ggplot(data = dat1, 
             aes(x = logFC, 
                 y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=5, 
             aes(color=change_lncr)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 20))
p1

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
  x1 = dat1 %>% 
    filter(change_lncr == "up") %>% 
    head(3)
  x2 = dat1 %>% 
    filter(change_lncr == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}

volcano_plot <- p1 +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = ID),
    data = for_label,
    color="black"
  )
volcano_plot

if(F) {
  #p值最小的前3下调和前3上调
  x1 = dat %>% 
    filter(change_circ == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change_circ == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
  volcano_plot <- p +
    geom_point(size = 3, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = ID),
      data = for_label,
      color="black"
    )
  volcano_plot
  }



ggsave(plot = volcano_plot,filename = paste0(gse,"volcano.png"))

#2.差异基因热图----
rm(list = ls())
load(file = 'step1_output.Rdata')
load(file = 'step2_output.Rdata')
load(file = 'step4_output.Rdata')
if(F){
  #全部差异基因
  cg = deg_circ$probe_id[deg_circ$change_circ !="stable"]
  length(cg)
}else{
  #取前30上调和前30下调
  x=deg_circ$logFC[deg_circ$change_circ !="stable"] 
  names(x)=deg_circ$probe_id[deg_circ$change_circ !="stable"] 
  cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))
  length(cg)
}
n=exp[cg,]
dim(n)

if(T){
  #全部差异基因
  cg1 = deg_lncr$probe_id[deg_lncr$change_lncr !="stable"]
  length(cg1)
}else{
  #取前30上调和前30下调
  x=deg_lncr$logFC[deg_lncr$change_lncr !="stable"] 
  names(x)=deg_lncr$probe_id[deg_lncr$change_lncr !="stable"] 
  cg1=c(names(head(sort(x),30)),names(tail(sort(x),30)))
  length(cg1)
}
n=exp[cg1,]
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
