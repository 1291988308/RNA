rm(list = ls())
#BiocManager::install("preprocessCore")
library("preprocessCore")
source("CIBERSORT.R")
if(T){result1 <- CIBERSORT('LM22.txt','outTab.txt', perm = 1000, QN = T)
result1 <- write.csv(result1,file = "result1.csv")}  ####已经运行过了，结果如下

library(dplyr)
library(tidyr)
library(tibble)
obj <- read.csv(file = 'result1.csv',head=T)
class(obj)
rownames(obj) <- obj$X
obj=obj[,-1]
dd1 <- obj %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  pivot_longer(cols=2:23,
               names_to= 'celltype',
               values_to = 'Proportion')


library(ggplot2)
ggplot(dd1,aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = 'stack',stat = 'identity')+
  theme_bw()+
  guides(fill=guide_legend(ncol=1))


dd1$Group <- c(rep("normal",110),rep("calcification",110),rep("normal",110),rep("calcification",110),rep("normal",110),rep("calcification",110),rep("calcification",66),rep("normal",66))
dd1 <- dd1[,c(7,1,5,6)]
TME_New <- dd1
plot_order = TME_New[TME_New$Group=="normal",] %>% 
  group_by(celltype) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(celltype)

## `summarise()` ungrouping output (override with `.groups` argument)

TME_New$celltype = factor(TME_New$celltype,levels = plot_order)


if(T){
  mytheme <- theme(plot.title = element_text(size = 20,color="black",hjust = 0.5),
                   axis.title = element_text(size = 20,color ="black"), 
                   axis.text = element_text(size= 20,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 20),
                   legend.title= element_text(size= 20)
  ) }
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
box_TME <- ggplot(TME_New, aes(x = celltype, y = Proportion))+ 
  labs(y="Cell composition",x= NULL,title = "calcification Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#EB7369", "#1CB4B8"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_TME

ggsave("免疫浸润.pdf",box_TME,height=20,width=30,unit="cm")

############################################################################
################################################################################
##相关性热图
rm(list = ls())
library(ggcorrplot)
obj <- read.csv(file = 'result1.csv',head=T)
class(obj)
rownames(obj) <- obj$X
immu_data=obj[,-1]
a <- immu_data
a <- a[,-c(23:25)]
corr <- round(cor(a), 2)
corr <- corr[-c(18,21),-c(18,21)]   ###删去  "NA"
head(corr[, 1:6])
p.mat <- cor_pmat(a)



head(p.mat[, 1:4])
ggcorrplot(corr)#method默认为square
ggcorrplot(corr, hc.order = F, outline.color = "white", lab = TRUE)
ggcorrplot(corr, hc.order = F, type = "lower", outline.color = "white")#下三角形
ggcorrplot(corr, hc.order = F, type = "lower", lab = TRUE)
ggcorrplot(corr, hc.order = F, type = "lower", p.mat = p.mat)
ggcorrplot(corr, p.mat = p.mat, hc.order=F, type = "lower", insig = "blank")

#ggcorrplot(corr, hc.order = TRUE, outline.color = "white", lab = TRUE)
#ggcorrplot(corr, hc.order = TRUE, type = "lower", outline.color = "white")#下三角形
#ggcorrplot(corr, hc.order = TRUE, type = "lower", lab = TRUE)
#ggcorrplot(corr, hc.order = TRUE, type = "lower", p.mat = p.mat)
#ggcorrplot(corr, p.mat = p.mat, hc.order=TRUE, type = "lower", insig = "blank")


#####################################################
########################################################
#基因与免疫细胞相关性分析
rm(list = ls())
library(xlsx)
##expr_data <- read.xlsx(file = 'outTab-exp-gene.xlsx',header = T,sheetName = "Sheet1",encoding = "UTF-8")
library(readxl)
outTab_exp_gene <- read_excel("outTab-exp-gene.xlsx")
View(outTab_exp_gene)
expr_data <- outTab_exp_gene
gene <- "SPP1"
y <- as.numeric(unlist(expr_data[,gene]))

### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation <- data.frame()
immu_data <- read.csv(file = 'result1.csv',head=T)
rownames(immu_data) <- immu_data$X
immu_data <- immu_data[,c(2:23)]
##immu_data <- immu_data [,-c(1,12,13,16,17)]    ##去掉标准差为0的列？？
### 2.批量把数据导出到容器
for(i in 1:length(colnames(immu_data))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(immu_data[,i]),y,method="spearman")
  ## 3.填充
  correlation[i,1] = colnames(immu_data)[i]
  correlation[i,2] = dd$estimate
  correlation[i,3] = dd$p.value
}
### 修改名称
colnames(correlation) <- c("cell","cor","p.value")




sig_gene <- c("SPP1","HMOX1","CD28")
#sig_gene <- c("PTEN","PTGS2")
library(psych)
x <- expr_data[,sig_gene]
y <- immu_data

y <- y[,-c(3,4,5,12,13,15,20,22)]######去掉标准差为0的列
#y=y[,-c(3,4,5,8,9)]######去掉标准差为0的列？？
library(psych)
d <- corr.test(x,y,use="complete",method = 'spearman')

r <- d$r
p <- d$p

library(ggcorrplot)
ggcorrplot(t(d$r), show.legend = T, 
           p.mat = t(d$p), digits = 2,lab_size = 2, sig.level = 0.05,insig = 'blank',lab = T)
f <- rbind(r,p)
write.csv(f,"CD28_HMOX1_SPP1_gene-cor_.csv")
c <- cbind(y,x)
write.csv(c,"CD28_HMOX1_SPP1_gene-immun_.csv")

