rm(list = ls())
options(stringsAsFactors = F)
if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray")
if(!require(AnnoProbe))devtools::install_github("jmzeng1314/AnnoProbe")
library(GEOquery)
library(stringr)
library(AnnoProbe)
library(ggplot2)
library(tinyarray)
gse = "GSE26953"

if(require(AnnoProbe)){
  if(!file.exists(paste0(gse,"_eSet.Rdata"))) geoChina(gse)
  load(paste0(gse,"_eSet.Rdata"))
  eSet <- gset
  rm(gset)
}else{
  eSet <- getGEO(gse, 
                 destdir = '.', 
                 getGPL = F)
}


#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
exp = log2(exp+1)
boxplot(exp)
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation

p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]


group_list=ifelse(str_detect(pd$title,"Fibrosa_Laminar"),"Fibrosa_Laminar",
                  ifelse(str_detect(pd$title,"Ventricularis_Laminar"),"Ventricularis_Laminar",
                  ifelse(str_detect(pd$title,"Fibrosa_Oscillatory"),"Fibrosa_Oscillatory","Ventricularis_Oscillatory")))
group_list=factor(group_list,levels = c("Fibrosa_Laminar","Ventricularis_Laminar","Fibrosa_Oscillatory","Ventricularis_Oscillatory"))


library(limma)
design=model.matrix(~0+group_list)
colnames(design)=levels(group_list)
px = levels(group_list)
contrast.matrix <- makeContrasts(paste0(px[2],"-",px[1]), 
                                 paste0(px[3],"-",px[2]), 
                                 paste0(px[4],"-",px[3]), 
                                 paste0(px[3],"-",px[1]),
                                 paste0(px[4],"-",px[1]),
                                 paste0(px[4],"-",px[2]),
                                 levels=design)
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg = lapply(1:ncol(contrast.matrix), function(i){    ###多组比较关键：用contrast.matrix替换了design
  topTable(fit2,coef = i,number = Inf)
})
names(deg) = colnames(contrast.matrix)



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


library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

for(i in 1:length(deg)){
  #1.加probe_id列
  deg[[i]] <- mutate(deg[[i]],probe_id=rownames(deg[[i]]))
  #2.id转换
  deg[[i]] <- inner_join(deg[[i]],ids,by="probe_id")
  deg[[i]] <- deg[[i]][!duplicated(deg[[i]]$symbol),]
  #logFC_t=with(deg[[i]],mean(abs(logFC)))
  logFC_t=0.58
  P.Value_t=0.05
  #3.change
  {
    change=ifelse(deg[[i]] $P.Value>P.Value_t,'stable', 
                  ifelse( deg[[i]]$logFC >logFC_t,'up', 
                          ifelse( deg[[i]]$logFC < -logFC_t,'down','stable') )
    )
    deg[[i]] <- mutate(deg[[i]],change)
    print(table(deg[[i]]$change))
    deg[[i]] <- mutate(deg[[i]],v = -log10(P.Value))
  }
  #4.加ENTREZID列，后面富集分析要用
  s2e <- bitr(unique(deg[[i]]$symbol), fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Hs.eg.db)
  deg[[i]] <- inner_join(deg[[i]],s2e,by=c("symbol"="SYMBOL"))
}
write.csv((deg$`Ventricularis_Laminar-Fibrosa_Laminar`),file = "Ventricularis_Laminar-Fibrosa_Laminar.csv")
write.csv((deg$`Fibrosa_Oscillatory-Ventricularis_Laminar`),file = "Fibrosa_Oscillatory-Ventricularis_Laminar.csv")
write.csv((deg$`Ventricularis_Oscillatory-Fibrosa_Oscillatory`),file = "Ventricularis_Oscillatory-Fibrosa_Oscillatory.csv")
write.csv((deg$`Fibrosa_Oscillatory-Fibrosa_Laminar`),file = "Fibrosa_Oscillatory-Fibrosa_Laminar.csv")
write.csv((deg$`Ventricularis_Oscillatory-Fibrosa_Laminar`),file = "Ventricularis_Oscillatory-Fibrosa_Laminar.csv")
write.csv((deg$`Ventricularis_Oscillatory-Ventricularis_Laminar`),file = "Ventricularis_Oscillatory-Ventricularis_Laminar.csv")
pca_plot <- draw_pca(exp,group_list)
vo = lapply(1:length(deg),function(k){
  draw_volcano(deg[[k]],
               lab = colnames(contrast.matrix)[k],
               pkg = 4,
               logFC_cutoff = logFC_t,
               pvalue_cutoff = P.Value_t,
               symmetry = T)
})
