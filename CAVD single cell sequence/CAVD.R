set.seed(123)  #设置随机数种子，使结果可重复
# step0:数据载入
rm(list=ls())
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)

dir = c('fsdownload/SRR10134385/filtered_feature_bc_matrix',
        'fsdownload/SRR10134386/filtered_feature_bc_matrix',
        'fsdownload/SRR10134387/filtered_feature_bc_matrix',
        'fsdownload/SRR10134388/filtered_feature_bc_matrix',
        'fsdownload/SRR10134389/filtered_feature_bc_matrix',
        'fsdownload/SRR10134390/filtered_feature_bc_matrix' )
names(dir) = c('CAVD4',"CAVD3","CAVD2","CAVD1","health1","health2")

obj.list <- list()
#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir))
{
  counts <- Read10X(data.dir = dir[i])
  obj.list[[i]] <- CreateSeuratObject(counts, min.cells=3,min.features = 200)
}     #保留至少在三个细胞里面表达的基因 ;#保留至少表达200个基因的细胞
names(obj.list) <- names(dir)
obj.list
obj.list$CAVD4@meta.data$group <- 'CAVD4'
obj.list$CAVD3@meta.data$group <- 'CAVD3'
obj.list$CAVD2@meta.data$group <- 'CAVD2'
obj.list$CAVD1@meta.data$group <- 'CAVD1'
obj.list$health1@meta.data$group <- 'health1'
obj.list$health2@meta.data$group <- 'health2'





#质控QC，计算每种特征在每个细胞里面的百分比
dir.create("QC")
setwd("QC")
for (x in names(obj.list)){
  obj.list[[x]][["percent.MT"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^MT-")
  obj.list[[x]][["percent.RP"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^RP[SL]")
  obj.list[[x]][["percent.HB"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^HB[^(p)]")
}

#循环绘制每个样本的QC图
qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.HB", "percent.MT",  "percent.RP")
for(x in names(obj.list)){
  pdf(file=paste0("1_",x,"_quality_control.pdf"),width = 15,height=7)
  print(VlnPlot(obj.list[[x]], features = qc_feature, ncol = 5, pt.size = 0.5))
  dev.off()
}

#循环绘制特征之间的相互关系图
for(x in names(obj.list)){
  pdf(file=paste0("1_",x,"_feature_relationship.pdf"),width = 12,height=10)
  p1 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.HB")
  p3 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.MT")
  p4 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.RP")
  print(p1 + p2 + p3 + p4)
  dev.off()
}

#过滤细胞
obj.list <- lapply(obj.list, function(x) {
  subset(x, subset = nFeature_RNA > 500 & 
           nFeature_RNA < 4000 &
           percent.MT < 10)})
obj.list







#多样本整合
#归一化
obj.list <- lapply(obj.list,function(x) {
  NormalizeData(x)
})

#寻找高变异度基因
obj.list <- lapply(obj.list, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
obj.list
#整合成一个对象
#寻找锚点
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20)
#根据锚点来整合
obj.combined <- IntegrateData(anchorset = obj.anchors, dims = 1:20)
DefaultAssay(obj.combined) <- "RNA"   #更改默认数组
#对RNA的数据进行scale标准化
all.genes <- rownames(obj.combined[["RNA"]]@data)    #针对所有基因
length(all.genes)
obj.combined <- ScaleData(obj.combined, features = all.genes)

DefaultAssay(obj.combined) <- "integrated"   #更改默认数组
#对整合后的数据进行尺度变换
all.genes <- rownames(obj.combined[["RNA"]]@data)    #针对所有基因
length(all.genes)
obj.combined <- ScaleData(obj.combined, features = all.genes)






#对整合之后的数据进行降维
#主成分分析

#save(obj.combined,file = "obj.combined.RData")
saveRDS(obj.combined, "obj.combined.rds")
dir.create("PCA")
setwd("PCA")

obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = FALSE)
#可视化PCA降维之后的结果
pdf('2_PCA.pdf', width = 10, height = 6)
DimPlot(obj.combined, reduction = "pca")
dev.off()
#热图展示前5个主成分
pdf(file="2_PCtop5_Heatmap.pdf",width=14,height=20)
DimHeatmap(obj.combined, dims = 1:5, balanced = TRUE)
dev.off()
#热图展示第6-10个主成分
pdf(file="2_PC6-10_Heatmap.pdf",width=14,height=20)
DimHeatmap(obj.combined, dims = 1:5, balanced = TRUE)
dev.off()
#确定用几个主成分做后续分析
#Elbowplot
pdf('2_ElbowPlot.pdf', width = 10, height = 6)
ElbowPlot(obj.combined, ndims = 50)
dev.off()
#对整合之后的数据进行聚类
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:20)
obj.combined <- RunTSNE(obj.combined, reduction = "pca", dims = 1:20)
#SWNE Visualizing single-cell RNA-seq datasets with Similarity Weighted Nonnegative Embedding (tSNE可准确的捕获数据集的局部结构，但是会扭曲数据集的全局结构，比如簇之间的距离)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  obj.combined <- FindClusters(obj.combined,resolution = res)
} 
#设置不同的分辨率，观察分群效果(选择哪一个？)
#可视化不同分辨率，分群效果
apply(obj.combined@meta.data[,grep("integrated_snn",colnames(obj.combined@meta.data))],2,table)
pdf('2_umap_resolution_high.pdf', width = 18)
plot_grid(ncol = 3, DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.8") + 
            ggtitle("louvain_0.8"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.1") + 
            ggtitle("louvain_1"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.3") + 
            ggtitle("louvain_0.3"))
dev.off()
pdf('2_umap_resolution_low.pdf', width = 18)
plot_grid(ncol = 3, DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.01") + 
            ggtitle("louvain_0.01"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.1") + 
            ggtitle("louvain_0.1"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.2") + 
            ggtitle("louvain_0.2"))
dev.off()
pdf('2_Tree_diff_resolution.pdf', width = 10,height = 10)
clustree(obj.combined@meta.data, prefix = "integrated_snn_res.",layout = "sugiyama")
dev.off()
#接下来分析，按照分辨率为1进行
sel.clust = "integrated_snn_res.1"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
table(obj.combined@active.ident)
#绘制UMAP聚类图
pdf('3_UMAP_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot(obj.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(obj.combined, reduction = "umap", label = TRUE)
p1 + p2
dev.off()
pdf('3_tsne_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot(obj.combined, reduction = "tsne", group.by = "group")
p2 <- DimPlot(obj.combined, reduction = "tsne", label = TRUE)
p1 + p2
dev.off()
#分别展示两组的UMAP聚类图
pdf('3_tsne_multi_samples_split.pdf',width = 10,height = 6)
DimPlot(obj.combined, reduction = "tsne", split.by = "group",label = TRUE)
dev.off()
#查看细胞周期影响
pdf(file="3_Pcna_Mki67_expression.pdf",width=12)
FeaturePlot(obj.combined,features = c("Pcna","Mki67"),reduction = "tsne")  #s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
dev.off()
#计算细胞周期分值，判断这些细胞分别处于什么期
cc_gene=unlist(cc.genes)  #cc.genes 为细胞周期相关的基因，为list，有两个元素s.genes，g2m.genes
s.genes=cc.genes$s.genes
#s.genes<-tolower(s.genes)
s.genes<-capitalize(s.genes)
g2m.genes=cc.genes$g2m.genes
#g2m.genes<-tolower(g2m.genes)
g2m.genes<-capitalize(g2m.genes)
obj.combined <- CellCycleScoring(obj.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pdf(file="3_cellcycle_score.pdf",width=8)
obj.combined@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
dev.off()
#去除细胞周期影响
obj.combined1 <- ScaleData(obj.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj.combined))
save(obj.combined,obj.combined1,file = "PCA_cellcycle.RData")




#细胞鉴定
#鉴定marker基因
sel.clust = "integrated_snn_res.1"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
sample.markers <- FindAllMarkers(obj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers_cluster4 <- FindMarkers(object = obj.combined, ident.1 = 4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster4_top10 <- markers_cluster4 %>% top_n(n = 10, wt = avg_log2FC)
#保存每个cluster top10的marker基因
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))
write.csv(file="4_top10_marker_genes.csv",top10_table,row.names=F)
write.csv(file="sample_marker_genes.csv",sample.markers,row.names=F)
save(obj.combined,sample.markers,top10,top10_table,file = "cluster.RData")
#细胞亚群注释——命名
if(T){
  new.cluster.ids <- c("VDSC","VDSC","VDSC","VDSC",
                       "VDSC","VDSC",  "VDSC", 
                       "VIC", "MACrophage", "monocyte","VEC","LYmphocyte","VIC","VIC","VEC")
  names(new.cluster.ids) <- levels(obj.combined)
  obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
  table(obj.combined@active.ident)
}
#细胞亚群注释后的UMAP
pdf('cluster_annotations.pdf',width = 9,height = 7)
DimPlot(obj.combined, reduction = "tsne", label = TRUE,split.by = "group", pt.size = 0.5) + NoLegend()
dev.off()
save(obj.combined,file = "rename_cluster.RData")






cell.num <- table(Idents(obj.combined))
cell.freq <- round(prop.table(table(Idents(obj.combined)))*100,2)
cell.combined <- rbind(cell.num, cell.freq)
write.csv(file="combined_cell_counts_freq.csv",cell.combined)
#分组计算细胞数目和比例
cell.num.group <- table(Idents(obj.combined), obj.combined$group) 
colnames(cell.num.group)<-paste0(colnames(cell.num.group),'_cell_counts')
cell.freq.group <- round(prop.table(table(Idents(obj.combined), obj.combined$group), margin = 2) *100,2)
colnames(cell.freq.group)<-paste0(colnames(cell.freq.group),'_cell_Freq')
cell.group <- cbind(cell.num.group, cell.freq.group)
write.csv(file="group_cell_counts_freq.csv",cell.group)
#表格和UMAP展示细胞数目变化
pdf('tsne_multi_samples_split_anno.pdf',width = 25,height = 8)
p<-DimPlot(obj.combined, reduction = "tsne", split.by = "group",label=T,repel=T) + NoLegend()
tb <- tableGrob(cell.group)
plot_grid(p, tb,nrow=2,rel_widths=c(0.6,0.4))
dev.off()    

#火山图展示 文献：Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis
#两组间差异基因
#####将CAVD1，2，3，4合并为CAVD，将health1，2合并为health
obj.combined@meta.data$group2 <- ifelse(str_detect(obj.combined@meta.data[["group"]],"CAVD"),"CAVD","health")

pdf('CAVD_health_tsne.pdf',width = 12,height = 9)
p<-DimPlot(obj.combined, reduction = "tsne", split.by = "group2",label=T,repel=T,label.size = 10) + 
  theme(axis.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 20),
        title = element_text(size = 20))
print(p)
dev.off()  
####将每个细胞的identity转换成celltype.group
DefaultAssay(obj.combined) <- "RNA"
obj.combined$celltype <- Idents(obj.combined) ##保存每个细胞的细胞类型
obj.combined$celltype.group <- paste(Idents(obj.combined), obj.combined$group2, sep = "_")   ##细胞类型前面加上分组信息
Idents(obj.combined) <- "celltype.group"
##鉴定VIC和VEC细胞在MI和Health组中差异表达的基因
VIC.diff <- FindMarkers(obj.combined, ident.1 = "VIC_CAVD", ident.2 = "VIC_health", verbose = FALSE)
write.csv(file="VIC_Diff_exp_genes.csv",VIC.diff)
VEC.diff <- FindMarkers(obj.combined, ident.1 = "VEC_CAVD", ident.2 = "VEC_health", verbose = FALSE)
write.csv(file="VEC_Diff_exp_genes.csv",VEC.diff)
##提取出VIC细胞
Idents(obj.combined) <- 'celltype'
VIC <- subset(obj.combined, idents = "VIC")
Idents(VIC) <- "group2"
avg.VIC <- log1p(AverageExpression(VIC, verbose = FALSE)$RNA)
avg.VIC <- data.frame(avg.VIC ,gene=rownames(avg.VIC))
write.csv(file = "avg.VIC_CAVD_health.csv",avg.VIC)
##提取出VEC细胞
Idents(obj.combined) <- 'celltype'
VEC <- subset(obj.combined, idents = "VEC")
Idents(VEC) <- "group2"
avg.VEC <- log1p(AverageExpression(VEC, verbose = FALSE)$RNA)
avg.VEC <- data.frame(avg.VEC ,gene=rownames(avg.VEC))
write.csv(file = "avg.VEC_CAVD_health.csv",avg.VEC)


#######################################################################

##VIC火山图可视化差异基因
low<-floor(range(VIC.diff$avg_log2FC)[1])
high<-ceiling(range(VIC.diff$avg_log2FC)[2])
pdf('DEG_VIC.pdf',width = 12,height = 9)
print(EnhancedVolcano(VIC.diff,
                      title = 'VIC CAVD versus health',
                      lab = rownames(VIC.diff),
                      selectLab = "",    ###不标识基因
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      xlim = c(-2.5,2)))
dev.off()

##VEC火山图可视化差异基因
low<-floor(range(VEC.diff$avg_log2FC)[1])
high<-ceiling(range(VEC.diff$avg_log2FC)[2])
pdf('DEG_VEC.pdf',width = 12,height = 9)
print(EnhancedVolcano(VEC.diff,
                      title = 'VEC CAVD versus health',
                      lab = rownames(VEC.diff),
                      selectLab = "",    ###不标识基因
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      xlim = c(-2.5, 2)))
dev.off()
