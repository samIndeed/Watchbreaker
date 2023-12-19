pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")

Idents(object = pbmc)=pbmc[["CMOassignment"]] 
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
# Filter out submerged conditions with few cells
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309","CMO311"))

pbmc <- subset(pbmc, cells = selected_c)
# Cluster and UMAP

pbmc <- RunPCA(pbmc, verbose = FALSE)

set.seed(100)

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
# Filter out mTO-FTM condition

Idents(object = pbmc)=pbmc[["CMOassignment"]] 

selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309"))

pbmc <- subset(pbmc, cells = selected_c)

pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO"))
pbmc$condition=Idents(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
DimPlot(pbmc, label = F)

# We can filter out the cluster with Ddit3 and Herpud1 as these are stressed cells
selected_c <- WhichCells(pbmc, idents=c("6"),invert = T)

pbmc <- subset(pbmc, cells = selected_c)

pbmc=RenameIdents(pbmc,c("0"="cTECI","1"="cTECIII","2"="mTECIII","3"="mTECII","4"="mcTEC-Prolif","5"="mcTEC","7"="TEC-Tuft"))




# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

pbmc$seurat_clusters=Idents(pbmc)

saveRDS(pbmc,"~/Projects/WatchBreaker/TECFig3.rds")
pbmc=readRDS("~/Projects/WatchBreaker/TECFig3.rds")

DimPlot(pbmc, label = F)
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Dark Red","Light Green","Dark Green","Purple"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Dark Red","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Green","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Dark Green","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","Purple"))

pdf(file = "~/Projects/WatchBreaker/Fig3A.pdf")
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Dark Red","Light Green","Dark Green","Purple"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Dark Red","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Green","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Dark Green","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","Purple"))

dev.off()
library(ggplot2)

pdf(file = "~/Projects/WatchBreaker/Fig3B.pdf")
ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
dev.off()
pdf(file = "~/Projects/WatchBreaker/Fig3c.pdf")
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
dev.off()

# pbmc=RenameIdents(pbmc,c("0"="cTEC high","1"="cTEC low","2"="mTEC III","3"="mTEC I","4"="mcTEC II","5"="mcTEC","6"="cTEC II","7"="TEC-tuft"))

FeaturePlot(pbmc, features = c("Ccl21a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt19"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt14"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Slpi"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt8"), min.cutoff = "q05", max.cutoff = "q95")


ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Dll4',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Psmb11',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Krt14',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Krt19',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Aire',group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Fezf2',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Slpi',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Pou2f3',group.by ="condition",log = T)


DimPlot(pbmc, label = F)

markers <- FindMarkers(object = pbmc, ident.1 = 2,group.by = "seurat_clusters")
head(x = markers,30)


# pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="RFTOC (FTOC proportions)","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
