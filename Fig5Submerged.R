pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")

Idents(pbmc) <- pbmc$CMOassignment
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
# Filter out submerged conditions with few cells
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO302","CMO304","CMO303","CMO309"))

pbmc <- subset(pbmc, cells = selected_c)
# Cluster and UMAP

pbmc <- RunPCA(pbmc, verbose = FALSE)

set.seed(100)

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)


pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO302"="Submerged FTOC","CMO303"="Submerged RFTOC","CMO304"="RFTOC","CMO309"="mTO"))
pbmc$condition=Idents(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.29)
DimPlot(pbmc, label = F)

FeaturePlot(pbmc, features = c("Foxn1"), min.cutoff = "q05", max.cutoff = "q95")
DimPlot(pbmc, label = F,group.by = "condition")
# We can filter out the cluster with Ddit3 and Herpud1 as these are stressed cells
selected_c <- WhichCells(pbmc, idents=c("7"),invert = T)

pbmc <- subset(pbmc, cells = selected_c)

pbmc=RenameIdents(pbmc,c("0"="cTECIII","1"="cTECI","2"="mTECI","3"="mTECIII","4"="mTECII","5"="mTEC-Prolif"))



# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

pbmc$seurat_clusters=Idents(pbmc)

# saveRDS(pbmc,"~/Projects/WatchBreaker/TECFig3.rds")
# pbmc=readRDS("~/Projects/WatchBreaker/Figures/TECFig3.rds")

DimPlot(pbmc, label = F)
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","#58a500","#01b6b9","#c850fe"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","#58a500","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","#01b6b9","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","#c850fe"))

pdf(file = "~/Projects/WatchBreaker/Figures/Fig5subA.pdf",width = 5,height=5)
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","#58a500","#01b6b9","#c850fe","Orange"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","Light Grey","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","#58a500","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","#01b6b9","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","#c850fe","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","Light Grey","Orange"))

dev.off()
library(ggplot2)

pdf(file = "~/Projects/WatchBreaker/Figures/Fig5subB.pdf",width = 9,height=5)
ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
dev.off()


# pdf(file = "~/Projects/WatchBreaker/Fig3subC.pdf")
# VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
# dev.off()
pdf(file = "~/Projects/WatchBreaker/Fig3subC.pdf")
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T)
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/Fig3SubD.pdf")
VlnPlot(object = pbmc,ncol=3, features = c('Foxn1','Psmb11','Dll4','Stat1','H2-K1','H2-Aa'),group.by ="condition",log = T)
dev.off()

VlnPlot(object = pbmc, features = c('Slpi'),group.by ="condition",split.by="seurat_clusters",log = F)



pdf(file = "~/Projects/WatchBreaker/Figures/FigSubS3.pdf",width=15,height=10)
FeaturePlot(pbmc,ncol=5, features = c("Foxn1","Psmb11","Dll4","Isg15",'H2-K1',"H2-Aa","Mki67","Krt14","Krt19","Aire","Fezf2","Slpi","Pou2f3","Krt5","Krt8"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/FigSubS4.pdf")
FeaturePlot(pbmc, features = c("Vim"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()
# pbmc=RenameIdents(pbmc,c("0"="cTEC high","1"="cTEC low","2"="mTEC III","3"="mTEC I","4"="mcTEC II","5"="mcTEC","6"="cTEC II","7"="TEC-tuft"))

FeaturePlot(pbmc, features = c("Hey1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt19"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt14"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Slpi"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt8"), min.cutoff = "q05", max.cutoff = "q95")


ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Sox9',group.by ="condition",log = T)
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


# Let's get percentages of Aire+ cells
Idents(pbmc) <- pbmc$condition


selected_c <- WhichCells(pbmc,  idents="FTOC")
pbmc1=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmc1, expression = Aire>1)
length(selected_c)/dim(pbmc1)[[2]]

selected_c <- WhichCells(pbmc,  idents="Sorted RFTOC")
pbmc1=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmc1, expression = Aire>1)
length(selected_c)/dim(pbmc1)[[2]]

selected_c <- WhichCells(pbmc,  idents="RFTOC")
pbmc1=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmc1, expression = Aire>1)
length(selected_c)/dim(pbmc1)[[2]]

selected_c <- WhichCells(pbmc,  idents="mTO")
pbmc1=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmc1, expression = Aire>1)
length(selected_c)/dim(pbmc1)[[2]]

# Make plots with FTOC, RFTOC and mTO and one with RFTOC and sorted RFTOC
Idents(pbmc) <- pbmc$condition

selected_c <- WhichCells(pbmc,  idents=c("FTOC","RFTOC","mTO"))
pbmc1=subset(pbmc, cells = selected_c)

pdf(file = "~/Projects/WatchBreaker/Figures/Fig3CA.pdf",width = 5,height=5)
# DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("#ff5053","#01b6b9","#c850fe"))
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("#ff5053","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("Light Grey","#01b6b9","Light Grey"))
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","#c850fe"))

dev.off()

selected_c <- WhichCells(pbmc,  idents=c("RFTOC","Sorted RFTOC"))
pbmc1=subset(pbmc, cells = selected_c)

pdf(file = "~/Projects/WatchBreaker/Figures/Fig3CB.pdf",width = 5,height=5)
# DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("#01b6b9","#58a500"))
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("Light Grey","#01b6b9","Light Grey"))
DimPlot(pbmc1, label = F,group.by = "condition",cols=c("Light Grey","#58a500"))

dev.off()

# Make Bar plot

df <- data.frame(condition=c("FTOC", "Sorted RFTOC", "RFTOC","mTO"),
                 Percentage=c(11,17,21,5))

df$condition <- as.character(df$condition)
df$condition <- factor(df$condition, levels=unique(df$condition))

p <- ggplot(data=df, aes(x=condition, y=Percentage,fill=condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic()

p+ scale_fill_manual(values=c("#ff5053","#58a500","#01b6b9","#c850fe"))+ theme(legend.position="none")+ labs( 
                                                                                                             x="Condition", y = "Percentage of Aire+ cells")
