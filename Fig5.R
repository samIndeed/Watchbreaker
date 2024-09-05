library(Seurat)

pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")

# What is the median number of unique genes per cell

median(pbmc$nFeature_RNA)
median(pbmc$nFeature_SCT)
mean(pbmc$nFeature_SCT)

# And including mesenchymal cells?

pbmc2=readRDS("~/Projects/WatchBreaker/WB.rds")

median(pbmc2$nFeature_RNA)


# Merge 312 and 311
Idents(object = pbmc)=pbmc[["CMOassignment"]] 

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

pbmc=RenameIdents(pbmc,c("0"="cTECI","1"="cTECIII","2"="mTECIII","3"="mTECII","4"="mTECI-Prolif","5"="mTECI","7"="TEC-Tuft"))




# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

pbmc$seurat_clusters=Idents(pbmc)

saveRDS(pbmc,"~/Projects/WatchBreaker/TECFig5.rds")
pbmc=readRDS("~/Projects/WatchBreaker/Figures/TECFig5.rds")

DimPlot(pbmc, label = F)
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","#58a500","#01b6b9","#c850fe"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","#58a500","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","#01b6b9","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","#c850fe"))

pdf(file = "~/Projects/WatchBreaker/Fig3A.pdf",width = 5,height=5)
# DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","#58a500","#01b6b9","#c850fe"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("#ff5053","Light Grey","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","#58a500","Light Grey","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","#01b6b9","Light Grey"))
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","#c850fe"))

dev.off()
library(ggplot2)

pdf(file = "~/Projects/WatchBreaker/Fig3B.pdf")
ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
dev.off()
pdf(file = "~/Projects/WatchBreaker/Figures/Fig5F.pdf",width = 16,height=10)
VlnPlot(object = pbmc,ncol=5, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3',"Ctsl","Prss16","Cxcl12","Flt3l","Il7","H2-K1","H2-Aa"),group.by ="condition",log = T,cols=c("#ff5053","#58a500","#01b6b9","#c850fe"))
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/Fig5Fextra.pdf",width = 18,height=10)
VlnPlot(object = pbmc,ncol=6, features = c('Foxn1','Psmb11','Dll4',"Cd83",'Krt14','Aire','Fezf2','Slpi',"Spink5",'Pou2f3',"Ctsl","Prss16","Cxcl12","Flt3l","Il7","H2-K1","H2-Aa","Plet1"),group.by ="condition",log = T,cols=c("#ff5053","#58a500","#01b6b9","#c850fe"))
dev.off()



# pbmc=RenameIdents(pbmc,c("0"="cTEC high","1"="cTEC low","2"="mTEC III","3"="mTEC I","4"="mcTEC II","5"="mcTEC","6"="cTEC II","7"="TEC-tuft"))

FeaturePlot(pbmc, features = c("Hey1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt19"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt14"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Slpi"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Notch1"), min.cutoff = "q05", max.cutoff = "q95")


ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Sox9',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Psmb11',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Krt14',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Krt19',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Aire',group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Fezf2',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Slpi',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Ifit1',group.by ="condition",log = T)


DimPlot(pbmc, label = F)

markers <- FindMarkers(object = pbmc, ident.1 = 2,group.by = "seurat_clusters")
head(x = markers,30)


# pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="RFTOC (FTOC proportions)","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))


# Let's compare conditions
Idents(pbmc) <- pbmc$condition


markersa <- FindMarkers(object = pbmc, ident.1 = "mTO")
head(x = markersa,30)

markersa["Foxn1",]

markers <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="FTOC")
head(x = markers,30)

markers2 <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="RFTOC")
head(x = markers2,30)


markers3 <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="Sorted RFTOC")
head(x = markers3,30)

markers["Fezf2",]
markers2["Aire",]

markers["Dll4",]
markersa["Cd83",]

markers["H2-T23",]
markers2["H2-T23",]
markers3["H2-T23",]
markers["H2-K1",]
markers2["H2-K1",]
markersa["H2-K1",]
markers["Krt14",]
markers2["Krt14",]
markers3["Krt14",]
markers2["Psmb11",]
markers["Eno1",]
markers2["Eno1",]
markers3["Eno1",]
markersa["Eno1",]
markers["Higd1a",]
markers2["Eno1",]
markers3["Eno1",]
markersa["Eno1",]

intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100)))
intersect(intersect(rownames(head(x = markers,200)),rownames(head(x = markers2,200))),rownames(head(x = markers3,200)))

markers[intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100))),]
markers2[intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100))),]
markers3[intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100))),]


VlnPlot(object = pbmc, features =intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100))),group.by ="condition",log = T)
VlnPlot(object = pbmc,ncol=7,features =intersect(intersect(rownames(head(x = markers,100)),rownames(head(x = markers2,100))),rownames(head(x = markers3,100))),group.by ="condition",log = T)+ stat_compare_means(comparisons = list(c(1,4),c(2,4),c(3,4)),label = "p.signif")
VlnPlot(object = pbmc,ncol=6,features =intersect(intersect(rownames(head(x = markers,200)),rownames(head(x = markers2,200))),rownames(head(x = markers3,200))),group.by ="condition",log = T)+ stat_compare_means(comparisons = list(c(1,4),c(2,4),c(3,4)),label = "p.signif")
VlnPlot(object = pbmc,features ="H2-T23",group.by ="condition",log = T)+ stat_compare_means(comparisons = list(c(1,4),c(2,4),c(3,4)),label = "p.signif")
# Let's focus on mTECII cells
Idents(pbmc) <- pbmc$seurat_clusters


selected_c <- WhichCells(pbmc,  idents="mTECII")
pbmc=subset(pbmc, cells = selected_c)

Idents(pbmc) <- pbmc$condition

mmarkersa <- FindMarkers(object = pbmc, ident.1 = "mTO")


mmarkers <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="FTOC")
head(x = mmarkers,30)
head(x = mmarkers[mmarkers$avg_log2FC<0,],30)
head(x = mmarkers2[mmarkers2$avg_log2FC<0,],30)
head(x = mmarkers3[mmarkers3$avg_log2FC<0,],30)

mmarkers2 <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="RFTOC")
head(x = mmarkers2,30)


mmarkers3 <- FindMarkers(object = pbmc, ident.1 = "mTO",ident.2="Sorted RFTOC")
head(x = mmarkers3,30)

VlnPlot(object = pbmc, features =intersect(rownames(head(x = mmarkers,50)),rownames(head(x = mmarkers2,50))),group.by ="condition",log = T)

write.csv(mmarkers,"~/Projects/WatchBreaker/mTECIImTOvsFTOC.csv")
write.csv(mmarkers2,"~/Projects/WatchBreaker/mTECIImTOvsRFTOC.csv")
write.csv(mmarkers3,"~/Projects/WatchBreaker/mTECIImTOvsSortedRFTOC.csv")
write.csv(mmarkersa,"~/Projects/WatchBreaker/mTECIImTOvsAll3.csv")
write.csv(markers,"~/Projects/WatchBreaker/mTOvsFTOC.csv")
write.csv(markers2,"~/Projects/WatchBreaker/mTOvsRFTOC.csv")
write.csv(markers3,"~/Projects/WatchBreaker/mTOvsSortedRFTOC.csv")
write.csv(markersa,"~/Projects/WatchBreaker/mTOvsAll3.csv")
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
# Make volcano plots of mTECII mTO vs FTOC

