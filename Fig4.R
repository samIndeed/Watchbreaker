library(Seurat)

pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")

Idents(object = pbmc)=pbmc[["CMOassignment"]] 
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
# Filter out submerged conditions with few cells
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309","CMO311"))

pbmc <- subset(pbmc, cells = selected_c)
pbmc$condition=Idents(pbmc)

# Cluster and UMAP

pbmc <- RunPCA(pbmc, verbose = FALSE)

set.seed(100)

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
# Filter out mTO-FTM condition
# 
# Idents(object = pbmc)=pbmc[["CMOassignment"]] 
# 
# selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309"))
# 
# pbmc <- subset(pbmc, cells = selected_c)
# 
pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
pbmc$condition=Idents(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
DimPlot(pbmc, label = F)

# We can filter out the cluster with Ddit3 and Herpud1 as these are stressed cells
selected_c <- WhichCells(pbmc, idents=c("7"),invert = T)

pbmc <- subset(pbmc, cells = selected_c)

pbmc=RenameIdents(pbmc,c("0"="cTECIII","1"="cTECI","2"="mTECIII","3"="mcTEC","4"="mTECII","5"="cTECneg","6"="mcTEC-Prolif","8"="TEC-Tuft"))




# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

pbmc$seurat_clusters=Idents(pbmc)

saveRDS(pbmc,"~/Projects/WatchBreaker/TECFig4.rds")
pbmc=readRDS("~/Projects/WatchBreaker/TECFig4.rds")

DimPlot(pbmc, label = F)
# DimPlot(pbmc, label = F,group.by = "seurat_clusters")


DimPlot(pbmc, label = F,group.by = "condition",cols=c("Purple","Light Grey","Light Grey","Light Grey","Light Grey","Light Grey"))

pdf(file = "~/Projects/WatchBreaker/Figures/Fig4A.pdf")
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Purple","Light Grey","Light Grey","Light Grey","Light Grey","Light Grey"))

dev.off()
library(ggplot2)

# pdf(file = "~/Projects/WatchBreaker/Fig3B.pdf")
# ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
# dev.off()
pdf(file = "~/Projects/WatchBreaker/Fig4c.pdf")
VlnPlot(object = pbmc,ncol=3, features = c('Foxn1','Psmb11','Dll4','Isg15','H2-K1','H2-Aa'),group.by ="condition",log = T)
dev.off()


pdf(file = "~/Projects/WatchBreaker/Figures/FigS3.pdf",width=15,height=10)
FeaturePlot(pbmc,ncol=5, features = c("Foxn1","Psmb11","Dll4","Isg15",'H2-K1',"H2-Aa","Mki67","Krt14","Krt19","Aire","Fezf2","Slpi","Pou2f3","Krt5","Krt8"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()

VlnPlot(object = pbmc, features = c('Foxn1'),group.by ="condition",log = F)
VlnPlot(object = pbmc, features = c('Icam1'),group.by ="condition",log = F)
VlnPlot(object = pbmc, features = c('H2-K1'),group.by ="condition",log = F)
VlnPlot(object = pbmc, features = c('H2-Aa'),group.by ="condition",log = F)
VlnPlot(object = pbmc, features = c('Isg15'),group.by ="condition",log = F)
VlnPlot(object = pbmc, features = c('Krt14'),group.by ="condition",log = F)

# pbmc=RenameIdents(pbmc,c("0"="cTEC high","1"="cTEC low","2"="mTEC III","3"="mTEC I","4"="mcTEC II","5"="mcTEC","6"="cTEC II","7"="TEC-tuft"))

FeaturePlot(pbmc, features = c("Mki67","H2-Aa"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Stat2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("H2-K1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cd83"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Bmp7"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Il2ra"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Gli2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ephb2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Jag1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Enpep"), min.cutoff = "q05", max.cutoff = "q95")


ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
VlnPlot(object = pbmc,ncol=4, features = c('Foxn1','Psmb11','Dll4','Krt14','Aire','Fezf2','Slpi','Pou2f3'),group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Zeb2',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Tbata',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Nfkbia',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Tgfbr2',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Aire',group.by ="condition",log = T,cols=c("Dark Red","Light Green","Dark Green","Purple"))
VlnPlot(object = pbmc, features = 'Fezf2',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Slpi',group.by ="condition",log = T)
VlnPlot(object = pbmc, features = 'Pou2f3',group.by ="condition",log = T)


DimPlot(pbmc, label = F)

markers <- FindMarkers(object = pbmc, ident.1 = 2,group.by = "seurat_clusters")
head(x = markers,30)


# pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="RFTOC (FTOC proportions)","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))

# Try to find correlations with Foxn1 and Isg15

Idents(object = pbmc)=pbmc[["CMOassignment"]] 
selected_c <- WhichCells(pbmc, idents=c("mTO","mTO without FTM"),invert = F)


selected_c <- WhichCells(pbmc, idents=c("cTECI","cTECIII","cTECneg","mcTEC","mcTEC-Prolif"),invert = F)
selected_c <- WhichCells(pbmc, idents=c("cTECI","cTECIII","cTECneg"),invert = F)

pbmc <- subset(pbmc, cells = selected_c)

matrix<-pbmc[["SCT"]]@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Foxn1",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})


gene2<-as.numeric(matrix_mod["Stat1",])
correlations2<-apply(matrix_mod,1,function(x){cor(gene2,x)})

df=data.frame(cbind(correlations,correlations2))

keep=df[df$correlations>0.3|df$correlations2>0.3,]
keep=df[c("Trp63","Tbata","H2-D1","Cd74","Ciita","Stat1","Foxn1","Isg15","Psmb11","Dll4","Cd83","H2-Aa","H2-K1","Prss16","Cxcl12","Ccl25","B2m","Ifi35","Psmb8","Ifitm3","Tap1"),]
keep$pathway="Antigen Presentation Class I"
keep[c("Cd74","Ciita","H2-Aa"),]$pathway="Antigen Presentation Class II"
keep[c("Stat1","Isg15","Ifi35","Psmb8","Ifitm3","Tap1"),]$pathway="Interferon Signalling"
keep[c("Trp63","Tbata","Psmb11","Dll4","Cd83","Prss16","Cxcl12","Ccl25","Foxn1"),]$pathway="Foxn1 Targets"
library(ggplot2)
library(ggrepel)
p=ggplot(keep, aes(x=correlations, y=correlations2,colour=pathway)) +
  geom_point() + 
  geom_text_repel(label=rownames(keep))+theme_classic()+ylim(c(0,1))+xlab("Correlation with Foxn1")+ylab("Correlation with Stat1")+ggtitle("mTO Cells")
p


# Let's get a volcano plot of mTO vs mTO without FTM

pbmc=readRDS("~/Projects/WatchBreaker/TECFig4.rds")


# Idents(object = pbmc)=pbmc[["CMOassignment"]] 
# # Merge 312 and 311
# pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
# 
# pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
# pbmc$condition=Idents(pbmc)


markers <- FindMarkers(object = pbmc, ident.1 = "mTO", ident.2 = "mTO without FTM")
head(x = markers,30)
library(EnhancedVolcano)

EnhancedVolcano(markers , 
                rownames(markers ),
                pCutoff = 10e-11,
                FCcutoff = 1.5,
                title = 'mTO vs mTO without FTM',
                subtitle=NULL,
                legendPosition=NULL,
                drawConnectors = TRUE,
                gridlines.major = F,
                gridlines.minor = F,
                x ="avg_log2FC", 
                y ="p_val_adj",selectLab="Trp63")


markers["Trp63",]
# See if Dll4+Icam+ or Dll4+Cd83+ cells correlate with DP cells


pbmc$Icam1Dll4="Icam1-Dll4-"
selected_c <- WhichCells(pbmc, expression = Icam1>1.75)
pbmc@meta.data[selected_c,]$Icam1Dll4="Icam1+Dll4-"
selected_c <- WhichCells(pbmc, expression = Dll4>1)
pbmc@meta.data[selected_c,]$Icam1Dll4="Icam1-Dll4+"
pbmcDl=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmcDl, expression = Icam1>1.75)
pbmc@meta.data[selected_c,]$Icam1Dll4="Icam1+Dll4+"

ggplot(pbmc@meta.data, aes(x=condition,fill=Icam1Dll4)) + geom_bar(position = "dodge")+theme_classic()+scale_y_log10()


pbmc$Cd83Dll4="Cd83-Dll4-"
selected_c <- WhichCells(pbmc, expression = Cd83>2)
pbmc@meta.data[selected_c,]$Cd83Dll4="Cd83+Dll4-"
selected_c <- WhichCells(pbmc, expression = Dll4>1)
pbmc@meta.data[selected_c,]$Cd83Dll4="Cd83-Dll4+"
pbmcDl=subset(pbmc, cells = selected_c)
selected_c <- WhichCells(pbmcDl, expression = Cd83>2)
pbmc@meta.data[selected_c,]$Cd83Dll4="Cd83+Dll4+"

ggplot(pbmc@meta.data, aes(x=condition,fill=Cd83Dll4)) + geom_bar(position = "dodge")+theme_classic()+scale_y_log10()

# Make violin plots split by cluster

VlnPlot(object = pbmc, features = c('Mki67'),group.by ="condition",split.by="seurat_clusters",log = F)
VlnPlot(object = pbmc, features = c('Mki67'),group.by ="seurat_clusters",split.by="condition",log = F)
VlnPlot(object = pbmc, features = c('Foxn1'),group.by ="seurat_clusters",split.by="condition",log = F)
VlnPlot(object = pbmc, features = c('Isg15'),group.by ="seurat_clusters",split.by="condition",log = F)
VlnPlot(object = pbmc, features = c('Foxn1'),group.by ="condition",split.by="seurat_clusters",log = F)
VlnPlot(object = pbmc, features = c('Mki67'),group.by ="condition",split.by="seurat_clusters",log = F)
