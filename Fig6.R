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

pbmc=RenameIdents(pbmc,c("0"="cTECIII","1"="cTECI","2"="mTECIII","3"="mTECI","4"="mTECII","5"="cTECneg","6"="mTECI-Prolif","8"="TEC-Tuft"))




# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

pbmc$seurat_clusters=Idents(pbmc)

saveRDS(pbmc,"~/Projects/WatchBreaker/Figures/TECFig6.rds")
pbmc=readRDS("~/Projects/WatchBreaker/Figures/TECFig6.rds")


# Make pdf

pdf(file = "~/Projects/WatchBreaker/Figures/Fig6D.pdf")
VlnPlot(object = pbmc,ncol=1, features = c("Mex3a"),group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = c("Stmn2"),group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = c("1700125H20Rik"),group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = c("Slc5a10"),group.by ="condition",log = T,split.by = "seurat_clusters")
dev.off()


DimPlot(pbmc, label = F)
# DimPlot(pbmc, label = F,group.by = "seurat_clusters")
FeaturePlot(pbmc, features = c("Mex3a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Bhlhe40"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Bnip3"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Zeb2"), min.cutoff = "q05", max.cutoff = "q95")
VlnPlot(object = pbmc, features = c('Rfng'),group.by ="condition",log = T)


DimPlot(pbmc, label = F,group.by = "condition",cols=c("Purple","Light Grey","Light Grey","Light Grey","Light Grey","Light Grey"))

pdf(file = "~/Projects/WatchBreaker/Figures/Fig4A.pdf",width = 5,height=5)
DimPlot(pbmc, label = F,group.by = "seurat_clusters")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("Light Grey","Light Grey","Light Grey","Light Grey","Purple"))
FeaturePlot(pbmc, features = c('Higd1a','Hif1a','Krt14','Dll4'), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Plet1"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()
library(ggplot2)

# pdf(file = "~/Projects/WatchBreaker/Fig3B.pdf")
# ggplot(pbmc@meta.data, aes(x=condition, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
# dev.off()
pdf(file = "~/Projects/WatchBreaker/Fig4c.pdf")
VlnPlot(object = pbmc,ncol=3, features = c('Bpifa1','Psmb11','Dll4','Stat1','H2-K1','H2-Aa'),group.by ="condition",log = T)
dev.off()


pdf(file = "~/Projects/WatchBreaker/Figures/FigS3.pdf",width=15,height=10)
FeaturePlot(pbmc,ncol=5, features = c("Foxn1","Psmb11","Dll4","Isg15",'H2-K1',"H2-Aa","Mki67","Krt14","Krt19","Aire","Fezf2","Slpi","Pou2f3","Krt5","Krt8"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()


pdf(file = "~/Projects/WatchBreaker/Figures/Fig6Bextra.pdf",width = 18,height=10)
VlnPlot(object = pbmc,ncol=6, features = c('Foxn1','Psmb11','Dll4',"Cd83",'Krt14','Aire','Fezf2','Slpi',"Spink5",'Pou2f3',"Ctsl","Prss16","Cxcl12","Flt3l","Il7","H2-K1","H2-Aa","Plet1"),group.by ="condition",log = T)
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/Fig6Dextra.pdf",width = 18,height=10)
VlnPlot(object = pbmc,ncol=3, features = c('Higd1a','Hif1a','Krt14'),group.by ="condition",log = T)
dev.off()

ifns=rownames(pbmc)[grep("Ifn",rownames(pbmc))]

# Show no IFns
pdf(file = "~/Projects/WatchBreaker/Figures/Fig5NoIFNs.pdf",width = 18,height=10)
VlnPlot(object = pbmc,ncol=6, features =ifns,group.by ="condition",log = T)
dev.off()

# # And including mesenchymal cells?
# 
# pbmc=readRDS("~/Projects/WatchBreaker/WB.rds")
# ifns=rownames(pbmc)[grep("Ifn",rownames(pbmc))]
# pdf(file = "~/Projects/WatchBreaker/Figures/AllCellsNoIFNs.pdf",width = 18,height=10)
# VlnPlot(object = pbmc,ncol=6, features =ifns,group.by ="CMOassignment",log = T)
# FeaturePlot(pbmc, features = c("Ifnz","Lama5"), min.cutoff = "q05", max.cutoff = "q95")
# dev.off()

# Just mTO

Idents(pbmc)=pbmc$condition
DimPlot(pbmc, label = F)
selected_c <- WhichCells(pbmc, idents=c("mTO","mTO without FTM"))

pbmc <- subset(pbmc, cells = selected_c)

pdf(file = "~/Projects/WatchBreaker/Figures/Fig6BextraextraJustmTO.pdf",width = 8,height=14)
VlnPlot(object = pbmc,ncol=5, features = c('Foxn1','Psmb11','Dll4',"Cd83",'Krt14','Aire','Fezf2','Slpi',"Spink5",'Pou2f3',"Ctsl","Prss16","Cxcl12","Flt3l","Il7","H2-K1","H2-Aa","Plet1","Ciita","Krt19"),group.by ="condition",log = T)
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/Fig6DextraJustmTO.pdf",width = 8,height=10)
VlnPlot(object = pbmc,ncol=6, features = c('Foxn1','Psmb11','Dll4',"Cd83",'Krt14','Aire','Fezf2','Slpi',"Spink5",'Pou2f3',"Ctsl","Prss16","Cxcl12","Flt3l","Il7","H2-K1","Stat1","Isg15"),group.by ="condition",log = T)
dev.off()




# Try to find correlations with Foxn1 and Stat1

# Idents(object = pbmc)=pbmc[["CMOassignment"]] 
# selected_c <- WhichCells(pbmc, idents=c("mTO","mTO without FTM"),invert = F)
# 
# 
# selected_c <- WhichCells(pbmc, idents=c("cTECI","cTECIII","cTECneg","mcTEC","mcTEC-Prolif"),invert = F)
# selected_c <- WhichCells(pbmc, idents=c("cTECI","cTECIII","cTECneg"),invert = F)

# pbmc <- subset(pbmc, cells = selected_c)

matrix<-pbmc[["SCT"]]@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Foxn1",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})


gene2<-as.numeric(matrix_mod["Stat1",])
correlations2<-apply(matrix_mod,1,function(x){cor(gene2,x)})

df=data.frame(cbind(correlations,correlations2))

# keep=df[df$correlations>0.3|df$correlations2>0.3,]
keep=df[c("Trp63","Stat2","Tbata","H2-D1","Cd74","Ciita","Stat1","Foxn1","Isg15","Psmb11","Dll4","Cd83","H2-Aa","H2-K1","Prss16","Cxcl12","Ccl25","B2m","Ifi35","Psmb8","Ifitm3","Tap1"),]
keep$pathway="Antigen Presentation Class I"
keep[c("Cd74","Ciita","H2-Aa"),]$pathway="Antigen Presentation Class II"
keep[c("Stat2","Stat1","Isg15","Ifi35","Psmb8","Ifitm3","Tap1"),]$pathway="Interferon Signalling"
keep[c("Trp63","Tbata","Psmb11","Dll4","Cd83","Prss16","Cxcl12","Ccl25","Foxn1"),]$pathway="Foxn1 Targets"
library(ggplot2)
library(ggrepel)
pdf(file = "~/Projects/WatchBreaker/Fig4d.pdf",width=7,height=4)
p=ggplot(keep, aes(x=correlations, y=correlations2,colour=pathway)) +
  geom_point() + 
  geom_text_repel(label=rownames(keep))+theme_classic()+ylim(c(0,1))+xlab("Correlation with Foxn1")+ylab("Correlation with Stat1")+ggtitle("All Cells")
p
dev.off()




# Let's get a volcano plot of mTO vs mTO without FTM

pbmc=readRDS("~/Projects/WatchBreaker/Figures/TECFig6.rds")

selected_c <- WhichCells(pbmc, idents=c("cTECneg","cTECI","mTECIII"))

pbmc <- subset(pbmc, cells = selected_c)

Idents(object = pbmc)=pbmc$CMOassignment
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))

pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
pbmc$condition=Idents(pbmc)


markers <- FindMarkers(object = pbmc, ident.1 = "mTO", ident.2 = "mTO without FTM")
head(x = markers,30)

write.csv(markers,"~/Projects/WatchBreaker/mTOvsmTO-FTMcTECnegcTECImTECIII.csv")
write.csv(markers,"~/Projects/WatchBreaker/mTOvsmTO-FTM.csv")
library(EnhancedVolcano)


markers$pathway="Other"
markers[c("H2-K1","B2m","H2-D1","H2-Q7"),]$pathway="Antigen Presentation Class I"
markers[c("Cd74","Ciita","H2-Aa"),]$pathway="Antigen Presentation Class II"
markers[c("Ifi27l2a","Stat2","Stat1","Isg15","Ifi35","Psmb8","Ifitm3","Tap1"),]$pathway="Interferon Signalling"
markers[c("Trp63","Tbata","Psmb11","Dll4","Cd83","Prss16","Cxcl12","Ccl25","Foxn1"),]$pathway="Foxn1 Targets"

markers$col="grey10"
markers[c("H2-K1","B2m","H2-D1","H2-Q7"),]$col="#ff5053"
markers[c("Cd74","Ciita","H2-Aa"),]$col="#58a500"
markers[c("Ifi27l2a","Stat2","Stat1","Isg15","Ifi35","Psmb8","Ifitm3","Tap1"),]$col="#01b6b9"
markers[c("Trp63","Tbata","Psmb11","Dll4","Cd83","Prss16","Cxcl12","Ccl25","Foxn1"),]$col="#c850fe"

keyvals.colour=markers$col
names(keyvals.colour)=markers$pathway

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
                max.overlaps =50,
                colCustom=keyvals.colour,
                shadeAlpha=1,
                colAlpha=3/4,
                boxedLabels=T,
                colConnectors="Grey",
                y ="p_val_adj",selectLab=c("H2-Q7","Ifi27l2a","Trp63","Stat2","Tbata","H2-D1","Cd74","Ciita","Stat1","Foxn1","Isg15","Psmb11","Dll4","Cd83","H2-Aa","H2-K1","Prss16","Cxcl12","Ccl25","B2m","Ifi35","Ifitm3","Tap1"))


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

# Make bar graph of ageing data proportions

props=matrix(rev(c(36,86,86,71,138,8,6,2)))
props=props*100/sum(props)
labs=rev(c("Mature cTEC","Perinatal cTEC","Intertypical TEC","Proliferating TEC","Mature mTEC","Post-Aire mTEC","Tuft mTEC","nTEC"))

rownames(props)=labs
colnames(props)="Week 1"
library(RColorBrewer)
coul <- brewer.pal(length(props), "Pastel2") 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

coul=rev(gg_color_hue(length(props)))
barplot(props,col=coul, xlab="Postnatal", 
        legend=labs)
barplot(props,col=coul, xlab="Postnatal")

props=matrix(rev(c(36,86,86,71,138,8,6)))
props=props*100/sum(props)
labs=rev(c("Mature cTEC","Perinatal cTEC","Intertypical TEC","Proliferating TEC","Mature mTEC","Post-Aire mTEC","Tuft mTEC"))

rownames(props)=labs
colnames(props)="Week 1"
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

coul=rev(gg_color_hue(length(props)))
barplot(props,col=coul, xlab="Postnatal", 
        legend=labs)
barplot(props,col=coul, xlab="Postnatal")


# GO analysis dot plot by cluster and condition

# CytoTRACE2

cytotrace2_result=cytotrace2(pbmc,is_seurat = T, slot_type = "data")
annotation <- data.frame(phenotype = pbmc@meta.data$seurat_clusters) %>% set_rownames(., colnames(pbmc))
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

DimPlot(cytotrace2_result, group.by = "CytoTRACE2_Potency")
plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Boxplot_byPheno

# Focus on mTECIII

DimPlot(pbmc)

selected_c <- WhichCells(pbmc, idents=c("mTECIII"))

pbmc <- subset(pbmc, cells = selected_c)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.2)
DimPlot(pbmc, label = F)


markers <- FindMarkers(object = pbmc, ident.1 = 1)
head(x = markers,30)
write.csv(markers,"~/Projects/WatchBreaker/mTECIIIHivsLo.csv")

pbmc=RenameIdents(pbmc,c("0"="mTECIII-Hi","1"="mTECIII-Lo"))

pbmc$seurat_clusters2=Idents(pbmc)
VlnPlot(object = pbmc, features = c('Spink5'),group.by ="condition",split.by="seurat_clusters2",log = F,split.plot = TRUE)
VlnPlot(object = pbmc, features = c('Slpi'),group.by="seurat_clusters2",log = F)
VlnPlot(object = pbmc, features = c('Notch1','Jag1','Krt19','Cldn3','Cldn4','Plet1','Wnt7b'),group.by ="seurat_clusters2",log = F,stack=T,flip=T,pt.size=20)
VlnPlot(object = pbmc, features = c('Wnt7b','Wnt11','Wnt9a','Wnt5b'),group.by ="seurat_clusters2",log = F,stack=T,flip=T,pt.size=20)
VlnPlot(object = pbmc, features = rownames(head(x = markers,15)),group.by ="seurat_clusters2",log = F,stack=T,flip=T,pt.size=20)

ggplot(pbmc@meta.data, aes(x=condition, fill=Idents(pbmc))) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())

# Make PCA of mTECII

# Subset to mTECII
selected_c <- WhichCells(pbmc, idents=c("mTECII"))

pbmc <- subset(pbmc, cells = selected_c)
# pbmc$condition=Idents(pbmc)

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
# pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
# pbmc$condition=Idents(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.32)
DimPlot(pbmc, label = F,reduction = "pca")
DimPlot(pbmc, label = F,group.by = "condition",reduction = "pca")
DimPlot(pbmc, label = F,group.by = "condition",cols=c("gray","gray","gray","purple"),reduction = "pca")
FeaturePlot(pbmc, features = c("Mex3a"), min.cutoff = "q05", max.cutoff = "q95",reduction = "pca")
FeaturePlot(pbmc, features = c("Fezf2"), min.cutoff = "q05", max.cutoff = "q95")

# Make pdf

pdf(file = "~/Projects/WatchBreaker/Figures/mTECIIPCA.pdf",height=4,width=5)
DimPlot(pbmc, label = F,group.by = "condition",reduction = "pca",pt.size = 1.5)
DimPlot(pbmc, label = F,group.by = "condition",cols=c("gray","gray","gray","purple"),reduction = "pca",pt.size = 1.5)
dev.off()

# Make volcano plot

library(EnhancedVolcano)

markers=read.csv("~/Projects/WatchBreaker/mTECIImTOvsRFTOC.csv",row.names = 1)


sig=row.names(markers[markers$p_val_adj<0.05,])
markers$pathway="Other"
markers[sig,]$pathway="Significant (Adjusted p-val<0.05)"
markers[c("Aire","Fezf2"),]$pathway="Aire"

markers$col="grey10"
markers[sig,]$col="#ff5053"
markers[c("Aire","Fezf2"),]$col="#01b6b9"

keyvals.colour=markers$col
names(keyvals.colour)=markers$pathway

pdf(file = "~/Projects/WatchBreaker/Figures/mTECIImTOvsRFTOC.pdf",height=7,width=10)
EnhancedVolcano(markers , 
                rownames(markers ),
                pCutoff = 3e-6,
                FCcutoff = 2,
                title = 'mTO vs RFTOC (mTECII)',
                subtitle=NULL,
                legendPosition=NULL,
                drawConnectors = T,
                gridlines.major = F,
                gridlines.minor = F,
                x ="avg_log2FC", 
                max.overlaps =500,
                colCustom=keyvals.colour,
                boxedLabels=T,
                colConnectors="Grey",
                y ="p_val",selectLab=c(sig,"Aire","Fezf2"))
dev.off()

markers=read.csv("~/Projects/WatchBreaker/mTOvsAll3.csv",row.names = 1)

markers$pathway="Other"
markers[sig,]$pathway="Significant (Adjusted p-val<0.05)"
markers[c("Aire","Fezf2"),]$pathway="TRA regulators"

markers$col="grey10"
markers[sig,]$col="#ff5053"
markers[c("Aire","Fezf2"),]$col="#01b6b9"

keyvals.colour=markers$col
names(keyvals.colour)=markers$pathway


EnhancedVolcano(markers , 
                rownames(markers ),
                pCutoff = 3e-6,
                FCcutoff = 2.5,
                title = 'mTO vs FTOC',
                subtitle=NULL,
                legendPosition=NULL,
                drawConnectors = T,
                gridlines.major = F,
                gridlines.minor = F,
                x ="avg_log2FC", 
                max.overlaps =500,
                colCustom=keyvals.colour,
                boxedLabels=T,
                colConnectors="Grey",
                y ="p_val",selectLab=c(sig,"Aire","Fezf2"))


markers[c("Aire","Fezf2"),]

pdf(file = "~/Projects/WatchBreaker/Figures/AireFezf2Stmn2.pdf",height=4,width=10)
VlnPlot(object = pbmc,ncol=1, features = "Aire",group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = "Fezf2",group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = "Stmn2",group.by ="condition",log = T,split.by = "seurat_clusters")
dev.off()


VlnPlot(object = pbmc,ncol=1, features = sig[2],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[3],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[4],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[5],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[6],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[7],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[8],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[9],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[10],group.by ="condition",log = T,split.by = "seurat_clusters")
VlnPlot(object = pbmc,ncol=1, features = sig[11],group.by ="condition",log = T,split.by = "seurat_clusters")


suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(org.Mm.eg.db))


geneList <- markers[markers$p_val_adj<0.05,]$p_val_adj
geneList <- markers$p_val_adj
xx <- as.list(org.Mm.egSYMBOL2EG)
names(geneList) <- xx[sig] # Convert to entrezgene IDs

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org , mapping = "org.Mm.eg.db")
resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")

tab2 <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)

head(tab2, 25)

allGO = genesInTerm(GOdata)
head(allGO[tab2$GO.ID],10)
yy <- as.list(org.Mm.egSYMBOL)
head(yy[allGO[tab2$GO.ID][[1]]])

