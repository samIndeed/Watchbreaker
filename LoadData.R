
library(Seurat)
data_dir <- "~/Projects/WatchBreaker/raw_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir)

seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['CMO']] = CreateAssayObject(counts = data$`Multiplexing Capture`)

# CMOs
library(data.table)

cellCMOs <- fread("~/Projects/WatchBreaker/assignment_confidence_table.csv",select = c("Barcode","Assignment"))
table(cellCMOs$Assignment)

seurat_object_use <- subset(seurat_object, cells = cellCMOs$Barcode)


pbmc <- seurat_object_use

pbmc$CMOassignment=cellCMOs$Assignment


selected_c <- WhichCells(pbmc, expression = nFeature_RNA > 2000)

pbmc <- subset(pbmc, cells = selected_c)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt")


# selected_c <- WhichCells(pbmc, expression = percent.mt < 8)
# 
# pbmc <- subset(pbmc, cells = selected_c)

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)



pbmc <- RunPCA(pbmc, verbose = FALSE)

set.seed(100)

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3)

# pbmc1=pbmc
# pbmc=pbmc1

DimPlot(pbmc, label = TRUE) + NoLegend()
DimPlot(pbmc, label = TRUE,group.by = "CMOassignment")

saveRDS(pbmc,"~/Projects/WatchBreaker/WB.rds")



markers <- FindMarkers(object = pbmc, group.by="CMOassignment",ident.1 = "CMO309",ident.2 = "CMO311")
markers <- FindMarkers(object = pbmc, ident.1 = 13)
head(x = markers,20)
write.csv(markers,file="~/Projects/WatchBreaker/WBFTM+vsFTM-.csv")


library(enrichR)

plt=enrichr(row.names(markers[markers$avg_log2FC>0,]),databases="GO_Biological_Process_2021")
pltdown=enrichr(row.names(markers[markers$avg_log2FC<0,]),databases="GO_Biological_Process_2021")

plt$GO_Biological_Process_2021$Term[1:10]
pltdown$GO_Biological_Process_2021$Term[1:10]

write.csv(plt,file="~/Projects/WatchBreaker/WBFTM+vsFTM-GO.csv")
write.csv(pltdown,file="~/Projects/WatchBreaker/WBFTM-vsFTM+GO.csv")


FeaturePlot(pbmc, features = c("percent.mt"), min.cutoff = "q01", max.cutoff = "q99")
FeaturePlot(pbmc, features = c("nFeature_RNA"), min.cutoff = "q01", max.cutoff = "q95")
# FeaturePlot(pbmc, features = c("Tcrg-C1"), min.cutoff = "q05", max.cutoff = "q95")
# FeaturePlot(pbmc, features = c("ab_MHCII"), min.cutoff = "q05", max.cutoff = "q95")
# FeaturePlot(pbmc, features = c("ab_CD40"), min.cutoff = "q05", max.cutoff = "q95")
# FeaturePlot(pbmc, features = c("ab_CD80"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Epcam"), min.cutoff = "q02", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Psmb11"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Pou2f3"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Avil"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO303"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO304"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO305"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO306"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO307"), min.cutoff = "q30", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO308"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO309"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO310"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO311"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("cmo_CMO312"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("mt-Nd6"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("mt-Co1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Vim"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Pdgfra"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Pdgfrb"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Dpp4"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cd3g"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cd4"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cd8a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Acta2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ccn3"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Mfap4"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Rrm2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Tm4sf1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Tubb6"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ccl7"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Agtr2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Dcn"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Enpp2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Lum"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Gdf10"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Itm2a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Mest"), min.cutoff = "q05", max.cutoff = "q95")



FeaturePlot(pbmc, features = c("Krt19"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ccl21a"), min.cutoff = "q05", max.cutoff = "q95")

# Cluster 4 is an mTEPC cluster. When FTM is absent, cells are stuck there and with FTM they move to the right.


library(ggplot2)

ggplot(pbmc@meta.data, aes(x=CMOassignment, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
VlnPlot(object = pbmc, features = 'Ciita',group.by ="CMOassignment",log = F)

# Export to Spring

write.table(as.matrix(GetAssayData(object = pbmc, slot = "counts")), 
            '~/Projects/WatchBreaker/counts.csv', 
            sep = ',', row.names = F, col.names = F, quote = F)

genenamess=row.names(pbmc)
write.table(genenamess,"~/Projects/WatchBreaker/genenames.csv", sep=",",row.names=FALSE,col.names = FALSE,quote=F)
# write.table(c("Day",pbmc@meta.data$Time),"~/Projects/Dedifferentiation/Times.csv", eol = ",",sep=",",row.names=F,col.names = FALSE,quote = F)

groupings=pbmc$CMOassignment

write.table(t(groupings),"~/Projects/WatchBreaker/Groupings.csv", sep=",",row.names=T,col.names = FALSE,quote = F)

# Subset to TEC


selected_c <- WhichCells(pbmc, idents=c(2,3,7), invert = T)

pbmc <- subset(pbmc, cells = selected_c)


# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)



pbmc <- RunPCA(pbmc, verbose = FALSE)

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)


pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.23)

# pbmc1=pbmc
# pbmc=pbmc1

DimPlot(pbmc, label = T) + NoLegend()


FeaturePlot(pbmc, features = c("Aire"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ccl21a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Spink5"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Pou2f3"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Prss16"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Myl9"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("H2-D1"), min.cutoff = "q05", max.cutoff = "q95")


# Compare FTM



markers <- FindMarkers(object = pbmc, group.by="CMOassignment",ident.1 = "CMO309",ident.2 = "CMO311")
markers <- FindMarkers(object = pbmc,ident.1 = 0,ident.2 = 3)
head(x = markers,20)
write.csv(markers,file="~/Projects/WatchBreaker/TECWBFTM+vsFTM-.csv")


library(enrichR)

plt=enrichr(row.names(markers[markers$avg_log2FC>0,]),databases="GO_Biological_Process_2021")
pltdown=enrichr(row.names(markers[markers$avg_log2FC<0,]),databases="GO_Biological_Process_2021")

plt$GO_Biological_Process_2021$Term[1:10]
pltdown$GO_Biological_Process_2021$Term[1:10]

write.csv(plt,file="~/Projects/WatchBreaker/cTECWBFTM+vsFTM-GO.csv")
write.csv(pltdown,file="~/Projects/WatchBreaker/cTECWBFTM-vsFTM+GO.csv")

saveRDS(pbmc,"~/Projects/WatchBreaker/TECSeurat2.rds")

pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")
pbmc=readRDS("~/Projects/WatchBreaker/WB.rds")

DimPlot(pbmc, label = F,group.by = "CMOassignment")
DimPlot(pbmc, label = TRUE,group.by = "seurat_clusters")
FeaturePlot(pbmc, features = c("Vegfa"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Higd1a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Usp18"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt5"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt14"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt15"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Krt8"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Tbata"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Dlk2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Enpep"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Oas1a"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Dll4"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Sphk1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Vegfa"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Zfp36"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cdh1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Itgb2"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ly51"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Foxn1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Prss16"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cxcl12"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("H2-K1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("H2-D1"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("B2m"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("H2-Aa"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Cd74"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Tslp"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Ifnz"), min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(pbmc, features = c("Il6"), min.cutoff = "q05", max.cutoff = "q95")

row.names(pbmc)[grep("Usp18",row.names(pbmc))]

# FTM produce IL6 and Tslp which might be causing interferon signalling in TEC
library(ggplot2)

ggplot(pbmc@meta.data, aes(x=CMOassignment, fill=seurat_clusters)) + geom_bar(position = "fill")+theme_classic()+scale_y_continuous(labels = scales::percent_format())
ggplot(pbmc@meta.data, aes(x=CMOassignment, fill=seurat_clusters)) + geom_bar()+theme_classic()+scale_y_continuous()
VlnPlot(object = pbmc, features = 'H2-K1',group.by ="seurat_clusters")
VlnPlot(object = pbmc, features = 'Avil',group.by ="CMOassignment")

Idents(object = pbmc)=pbmc[["CMOassignment"]] 
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309","CMO311"))
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309"))
selected_c <- WhichCells(pbmc, idents=c("CMO302","CMO307"))
selected_c <- WhichCells(pbmc, idents=c("CMO311","CMO312"))

pbmc <- subset(pbmc, cells = selected_c)
