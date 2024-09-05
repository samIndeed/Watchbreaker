pbmc=readRDS("~/Projects/WatchBreaker/WB.rds")


pdf(file = "~/Projects/WatchBreaker/Figures/FigS10B.pdf",width=8,height=4)
FeaturePlot(pbmc,ncol=2, features = c("Foxn1","Epcam"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()


pdf(file = "~/Projects/WatchBreaker/Figures/FigS10CD.pdf",width=16,height=8)
FeaturePlot(pbmc,ncol=4, features = c("Vim","Pdgfra","Pdgfrb","Zeb2","Cd4","Cd8a","Ptprc","Cd3g"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()



Idents(object = pbmc)=pbmc$CMOassignment
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
pbmc=RenameIdents(pbmc,c("CMO307"="Unsorted mTO","CMO301"="FTOC","CMO302"="Submerged FTOC","CMO303"="Submerged RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
pbmc=RenameIdents(pbmc,c("Blank"="Blank/Unassigned/Multiplet","Unassigned"="Blank/Unassigned/Multiplet","Multiplet"="Blank/Unassigned/Multiplet"))
Idents(pbmc) <- factor(x = Idents(pbmc), levels = sort(levels(pbmc)))

DimPlot(pbmc)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue_grey <- function(n,k) {
  hues = seq(15, 375, length = n + 1)
  greyss=rep("grey",n)
  greyss[k]=hcl(h = hues, l = 65, c = 100)[k]
  greyss
}


pdf(file = "~/Projects/WatchBreaker/Figures/FigS10A.pdf",width=6,height=4)
DimPlot(pbmc)
coul=gg_color_hue_grey(9,1)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,2)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,3)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,4)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,5)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,6)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,7)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,8)
DimPlot(pbmc,cols=coul)
coul=gg_color_hue_grey(9,9)
DimPlot(pbmc,cols=coul)
dev.off()

