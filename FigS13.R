pbmc=readRDS("~/Projects/WatchBreaker/WB.rds")





Idents(object = pbmc)=pbmc$CMOassignment
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO306","CMO304","CMO309","CMO311"))
pbmc <- subset(pbmc, cells = selected_c)

pbmc=RenameIdents(pbmc,c("CMO307"="Unsorted mTO","CMO301"="FTOC","CMO302"="Submerged FTOC","CMO303"="Submerged RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))



pbmc$condition=Idents(pbmc)


ifns=rownames(pbmc)[grep("Ifn",rownames(pbmc))]

# Show no IFns
pdf(file = "~/Projects/WatchBreaker/Figures/Fig513B.pdf",width = 18,height=10)
VlnPlot(object = pbmc,ncol=6, features =ifns,group.by ="condition",log = T)
dev.off()

pdf(file = "~/Projects/WatchBreaker/Figures/Fig513C.pdf",width = 6,height=4)
FeaturePlot(pbmc, features = c("Ifnz"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()

