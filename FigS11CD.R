
pbmc=readRDS("~/Projects/WatchBreaker/TECSeurat2.rds")

Idents(pbmc) <- pbmc$CMOassignment
# Merge 312 and 311
pbmc=RenameIdents(pbmc,c("CMO312"="CMO311"))
selected_c <- WhichCells(pbmc, idents=c("CMO301","CMO302","CMO304","CMO303","CMO309","CMO306","CMO311"))

pbmc <- subset(pbmc, cells = selected_c)
pbmc=RenameIdents(pbmc,c("CMO301"="FTOC","CMO302"="Submerged FTOC","CMO303"="Submerged RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO301"="FTOC","CMO306"="Sorted RFTOC","CMO304"="RFTOC","CMO309"="mTO","CMO311"="mTO without FTM"))
pbmc$condition=Idents(pbmc)

VlnPlot(object = pbmc, features = c('Higd1a',"Hif1a","Krt14","Eno1","Bnip3","Bhlhe40"),group.by ="condition",log = T)



pdf(file = "~/Projects/WatchBreaker/Figures/FigS11CD.pdf",width = 14,height=8)

VlnPlot(object = pbmc,ncol=4, features = c("Foxn1","Psmb11","Dll4","Krt14","Aire","Fezf2","Slpi","Pou2f3"),group.by ="condition",log = T)
VlnPlot(object = pbmc, features = c('Higd1a',"Hif1a","Krt14","Eno1","Bnip3","Bhlhe40"),group.by ="condition",log = T)

dev.off()
