pbmc=readRDS("~/Projects/WatchBreaker/Figures/TECFig6.rds")


pdf(file = "~/Projects/WatchBreaker/Figures/FigS11.pdf",width=15,height=16)
FeaturePlot(pbmc,ncol=5, features = c("Foxn1","Psmb11","Dll4","Isg15",'H2-K1',"H2-Aa","Mki67","Krt14","Krt19","Aire","Fezf2","Slpi","Pou2f3","Krt5","Krt8","Cldn3","Cldn4","Plet1","Sox9","Tnfrsf11a","Prss16","Ccl21a","Ccl25","Cxcl12","Ly75"), min.cutoff = "q05", max.cutoff = "q95")
dev.off()
