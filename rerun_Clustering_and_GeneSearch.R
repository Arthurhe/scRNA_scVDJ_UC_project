suppressMessages(require(data.table))
suppressMessages(require(matrixStats))
suppressMessages(require(Seurat))

load("180521_RUVscale_800.Rdata")
Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:25,resolution = 0.6, print.output = F, save.SNN = TRUE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Clustering done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#find marker gene, Seurat
marker.genes=FindAllMarkers(object = Exp_Seurat, only.pos = TRUE, min.pct = 0.25)

save(Exp_Seurat,marker.genes,file="reClustering_out.Rdata")