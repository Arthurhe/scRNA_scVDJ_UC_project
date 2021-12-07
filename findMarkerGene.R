require(Seurat)

load("180607_RUVscale.Rdata")
marker_genes=FindAllMarkers(object = Exp_Seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
saveRDS(marker_genes, "markerGenes.rds")