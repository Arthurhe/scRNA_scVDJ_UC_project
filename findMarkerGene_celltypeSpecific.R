require(Seurat)

tag_celltype=c("macrophage","NKcell","Bcell","Tcell")
base_Rdata_name="180607_RUVscale"

for(i in 1:length(tag_celltype)){
    Exp_Seurat=readRDS(paste0("Exp_Seurat_",tag_celltype[i],".Rds"))
    marker_genes=FindAllMarkers(object = Exp_Seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
    saveRDS(marker_genes, paste0("markerGenes_",base_Rdata_name,"_",tag_celltype[i],".rds"))
}
