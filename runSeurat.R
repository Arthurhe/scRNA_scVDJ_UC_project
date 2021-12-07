require("data.table")
require("Seurat")
devtools::load_all("/home/ahe/tools/Lightbulb")

load("ExpMatRaw_VDJ.Rdata")
batch=sample_gp
for(i in 1:max(batch)){
    tag=which(batch==i)
    rownames(ExpressionMat)[tag]=gsub("1$", i, rownames(ExpressionMat)[tag])
}

ptm <- proc.time()

Exp_Seurat <- CreateSeuratObject(raw.data = t(ExpressionMat), project = "VDJ", min.cells = 5)
Exp_Seurat@meta.data$batch <- batch
Exp_Seurat@meta.data$cell_assignment <- cell_assignment
Exp_Seurat <- FilterCells(Exp_Seurat, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
Exp_Seurat <- NormalizeData(Exp_Seurat)
Exp_Seurat <- ScaleData(Exp_Seurat, display.progress = F)
Exp_Seurat <- FindVariableGenes(Exp_Seurat, do.plot = F)
genes.use <- head(rownames(Exp_Seurat@hvg.info), 1000)

SeuratList=list()
for(i in 1:max(batch)){
    tag=which(Exp_Seurat@meta.data$batch==i)
    SeuratList[[i]]=SubsetData(Exp_Seurat,cells.use =Exp_Seurat@cell.names[tag])
    SeuratList[[i]]@meta.data$group <- i
}
names(SeuratList)=paste0("batch",1:max(batch))

message("preprocess done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

Exp_Seurat_CCA <- RunMultiCCA(SeuratList, genes.use = genes.use, num.cc = 20)

message("RunMultiCCA done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

Exp_Seurat_CCA <- AlignSubspace(Exp_Seurat_CCA, reduction.type = "cca", grouping.var = "batch",dims.align = 1:20,verbose=F)

message("AlignSubspace done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

# t-SNE
Exp_Seurat_CCA <- RunTSNE(Exp_Seurat_CCA, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)

message("t-SNE done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

# Clustering
Exp_Seurat_CCA <- FindClusters(Exp_Seurat_CCA, reduction.type = "cca.aligned",resolution = 0.6, dims.use = 1:20,print.output = F)

message("Clustering done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

# Identify conserved cell type markers
#nk.markers <- FindConservedMarkers(Exp_Seurat_CCA, ident.1 = 7, grouping.var = "batch", print.bar = FALSE)

#message("nk.markers done")
#t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
#message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

save(Exp_Seurat,genes.use,Exp_Seurat_CCA,file="Seurat_out_VDJ.Rdata")

message("save done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))