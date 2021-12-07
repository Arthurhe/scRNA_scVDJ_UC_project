#source active R35
suppressMessages(require(data.table))
suppressMessages(require(Seurat))
suppressMessages(require(matrixStats))
require(RUVnormalize)
#require(umap)
require(Matrix)
source("VDJ_function_pack.R")

housekeeping_gene=fread("/home/zhh033/genomeFiles/housekeeping_gene_human.txt",data.table=F)
housekeeping_gene=housekeeping_gene[,1]
t2g <- readRDS("/home/zhh033/genomeFiles/hg38_t2g_R35.Rds")


Exp_Seurat=readRDS(file="190609/Exp_Seurat_macrophage_R2.Rds")
tag_celltype="macrophage"
ptm <- proc.time()

#find top variant gene / use mt gene and house keeping gene as control
taggene_mt=t2g$ext_gene[t2g$chromosome_name=="MT"]
tagmat=t(Exp_Seurat@data)[,rownames(Exp_Seurat@data) %in% c(housekeeping_gene)]
taggene_hk=colnames(tagmat)[order(colSds(as.matrix(tagmat))/colMeans(as.matrix(tagmat)),decreasing = F)[1:40]]
taggene=c(taggene_hk,taggene_mt) #taggene_mt
taggene_posi=which(rownames(Exp_Seurat@data) %in% taggene)
message(paste0(length(taggene_posi)," gene selected"))

ExpressionNormed=scale(t(Exp_Seurat@raw.data),scale=F)
ExpressionNormed=naiveRandRUV(ExpressionNormed, which(colnames(ExpressionNormed) %in% taggene), nu.coeff=1e-3, k=10)
ExpressionNormed=scale(ExpressionNormed) #[,-which(colnames(ExpressionNormed) %in% taggene)])
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Normed done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

Exp_Seurat@scale.data=t(ExpressionNormed)



Exp_Seurat <- FindVariableGenes(Exp_Seurat, do.plot = F)
HVG <- head(rownames(Exp_Seurat@hvg.info), 5000)
#save(list=ls(),file="loading_intermediates.Rdata")


ptm <- proc.time()
Exp_Seurat <- RunPCA(object = Exp_Seurat, pc.genes = HVG,pcs.compute = 50,do.print = F)

ptm <- proc.time()
Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:20,resolution = 2, print.output = F, save.SNN = TRUE)

ptm <- proc.time()
Exp_Seurat <- RunTSNE(object = Exp_Seurat, dims.use = 1:50, do.fast = TRUE)

ptm <- proc.time()
Exp_Seurat <- umap_seurat(Exp_Seurat,pca_dim=20)

saveRDS(Exp_Seurat,file="190609/Exp_Seurat_macrophage_normed.Rds")