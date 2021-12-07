require(data.table)
require(matrixStats)
require(Seurat)
require(Matrix)
require(Rmagic)

second_to_humanReadableTime=function(t){
  #change second to a vector of hour,min,second
  h=floor(t/3600)
  t=t-h*3600
  m=floor(t/60)
  t=t-m*60
  s=t
  t=c(h,m,s)
  return(t)
}

ptm <- proc.time()
load("190609/190609_RUVscale_misc.Rdata")
Exp_Seurat=readRDS("190609/190609_RUVscale_Exp_Seurat_R2.Rds")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Loading done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

ptm <- proc.time()
ExpressionMat_magic=magic(as.matrix(t(Exp_Seurat@data)))
ExpressionMat_magic=as.matrix(ExpressionMat_magic)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

saveRDS(ExpressionMat_magic,file="190609/190609_logTPM_magic.rds")

##########################
Exp_Seurat=readRDS(file="190609/Exp_Seurat_Tcell_R2.Rds")
tag_celltype="Tcell"
ptm <- proc.time()
#MAGIC

ptm <- proc.time()
ExpressionMat_magic=magic(as.matrix(t(Exp_Seurat@data)))
ExpressionMat_magic=as.matrix(ExpressionMat_magic)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

saveRDS(ExpressionMat_magic,file="190609/Tcell_logTPM_magic.rds")

##########################
Exp_Seurat=readRDS(file="190609/Exp_Seurat_Bcell_R2.Rds")
tag_celltype="Bcell"
ptm <- proc.time()
#MAGIC

ptm <- proc.time()
ExpressionMat_magic=magic(as.matrix(t(Exp_Seurat@data)))
ExpressionMat_magic=as.matrix(ExpressionMat_magic)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

saveRDS(ExpressionMat_magic,file="190609/Bcell_logTPM_magic.rds")

##########################
Exp_Seurat=readRDS(file="190609/Exp_Seurat_NKcell_R2.Rds")
tag_celltype="NKcell"
ptm <- proc.time()
#MAGIC

ptm <- proc.time()
ExpressionMat_magic=magic(as.matrix(t(Exp_Seurat@data)))
ExpressionMat_magic=as.matrix(ExpressionMat_magic)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

saveRDS(ExpressionMat_magic,file="190609/NKcell_logTPM_magic.rds")

##########################
Exp_Seurat=readRDS(file="190609/Exp_Seurat_macrophage_R2.Rds")
tag_celltype="macrophage"
ptm <- proc.time()
#MAGIC

ptm <- proc.time()
ExpressionMat_magic=magic(as.matrix(t(Exp_Seurat@data)))
ExpressionMat_magic=as.matrix(ExpressionMat_magic)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

saveRDS(ExpressionMat_magic,file="190609/macrophage_logTPM_magic.rds")