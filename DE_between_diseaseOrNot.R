suppressMessages(require(Seurat))
require(xlsx)
source("VDJ_function_pack.R")
load("180607_RUVscale.Rdata")

#cell type assignment
system("rm DE_healthyVSdisease_within_cluster_all.xlsx")
celltype_assignment_df=data.frame(cluster=c("17,16","8,2","13,5","6,0","7,10","11,1","3","18,12,9,19","21,20","4,14,15"),
                                  celltype=c("Naive_B","Memory_B","Plasma","NKT","CD8","CD4","gdT","NK","unknown","M_phi"),
                                  stringsAsFactors=F)

celltype_assignment_expanded=cellTypeDF_processing(celltype_assignment_df)

celltype_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),celltype_assignment_expanded$cluster,celltype_assignment_expanded$celltype)
subcelltype_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),celltype_assignment_expanded$cluster,celltype_assignment_expanded$subcelltype)

DEG=list()
lvl= levels(Exp_Seurat@ident)
for(i in 1:length(lvl)){
    tag_exp=t(2^Exp_Seurat@data[,Exp_Seurat@ident==lvl[i]]-1)
    tag_meta=Exp_Seurat@meta.data$disease_assignment[Exp_Seurat@ident==lvl[i]]
    gp_healthy=tag_exp[tag_meta==1,]
    gp_disease=tag_exp[tag_meta==2,]
    DEG[[i]]=DEG_wilcox_UMI(gp_healthy,gp_disease,p_threshold=0.05,log2fold_threshold=1)
    names(DEG)[i]=celltype_assignment_expanded$subcelltype[celltype_assignment_expanded$cluster==lvl[i]]
    if(sum(DEG[[i]]$DE)==0){
        write.xlsx(data.frame(out="no DE"), file="DE_healthyVSdisease_within_cluster_all.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }else{
        write.xlsx(DEG[[i]][DEG[[i]]$DE,], file="DE_healthyVSdisease_within_cluster_all.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }
}
saveRDS(DEG,file="DE_healthyVSdisease_within_cluster_all.Rds")

###############
Exp_Seurat=readRDS(file="Exp_Seurat_Tcell.Rds")
tag_celltype="Tcell"

#cell type assignment
Tcell_assignment_df=data.frame(cluster=c("12,9,5,15","6,16","11","1,0,10,7","14","4,17,3,8,2,13"),
                                  celltype=c("CD8_TEM_effector","CD8","CD8_TCM","NKT","gdT","CD4"),
                                  stringsAsFactors=F)
Tcell_assignment_expanded=cellTypeDF_processing(Tcell_assignment_df)
Tcell_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),Tcell_assignment_expanded$cluster,Tcell_assignment_expanded$celltype)
Tcellsubcelltype_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),Tcell_assignment_expanded$cluster,Tcell_assignment_expanded$subcelltype)

system("rm DE_healthyVSdisease_within_cluster_Tcell.xlsx")
DEG=list()
lvl= levels(Exp_Seurat@ident)
for(i in 1:length(lvl)){
    tag_exp=t(2^Exp_Seurat@data[,Exp_Seurat@ident==lvl[i]]-1)
    tag_meta=Exp_Seurat@meta.data$disease_assignment[Exp_Seurat@ident==lvl[i]]
    gp_healthy=tag_exp[tag_meta==1,]
    gp_disease=tag_exp[tag_meta==2,]
    DEG[[i]]=DEG_wilcox_UMI(gp_healthy,gp_disease,p_threshold=0.05,log2fold_threshold=1)
    names(DEG)[i]=Tcell_assignment_expanded$subcelltype[Tcell_assignment_expanded$cluster==lvl[i]]
    if(sum(DEG[[i]]$DE)==0){
        write.xlsx(data.frame(out="no DE"), file="DE_healthyVSdisease_within_cluster_Tcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }else{
        write.xlsx(DEG[[i]][DEG[[i]]$DE,], file="DE_healthyVSdisease_within_cluster_Tcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }
}
saveRDS(DEG,file="DE_healthyVSdisease_within_cluster_Tcell.Rds")

###############
Exp_Seurat=readRDS(file="Exp_Seurat_Bcell.Rds")
tag_celltype="Bcell"

#cell type assignment
Bcell_assignment_df=data.frame(cluster=c("6,14,5,2,13,10,12","3,0,1,8","11","4,7,9"),
                                  celltype=c("Plasma","Memory_B","Breg","Naive_B"),
                                  stringsAsFactors=F)
Bcell_assignment_expanded=cellTypeDF_processing(Bcell_assignment_df)
Bcell_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),Bcell_assignment_expanded$cluster,Bcell_assignment_expanded$celltype)
Bcellsubcelltype_assignment=gp_name_replacing(as.character(Exp_Seurat@ident),Bcell_assignment_expanded$cluster,Bcell_assignment_expanded$subcelltype)

system("rm DE_healthyVSdisease_within_cluster_Bcell.xlsx")
DEG=list()
lvl= levels(Exp_Seurat@ident)
for(i in 1:length(lvl)){
    tag_exp=t(2^Exp_Seurat@data[,Exp_Seurat@ident==lvl[i]]-1)
    tag_meta=Exp_Seurat@meta.data$disease_assignment[Exp_Seurat@ident==lvl[i]]
    gp_healthy=tag_exp[tag_meta==1,]
    gp_disease=tag_exp[tag_meta==2,]
    DEG[[i]]=DEG_wilcox_UMI(gp_healthy,gp_disease,p_threshold=0.05,log2fold_threshold=1)
    names(DEG)[i]=Bcell_assignment_expanded$subcelltype[Bcell_assignment_expanded$cluster==lvl[i]]
    if(sum(DEG[[i]]$DE)==0){
        write.xlsx(data.frame(out="no DE"), file="DE_healthyVSdisease_within_cluster_Bcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }else{
        write.xlsx(DEG[[i]][DEG[[i]]$DE,], file="DE_healthyVSdisease_within_cluster_Bcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }
}
saveRDS(DEG,file="DE_healthyVSdisease_within_cluster_Bcell.Rds")

###############
Exp_Seurat=readRDS(file="Exp_Seurat_NKcell.Rds")
tag_celltype="NKcell"

system("rm DE_healthyVSdisease_within_cluster_NKcell.xlsx")
DEG=list()
lvl= levels(Exp_Seurat@ident)
for(i in 1:length(lvl)){
    tag_exp=t(2^Exp_Seurat@data[,Exp_Seurat@ident==lvl[i]]-1)
    tag_meta=Exp_Seurat@meta.data$disease_assignment[Exp_Seurat@ident==lvl[i]]
    gp_healthy=tag_exp[tag_meta==1,]
    gp_disease=tag_exp[tag_meta==2,]
    DEG[[i]]=DEG_wilcox_UMI(gp_healthy,gp_disease,p_threshold=0.05,log2fold_threshold=1)
    names(DEG)[i]=paste0("subcluster_",lvl[i])
    if(sum(DEG[[i]]$DE)==0){
        write.xlsx(data.frame(out="no DE"), file="DE_healthyVSdisease_within_cluster_NKcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }else{
        write.xlsx(DEG[[i]][DEG[[i]]$DE,], file="DE_healthyVSdisease_within_cluster_NKcell.xlsx", sheetName=names(DEG)[i], append=TRUE, row.names=T)
    }
}
saveRDS(DEG,file="DE_healthyVSdisease_within_cluster_NKcell.Rds")