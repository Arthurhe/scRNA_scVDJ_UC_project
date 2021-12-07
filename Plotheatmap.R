suppressMessages(require(ggplot2))
suppressMessages(require(data.table))
suppressMessages(require(gplots))
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
source("VDJ_ploting_suite.R")

Fill_Seurat_DR=function(SeuratOBJ,DRmatrix,DRtype="pca"){
    DR=new("dim.reduction", cell.embeddings = DRmatrix)
    SeuratOBJ@dr[[length(SeuratOBJ@dr)+1]]=DR
    names(SeuratOBJ@dr)[length(SeuratOBJ@dr)]=DRtype
    return(SeuratOBJ)
}

load("180521_RUVscale_800.Rdata")
#load("180521_RUVscale_800_magic.Rdata")
ls()

pdf("exp_heatmap.pdf",width=12,height=12)
tmp=t(Exp_Seurat@scale.data)
tmp[tmp>5]=5
tmp[tmp<(-5)]=-5
col2=rev(RColorBrewer::brewer.pal(10,"RdBu"))
cluster_col=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(levels(Exp_Seurat@ident)))
tagcell=sample(1:ncol(Exp_Seurat@scale.data),4000)
#heatmap.2(tmp[tagcell,HVG],col=colorRampPalette(col2),labRow=NULL,labCol=NULL,
#          lhei=c(1,5),lwid=c(0.1,5),margins = c(0.5,2),trace="none",density.info="none",dendrogram='none',Colv=T,Rowv=T,
#          notecol="black",key=F,main="gene_expression",RowSideColors=cluster_col[Exp_Seurat@ident[tagcell]])

heatmap.2(tmp[order(Exp_Seurat@ident)[tagcell],HVG],col=colorRampPalette(col2),labRow=NULL,labCol=NULL,
          lhei=c(1,5),lwid=c(0.1,5),margins = c(0.5,2),trace="none",density.info="none",dendrogram='none',Colv=T,Rowv=F,
          notecol="black",key=F,main="gene_expression",RowSideColors=cluster_col[sort(Exp_Seurat@ident)[tagcell]])

dev.off()