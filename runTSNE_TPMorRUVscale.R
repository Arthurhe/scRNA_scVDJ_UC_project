suppressMessages(require(Rtsne))
suppressMessages(require(cellrangerRkit))
suppressMessages(require(ggplot2))
suppressMessages(require(data.table))
suppressMessages(require(gplots))
suppressMessages(require(matrixStats))
suppressMessages(require(RColorBrewer))
require(RUVnormalize)
venn_matrix=function(inlist){
    #mk matrix
    out_matrix=matrix(NA,length(inlist),length(inlist))
    colnames(out_matrix)=names(inlist)
    rownames(out_matrix)=names(inlist)
    for(i in 2:length(inlist)){
        for(j in 1:(i-1)){
            out_matrix[i,j]=length(intersect(inlist[[i]],inlist[[j]]))
        }
    }
    diag(out_matrix)=sapply(inlist,length)
    return(out_matrix)
}

center_by_batch=function(inMatrix){
    inMatrix=scale(inMatrix,scale = F)
    batch_ids=as.numeric(stringi::stri_sub(rownames(inMatrix),-1, -1))
    unique_batch_ids=unique(batch_ids)
    for(i in 1:length(unique_batch_ids)){
        tagrow=batch_ids == unique_batch_ids[i]
        inMatrix[tagrow,]=scale(inMatrix[tagrow,],scale = F)
    }
    return(inMatrix)
}

genome="GRCh38"
sample_list=c("C12pBMC","C12Rectum","C9pBMC","C9Rectum") #
cell_assignment_gps=c("failed BCR/TCR","B cell","T cell","dual-label","unknown")
SC_matrix=list()
ercc_matrix=list()
cell_assignment=list()
sample_gp=c()
BCR_annoation_list=list()
TCR_annoation_list=list()
BCR_annoation_list_weak=list()
TCR_annoation_list_weak=list()

for(i in 1:length(sample_list)){
    #sc
    valid_sc_barcode=load_cellranger_matrix_h5(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/",sample_list[i],"SC"), genome=genome,barcode_filtered =T)
    valid_sc_barcode=valid_sc_barcode$barcode
    SC_matrix[[i]]=load_cellranger_matrix_h5(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/",sample_list[i],"SC"), genome=genome,barcode_filtered =F)
    SC_matrix[[i]]=exprs(SC_matrix[[i]])
    #vdj
    BCR_annoation_list[[i]]=fread(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/",sample_list[i],"BCR/outs/all_contig_annotations.csv"),sep=",",data.table=F,na.strings = "None")
    TCR_annoation_list[[i]]=fread(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/",sample_list[i],"TCR/outs/all_contig_annotations.csv"),sep=",",data.table=F,na.strings = "None")   

    #continue BCR TCR processing
    valid_BCR_list=which(BCR_annoation_list[[i]]$is_cell & BCR_annoation_list[[i]]$high_confidence & BCR_annoation_list[[i]]$productive)
    valid_TCR_list=which(TCR_annoation_list[[i]]$is_cell & TCR_annoation_list[[i]]$high_confidence & TCR_annoation_list[[i]]$productive)
    valid_B_list=unique(BCR_annoation_list[[i]]$barcode[valid_BCR_list])
    valid_T_list=unique(TCR_annoation_list[[i]]$barcode[valid_TCR_list])
    tag_B=BCR_annoation_list[[i]]$barcode %in% valid_B_list
    tag_T=TCR_annoation_list[[i]]$barcode %in% valid_T_list
    BCR_annoation_list_weak[[i]]=BCR_annoation_list[[i]][-tag_B,]
    TCR_annoation_list_weak[[i]]=TCR_annoation_list[[i]][-tag_T,]
    BCR_annoation_list_weak[[i]]=BCR_annoation_list_weak[[i]][BCR_annoation_list_weak[[i]]$is_cell | BCR_annoation_list_weak[[i]]$high_confidence,]
    TCR_annoation_list_weak[[i]]=TCR_annoation_list_weak[[i]][TCR_annoation_list_weak[[i]]$is_cell | TCR_annoation_list_weak[[i]]$high_confidence,]
    BCR_annoation_list[[i]]=BCR_annoation_list[[i]][tag_B,]
    TCR_annoation_list[[i]]=TCR_annoation_list[[i]][tag_T,]

    #the cell selection
    valid_cell_barcodes=unique(c(valid_sc_barcode,BCR_annoation_list[[i]]$barcode,TCR_annoation_list[[i]]$barcode))
    SC_matrix[[i]]=SC_matrix[[i]][,colnames(SC_matrix[[i]]) %in% valid_cell_barcodes]
    BCR_annoation_list[[i]]=BCR_annoation_list[[i]][BCR_annoation_list[[i]]$barcode %in% valid_cell_barcodes,]
    TCR_annoation_list[[i]]=TCR_annoation_list[[i]][TCR_annoation_list[[i]]$barcode %in% valid_cell_barcodes,]
    #continue SC MATRIX PROCESSION
    SC_matrix[[i]]=t(as.matrix(SC_matrix[[i]]))
    SC_matrix[[i]]=SC_matrix[[i]][rowSums(SC_matrix[[i]])>100,]
    #ercc
    ercc_matrix[[i]]=load_cellranger_matrix_h5(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/",sample_list[i],"SC_ercc"), genome="ercc92",barcode_filtered =F)
    ercc_matrix[[i]]=exprs(ercc_matrix[[i]])
    ercc_matrix[[i]]=ercc_matrix[[i]][,match(rownames(SC_matrix[[i]]),colnames(ercc_matrix[[i]]))]
    ercc_matrix[[i]]=t(as.matrix(ercc_matrix[[i]]))
    SC_matrix[[i]]=cbind(ercc_matrix[[i]],SC_matrix[[i]])
    sample_gp=c(sample_gp,rep(i,nrow(SC_matrix[[i]])))
    #renaming
    rownames(SC_matrix[[i]])=paste0(gsub("1$", "", rownames(SC_matrix[[i]])),i)
    TCR_annoation_list[[i]]$barcode=paste0(gsub("1$", "", TCR_annoation_list[[i]]$barcode),i)
    BCR_annoation_list[[i]]$barcode=paste0(gsub("1$", "", BCR_annoation_list[[i]]$barcode),i)
    #cell labeling
    cell_assignment[[i]]=rep("unknown",nrow(SC_matrix[[i]]))
    cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% BCR_annoation_list_weak[[i]]$barcode |
                         rownames(SC_matrix[[i]]) %in% TCR_annoation_list_weak[[i]]$barcode]="failed BCR/TCR"
    cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% BCR_annoation_list[[i]]$barcode]="B cell"
    cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% TCR_annoation_list[[i]]$barcode]="T cell"
    cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% TCR_annoation_list[[i]]$barcode & 
                         rownames(SC_matrix[[i]]) %in% BCR_annoation_list[[i]]$barcode]="dual-label"
    
    
}

ExpressionMat=do.call(rbind,SC_matrix)
ExpressionMat=ExpressionMat[,colSums(ExpressionMat)>2]
ExpressionBinaryMat=ExpressionMat>=1
ExpressionMat=ExpressionMat[,colSums(ExpressionBinaryMat)>nrow(ExpressionMat)/100]
#ExpressionMat=t(apply(ExpressionMat,1,function(x){x/sum(x)*10000000}))

cell_assignment=unlist(cell_assignment)
BCR_annoation=do.call(rbind,BCR_annoation_list)
TCR_annoation=do.call(rbind,TCR_annoation_list)
BCR_annoation_w=do.call(rbind,BCR_annoation_list_weak)
TCR_annoation_w=do.call(rbind,TCR_annoation_list_weak)
rm(SC_matrix,BCR_annoation_list,TCR_annoation_list)

#find top variant gene
sdY <- apply(ExpressionMat[,-grep("ERCC",colnames(ExpressionMat))], 2, sd)
ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
ExpressionFiltered=ExpressionMat[,c(grep("ERCC",colnames(ExpressionMat)),ssd[1:3000]+max(grep("ERCC",colnames(ExpressionMat))))]

#TPMscale
ExpressionNormed=t(apply(ExpressionFiltered,1,function(x){x/sum(x)*10000000}))
ExpressionNormed=scale(ExpressionNormed[,-grep("ERCC",colnames(ExpressionMat))])

rtsne_allcell=Rtsne(ExpressionNormed,dims=2,max_iter = 1000)
pca_allcell=prcomp(ExpressionNormed)
#save(pca_allcell,rtsne_allcell,file="TSNE_PCA_normed_out_perbatchscale_RUV2shrink.Rdata")
save(ExpressionNormed,ExpressionMat,pca_allcell,rtsne_allcell,BCR_annoation,cell_assignment,TCR_annoation,BCR_annoation_w,TCR_annoation_w,file="testing_TPMscale.Rdata")

require(repr)
pdf("TSNE_PCA scRNAseq_detectedTCRBCR_TPMscale.pdf",width=10,height=10)
par(mfcol=c(2,2))

coltouse=cell_assignment
tagcol=RColorBrewer::brewer.pal(length(cell_assignment_gps),"Set1")
for (i in 1:length(cell_assignment_gps)){
    coltouse[coltouse==cell_assignment_gps[i]]=tagcol[i]
}
plot(rtsne_allcell$Y,pch=19,col=coltouse,cex=0.5,main="TSNE B/T cell label",xlim=c(min(rtsne_allcell$Y[,1]),max(rtsne_allcell$Y[,1])+20))
legend("topright",cell_assignment_gps,bty = "n",lty=0,pch=19,col=tagcol)
plot(pca_allcell$x[,1:2],pch=19,col=coltouse,cex=0.5,main="PCA B/T cell label",xlim=c(min(pca_allcell$x[,1]),max(pca_allcell$x[,1])+20))
legend("topright",cell_assignment_gps,bty = "n",lty=0,pch=19,col=tagcol)

tagcol=RColorBrewer::brewer.pal(length(sample_list),"Set1")
plot(rtsne_allcell$Y,pch=19,col=tagcol[sample_gp],cex=0.5,main="TSNE PBMC/Rectum label",xlim=c(min(rtsne_allcell$Y[,1]),max(rtsne_allcell$Y[,1])+10))
legend("topright",sample_list,bty = "n",lty=0,pch=19,col=tagcol)
plot(pca_allcell$x[,1:2],pch=19,col=tagcol[sample_gp],cex=0.5,main="PCA PBMC/Rectum label",xlim=c(min(pca_allcell$x[,1]),max(pca_allcell$x[,1])+10))
legend("topright",sample_list,bty = "n",lty=0,pch=19,col=tagcol)
dev.off()


###########################################
#RUVscale
ExpressionNormed=scale(ExpressionFiltered,scale=F)
ExpressionNormed=naiveRandRUV(ExpressionNormed, grep("ERCC",colnames(ExpressionNormed)), nu.coeff=1e-3, k=20)
ExpressionNormed=scale(ExpressionNormed[,-grep("ERCC",colnames(ExpressionNormed))])

rtsne_allcell=Rtsne(ExpressionNormed,dims=2,max_iter = 1000)
pca_allcell=prcomp(ExpressionNormed)
#save(pca_allcell,rtsne_allcell,file="TSNE_PCA_normed_out_perbatchscale_RUV2shrink.Rdata")
save(ExpressionNormed,ExpressionMat,pca_allcell,rtsne_allcell,BCR_annoation,cell_assignment,TCR_annoation,BCR_annoation_w,TCR_annoation_w,file="testing_RUVscale.Rdata")

require(repr)
pdf("TSNE_PCA scRNAseq_detectedTCRBCR_RUVscale.pdf",width=10,height=10)
par(mfcol=c(2,2))

coltouse=cell_assignment
tagcol=RColorBrewer::brewer.pal(length(cell_assignment_gps),"Set1")
for (i in 1:length(cell_assignment_gps)){
    coltouse[coltouse==cell_assignment_gps[i]]=tagcol[i]
}
plot(rtsne_allcell$Y,pch=19,col=coltouse,cex=0.5,main="TSNE B/T cell label",xlim=c(min(rtsne_allcell$Y[,1]),max(rtsne_allcell$Y[,1])+20))
legend("topright",cell_assignment_gps,bty = "n",lty=0,pch=19,col=tagcol)
plot(pca_allcell$x[,1:2],pch=19,col=coltouse,cex=0.5,main="PCA B/T cell label",xlim=c(min(pca_allcell$x[,1]),max(pca_allcell$x[,1])+20))
legend("topright",cell_assignment_gps,bty = "n",lty=0,pch=19,col=tagcol)

tagcol=RColorBrewer::brewer.pal(length(sample_list),"Set1")
plot(rtsne_allcell$Y,pch=19,col=tagcol[sample_gp],cex=0.5,main="TSNE PBMC/Rectum label",xlim=c(min(rtsne_allcell$Y[,1]),max(rtsne_allcell$Y[,1])+10))
legend("topright",sample_list,bty = "n",lty=0,pch=19,col=tagcol)
plot(pca_allcell$x[,1:2],pch=19,col=tagcol[sample_gp],cex=0.5,main="PCA PBMC/Rectum label",xlim=c(min(pca_allcell$x[,1]),max(pca_allcell$x[,1])+10))
legend("topright",sample_list,bty = "n",lty=0,pch=19,col=tagcol)
dev.off()
