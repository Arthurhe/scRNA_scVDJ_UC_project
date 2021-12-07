##source activate R35

suppressMessages(require(Rtsne))
suppressMessages(require(data.table))
suppressMessages(require(matrixStats))
suppressMessages(require(Seurat))
require(RUVnormalize)
#require(umap)
require(Matrix)
source("VDJ_function_pack.R")

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

cross_combining=function(string1,string2){
    o=rep("",length(string1)*length(string2))
    k=0
    for(i in 1:length(string1)){
        for(j in 1:length(string2)){
            k=k+1
            o[k]=paste0(string1[i],"_",string2[j])
        } 
    }
    return(o)
}

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

Fill_Seurat_DR=function(SeuratOBJ,DRmatrix,DRtype="pca"){
    DR=new("dim.reduction", cell.embeddings = DRmatrix)
    SeuratOBJ@dr[[length(SeuratOBJ@dr)+1]]=DR
    names(SeuratOBJ@dr)[length(SeuratOBJ@dr)]=DRtype
    return(SeuratOBJ)
}

ptm <- proc.time()

genome="GRCh38"
tag_dir="/home/ahe/Analysis/201801_JohnVDJ/data"
sample_list=c(cross_combining(c("C9"),c("R")),
              cross_combining(c("C12","C16","U4","U5","U34","U35","U41","U44","U45"),c("PBMC","R")),
              cross_combining(c("C17","C18","C19","C21","C30"),c("PBMC","R")), #,"I"
              cross_combining(c("C33"),c("PBMC"))) #,"I""U3","C10","C22","C23",
outprefix="190609_RUVscale" #180919_RUVscale

num_of_BT_cell_at_most_per_lirary=2000
num_of_failed_BT_cell_at_most_per_lirary=1000
num_of_unknown_cell_at_most_per_lirary=1000
nGeneThreshold=400
mitoUpperThreshold=0.1 #smaller than 1
mitoLowerThreshold=0.005
upper_ranking_threshold_pecentage=0.1/100
set.seed(123)

housekeeping_gene=fread("/home/zhh033/genomeFiles/housekeeping_gene_human.txt",data.table=F)
housekeeping_gene=housekeeping_gene[,1]

cell_assignment_gps=c("failed BCR/TCR","B cell","T cell","dual-label","unknown")

t2g <- readRDS("/home/zhh033/genomeFiles/hg38_t2g_R35.Rds")


#loading the real stuff
load(paste0(outprefix,"_tmp.Rdata"))

#seurat construction prt.I
ptm <- proc.time()
Exp_Seurat <- CreateSeuratObject(raw.data = t(ExpressionMat))
Exp_Seurat@data=log2(apply(ExpressionMat,1,function(x){x/sum(x)*1000000})+1)
rm(ExpressionMat)
Exp_Seurat@raw.data=as(Exp_Seurat@raw.data, "sparseMatrix")
Exp_Seurat@data=as(Exp_Seurat@data, "sparseMatrix")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("seurat construction prt1 done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

ptm <- proc.time()
#find top variant gene / use mt gene and house keeping gene as control
taggene_mt=t2g$ext_gene[t2g$chromosome_name=="MT"]
tagmat=t(Exp_Seurat@data)[,rownames(Exp_Seurat@data) %in% c(housekeeping_gene)]
taggene_hk=colnames(tagmat)[order(colSds(as.matrix(tagmat))/colMeans(as.matrix(tagmat)),decreasing = F)[1:40]]
taggene=c(taggene_hk,taggene_mt) #taggene_mt
taggene_posi=which(rownames(Exp_Seurat@data) %in% taggene)
message(paste0(length(taggene_posi)," gene selected"))
    
#TPM + RUV scale
#ExpressionNormed=t(apply(ExpressionFiltered,1,function(x){x/sum(x)*1000000}))
ExpressionNormed=scale(t(Exp_Seurat@raw.data),scale=F)
ExpressionNormed=naiveRandRUV(ExpressionNormed, which(colnames(ExpressionNormed) %in% taggene), nu.coeff=1e-3, k=10)
ExpressionNormed=scale(ExpressionNormed) #[,-which(colnames(ExpressionNormed) %in% taggene)])
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Normed done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))


#seurat construction prt.II
ptm <- proc.time()
patient_type=unique(gsub("_.*","",sample_list))
patient_assignment=rep(0,length(sample_gp))
for(i in 1:length(patient_type)){
    patient_assignment[sample_gp %in% grep(patient_type[i],sample_list)]=i
}

tissue_type=c("PBMC","R") #"I",
tissue_assignment=rep(0,length(sample_gp))
for(i in 1:length(tissue_type)){
    tissue_assignment[sample_gp %in% grep(tissue_type[i],sample_list)]=i
}

diseaseOrNot=c("healthy","disease")
disease_assignment=rep(1,length(patient_assignment))
disease_assignment[patient_assignment %in% grep("U",patient_type)]=2

Exp_Seurat@meta.data$sample_gp <- sample_gp
Exp_Seurat@meta.data$cell_assignment <- cell_assignment
Exp_Seurat@meta.data$patient_assignment <- patient_assignment
Exp_Seurat@meta.data$tissue_assignment <- tissue_assignment
Exp_Seurat@meta.data$disease_assignment <- disease_assignment
Exp_Seurat@scale.data=t(ExpressionNormed)
rm(ExpressionNormed)

Exp_Seurat <- FindVariableGenes(Exp_Seurat, do.plot = F)
HVG <- head(rownames(Exp_Seurat@hvg.info), 5000)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("seurat construction done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))
#save(list=ls(),file="loading_intermediates.Rdata")


ptm <- proc.time()
Exp_Seurat <- RunPCA(object = Exp_Seurat, pc.genes = HVG,pcs.compute = 50,do.print = F)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("PCA done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

ptm <- proc.time()
Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:50,resolution = 2, print.output = F, save.SNN = TRUE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Clustering done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

ptm <- proc.time()
Exp_Seurat <- RunTSNE(object = Exp_Seurat, dims.use = 1:50, do.fast = TRUE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("TSNE done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

ptm <- proc.time()
Exp_Seurat <- umap_seurat(Exp_Seurat,pca_dim=50)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("umap done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

save(sample_list,cell_assignment_gps,patient_type,tissue_type,t2g,taggene_mt,taggene_hk,HVG,file=paste0(outprefix,"_misc.Rdata"))
saveRDS(Exp_Seurat,file=paste0(outprefix,"_Exp_Seurat_R2.Rds"))

ptm <- proc.time()
#contrust matrix, each row is a state
ExpressionMat_by_gp=cbind(Exp_Seurat@ident,2^t(Exp_Seurat@data))
colnames(ExpressionMat_by_gp)[1]="states"
tagV=c("TRAV","TRAC","TRBV","TRBC","TRGV","TRGC","TRDV","TRDC",
       "IGHV","IGHG","IGHA","IGHD","IGHM","IGHE","IGLV","IGLL","IGKV","IGLC","IGKC")
togetrid=c("IGHMBP2")
for(i in 1:length(tagV)){
    tagG=grep(paste0("^",tagV[i]),colnames(ExpressionMat_by_gp))
    tagG=setdiff(tagG,togetrid)
    if(length(tagG)==0){
        next
    }
    tmp=ExpressionMat_by_gp[,tagG]
    if(!is.null(ncol(tmp))){
        tmp=rowSums(tmp)
    }
    ExpressionMat_by_gp=cbind(ExpressionMat_by_gp,tmp)
    colnames(ExpressionMat_by_gp)[ncol(ExpressionMat_by_gp)]=tagV[i]
}
ExpressionMat_by_gp=data.table(as.matrix(ExpressionMat_by_gp))
ExpressionMat_by_gp=ExpressionMat_by_gp[,lapply(.SD, mean), by=states]
ExpressionMat_by_gp=ExpressionMat_by_gp[order(ExpressionMat_by_gp$states)]
ExpressionMat_by_gp[,states:=NULL]
ExpressionMat_by_gp=data.frame(ExpressionMat_by_gp)
colnames(ExpressionMat_by_gp)=gsub("\\.","-",colnames(ExpressionMat_by_gp))
ExpressionMat_by_gp=log2(ExpressionMat_by_gp)
rownames(ExpressionMat_by_gp)=paste0("cluster_",levels(Exp_Seurat@ident))
write.table(t(ExpressionMat_by_gp),"mediam_gene_expression.tsv",sep="\t",col.names=T,row.names=T,quote=F)
saveRDS(ExpressionMat_by_gp,file=paste0(outprefix,"_ExpressionMat_by_gp.Rds"))
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("merging done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

rm(list=setdiff(ls(),c("outprefix","Exp_Seurat","second_to_humanReadableTime")))
