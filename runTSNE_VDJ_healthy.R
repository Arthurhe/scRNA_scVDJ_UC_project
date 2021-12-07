suppressMessages(require(Rtsne))
suppressMessages(require(cellrangerRkit))
suppressMessages(require(data.table))
suppressMessages(require(matrixStats))
suppressMessages(require(Seurat))
require(RUVnormalize)
require(Matrix)

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
sample_list=c(cross_combining(c("C9","C12","C16"),c("pBMC","R")),
              cross_combining(c("C10","C17","C18","C19","C21","C22","C23","C24","C30","C33"),c("pBMC","R","I"))) #
#sample_list=cross_combining(c("C9"),c("pBMC","R"))
num_of_unknown_cell_at_most_per_lirary=800
nGeneThreshold=300
mitoUpperThreshold=0.2
upper_ranking_threshold_pecentage=0.1/100

#check all sample presenting
existing10x=list.dirs(paste0(tag_dir,"/10x"),recursive = F,full.names = F)
existing10x=gsub("_10x","",existing10x)
if(any(!sample_list %in% existing10x)){
    stop(paste0(paste0(sample_list[!sample_list %in% existing10x],collapse = ", ")," don't exist"))
}

housekeeping_gene=fread("~/Analysis/genomeFiles/housekeeping_gene_human.txt",data.table=F)
housekeeping_gene=housekeeping_gene[,1]

cell_assignment_gps=c("failed BCR/TCR","B cell","T cell","dual-label","unknown")

t2g <- readRDS("/home/ahe/Analysis/genomeFiles/hg38_t2g_R35.Rds")

MTlist=t2g$ens_gene[t2g$chromosome_name=="MT"]

SC_matrix=list()
ercc_matrix=list()
cell_assignment=list()
sample_gp=c()
BCR_annoation_list=list()
TCR_annoation_list=list()
BCR_annoation_list_weak=list()
TCR_annoation_list_weak=list()

for(i in 1:length(sample_list)){
    #check file exist (file.exists(destfile))
    SC2read=paste0(tag_dir,"/10x/",sample_list[i],"_10x")
    BCR2read=paste0(tag_dir,"/VDJ/",sample_list[i],"_BCR/outs/all_contig_annotations.csv")
    TCR2read=paste0(tag_dir,"/VDJ/",sample_list[i],"_TCR/outs/all_contig_annotations.csv")
    #read sc
    if(file.exists(BCR2read) | file.exists(TCR2read)){
        valid_cell_barcodes=load_cellranger_matrix_h5(SC2read, genome=genome,barcode_filtered =T)
        valid_cell_barcodes=valid_cell_barcodes$barcode
        SC_matrix[[i]]=load_cellranger_matrix_h5(SC2read, genome=genome,barcode_filtered =F)
        SC_matrix[[i]]=exprs(SC_matrix[[i]])
    }else{
        SC_matrix[[i]]=load_cellranger_matrix_h5(SC2read, genome=genome,barcode_filtered =T)
        valid_cell_barcodes=SC_matrix[[i]]$barcode
        SC_matrix[[i]]=exprs(SC_matrix[[i]])
    }
    
    #if BCR exist, read it
    if(file.exists(BCR2read)){
        BCR_annoation_list[[i]]=fread(BCR2read,sep=",",data.table=F,na.strings = "None")
        #filter the error signal
        BCR_annoation_list[[i]]=BCR_annoation_list[[i]][-grep("TR",BCR_annoation_list[[i]]$chain),]
        #continue BCR TCR processing
        valid_BCR_list=which(BCR_annoation_list[[i]]$is_cell & BCR_annoation_list[[i]]$high_confidence & BCR_annoation_list[[i]]$productive)
        valid_B_list=unique(BCR_annoation_list[[i]]$barcode[valid_BCR_list])
        tag_B=BCR_annoation_list[[i]]$barcode %in% valid_B_list
        BCR_annoation_list_weak[[i]]=BCR_annoation_list[[i]][!tag_B,]
        BCR_annoation_list_weak[[i]]=BCR_annoation_list_weak[[i]][BCR_annoation_list_weak[[i]]$is_cell & BCR_annoation_list_weak[[i]]$high_confidence,]
        BCR_annoation_list[[i]]=BCR_annoation_list[[i]][tag_B,]
        valid_cell_barcodes=unique(c(valid_cell_barcodes,valid_B_list))
    }
    #if TCR exist, read it
    if(file.exists(TCR2read)){
        TCR_annoation_list[[i]]=fread(TCR2read,sep=",",data.table=F,na.strings = "None")   
        #filter the error signal
        TCR_annoation_list[[i]]=TCR_annoation_list[[i]][-grep("IG",TCR_annoation_list[[i]]$chain),]
        #continue BCR TCR processing 
        valid_TCR_list=which(TCR_annoation_list[[i]]$is_cell & TCR_annoation_list[[i]]$high_confidence & TCR_annoation_list[[i]]$productive) 
        valid_T_list=unique(TCR_annoation_list[[i]]$barcode[valid_TCR_list])
        tag_T=TCR_annoation_list[[i]]$barcode %in% valid_T_list
        TCR_annoation_list_weak[[i]]=TCR_annoation_list[[i]][!tag_T,]
        TCR_annoation_list_weak[[i]]=TCR_annoation_list_weak[[i]][TCR_annoation_list_weak[[i]]$is_cell & TCR_annoation_list_weak[[i]]$high_confidence,]
        TCR_annoation_list[[i]]=TCR_annoation_list[[i]][tag_T,]
        valid_cell_barcodes=unique(c(valid_cell_barcodes,valid_T_list))
    }
    
    #the cell selection
    if(file.exists(BCR2read) | file.exists(TCR2read)){
        SC_matrix[[i]]=SC_matrix[[i]][,colnames(SC_matrix[[i]]) %in% valid_cell_barcodes]
    }
    if(file.exists(BCR2read)){
        BCR_annoation_list[[i]]=BCR_annoation_list[[i]][BCR_annoation_list[[i]]$barcode %in% valid_cell_barcodes,]
        BCR_annoation_list[[i]]$barcode=gsub("1$", i, BCR_annoation_list[[i]]$barcode)
        BCR_annoation_list_weak[[i]]$barcode=gsub("1$", i, BCR_annoation_list_weak[[i]]$barcode)  
    }
    if(file.exists(TCR2read)){
        TCR_annoation_list[[i]]=TCR_annoation_list[[i]][TCR_annoation_list[[i]]$barcode %in% valid_cell_barcodes,]
        TCR_annoation_list[[i]]$barcode=gsub("1$", i, TCR_annoation_list[[i]]$barcode)
        TCR_annoation_list_weak[[i]]$barcode=gsub("1$", i, TCR_annoation_list_weak[[i]]$barcode)
    }
    
    #continue SC MATRIX PROCESSION
    SC_matrix[[i]]=t(as.matrix(SC_matrix[[i]]))
    
    #filtering
    enoughGene=rowSums(SC_matrix[[i]]>0)>=nGeneThreshold
    SC_matrix[[i]]=SC_matrix[[i]][enoughGene,]
    MTsum=rowSums(SC_matrix[[i]][,colnames(SC_matrix[[i]]) %in% MTlist])
    toomuchmito=MTsum/rowSums(SC_matrix[[i]])>mitoUpperThreshold
    nomito=MTsum==0
    SC_matrix[[i]]=SC_matrix[[i]][(!toomuchmito) & (!nomito),]
    
    #renaming
    rownames(SC_matrix[[i]])=gsub("1$", i, rownames(SC_matrix[[i]]))
    
    #cell labeling
    cell_assignment[[i]]=rep("unknown",nrow(SC_matrix[[i]]))    
    if(file.exists(BCR2read) & file.exists(TCR2read)){
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% BCR_annoation_list_weak[[i]]$barcode |
                             rownames(SC_matrix[[i]]) %in% TCR_annoation_list_weak[[i]]$barcode]="failed BCR/TCR"
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% TCR_annoation_list[[i]]$barcode & 
                             rownames(SC_matrix[[i]]) %in% BCR_annoation_list[[i]]$barcode]="dual-label"
    }
    if(file.exists(BCR2read)){
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% BCR_annoation_list_weak[[i]]$barcode]="failed BCR/TCR"
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% BCR_annoation_list[[i]]$barcode]="B cell"
    }
    if(file.exists(TCR2read)){
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% TCR_annoation_list_weak[[i]]$barcode]="failed BCR/TCR"
        cell_assignment[[i]][rownames(SC_matrix[[i]]) %in% TCR_annoation_list[[i]]$barcode]="T cell"
    }
    
    
    #remove the redundant unknown cells
    if(sum(cell_assignment[[i]]=="unknown")>num_of_unknown_cell_at_most_per_lirary){
        tag_unknown=which(cell_assignment[[i]]=="unknown")
        discard_unknown=sample(tag_unknown,length(tag_unknown)-num_of_unknown_cell_at_most_per_lirary,replace = F)
        SC_matrix[[i]]=SC_matrix[[i]][-discard_unknown,]
        cell_assignment[[i]]=cell_assignment[[i]][-discard_unknown]
    }    
    
    sample_gp=c(sample_gp,rep(i,nrow(SC_matrix[[i]])))
    #ercc
    #ercc_matrix[[i]]=load_cellranger_matrix_h5(paste0("/home/ahe/Analysis/201801_JohnVDJ/data/ercc/",sample_list[i],"_10x_ercc"), genome="ercc92",barcode_filtered =F)
    #ercc_matrix[[i]]=exprs(ercc_matrix[[i]])
    #ercc_matrix[[i]]=ercc_matrix[[i]][,match(rownames(SC_matrix[[i]]),colnames(ercc_matrix[[i]]))]
    #ercc_matrix[[i]]=t(as.matrix(ercc_matrix[[i]]))
    #SC_matrix[[i]]=cbind(ercc_matrix[[i]],SC_matrix[[i]])
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
    message(paste0(i,"/",length(sample_list)," loading done: ",t[1],"h-",t[2],"m-",round(t[3]),"s"))
    message(paste0("droped by too few gene:",sum(!enoughGene),"| mito:",sum(toomuchmito),"|",sum(nomito)))
}

ExpressionMat=do.call(rbind,SC_matrix)
ExpressionMat=ExpressionMat[,colSums(ExpressionMat)>2]
ExpressionBinaryMat=ExpressionMat>=1
ExpressionMat=ExpressionMat[,colSums(ExpressionBinaryMat)>nrow(ExpressionMat)/100]
#ExpressionMat=t(apply(ExpressionMat,1,function(x){x/sum(x)*10000000}))

#filter the outlier for each gene:
upper_ranking_threshold=ceiling(upper_ranking_threshold_pecentage*nrow(ExpressionMat))
ExpressionMat=apply(ExpressionMat,2,function(x){
    ordered_x=sort(x,decreasing=T)
    max_threshold=ordered_x[upper_ranking_threshold]
    x[x>max_threshold]=max_threshold
    return(x)
})

#change gene name
ExpressionMat=ExpressionMat[,which(colnames(ExpressionMat) %in% t2g$ens_gene)]
colnames(ExpressionMat)=t2g$ext_gene[match(colnames(ExpressionMat),t2g$ens_gene)]
colnames(ExpressionMat)[which(duplicated(colnames(ExpressionMat)))]=paste0(colnames(ExpressionMat)[which(duplicated(colnames(ExpressionMat)))],"_1")

cell_assignment=unlist(cell_assignment)
BCR_annoation=do.call(rbind,BCR_annoation_list)
TCR_annoation=do.call(rbind,TCR_annoation_list)
BCR_annoation_w=do.call(rbind,BCR_annoation_list_weak)
TCR_annoation_w=do.call(rbind,TCR_annoation_list_weak)
rm(SC_matrix,BCR_annoation_list,TCR_annoation_list,TCR_annoation_list_weak,BCR_annoation_list_weak)

#find top variant gene / use mt gene and house keeping gene as control
taggene_mt=t2g$ext_gene[t2g$chromosome_name=="MT"]
tagmat=t(apply(ExpressionMat,1,function(x){x/sum(x)*1000000}))[,colnames(ExpressionMat) %in% housekeeping_gene]
taggene_hk=colnames(tagmat)[order(colSds(tagmat),decreasing = T)[1:30]]
taggene=c(taggene_mt,housekeeping_gene) #
taggene_posi=which(colnames(ExpressionMat) %in% taggene)

#TPM + RUV scale
#ExpressionNormed=t(apply(ExpressionFiltered,1,function(x){x/sum(x)*1000000}))
ExpressionNormed=scale(ExpressionMat,scale=F)
ExpressionNormed=naiveRandRUV(ExpressionNormed, which(colnames(ExpressionNormed) %in% taggene), nu.coeff=1e-3, k=10)
ExpressionNormed=scale(ExpressionNormed) #[,-which(colnames(ExpressionNormed) %in% taggene)])

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Normed done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

patient_type=unique(gsub("_.*","",sample_list))
patient_assignment=rep(0,length(sample_gp))
for(i in 1:length(patient_type)){
    patient_assignment[sample_gp %in% grep(patient_type[i],sample_list)]=i
}

tissue_type=c("pBMC","I","R")
tissue_assignment=rep(0,length(sample_gp))
for(i in 1:length(tissue_type)){
    tissue_assignment[sample_gp %in% grep(tissue_type[i],sample_list)]=i
}

diseaseOrNot=c("healthy","disease")
disease_assignment=rep(1,length(patient_assignment))
disease_assignment[patient_assignment %in% grep("U",patient_type)]=2

Exp_Seurat <- CreateSeuratObject(raw.data = t(ExpressionMat))
Exp_Seurat@raw.data=as(t(ExpressionMat), "sparseMatrix")
Exp_Seurat@meta.data$sample_gp <- sample_gp
Exp_Seurat@meta.data$cell_assignment <- cell_assignment
Exp_Seurat@meta.data$patient_assignment <- patient_assignment
Exp_Seurat@meta.data$tissue_assignment <- tissue_assignment
Exp_Seurat@meta.data$disease_assignment <- disease_assignment
Exp_Seurat@data=log2(apply(ExpressionMat,1,function(x){x/sum(x)*1000000})+1)
Exp_Seurat@scale.data=t(ExpressionNormed)
Exp_Seurat <- FindVariableGenes(Exp_Seurat, do.plot = F)
HVG <- head(rownames(Exp_Seurat@hvg.info), 5000)

#save(list=ls(),file="loading_intermediates.Rdata")
rm(ExpressionMat,ExpressionNormed)

ptm <- proc.time()
Exp_Seurat <- RunPCA(object = Exp_Seurat, pc.genes = HVG[1:3000],pcs.compute = 50,do.print = F)
#pca_allcell=prcomp(ExpressionNormed[,HVG[1:3000]])
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("PCA done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))


Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:25,resolution = 0.6, print.output = F, save.SNN = TRUE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("Clustering done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#Exp_Seurat <- RunTSNE(object = Exp_Seurat, reduction.use = "pca", dims.use = 1:25, tsne.method = "FIt-SNE", nthreads = 4,reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "/home/ahe/tools/FIt-SNE/bin/fast_tsne",max_iter = 2000)
Exp_Seurat <- RunTSNE(object = Exp_Seurat, dims.use = 1:25, do.fast = TRUE)
#rtsne_allcell=Rtsne(EExp_Seurat@dr$pca@cell.embeddings[,1:min(25,ncol(pca_allcell$x))],dims=2,max_iter = 1500, pca=F)
#Exp_Seurat=Fill_Seurat_DR(Exp_Seurat,rtsne_allcell$Y[1:1000,],DRtype = "tsne")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("TSNE done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#Exp_Seurat <- RunUMAP(object = Exp_Seurat, reduction.use = "pca", dims.use = 1:25, min_dist = 0.75)
#t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
#message(paste("UMAP done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#contrust matrix, each row is a state
ExpressionMat_by_gp=cbind(Exp_Seurat@ident,t(Exp_Seurat@data))
colnames(ExpressionMat_by_gp)[1]="states"
ExpressionMat_by_gp=data.table(ExpressionMat_by_gp)
ExpressionMat_by_gp=ExpressionMat_by_gp[,lapply(.SD, median), by=states]
ExpressionMat_by_gp=ExpressionMat_by_gp[order(ExpressionMat_by_gp$states)]
ExpressionMat_by_gp[,states:=NULL]
ExpressionMat_by_gp=data.frame(ExpressionMat_by_gp)
#fwrite(t(ExpressionMat_by_gp),"mediam_gene_expression.tsv",sep="\t",col.names=T,row.names=T,quote=F)

save(ExpressionMat_by_gp,Exp_Seurat,sample_list,cell_assignment_gps,patient_type,tissue_type,BCR_annoation,TCR_annoation,BCR_annoation_w,TCR_annoation_w,t2g,taggene_mt,taggene_hk,HVG,file="180529_RUVscale_800_healthy.Rdata")

library(pryr)
message(paste0(mem_used()/1000/1000/1000,"Gb"))

rm(list=setdiff(ls(),c("Exp_Seurat",second_to_humanReadableTime)))

suppressMessages(require(feather))

#MAGIC
forMagic="For_magic.feather"
afterMagic="Magic_out.feather"

temp_Exp=as.data.frame(t(Exp_Seurat@scale.data))
#magic
write_feather(temp_Exp, forMagic)
system(paste0("python /home/ahe/tools/Lightbulb/src/Lightbulb_MagicWrapper.py --matx ",forMagic," --out ",afterMagic))
ExpressionMat_magic = data.matrix(read_feather(afterMagic))
colnames(ExpressionMat_magic)=gsub("MAGIC ", "", colnames(ExpressionMat_magic))

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

save(ExpressionMat_magic,file="180529_RUVscale_800_magic_healthy.Rdata")

#save(ExpressionMat_magic,sample_list,ExpressionNormed,ExpressionMat,pca_allcell,rtsne_allcell,BCR_annoation,cell_assignment,cell_assignment_gps,sample_gp,TCR_annoation,BCR_annoation_w,TCR_annoation_w,patient_type,tissue_type,patient_assignment,t2g,tissue_assignment,taggene_mt,taggene_hk,HVG,file="180521_RUVscale_800.Rdata")




library(pryr)
message(paste0(mem_used()/1000/1000/1000,"Gb"))
