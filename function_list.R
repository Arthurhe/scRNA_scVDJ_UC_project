Sampling_merging_seurat=function(tagSeurat,k_filter=15,k_merge=100,n=5000,tag_cell=NULL,verbose=1,sampling_ref=rep(1,ncol(tagSeurat@data))){
    #input a seurat object (tagSeurat), get randomly selected n cells as super cell seed. For each seed, merge k_merge nearest neighbor cells to the seed cell and calculate average expression.
    # tagSeurat@scale.data were used for calculation.
    #setup
    ptm <- proc.time()
    pca_out=tagSeurat@dr$pca@cell.embeddings
    kmax=max(k_filter,k_merge)
    kmax_plus1=kmax+1
    kfilter_plus1=k_filter+1
    kmerge_plus1=k_merge+1
    
    #mnn construction
    nn_idx=nn2(pca_out,k=kmax_plus1)$nn.idx
    for(i in 1:nrow(nn_idx)){ #<10s for 20K cells and k=15
        posi=nn_idx[i,] %in% i      
        if(sum(posi)>0){
            posi=which(posi)
            nn_idx[i,posi:kmax]=nn_idx[i,(posi+1):kmax_plus1]
        }
    }
    nn_idx=nn_idx[,-kmax_plus1]
    nn_idx_plusSelf=cbind(1:nrow(nn_idx),nn_idx[,1:k_merge])
    
    edgelist=data.frame(q=rep(1:nrow(pca_out),k_filter),d=(1:nrow(pca_out))[nn_idx[,1:k_filter]])
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("MNN done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}

    ptm <- proc.time()
    colnames(edgelist)=c("from","to")
    edgelist$from=paste0("V",edgelist$from)
    edgelist$to=paste0("V",edgelist$to)

    vertex_mat=data.frame(name=paste0("V",1:nrow(pca_out)),
                          id=rownames(tagSeurat@meta.data),
                          #sample_gp=tagSeurat@meta.data$sample_gp,
                          #patient_assignment=tagSeurat@meta.data$patient_assignment,
                          #tissue_assignment=tagSeurat@meta.data$tissue_assignment,
                          #disease_assignment=tagSeurat@meta.data$disease_assignment,
                          stringsAsFactors=F)
    net <- graph_from_data_frame(edgelist, vertices=vertex_mat, directed=T)
    net <- delete_edges(net, which(!which_mutual(net)))
    net <- as.undirected(net,mode="collapse")
    #net <- simplify(net)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("network constructed: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}

    if(is.null(tag_cell)){
        #valid targets
        tags=which(degree(net)>0)

        #select sample
        if(length(tags)<=n){
            warning(paste0("n bigger than meaningful sample number, n set to ",length(tags)))
            tag_cell=tags
        }else{
            tag_cell=tags[Subsample_by_group(sampling_ref[tags],n)]
        }
    }

    #merging cells to pseudo cells according to target cells
    merging_table=data.frame(center=rep(tag_cell,each=kmerge_plus1),
                             surrounding=as.vector(t(nn_idx_plusSelf[tag_cell,])),
                             stringsAsFactors=F)
    
    #diving merging into multiple chunks
    chunk_cap=40000
    pesudocell_per_chunk=ceiling(chunk_cap/kmerge_plus1)
    row_per_chunk=pesudocell_per_chunk*kmerge_plus1
    number_of_chunk=ceiling(nrow(merging_table)/row_per_chunk)

    ptm <- proc.time()
    chunks=list()
    for(i in 1:number_of_chunk){
        tagrow=((i-1)*row_per_chunk+1):min((i*row_per_chunk),nrow(merging_table))
        chunks[[i]]=t(as.matrix(tagSeurat@scale.data))[merging_table$surrounding[tagrow],]
        chunks[[i]]=cbind(merging_table$center[tagrow],chunks[[i]])
        colnames(chunks[[i]])[1]="cell"
        chunks[[i]]=data.table(chunks[[i]])
        chunks[[i]]=chunks[[i]][,lapply(.SD, mean), by=cell]
        #chunks[[i]]=chunks[[i]][order(chunks[[i]]$cell)]
        rn=chunks[[i]]$cell
        chunks[[i]][,cell:=NULL]
        chunks[[i]]=data.frame(chunks[[i]])
        colnames(chunks[[i]])=rownames(as.matrix(tagSeurat@scale.data))
        rownames(chunks[[i]])=rn
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose==2){message(paste("chunk ",i,"/",number_of_chunk," finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }
    tagMat=do.call(rbind,chunks)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("merging finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    return(tagMat)
}
        
Sampling_merging=function(tagMat,k_filter=15,k_merge=100,n=5000,tag_cell=NULL,verbose=1){
    #input a cell-gene matrix (tagMat), get randomly selected n cells as super cell seed. For each seed, merge k_merge nearest neighbor cells to the seed cell and calculate average expression. 
    tagMat_filtered=as.matrix(tagMat[,colSums(as.matrix(tagMat)>0)>ceiling(nrow(tagMat)/100)])
    
    #setup
    ptm <- proc.time()
    #variable gene
    HVG_idx=order(matrixStats::colSds(tagMat_filtered)/colMeans(tagMat_filtered),decreasing = T)[1:3000]
    
    pca_out=rsvd::rpca(tagMat_filtered[,HVG_idx],30)$x
    kmax=max(k_filter,k_merge)
    kmax_plus1=kmax+1
    kfilter_plus1=k_filter+1
    kmerge_plus1=k_merge+1
    
    #mnn construction
    nn_idx=nn2(pca_out,k=kmax_plus1)$nn.idx
    for(i in 1:nrow(nn_idx)){ #<10s for 20K cells and k=15
        posi=nn_idx[i,] %in% i      
        if(sum(posi)>0){
            posi=which(posi)
            nn_idx[i,posi:kmax]=nn_idx[i,(posi+1):kmax_plus1]
        }
    }
    nn_idx=nn_idx[,-kmax_plus1]
    nn_idx_plusSelf=cbind(1:nrow(nn_idx),nn_idx[,1:k_merge])

    edgelist=data.frame(q=rep(1:nrow(pca_out),k_filter),d=(1:nrow(pca_out))[nn_idx[,1:k_filter]])
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("MNN done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}

    ptm <- proc.time()
    colnames(edgelist)=c("from","to")
    edgelist$from=paste0("V",edgelist$from)
    edgelist$to=paste0("V",edgelist$to)

    vertex_mat=data.frame(name=paste0("V",1:nrow(pca_out)),
                          id=rownames(tagMat),
                          stringsAsFactors=F)
    net <- graph_from_data_frame(edgelist, vertices=vertex_mat, directed=T)
    net <- delete_edges(net, which(!which_mutual(net)))
    net <- as.undirected(net,mode="collapse")
    #net <- simplify(net)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("network constructed: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    
    if(is.null(tag_cell)){
        #valid targets
        tags=which(degree(net)>0)

        #select sample
        if(length(tags)<=n){
            warning(paste0("n bigger than meaningful sample number, n set to ",length(tags)))
            tag_cell=tags
        }else{
            tag_cell=sample(tags,n)
        }
    }

    #merging cells to pseudo cells according to target cells
    merging_table=data.frame(center=rep(tag_cell,each=kmerge_plus1),
                             surrounding=as.vector(t(nn_idx_plusSelf[tag_cell,])),
                             stringsAsFactors=F)
    
    #diving merging into multiple chunks
    chunk_cap=40000
    pesudocell_per_chunk=ceiling(chunk_cap/kmerge_plus1)
    row_per_chunk=pesudocell_per_chunk*kmerge_plus1
    number_of_chunk=ceiling(nrow(merging_table)/row_per_chunk)

    ptm <- proc.time()
    chunks=list()
    for(i in 1:number_of_chunk){
        tagrow=((i-1)*row_per_chunk+1):min((i*row_per_chunk),nrow(merging_table))
        chunks[[i]]=as.matrix(tagMat)[merging_table$surrounding[tagrow],]
        chunks[[i]]=cbind(merging_table$center[tagrow],chunks[[i]])
        colnames(chunks[[i]])[1]="cell"
        chunks[[i]]=data.table(chunks[[i]])
        chunks[[i]]=chunks[[i]][,lapply(.SD, mean), by=cell]
        #chunks[[i]]=chunks[[i]][order(chunks[[i]]$cell)]
        rn=chunks[[i]]$cell
        chunks[[i]][,cell:=NULL]
        chunks[[i]]=data.frame(chunks[[i]])
        colnames(chunks[[i]])=colnames(as.matrix(tagMat))
        rownames(chunks[[i]])=rn
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose==2){message(paste("chunk ",i,"/",number_of_chunk," finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }
    tagMat=do.call(rbind,chunks)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("merging finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    return(tagMat)
}
