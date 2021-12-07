plotByGroup=function(Exp_Seurat,
                     Group_assignment,
                     Group_type=NULL, #provide the order of Group_type
                     random_order=T,
                     embedding_type="tsne",
                     plotGroupOnly=NULL, #only plot the identified group
                     backGroundGroup=NULL,
                     Addtext=F,
                     main_title="",
                     target_dot_number=10000,
                     tagcol=NULL){
    if(sum(is.na(Group_assignment))!=0){
        warning("NAs in Group_assignment")
    }
    #rm groups that's too small
    gp_size=table(Group_assignment)
    gp2rm=names(gp_size)[gp_size*target_dot_number*2/length(Group_assignment)<=20]
    if(length(gp2rm)>0){
        for(i in gp2rm){
            Group_assignment[Group_assignment==i]=NA
        }
    }
    
    Group_assignment=as.character(Group_assignment)
   
    if(is.null(Group_type)){
        Group_type=unique(Group_assignment)
        Group_type=Group_type[gtools::mixedorder(Group_type)]
    }else{
        Group_type=Group_type[Group_type %in% Group_assignment]
    }
    
    #color by group id
    if(is.null(tagcol)){
        if(length(Group_type)<=9){
            tagcol=RColorBrewer::brewer.pal(max(3,length(Group_type)),"Set1")
        }else{
            tagcol=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))
            tagcol=tagcol(length(Group_type))
        }
    }
    
    eval(parse(text=paste0("plotting_coordinates=Exp_Seurat@dr$",embedding_type,"@cell.embeddings[,1:2]")))
    #reduce dot number to target_dot_number
    tot_dot_num=nrow(plotting_coordinates)
    if(tot_dot_num>target_dot_number){
        get1w=Subsample_by_group(Group_assignment,target_dot_number)
        plotting_coordinates=plotting_coordinates[get1w,]
        cluster_assign=Group_assignment[get1w]
    }else{
        cluster_assign=Group_assignment
    }
    
    if(random_order){
        random_order=sample(1:length(cluster_assign),length(cluster_assign))
    }else{
        random_order=1:length(cluster_assign)
    }
    
    #calculate cluster center
    if(Addtext){
        cluster_center=matrix(0,length(Group_type),2)
        not_present=c()
        for(i in 1:length(Group_type)){
            tagcoord=plotting_coordinates[cluster_assign==Group_type[i],]
            if(sum(cluster_assign==Group_type[i],na.rm=T)==0){
                not_present=c(not_present,i)
                next
            }
            quick_kmean=kmeans(tagcoord,centers=2)
            center=quick_kmean$centers[which.max(table(quick_kmean$cluster)),]
            cluster_center[i,]=tagcoord[which.min(proxy::dist(tagcoord,t(center))),]
        }
    }
    
    maxx=max(plotting_coordinates[,1])
    minx=min(plotting_coordinates[,1])
    maxy=max(plotting_coordinates[,2])
    miny=min(plotting_coordinates[,2])  
    if(is.null(plotGroupOnly)){
        if(!is.null(backGroundGroup)){
            background_cell=which(cluster_assign %in% backGroundGroup)
            random_order_b=random_order[random_order %in% background_cell]
            if(length(background_cell)==0){stop("backGroundGroup not found in assignment")}
            plot(plotting_coordinates[random_order_b,],pch=19,col="gray",cex=0.5,
                 main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))
            foreground_cell=which(!cluster_assign %in% backGroundGroup)
            random_order=random_order[random_order %in% foreground_cell]
            points(plotting_coordinates[random_order,],pch=19,col=tagcol[match(cluster_assign,Group_type)][random_order],cex=0.5)
            tagcol[Group_type %in% backGroundGroup]="gray"
        }else{
            plot(plotting_coordinates[random_order,],pch=19,col=tagcol[match(cluster_assign,Group_type)][random_order],cex=0.5,
                 main=paste(main_title,embedding_type),xlim=c(minx,maxx+0.3*(maxx-minx)))        
        }
        if(Addtext){text(cluster_center[,1],cluster_center[,2],label=Group_type,font=2,cex=1.2)}
        legend("topright",Group_type,bty = "n",lty=0,pch=19,col=tagcol)
    }else{
        if(!plotGroupOnly %in% Group_type){
            stop("cell type to plot is missing in current data")
        }
        tag_cluster_num=match(plotGroupOnly,Group_type)
        plot(plotting_coordinates[!cluster_assign %in% plotGroupOnly,],pch=19,col="gray",cex=0.5,
             main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))
        for(i in 1:length(plotGroupOnly)){
            points(plotting_coordinates[cluster_assign==plotGroupOnly[i],],pch=19,col=tagcol[which(Group_type==plotGroupOnly[i])],cex=0.5)
        }
        if(Addtext){text(cluster_center[tag_cluster_num,1],cluster_center[tag_cluster_num,2]
                         ,label=Group_type[tag_cluster_num],font=2,cex=1.2)}      
        legend("topright",Group_type[tag_cluster_num],bty = "n",lty=0,pch=19,col=tagcol[tag_cluster_num])
    }
}

plotByGroup_1by1=function(Exp_Seurat,
                          Group_assignment,
                          Group_type=NULL, #the group ordering
                          embedding_type="tsne",
                          plotGroupOnly=Group_type, #the group to plot
                          main_title="",
                          target_dot_number=10000
                         ){
    for(i in 1:length(plotGroupOnly)){
        plotByGroup(Exp_Seurat,
                       Group_assignment=Group_assignment,
                       Group_type=Group_type,
                       random_order=F,
                       embedding_type=embedding_type,
                       plotGroupOnly=Group_type[i],
                       Addtext=T,
                       main_title=main_title,
                       target_dot_number=target_dot_number)
    }
}


plotByScore=function(Exp_Seurat,
                     Score_assignment,
                     random_order=T,
                     embedding_type="tsne",
                     main_title="",
                     colorSD=1,
                     target_dot_number=10000)
{   
    Score_assignment=as.numeric(Score_assignment)
    eval(parse(text=paste0("plotting_coordinates=Exp_Seurat@dr$",embedding_type,"@cell.embeddings[,1:2]")))
    #reduce dot number to target_dot_number
    tot_dot_num=nrow(plotting_coordinates)
    if(tot_dot_num>target_dot_number){
        get1w=Subsample_by_group(Exp_Seurat@ident,target_dot_number)
        plotting_coordinates=plotting_coordinates[get1w,]
        Score_assignment=Score_assignment[get1w]
    }
    
    if(random_order){
        random_order=sample(1:length(Score_assignment),length(Score_assignment))
    }else{
        random_order=1:length(Score_assignment)
    }
     
    rbPal <- colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdBu")))
    lowerbound=mean(Score_assignment)-colorSD*sd(Score_assignment)
    upperbound=mean(Score_assignment)+colorSD*sd(Score_assignment)
    breaks=c(-Inf,seq(lowerbound, upperbound, length.out=25),Inf)
    cols=rbPal(25)[as.numeric(cut(Score_assignment,breaks = breaks, include.lowest=TRUE))]

    maxx=max(plotting_coordinates[,1])
    minx=min(plotting_coordinates[,1])
    maxy=max(plotting_coordinates[,2])
    miny=min(plotting_coordinates[,2]) 
    plot(plotting_coordinates[random_order,],pch=19,cex=0.5,col=cols[random_order],
         main=paste(main_title,embedding_type),xlim=c(minx,maxx+0.3*(maxx-minx)))  
}


