suppressMessages(require(data.table))
suppressMessages(require(Rtsne))
devtools::load_all("/home/ahe/tools/Lightbulb")

ptm <- proc.time()

load("mnnReturn_VDJ.Rdata")

message("loading  done")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#PCA
PCAresult=prcomp(ExpressionMat_norm)
save(PCAresult,file="TSNE_loading_TPM_PCA_VDJ_MNNNormed.Rdata")
message("PCA out")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

#TSNE
rtsne_result3d=Rtsne(PCAresult$x[,1:23],dims=3,max_iter = 1500, pca=F)
rtsne_result2d=Rtsne(PCAresult$x[,1:23],dims=2,max_iter = 1500, pca=F)
save(PCAresult,rtsne_result2d,rtsne_result3d,ExpressionMat_norm,file="Rtsne_result_1st_VDJ_MNNNormed.Rdata")
message("1st result out")
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))