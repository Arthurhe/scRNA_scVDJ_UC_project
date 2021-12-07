suppressMessages(require(data.table))
suppressMessages(require(gplots))
suppressMessages(require(matrixStats))
suppressMessages(require(RColorBrewer))
require(NetworkToolbox)

load("180505_DistLouvain.Rdata")

louvain_out=louvain(as.matrix(dist_Exp), 1)
save(dist_Exp,louvain_out,file="180505_DistLouvain.Rdata")

hclust_out=hclust(dist_Exp, method="ward.D2")
save(dist_Exp,louvain_out,hclust_out,file="180505_DistLouvain.Rdata")
