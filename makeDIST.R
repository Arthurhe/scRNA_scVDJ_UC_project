require(NetworkToolbox)
load("180505_TPMandRUVscale_1k.Rdata")

dist_Exp=dist(ExpressionNormed)
save(dist_Exp,file="180505_DistLouvain.Rdata")
louvain_out=louvain(dist_Exp, 1)
save(louvain_out,dist_Exp,file="180505_DistLouvain.Rdata")