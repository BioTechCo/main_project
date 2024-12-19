# nolint start
library(ChAMP)
myLoad <- champ.load("data/TCGA_COAD_READ")#, SampleCutoff = 0.2)#, arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,plotBMIQ=FALSE,cores=5)#, arraytype = "EPIC")

write.csv(myNorm,file="all_beta_normalized_TCGA_COAD_READ.csv",quote=F,row.names=T)

myDMP <- champ.DMP(beta=myNorm, pheno=myLoad$pd$Sample_Group)#, arraytype = "EPIC")

write.csv(myDMP[[1]], file="DMP_result_TCGA_COAD_READ.csv", quote=F)

# nolint end
