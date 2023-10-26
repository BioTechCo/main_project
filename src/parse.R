# nolint start
library(ChAMP)
myLoad <- champ.load("train/stage_III")#, SampleCutoff = 0.2)#, arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,plotBMIQ=FALSE,cores=5)#, arraytype = "EPIC")

write.csv(myNorm,file="result/all_beta_normalized_stage_III.csv",quote=F,row.names=T)

myDMP <- champ.DMP(beta=myNorm, pheno=myLoad$pd$Sample_Group)#, arraytype = "EPIC")

write.csv(myDMP[[1]], file="result/DMP_result_stage_III.csv", quote=F)

# nolint end
