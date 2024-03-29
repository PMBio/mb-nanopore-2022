folder = "C:\\OnlineFolders\\BitSync\\CurrentWork\\Medulloblastoma\\"
#folder = "D:\\OnlineFolders\\BitSync\\CurrentWork\\Medulloblastoma\\"
setwd(folder)

##GTEx gtf.
gtf = read.delim("gencode.v26.GRCh38.genes.gtf.gz",as.is=T,skip = 6, header=F)

# gtf_exon = gtf[which(gtf[,3]=="exon"),]
# rm(gtf)
# 
# gene_id = gsub(pattern = "gene_id ",replacement = "",lapply((strsplit(gtf_exon$V9,split = ";")),'[[',1))
# gene_id = unlist(lapply(strsplit(gene_id,"\\."),'[[',1))
# gtf_exon = cbind(gtf_exon,gene_id)
# gtf_exon$gene_id = as.character(gtf_exon$gene_id)
# gtf_exon = gtf_exon[,c(10,1,4,5,7)]

#Stuttgart.
matrix1 <- read.delim("./expression/StuttgartSample.gene_tpm.gct.gz",as.is=T)
matrix1 = matrix1[order(matrix1[,1]),]

#Reference.
matrix2 <- read.delim("./expression/SRP151960.gene_tpm.gct.gz",as.is=T)
matrix2 = matrix2[order(matrix2[,1]),]

all(which(matrix1[,1] == matrix2[,1]))

combined = cbind(matrix1,matrix2[,c(3:ncol(matrix2))])


##GTEx meta
metaData = read.delim("./GTEx/GTEX_AE_METADATA.txt",as.is=T)
metaData = metaData[which(metaData$Comment.original.body.site.annotation.=="Brain - Cerebellum"),]
metaData[,1] = gsub(metaData[,1],pattern = "-",replacement = ".")
# GTExExp = read.delim("./GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",as.is=T,skip=2)
# saveRDS(GTExExp, "./GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.Rds")
GTExExp <- readRDS("./GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.Rds")

GTExExp = GTExExp[,c(1,2,which(colnames(GTExExp) %in% metaData[,1]))]
GTExExp = GTExExp[order(GTExExp[,1]),]

all(which(combined[,1] == GTExExp[,1]))

combined = cbind(combined,GTExExp[,c(3:ncol(GTExExp))])

##
annotation<- combined[,c(1,2)]
combined = combined[,-2]
rownames(combined)=combined[,1]
combined = combined[,-1]

nonExpressedGenes = which(rowSums(combined)==0)
#View(combined[nonExpressedGenes,])
combined = combined[-nonExpressedGenes,]
###

pcs = prcomp(t(combined))

write.table(combined,"./expression/GTEX_SRP151960_Medulloblastoma.exp_gene_tpm.gct.txt",quote=F, sep="\t",col.names=NA)

plot(pcs$x)

pcs = prcomp(t(combined[,-c(1,2)]),center = T)

combinedNorm = sweep(x = combined, MARGIN = 1, STATS = pcs$center,FUN = "-")
pcValues = t(as.matrix(combinedNorm)) %*% as.matrix(pcs$rotation)

write.table(pcValues,"./expression/GTEX_SRP151960_Medulloblastoma.exp_gene_tpm.gct.PCs.txt",quote=F, sep="\t",col.names=NA)

combinedNorm = sweep(x = combined, MARGIN = 1, STATS = pcs$center,FUN = "-")
zMatrix = combinedNorm
for(g in 1:nrow(combined)){
  combinedNorm[g,] = residuals(lm(as.numeric(combined[g,])~(as.numeric(pcValues[,1]))+(as.numeric(pcValues[,2]))+(as.numeric(pcValues[,3]))+(as.numeric(pcValues[,4]))+(as.numeric(pcValues[,5]))))
  meanV = mean(as.numeric(combinedNorm[g,]))
  sdV = sd(as.numeric(combinedNorm[g,]))
  zMatrix[g,] = ((combinedNorm[g,]-meanV)/sdV)
}

write.table(combinedNorm,"./expression/GTEX_SRP151960_Medulloblastoma.exp_gene_tpm.gct.PC.corrected.txt",quote=F, sep="\t",col.names=NA)
write.table(combinedNorm,"./expression/GTEX_SRP151960_Medulloblastoma.exp_gene_tpm.gct.PC.corrected.Zmatrix.txt",quote=F, sep="\t",col.names=NA)