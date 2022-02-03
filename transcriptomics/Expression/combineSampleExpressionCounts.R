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

#Medullo.
expFolder = "./expression/StuttgartSample/"
fileOut = "./expression/StuttgartSample.gene_tpm.gct"

#Reference.
expFolder = "./expression/SRP151960/"
fileOut = "./expression/SRP151960.gene_tpm.gct"


toParse <- list.files(expFolder,full.names =T )

expOut = NULL

for(f in toParse){
  file = t(read.delim(f,as.is=T,skip=2,sep="\t"))
  if(is.null(expOut)){
    expOut=file
    rownames(expOut)[nrow(expOut)] = strsplit(f,split = "/")[[1]][4]
  } else {
    if(all(colnames(expOut)==colnames(file))){
      expOut = rbind(expOut,file[nrow(file),])
      rownames(expOut)[nrow(expOut)] = strsplit(f,split = "/")[[1]][4]
    } else{
      stop();
    }
  }
}
expOut = t(expOut)

for(col in 3:ncol(expOut)){
  expOut[,col] = as.numeric(expOut[,col])
}

write.table(expOut, file=fileOut, quote=F, sep="\t",row.names=F)
