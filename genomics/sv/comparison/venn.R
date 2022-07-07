library(ggVennDiagram)
library(ggplot2)
library(eulerr)

ill = read.table("illumina.calls",header=F)
ont = read.table("ont.calls",header=F)
x = list(Illumina=as.vector(as.character(ill$V1)), ONT=as.vector(as.character(ont$V1)))

p = ggVennDiagram(x)
p

euler_plot = euler(x)
plot(euler_plot, quantities = TRUE, labels=list(font=4))

