library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)
colnames(x) = c("pos", "af", "gt", "h1", "h2")
x = melt(x, id.vars=c("pos", "af", "gt"))

png(paste0(args[1], ".png"), width=1200, height=600)
p = ggplot(data=x, aes(x=pos, y=value))
p = p + geom_point(aes(color=variable), size=0.3)
p = p + ylab("Phased het. variant allele frequency")
p = p + xlab("Position")
p = p + labs(color="Haplotype")
p = p + theme(legend.position = "none")
p
dev.off()
