library(reshape2)
library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.y=element_line(size=0.7, color="black"), axis.line.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize), axis.ticks.x=element_blank())


args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)
colnames(x) = c("pos", "af", "gt", "h1", "h2")
x = melt(x, id.vars=c("pos", "af", "gt"))

png(paste0(args[1], ".png"), width=1200, height=600)
p = ggplot(data=x, aes(x=pos, y=value))
p = p + geom_point(aes(color=variable), size=0.3)
p = p + ylab("Phased het. variant allele frequency")
p = p + xlab("")
p = p + scienceTheme
p = p + labs(color="Haplotype")
p = p + theme(legend.position = "none")
p
dev.off()
