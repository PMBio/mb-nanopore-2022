library(ggplot2)
library(scales)
library(gtable)
library(grid)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.y=element_line(size=0.7, color="black"), axis.line.x=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], sep="\t", header=T)

png("vaf.png")
p = ggplot(data=x, aes(x=tumor_vaf, y=relapse_vaf))
p = p + geom_jitter(width=0.01, height=0.01)
p = p + xlab("Tumor variant allele frequency")
p = p + ylab("Relapse variant allele frequency")
p = p + scienceTheme
p
dev.off()



