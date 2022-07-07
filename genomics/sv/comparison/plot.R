library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))


x=read.table("size.tsv", header=F)
colnames(x)=c('chrom', 'pos', 'end', 'id', 'type', 'size', 'tech')
small = x[x$size>-1000 & x$size<1000,]

png("small.png", width=1200, height=600)
p = ggplot(data=small, aes(x=size))
p = p + geom_freqpoly(aes(group=tech, color=tech), bins=100)
p = p + xlab("Deletion size") + ylab("Number of ascertained deletions")
p = p + scienceTheme + scale_x_continuous(breaks=((1:20)*100 - 1000))
p = p + theme(legend.position="top")
p = p + labs(color="Sequencing Technology")
p
dev.off()
