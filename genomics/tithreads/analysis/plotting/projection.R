library(ggplot2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
colnames(x) = c("chr", "start", "end", "read", "cstart", "cend", "direction", "mapq", "idchr", "idstart", "idend", "size")
x$idname = paste0(x$idchr, ":", x$idstart, "-", x$idend)
x$idname = factor(x$idname)

df = x[x$direction == "fwd",]
x = x[x$direction == "rev",]
x$tmp = x$start
x$start = x$end
x$end = x$tmp
x = x[,!(colnames(x) %in% c("tmp"))]
df = rbind(df, x)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))


png(paste0(args[1], ".png"), width=800, height=600)
#svg("alignments.svg", width=16, height=6)
p = ggplot(data=df)
p = p + geom_segment(aes(y=cstart, yend=cend, x=start, xend=end, color=direction), size=1.2)
p = p + facet_wrap(chr~idname, scales="free_x")
#p = p + facet_wrap(~chr, scales="free_x")
p = p + scale_x_continuous(labels=comma)
p = p + scale_y_continuous(labels=comma, limits=c(0, 60000))
p = p + ylab("ONT read") + xlab("Templated insertion source locus")
p = p + scienceTheme
p = p + scale_color_manual(name="Direction", values=c("#404040", "#ca0020"), labels=c("forward match", "reverse match"))
p = p + theme(legend.position="top")
p
dev.off()

