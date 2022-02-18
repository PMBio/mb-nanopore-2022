library(ggplot2)
library(scales)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
read = args[2]
x = x[x$read==read,]
print(summary(x))
#x$idname = paste0(x$chr, ":", comma(x$refstart, accuracy=1), "-", comma(x$refend, accuracy=1), " ")
x$idname = paste0(comma(x$refstart, accuracy=1), " - ", comma(x$refend, accuracy=1), " ")
x$idname = factor(x$idname)
x$align = '←'
x[x$direction=="fwd",]$align = '→'


png(paste0(read, ".png"), width=500, height=400)
p = ggplot(data=x)
p = p + geom_hline(aes(yintercept=idname), linetype="dashed", color="lightgray")
p = p + geom_segment(aes(x=readstart, xend=readend, y=idname, yend=idname, color=chr), size=10)
p = p + geom_text(aes(x=((readstart+readend)/2), y=idname, label=align), size=15)
p = p + xlab("ONT read")
p = p + ylab("Reference matches")
p = p + scienceTheme
p = p + ggtitle(args[3])
p = p + labs(color="")
p = p + theme(legend.position="top")
p
dev.off()

