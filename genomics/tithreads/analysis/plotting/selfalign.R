library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

# Dotplot
png(paste0(args[1], ".png"), width=800, height=600)
#svg("selfalign.svg", width=6, height=6)
p = ggplot(data=x)
p = p + geom_segment(aes(x=xstart, xend=xend, y=ystart, yend=yend, group=query, color=direction), size=1.2)
p = p + ylab("ONT read") + xlab("ONT read")
p = p + scienceTheme
p = p + theme(legend.position = c(0.9,0.5))
p = p + scale_color_manual(name="Direction", values=c("#404040", "#ca0020"), labels=c("forward match", "reverse match"))
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
p
dev.off()
