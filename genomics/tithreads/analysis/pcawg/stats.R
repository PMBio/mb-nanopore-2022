library(ggplot2)
library(dplyr)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))


args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

x = x[x$affected != "unclear",]
x$affected = factor(x$affected, levels=c("yes", "no"))
print(summary(x))


png("stats.png", width=1200, height=400)
p = ggplot(data=x, aes(x=project))
p = p + geom_bar(aes(fill=affected), stat="count", color="black")
p = p + xlab("Project") + ylab("Count")
p = p + scienceTheme
p = p + xlab("PCAWG project")
p
dev.off()

