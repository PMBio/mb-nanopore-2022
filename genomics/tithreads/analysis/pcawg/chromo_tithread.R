library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))

x=read.table("table.tsv", header=T, comment.char="$", sep="\t")
s=read.table("stats.tsv", header=T, sep="\t")
print(head(x))
print(head(s))
x=merge(x, s, by=c("tumour_specimen_aliquot_id"))
x$tithread = 0
x[x$affected == "yes",]$tithread = 1
x$PCAWG_chromo = factor(x$PCAWG_chromo, levels=c(0, 1))
x$tithread = factor(x$tithread, levels=c(0, 1))
x$combined=paste0(x$tithread, "-", x$PCAWG_chromo)

print("Tumor samples")
print(nrow(x))
print("Fraction TI threads")
print(mean(as.numeric(x$tithread) - 1))
print("Table by histology")
print(table(x$tithread, x$histology_abbreviation))
print("Fraction by histology")
print(aggregate(as.numeric(x$tithread) - 1, list(x$histology_abbreviation), mean))
print("Chromothripsis")
print(table(x$tithread, x$PCAWG_chromo))

model = glm(x$PCAWG_chromo ~ histology_abbreviation + gender + ancestry_primary + tithread, family=binomial, data=x)
#model = glm(x$PCAWG_chromo ~ SV.events + histology_abbreviation + gender + ancestry_primary + tithread, family=binomial, data=x)
print(summary(model))

#print(fisher.test(table(x$tithread, x$PCAWG_chromo)))


png("dist.png", height=800, width=1200)
p = ggplot(data=x, aes(x=histology_abbreviation))
p = p + geom_bar(aes(fill=combined))
p = p + xlab("Tumor histology")
p = p + ylab("Count")
p = p + scale_fill_manual(values=c("grey", "darkgreen", "darkblue", "darkred"), labels=c("No chromothripsis & no TI thread", "Chromothripsis", "TI Thread", "Chromothripsis & TI Thread")) + labs(fill="Rearrangements")
p = p + scienceTheme
p = p + theme(legend.position="top", axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
p
dev.off()

png("sv.ti.png", height=800, width=1200)
p = ggplot(data=x, aes(x=histology_abbreviation, y=SV.events))
p = p + geom_point(aes(color=tithread))
p = p + xlab("Tumor histology")
p = p + ylab("Total number of SVs")
p = p + scienceTheme
p = p + theme(axis.text.x = element_text(angle=90, hjust=1))
p
dev.off()
