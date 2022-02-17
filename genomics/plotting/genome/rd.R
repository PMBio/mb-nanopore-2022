library(ggplot2)
library(scales)
library(gtable)
library(grid)

txtFontSize=16
axisFontSize=32
axisTtlFontSize=32
stripTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=stripTtlFontSize), strip.text.y=element_text(size=stripTtlFontSize))


chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX")
chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
maxCN = 8
seg = data.frame()
if (length(args)>1) {
   seg = read.table(args[2], header=F, sep="\t")
   colnames(seg) = c("chr", "start", "end", "id", "cn")
}

# Fix chromosome ordering
if (sum(x$chr %in% chrNamesLong) > sum(x$chr %in% chrNamesShort)) { chrs = chrNamesLong; } else { chrs = chrNamesShort; }
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)
if (nrow(seg) > 0) {
 seg = seg[seg$chr %in% chrs,]
 seg$chr = factor(seg$chr, levels=chrs)
}

# Whole genome
png("plot.wholegenome.png", width=600, height=1800)
p = ggplot(data=x, aes(x=start, y=x[,6]))
p = p + geom_point(pch=21, color="black", fill="black", size=0.5)
p = p + xlab("Chromosome")
p = p + ylab("Copy-number")
p = p + scale_x_continuous(labels=comma)
if (nrow(seg)) { p = p + geom_segment(data=seg, aes(x=start, y=cn, xend=end, yend=cn), color="#31a354", size=1.2); }
p = p + facet_wrap(~ chr, nrow = 6, scales="free_x")
p = p + ylim(0, maxCN)
p = p + scienceTheme
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p
dev.off()

