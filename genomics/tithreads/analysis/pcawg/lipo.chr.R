library(ggplot2)

x = read.table("lipo.sarc-us.tmpl_ins.bed", comment.char="$", quote="$", sep="\t", header=F)
colnames(x) = c("sample", "chr", "start", "end", "nodeid", "selfdegree", "degree", "estcn", "clusterid", "edges")

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))


#chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX", "chrY")
chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X", "Y")

x$chr = factor(x$chr, levels=chrNamesShort)

png("chr.dist.png", width=500, height=400)
p = ggplot(data=x, aes(x=chr))
p = p + geom_histogram(stat="count", fill="darkgrey", color="black")
p = p + xlab("Chromosome")
p = p + ylab("Ascertained\ntemplated insertions")
p = p + scale_x_discrete(drop=F)
p = p + scienceTheme
p = p + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
p
dev.off()

x = x[x$chr == "12",]
png("chr12.png", width=200, height=150)
p = ggplot(data=x, aes(x=start, y=sample))
p = p + geom_point(size=2)
p = p + xlab("chr12")
p = p + ylab("Liposarcomas")
p = p + scienceTheme
p = p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
p = p + geom_hline(aes(yintercept=sample), color="lightgray")
p
dev.off()


