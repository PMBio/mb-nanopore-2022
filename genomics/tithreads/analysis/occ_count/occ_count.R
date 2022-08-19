library(ggplot2)

x=read.table("segment.counts.cov.bed", header=T)

print(cor(x$count, x$ill, method="pearson"))
print(cor(x$count, x$ill, method="spearman"))
png("occ_count.png")
p = ggplot(data=x, aes(x=count, y=ill))
p = p + geom_point()
p = p + xlab("Occurence count of source segment in assembled Shasta contig")
p = p + ylab("Read-depth estimated copy-number (short-read data)")
p = p + xlim(0,80)
p = p + ylim(0,80)
p = p + geom_hline(yintercept=2, linetype="dashed")
p
dev.off()
