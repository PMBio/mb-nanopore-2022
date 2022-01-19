library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
facetTtlFontSize=20
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=facetTtlFontSize))


args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
print(summary(x))

exp = x[x$obsexp=="expected",]
obs = x[x$obsexp=="observed",]

# Calculate p-values
print(obs)
mu = mean(exp$value)
sdval = sd(exp$value)
exp$zscore = (exp$value - mu) / sdval
obs$zscore = (obs$value - mu) / sdval
obs$p.value = 2*pnorm(q=abs(obs$zscore), lower.tail=FALSE)
print(obs)

# Plot
png(paste0(args[1], ".png"))
p = ggplot(data=exp, aes(x=value))
p = p + geom_histogram(aes(fill="lightblue"), color="black")
p = p + xlab("Number of overlaps") + ylab("Counts")
p = p + theme_classic()
p = p + geom_vline(xintercept=obs$value, linetype="dashed", color="red")
p = p + ggtitle(paste0("Enrichment\n", "p-value = ", round(obs$p.value, digits=8)))
p = p + scale_fill_manual(values=c("lightblue"), labels=c("Expected"), name="")
p = p + scienceTheme
p = p + theme(legend.position=c(0.2,1))
p
dev.off()
