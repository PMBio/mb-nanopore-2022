library(ggplot2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))


dat = read.table("summary.stats.tsv", header=T)
dat$coverage = factor(dat$coverage)
dat$sd = factor(dat$sd)

# F1 score lines
df = data.frame()
for (f1 in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  for(p in (1:1000)/1000) {
     r = f1 * p / (2 * p - f1)
     if ((r >= 0) && (r <= 1)) {
        df = rbind(df, data.frame(recall=r, precision=p, f1=f1))
     }
  }
}
df$f1 = factor(df$f1, levels=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
labels = df[df$precision == 1,]
labels$text = paste0("F1=", labels$f1)
labels$recall = labels$recall - 0.01

p = ggplot(data=dat, aes(x=precision, y=recall))
p = p + geom_point(aes(shape=coverage, color=sd), size=2)
p = p + xlim(0,1) + ylim(0,1)
p = p + geom_line(data=df, aes(x=precision, y=recall, group=f1), color="gray", linetype="dashed")
p = p + geom_text(data=labels, aes(x=precision, y=recall, label=text), color="darkgray")
p = p + xlab("Precision")
p = p + ylab("Recall")
p = p + labs(color="Sd", shape="Coverage")
p = p + scienceTheme
p






  