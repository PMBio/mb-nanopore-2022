library(DNAcopy)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)


txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))


fmt_decimals = function(decimals=2) {
	     function(x) format(x, nsmall=decimals, scientific=F)
}

gg_color_hue = function(n) {
  hues=seq(15,375,length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

args=commandArgs(trailingOnly=TRUE)
idname = gsub(".vaf.tsv", "", args[2])
cov = read.table(args[1], header=T)
cov = cov[,c(1,2,3,6)]
colnames(cov) = c("chr", "start", "end", "cn")
snv = read.table(args[2], header=F)
colnames(snv) = c("chr", "pos", "h1af")
allSV = read.table(args[3], header=F)
colnames(allSV) = c("chr", "start", "chr2", "end", "svtype")
allSV$svtype = factor(allSV$svtype, levels=c("DEL", "DUP", "INS", "INV", "BND"))
svPalette = gg_color_hue(5)
reg = args[4]

chrname=unlist(strsplit(reg, ':'))[1]
interval=unlist(strsplit(reg, ':'))[2]
startC = as.integer(unlist(strsplit(interval, '-'))[1])
endC = as.integer(unlist(strsplit(interval, '-'))[2])
svStart = startC
svEnd = endC

# Subset to region of interest
cov = cov[cov$chr==chrname,]
cov = cov[cov$start >= svStart,]
cov = cov[cov$end <= svEnd,]
snv = snv[snv$chr==chrname,]
if (nrow(snv) > 0) {
 snv = snv[snv$pos >= svStart,]
 snv = snv[snv$pos <= svEnd,]
 minX = min(min(snv$pos), min(cov$start))
 maxX = max(max(snv$pos), max(cov$end))
} else {
 minX = min(cov$start)
 maxX = max(cov$end)
}

for (chrlabel in unique(cov$chr)) {
    snvchr = snv[snv$chr==chrlabel,]
    sv = allSV[allSV$chr==chrlabel | allSV$chr2==chrlabel,]

    # Segment coverage
    reschr = cov[cov$chr==chrlabel,]$cn
    mid = (cov[cov$chr==chrlabel,2] + cov[cov$chr==chrlabel,3])/2
    cnaData = CNA(reschr, maploc=mid, chrom=rep(chrlabel, length(mid)), sampleid=idname, presorted=T)
    cnaDataSmooth=smooth.CNA(cnaData)
    cnaSegments=segment(cnaDataSmooth, undo.splits="sdundo", undo.SD=1)
    cnaSummary=segments.summary(cnaSegments)
    #write.table(cnaSummary, file=paste0(idname, ".", chrlabel, ".cov.segment"), sep="\t", quote=FALSE, row.names=FALSE)
    df = data.frame(pos=mid, signal=reschr)

    # Segment h1
   if (nrow(snvchr) > 0) {  
    reschr = snv[snv$chr==chrlabel,3]
    start = snv[snv$chr==chrlabel,2]
    cnaData = CNA(reschr, maploc=start, chrom=rep(chrlabel, length(start)), sampleid=idname, presorted=T)
    cnaDataSmooth=smooth.CNA(cnaData)
    cnaSegments=segment(cnaDataSmooth, undo.splits="sdundo", undo.SD=1)
    h1 = segments.summary(cnaSegments)
    #write.table(h1, file=paste0(idname, ".", chrlabel, ".h1.segment"), sep="\t", quote=FALSE, row.names=FALSE)
    reschr = 1 - snv[snv$chr==chrlabel,3]
    start = snv[snv$chr==chrlabel,2]
    cnaData = CNA(reschr, maploc=start, chrom=rep(chrlabel, length(start)), sampleid=idname, presorted=T)
    cnaDataSmooth=smooth.CNA(cnaData)
    cnaSegments=segment(cnaDataSmooth, undo.splits="sdundo", undo.SD=1)
    h2 = segments.summary(cnaSegments)
    #write.table(h2, file=paste0(idname, ".", chrlabel, ".h2.segment"), sep="\t", quote=FALSE, row.names=FALSE)
   }
    p1 = ggplot(data=df, aes(x=pos, y=signal)) + geom_point(size=0.5, pch=21, fill="black", color="black")
    p1 = p1 + xlab(chrlabel) + ylab("Copy-number") + scale_x_continuous(labels=comma, limits=c(minX, maxX))
    p1 = p1 + scienceTheme
    #p1 = p1 + geom_segment(data=cnaSummary, aes(x=loc.start, y=seg.median, xend=loc.end, yend=seg.median), colour="lightblue")
    p1 = p1 + scale_y_continuous(labels=fmt_decimals(2), breaks=c(0,1,2,3,4,5,6,7,8), limits=c(0,10)) + theme(legend.position="bottom")
    if (nrow(sv)) {
       svplotpos = 7
       svarc = sv[sv$svtype != "BND" & (sv$end - sv$start > 1000000),]
       if (nrow(svarc)) { p1=p1 + geom_curve(data=svarc, aes(x=start, xend=end, color=svtype), y=svplotpos, yend=svplotpos, curvature=-0.3, arrow=arrow(length=unit(0.25, "cm"), ends="both", type="closed")); }
       svarc = sv[sv$svtype != "BND" & (sv$end - sv$start <= 1000000),]
       if (nrow(svarc)) {
       	  p1=p1 + geom_curve(data=svarc, aes(x=start, xend=start, color=svtype), y=svplotpos + 1.5, yend=svplotpos, curvature=0.8, angle=60, arrow=arrow(length=unit(0.25, "cm"), type="closed"));
      	  p1=p1 + geom_curve(data=svarc, aes(x=start, xend=start, color=svtype), y=svplotpos + 1.5, yend=svplotpos, curvature=-0.8, angle=-60, arrow=arrow(length=unit(0.25, "cm"), type="closed"));
       }
       chr1sv = sv[sv$svtype == "BND" & sv$chr==chrlabel,]
       if (nrow(chr1sv)) { p1 = p1 + geom_segment(data=chr1sv, aes(x=start, xend=start, color=svtype), y=svplotpos + 1.5, yend=svplotpos, arrow=arrow(length=unit(0.25, "cm"), type="closed")); }
       chr2sv = sv[sv$svtype == "BND" & sv$chr2==chrlabel,]
       if (nrow(chr2sv)) { p1 = p1 + geom_segment(data=chr2sv, aes(x=end, xend=end, y=svplotpos + 1.5, yend=svplotpos, color=svtype), arrow=arrow(length=unit(0.25, "cm"), type="closed")); }
       p1 = p1 + scale_color_manual(values=svPalette, drop=F) + labs(colour="Structural variants") + theme(legend.key.size=unit(2, "line"))
    }
    p2 = ggplot(data=snvchr, aes(x=pos, y=h1af)) + geom_jitter(color=gg_color_hue(2)[1], size=0.25, width=0, height=0)
    p2 = p2 + geom_jitter(data=snvchr, aes(x=pos, y=1-h1af), color=gg_color_hue(2)[2], size=0.25, width=0, height=0)
    #p2 = p2 + geom_segment(data=h1, aes(x=loc.start, y=seg.median, xend=loc.end, yend=seg.median), colour=gg_color_hue(2)[1], size=0.8)
    #p2 = p2 + geom_segment(data=h2, aes(x=loc.start, y=seg.median, xend=loc.end, yend=seg.median), colour=gg_color_hue(2)[2], size=0.8)
    p2 = p2 + xlab(chrlabel) + ylab("Phased het.\nSNV VAF") + scale_x_continuous(labels=comma,limits=c(minX, maxX)) + ylim(0,1)
    p2 = p2 + scienceTheme
    p2 = p2 + theme(axis.title.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())
    plot_grid(p2, p1, align="v", nrow=2, rel_heights=c(1/3, 2/3))
    #ggsave(paste0(idname, ".", chrlabel, ".pdf"), width=16, height=8)
    ggsave(paste0(idname, ".", chrlabel, ".png"), width=16, height=6)
}
print(warnings())

