library(Gviz)

itrack = IdeogramTrack(genome="hg38", chromosome="4", name="")
png("chr4.png", width=600, height=100)
plotTracks(itrack, from=168398333, to=168398412)
dev.off()

itrack = IdeogramTrack(genome="hg38", chromosome="5", name="")
png("chr5.png", width=600, height=100)
plotTracks(itrack, from=18569089, to=18569221)
dev.off()

itrack = IdeogramTrack(genome="hg38", chromosome="7", name="")
png("chr7.reg1.png", width=600, height=100)
plotTracks(itrack, from=7805332, to=7805395)
dev.off()

itrack = IdeogramTrack(genome="hg38", chromosome="7", name="")
png("chr7.reg2.png", width=600, height=100)
plotTracks(itrack, from=82485397, to=82485523)
dev.off()



