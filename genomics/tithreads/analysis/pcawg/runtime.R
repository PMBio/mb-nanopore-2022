library(lubridate)

x=read.table("runtime.tsv", header=T, sep="\t")
x$begin = ymd_hms(x$start, tz = "Pacific/Auckland")
x$finish = ymd_hms(x$end, tz = "Pacific/Auckland")
print(median(with(x, difftime(finish,begin,units="secs"))))
