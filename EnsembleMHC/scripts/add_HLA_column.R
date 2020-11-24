a <- commandArgs(trailingOnly = T)

x <- read.csv(a[1])

x$HLA <- a[2]

write.csv(x, a[1], row.names = F)
