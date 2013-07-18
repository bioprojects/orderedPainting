#
# R --vanilla --quiet < lib/check_missingCnt_distStat.R --args  jejuni200_default.gap20.absPos2missingInd  jejuni200_default.gap20.bgl_orderedS1_rnd_1_10_both_dirs_results/results_siteStats.txt.gz
#

par(cex.axis=1.5)

c_args <- commandArgs(trailingOnly=T)
if (length(c_args) != 2) {
  stop("Error: num of command line args must be 2")
}

###################################################

d_pos2missingInd <- read.table(c_args[1],header=F)

d_pos2missingCnt <- data.frame()
d_pos2missingCnt <- as.data.frame(table(d_pos2missingInd[,1]))
names(d_pos2missingCnt) <- c("pos","missingCnt")

d_siteStats <- read.table(gzfile(c_args[2]),header=T)

d_merged    <- merge(d_siteStats, d_pos2missingCnt, all=T)
d_merged[is.na(d_merged$missingCnt), "missingCnt"] <- 0

stopifnot(nrow(d_merged)==nrow(d_siteStats))

#cor.value <- cor(d_merged$missingCnt, d_merged$sum_distScore,  method = c("spearman"))

#
# binning
#
bin_size  <- 10
png_fname <- sub(".txt.gz","_missingCnt.png", c_args[2])

d_merged_xsorted <- d_merged[order(d_merged$missingCnt),]

list_d_each_bin <- list()
d_bin_value <- data.frame()
for (i in 1:round(nrow(d_merged_xsorted)/bin_size)) {
  range <- ((i-1)*bin_size+1):(i*bin_size)
  d_each_bin <- d_merged_xsorted[ range, ]
  d_each_bin <- d_each_bin[!is.na(d_each_bin[,1]), ] # to deal with the last bin

  mean_per_bin <- mean(d_each_bin$sum_distScore)
  d_bin_value <- rbind(d_bin_value, data.frame(i=i, mean=mean_per_bin))

  list_d_each_bin[[i]] <- d_each_bin
}

png(png_fname,width=2000)
print(png_fname)
barplot(t(d_bin_value$mean), ylim=c(summary(d_bin_value$mean)[1], max(d_bin_value$mean))
  , main="SNP bins sorted by missing frequency (from low to high)"
  , ylab="Distance statistic per bin on average"
)
dev.off()
