#
# R --vanilla --quiet < lib/visualization1.R --args lib/plotHeatmap.R  results_siteStats.txt.gz  
#
library(MASS)

c_args <- commandArgs(trailingOnly=T)
if (length(c_args) != 2) {
  stop("Error: num of command line args must be 2")
}

source(c_args[1])

###################################################

gzfile_site_dist <- c_args[2]
if (!file.exists(gzfile_site_dist)) {
  stop(sprintf("Error: %s doesn't exist"))
}

d_site_dist <- read.table(gzfile(gzfile_site_dist),header=T)

# top 1%
dist_top_percentile <- round( d_site_dist[ round(nrow(d_site_dist)*0.01), 2] )

# histgram
out_histfile <- sub(".txt.gz","_hist.png", gzfile_site_dist)
png(out_histfile)
truehist( d_site_dist$sum_distScore ,h=1,col=F,prob=F
          , xlab="D_i (intensity of recombination)" 
          , main=sprintf("top percentile = %s", dist_top_percentile)
          ) 
dev.off()

# plot along the sequence
out_along_sequence_file <- sub(".txt.gz","_along_seq.png", gzfile_site_dist)
png(out_along_sequence_file)
plot( d_site_dist$pos, d_site_dist$sum_distScore, pch=16,cex=0.2
      ,  xlab="genomic position" 
      ,  ylab="D_i (intensity of recombination)" 
      )
dev.off()

