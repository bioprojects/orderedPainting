#
# R --vanilla --quiet < lib/visualization2.R --args lib/plotHeatmap.R  pos_matrix.txt  min  max [ visualStrainOrder.list ]
#
library(MASS)

c_args <- commandArgs(trailingOnly=T)
if (length(c_args) < 4) {
  stop("Error: num of command line args must be >= 4")
}

source(c_args[1])

###################################################

pos_matrixFile <- c_args[2]
if (!file.exists(pos_matrixFile)) {
  stop(sprintf("Error: %s doesn't exist", pos_matrixFile))
}

min_among_sites <- as.numeric(c_args[3])
max_among_sites <- as.numeric(c_args[4])

visual_orderFile <- ""
if (length(c_args) > 4) {
  visual_orderFile <- c_args[5]
}


#
# read
# type distRankDesc pos ... (name of each strain)
#
d_each_target_site <- read.table(pos_matrixFile,header=T)

#
# set min and max (in the same scale)
#
abs_min_among_sites <- abs(min_among_sites)

if (max_among_sites < abs_min_among_sites) {
  max_among_sites <- abs_min_among_sites
  min_among_sites <- min_among_sites
} else {
  max_among_sites <- max_among_sites
  min_among_sites <- -1*max_among_sites
}

#
# plot heatmaps
#
outprefix <- sub(".txt","", pos_matrixFile)

pos          <- d_each_target_site[1,]$pos
type         <- d_each_target_site[1,]$type
distRankDesc <- d_each_target_site[1,]$distRankDesc

d_each_for_plot <- d_each_target_site[,-c(1:3)]

names(d_each_for_plot) <- sub("^X","",names(d_each_for_plot))

if (visual_orderFile != "") {
  c_visual_order <- c()

  d_visual_order <- read.table(visual_orderFile, header=F)
  for (each in d_visual_order[,1]) {
    c_visual_order <- c(c_visual_order, which(names(d_each_for_plot) == each))
  }
  
  d_each_for_plot <- d_each_for_plot[c_visual_order,c_visual_order]
}

rownames(d_each_for_plot) <- names(d_each_for_plot)

png_fpath <- sprintf("%s_%s.png", outprefix, type) # "rank" is included in the outprefix

png(png_fpath)
print(png_fpath)
plotHeatmap(as.matrix(d_each_for_plot), cex.axis=0.7, min_colscale=min_among_sites, max_colscale=max_among_sites, ordering_flag=T)
dev.off()

