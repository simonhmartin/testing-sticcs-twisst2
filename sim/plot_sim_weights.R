source("~/Dropbox/Research/dev/twisst2/plot_twisst/plot_twisst.R")

dir <- "/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/"
setwd(dir)


# twisst_data <- import.twisst(topocounts_files="admix3_n8_123_l2e5_r1e8.true.topocounts.tsv",
#                              intervals_files="admix3_n8_123_l2e5_r1e8.true.intervals.tsv", ignore_extra_columns=FALSE)
# 
# weights_order = order(twisst_data$weights_mean, decreasing=T)
# 
# topocols = topo_cols[order(weights_order)]
# 
# plot.twisst(twisst_data, mode=3, include_region_names=TRUE, tree_x_lim = c(-3,8), rel_height=2, xlab="", margins=c(2,4,1,1), cols = topocols)


##################################################### ALL METHODS COMPARISON ############################################

sim_name <- "admix_n20_seed123_l200000.0_r1e-08"
subtrees = 512

sim_name <- "admix3_n20_seed6_l200000.0_r1e-08"
subtrees = 8000

true_topocounts_file <- paste0(dir, sim_name, ".true.topocounts.tsv")
true_intervals_file <- paste0(dir, sim_name, ".true.intervals.tsv")

mut_name <- "_mu1e-08"

topocounts_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.topocounts.tsv",
                                                      ".relate.topocounts.tsv",
                                                      ".singer.topocounts.tsv",
                                                      ".argweaver.topocounts.tsv",
                                                    paste0(".sticcs_ploidy1_subtre", subtrees, ".topocounts.tsv"),
                                                    paste0(".sticcs_ploidy2_subtre", subtrees, ".topocounts.tsv"),
                                                    paste0(".sticcs_ploidy4_subtre", subtrees, ".topocounts.tsv")))

intervals_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.intervals.tsv",
                                                     ".relate.intervals.tsv",
                                                     ".singer.intervals.tsv",
                                                     ".argweaver.intervals.tsv",
                                                        paste0(".sticcs_ploidy1_subtre", subtrees, ".intervals.tsv"),
                                                        paste0(".sticcs_ploidy2_subtre", subtrees, ".intervals.tsv"),
                                                        paste0(".sticcs_ploidy4_subtre", subtrees, ".intervals.tsv")))

topocounts_files <- c(true_topocounts_file, topocounts_files)
intervals_files <- c(true_intervals_file, intervals_files)

twisst_data <- import.twisst(topocounts_files=topocounts_files,
                             intervals_files=intervals_files, ignore_extra_columns=FALSE,
                             names = c("Simulated truth", "tsinfer (phased haplotypes)", "Relate (phased haplotypes)", "Singer  (phased haplotypes)", "Argweaver (phased haplotypes)", "sticcstack (phased haplotypes)", "sticcstack (unphased diploid)", "sticcstack (unphased tetraploid)"))

ntopos = length(twisst_data$topos)

weights_order = order(twisst_data$weights_mean[1,-(ntopos + 1)], decreasing=T)

topocols = palettes[["sashamaps"]][order(weights_order)]

# svg(paste0(dir, sim_name, mut_name, ".weights.svg"), width=5, height=8)
# plot.twisst(twisst_data, show_topos=FALSE, mode=5, include_region_names=TRUE, tree_x_lim = c(-3,8),
#             ncol_topos = ifelse(length(twisst_data$topos) == 3, 3, 5), rel_height=2, xlab="", margins=c(2,4,1,1), cols=topocols)
# dev.off()
# 
# twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 10000, new_positions = seq(1001, 200000, 1000))
# 
# svg(paste0(dir, sim_name, mut_name, ".weights_smooth20_trans.svg"), width=5, height=8)
# plot.twisst(twisst_data_smooth, show_topos=FALSE, mode=2, include_region_names=TRUE, tree_x_lim = c(-3,8),
#             ncol_topos = ifelse(length(twisst_data_smooth$topos) == 3, 3, 5), rel_height=2, margins=c(2,4,1,1), xlab="", cols=topocols)
# dev.off()


############################### line plots #######################################

focal_topo <- 1
focal_topo <- 5

single_line_cols <- single_gray_cols <- rep(NA, ntopos)

single_gray_cols[focal_topo] <- "gray80"
single_line_cols[focal_topo] <- topocols[focal_topo]


ntests = twisst_data$n_regions

svg(paste0(dir, sim_name, mut_name, ".weights_comparison.svg"), width=12, height=12)

par(mfrow = c(ntests, 2), mar = c(3.5,3.5,1.5,1))

for (i in 1:ntests){
    plot(0, cex=0, xlim = c(0, twisst_data$lengths[i]), ylim = c(0,1), yaxt="n", ylab="", xlab="", bty="l")
    axis(2, at = c(0,0.5,1))
    mtext(2, text="Topology weighting", cex=0.8, line=2.5)
    mtext(3, text=names(twisst_data$weights)[i], adj=0)
    if (i == ntests) mtext(1, text="Chromosome position", cex=0.8, line=2.5)

    plot.weights(twisst_data$weights[[i]], twisst_data$interval_data[[i]][,c(2,3)], add=TRUE, line_cols=topocols,fill_cols=topocols, stack=T)

    plot(0, cex=0, xlim = c(0, twisst_data$lengths[i]), ylim = c(0,1), yaxt="n", ylab="", xlab="", bty="l")
    axis(2, at = c(0,0.5,1))
    mtext(3, text=names(twisst_data$weights)[i], adj=0)

    plot.weights(twisst_data$weights[[1]], twisst_data$interval_data[[1]][,c(2,3)], add=TRUE, line_cols=single_gray_cols, lwd=1)
    plot.weights(twisst_data$weights[[i]], twisst_data$interval_data[[i]][,c(2,3)], add=TRUE, line_cols=single_line_cols)
    if (i == ntests) mtext(1, text="Chromosome position", cex=0.8, line=2.5)
    }

dev.off()


################################ weights only for concept figure #################################

sim_name <- "admix3_n20_6_l2e5_r1e8"

mut_name <- "_mu1e8_er0_mis0"

topocounts_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.topocounts.tsv",
                                                      ".sticcs_ploidy2_subtre8000.topocounts.tsv",
                                                      paste0(".sticcs_ploidy2_subtre8.topocounts.",1:3,".tsv")))

intervals_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.intervals.tsv",
                                                      ".sticcs_ploidy2_subtre8000.intervals.tsv",
                                                      paste0(".sticcs_ploidy2_subtre8.intervals.",1:3,".tsv")))


twisst_data <- import.twisst(topocounts_files=topocounts_files, intervals_files=intervals_files, ignore_extra_columns=TRUE,
                             names = c("tsinfer (phased haplotypes)", "sticcstack (unphased diploid)", "sticcstack itr1", "sticcstack itr2", "sticcstack itr3"))


ntopos = length(twisst_data$topos)

weights_order = order(twisst_data$weights_mean[1,-(ntopos + 1)], decreasing=T)

topocols = palettes[["sashamaps"]][order(weights_order)]

#tsinfer
svg(paste0(dir, sim_name, mut_name, ".tsinfer.weights_only.svg"), width=5, height=2)
par(mar=c(0,0,0,0))
plot.weights(twisst_data$weights[[1]], twisst_data$interval_data[[1]][,c("start","end")], stacked=T, fill_cols = topocols, line_cols = NA, yaxt="n", xaxt="n", ylab="", xlab="", xlim = c(0, 50000))
dev.off()

#sticcstack
svg(paste0(dir, sim_name, mut_name, ".sticcstack.weights_only.svg"), width=5, height=2)
par(mar=c(0,0,0,0))
plot.weights(twisst_data$weights[[2]], twisst_data$interval_data[[2]][,c("start","end")], stacked=T, fill_cols = topocols, line_cols = NA, yaxt="n", xaxt="n", ylab="", xlab="", xlim = c(0, 50000))
dev.off()

#sticcstack itr 1
svg(paste0(dir, sim_name, mut_name, ".sticcstack_1.weights_only.svg"), width=5, height=0.5)
par(mar=c(0,0,0,0))
plot.weights(twisst_data$weights[[3]], twisst_data$interval_data[[3]][,c("start","end")], stacked=T, fill_cols = topocols, line_cols = NA, yaxt="n", xaxt="n", ylab="", xlab="", xlim = c(0, 50000))
dev.off()

#sticcstack itr 2
svg(paste0(dir, sim_name, mut_name, ".sticcstack_2.weights_only.svg"), width=5, height=0.5)
par(mar=c(0,0,0,0))
plot.weights(twisst_data$weights[[4]], twisst_data$interval_data[[4]][,c("start","end")], stacked=T, fill_cols = topocols, line_cols = NA, yaxt="n", xaxt="n", ylab="", xlab="", xlim = c(0, 50000))
dev.off()

svg(paste0(dir, sim_name, mut_name, ".sticcstack_3.weights_only.svg"), width=5, height=0.5)
par(mar=c(0,0,0,0))
plot.weights(twisst_data$weights[[5]], twisst_data$interval_data[[5]][,c("start","end")], stacked=T, fill_cols = topocols, line_cols = NA, yaxt="n", xaxt="n", ylab="", xlab="", xlim = c(0, 50000))
dev.off()


#trees
svg(paste0(dir, sim_name, ".topos.svg"), width=5, height=1)
par(mar=c(1,0,0,0), xpd=NA)
plot(0,cex=0,xlim = c(0.2, 3), xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1), bty="n")

#now run the draw.tree function for each topology. You can set x_scale and y_scale to alter the tree width and height.
for (i in 1:3){
    draw.tree(twisst_data$topos[[i]], x=i-1, y=0, x_scale=0.3, y_scale=0.3,
            col=topocols[i], label_offset=0.05, cex=1, lwd=2, direction="down",
            label_alias=c(group1="grp1", group2="grp2", group3="grp3"))
    }
dev.off()

####################  comparison of averaged 5k  blocks #################################

sim_name <- "admix_n8_123_l2e5_r1e8"
sim_name <- "admix3_n8_123_l2e5_r1e8"

true_topocounts_file <- paste0(dir, sim_name, ".true.ave5k.topocounts.tsv")
true_intervals_file <- paste0(dir, sim_name, ".true.ave5k.intervals.tsv")

mut_name <- "_mu1e8_er0_mis0"


topocounts_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.ave5k.topocounts.tsv",
                                                    ".sticcs_ploidy1.ave5k.topocounts.tsv",
                                                    ".sticcs_ploidy2.ave5k.topocounts.tsv",
                                                    ".sticcs_ploidy4.ave5k.topocounts.tsv"))

intervals_files <- paste0(dir, sim_name, mut_name, c(".tsinfer.ave5k.intervals.tsv",
                                                        ".sticcs_ploidy1.ave5k.intervals.tsv",
                                                        ".sticcs_ploidy2.ave5k.intervals.tsv",
                                                        ".sticcs_ploidy4.ave5k.intervals.tsv"))

topocounts_files <- c(true_topocounts_file, topocounts_files)
intervals_files <- c(true_intervals_file, intervals_files)

twisst_data <- import.twisst(topocounts_files=topocounts_files, intervals_files=intervals_files, ignore_extra_columns=FALSE,
                             names = c("Simulated", "tsinfer", "Twisst2 (phased)", "Twisst2 (unphased diploid)", "Twisst2 (unphased tetraploid)"))

weights_order = order(twisst_data$weights_mean[1,], decreasing=T)

topocols = topo_cols[order(weights_order)]

svg(paste0(dir, sim_name, mut_name, ".weights.svg"), width=12, height=7)
plot.twisst(twisst_data, mode=3, include_region_names=TRUE, tree_x_lim = c(-3,8), rel_height=2, xlab="", margins=c(2,4,1,1), cols=topocols)
dev.off()


#######################################################################  distances based on smoothed weights (no longer using this for the paper ) ############################################

#double scaled euclidean distance - http://www.pbarrett.net/techpapers/euclid.pdf
#set A_is_target to FALSE to calculate maximum distance as distance from valuesA
#If independent is FALSE, will scale by n-1 rather than by n.
ds_euc_dist <- function(valuesA, valuesB, max_values, min_values, A_is_target = FALSE, independent=TRUE){
    #calculate maximum distance for each variable
    if (A_is_target == FALSE) md <- (max_values - min_values)**2 
    #if A is target, use maximum distance from A for md
    else md <- (apply(abs(cbind(valuesA-max_values,valuesA-min_values)), 1, max))**2 
    
    if (independent == TRUE) sqrt(sum(((valuesA - valuesB)**2)/md) / length(valuesA))
    else sqrt(sum(((valuesA - valuesB)**2)/md) /  (length(valuesA)-1))
    }

nfiles <- length(topocounts_files)

spans <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000)

# spans <- c(5000, 10000, 20000, 50000, 100000)

dists_by_span <- array(dim = c(length(spans), nfiles-1))

for (x in 1:length(spans)){
    print(spans[x])
    twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = spans[x], new_positions = seq(1001, 500000, 1000))
#     plot.twisst(twisst_data_smooth, mode=2) #mode 2 overlays polygons, mode 3 would stack them
    
    for (y in 2:nfiles){
        dist <- mean(sapply(1:nrow(twisst_data_smooth$weights[[1]]),
                            function(i){ds_euc_dist(twisst_data_smooth$weights[[1]][i,],
                                                    twisst_data_smooth$weights[[y]][i,],
                                                    rep(1,3), rep(0,3), A_is_target=T, independent=F)}), na.rm=T)

        dists_by_span[x, y-1] <- dist
        }
    }




svg(paste0(dir, sim_name, mut_name, ".weights_accuracy.svg"), width=6, height=6)
plot(0, cex=0, xlim= c(1000, 100000), ylim=c(.6,1), log="x", xlab = "Resolution (loess span)", ylab="Accuarcy (1 - error rate)")

grid(col="gray90", lty=1)

for (y in 1:(nfiles-1)){
    lines(spans, 1-dists_by_span[,y], type="b", lty=y, pch=y)
    }

legend("bottom", legend=c("tsinfer + Twisst (phased)", "Twigs (phased)", "Twigs (unphased diploid)", "Twisst2 (unphased diploid)", "Twisst2 (unphased tetraploid)"),
       lty = 1:(nfiles-1), pch = 1:(nfiles-1))

dev.off()


