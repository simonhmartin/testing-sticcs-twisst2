source("~/Dropbox/Research/dev/twisst2/plot_twisst/plot_twisst.R")

setwd("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/Heliconius/")


parse_gff_metadata = function(metadata){
    data_split <- strsplit(metadata, "[,=,;]", perl=T)[[1]]
    l <- length(data_split)
    output <- data_split[seq(2,length(data_split), 2)]
    names(output) <- data_split[seq(1,length(data_split), 2)]
    output
    }


prefix = "bar92.chr18.DP8MP4BIMAC2HET75.mpg_ama_txn_chi.chr18"
# prefix = "bar92.chr18.DP8MP4BIMAC2HET75.mpg_mal_flo_chi.chr18"

topocounts_files <- paste0(prefix,".topocounts.tsv.gz")
intervals_files <- paste0(prefix,".intervals.tsv.gz")


twisst_data <- import.twisst(topocounts_files=topocounts_files, intervals_files=intervals_files, ignore_extra_columns=TRUE)

twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 100000, new_positions = seq(2501, 18e6, 5000))


topocols = c(
'#4363d8', #azure blue
'#46f0f0', #cyan
'#bcf60c', #green banana
'#3cb44b', #grass green
'#aaffc3', #light sage
'#9a6324', #coconut bown
'#e6194b', #pomergranate red
'#fabebe', #pale blush
'#e6beff', #mauve
'#f58231', #chalk orange
'#ffe119', #custard yellow
'#f032e6', #hot pink
'#911eb4', #violet
'#008080', #teal
'#000075', #navy blue
'#800000', #red leather
'#808000', #olive
'#ffd8b1', #white skin
'#808080', #grey
'#fffac8') #beach sand


svg(paste0(prefix, ".boxplot.svg"), width=9, height=5, bg=NA)

plot.twisst.summary.boxplot(twisst_data, cols=topocols, order_by_weights=F, lwd=2, trees_below=FALSE, cex=0.8)

dev.off()







anno <- read.delim("Hmel2.5.chromosomes.fullNames.gff3", sep="\t", as.is=T, comment.char="#", header=FALSE)


genes_f <- subset(anno, anno[,1] == "chr18" & anno[,3] == "gene" & anno[,7] == "+")
genes_r <- subset(anno, anno[,1] == "chr18" & anno[,3] == "gene" & anno[,7] == "-")

CDS_f <- subset(anno, anno[,1] == "chr18" & anno[,3] == "CDS" & anno[,7] == "+")
CDS_r <- subset(anno, anno[,1] == "chr18" & anno[,3] == "CDS" & anno[,7] == "-")

genes_f_IDs <- sapply(genes_f[,9], function(x) parse_gff_metadata(x)["ID"])
genes_r_IDs <- sapply(genes_r[,9], function(x) parse_gff_metadata(x)["ID"])

CDS_f_parents <- sapply(CDS_f[,9], function(x) parse_gff_metadata(x)["Parent"], simplify=F)
CDS_r_parents <- sapply(CDS_r[,9], function(x) parse_gff_metadata(x)["Parent"], simplify=F)

key_gene_IDs <- c("HMEL001028", "HMEL013552g1", "HMEL034199g1")
key_transcript_IDs <- c("HMEL001028-RA", "HMEL013552g1.t1", "HMEL034199g1.t1")


svg(paste0(prefix, ".weights_smooth.svg"), width=12, height=8, bg=NA)

par(mfrow=c(3,1), mar=c(4,2,1,1))

plot.weights(twisst_data_smooth$weights[[1]][,1:15], twisst_data_smooth$pos[[1]],lwd=0,stacked=FALSE, xlim=c(0,17e6),
             fill_cols = topocols, line_cols=NA, xlab="", ylab="", xaxt="n", yaxt="n")

# plot.weights(twisst_data$weights[[1]][,1:15], twisst_data$interval_data[[1]][,c("start","end")],lwd=0,stacked=FALSE, xlim=c(0,17e6),
#              fill_cols = topocols, line_cols=NA, xlab="", ylab="", xaxt="n", yaxt="n")


axis(2, line=-2)
mtext(2, text="Weighting", line=.8)

axis(1, at = seq(0, 17e6, 0.5e6), labels = seq(0, 17, 0.5))

rect(0,0,1.5e6,1)



left_trim = 0
right_trim = 1500000

rows <- which(twisst_data$interval_data[[1]]$start >= left_trim & twisst_data$interval_data[[1]]$end <= right_trim)

genes_f_rows <- which(genes_f[,4] >= left_trim & genes_f[,5] <= right_trim)
genes_r_rows <- which(genes_r[,4] >= left_trim & genes_r[,5] <= right_trim)



par(xpd=NA, mar=c(7,2,1,1))

plot.weights(twisst_data$weights[[1]][rows,1:15], twisst_data$interval_data[[1]][rows,c("start","end")],lwd=0,stacked=FALSE, xlim=c(left_trim,right_trim),
             fill_cols = topocols, line_cols=NA, xlab="", ylab="", xaxt="n", yaxt="n")


rect(genes_f[genes_f_rows,4],-0.02,genes_f[genes_f_rows,5],-0.06,border=NA, col=ifelse(genes_f_IDs[genes_f_rows] %in% key_gene_IDs,  '#f58231', "black"))

rect(genes_r[genes_r_rows,4],-0.06,genes_r[genes_r_rows,5],-0.1,border=NA, col=ifelse(genes_r_IDs[genes_r_rows] %in% key_gene_IDs,  '#f58231', "black"))


axis(1, line=2, at = seq(left_trim, right_trim, 100000), labels = seq(left_trim, right_trim, 100000)/1e6)
mtext(1, text="Position on chr18 (Mb)", line=5)

axis(2, line=-2)
mtext(2, text="Weighting", line=.8)


regucalcin_pos <- 476883
optix_pos <- 1059180

text(regucalcin_pos, -0.15, labels="regucalcin", cex=0.8, col='#f58231')
text(optix_pos, -0.15, labels="optix", cex=0.8, col='#f58231')

dev.off()



########################### gene locations #####################
 
genes <- genes_f <- subset(anno, anno[,1] == "chr18" & anno[,3] == "gene")
 
#genes around the regucalcin introgression
genes[genes[,4] > 400000 & genes[,5] < 500000,]

#based on blasting these, I determined that regucalcin matches two genes on the reverse strand:
#"HMEL013552g1", "HMEL034199g1"
#These are two parts of the same gene, so I just approximated their location in the gff as the midpoint
regucalcin_pos <- 476883


#and optix is already known

genes[genes[,4] > 1050000 & genes[,5] < 1100000,]
# which confirmed optix is "HMEL001028"
optix_midpoint <- 1059180

key_gene_IDs <- c("HMEL001028", "HMEL013552g1", "HMEL034199g1")
