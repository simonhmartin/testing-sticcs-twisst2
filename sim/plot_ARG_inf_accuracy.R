


cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#4363d8', #azure blue
'#3cb44b', #grass green
'#f58231', #chalk orange
'#e6194b') #pomergranate red

methods <- c("sticcs_sc", "tsinfer", "relate", "singer", "argweaver")

method_names <- c("sticcs", "tsinfer", "Relate", "Singer", "Argweaver")

names(cols) <- methods

pch = 21:25
names(pch)<-methods

################################################### accuracy and speed #############################################################

params <- read.csv("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/params.csv", as.is=T)

inference_data <- read.csv("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/arg_inference_data.csv", as.is=T)

treatments <- list()
treatments[["infinite_sites"]] <- inference_data[1:5,]
treatments[["finite_sites"]] <- inference_data[6:10,]
treatments[["geno_error"]] <- inference_data[11:15,]
treatments[["pol_error"]] <- inference_data[16:20,]
treatments[["low_mu/rho"]] <- inference_data[21:25,]
treatments[["gene_conversion"]] <- inference_data[26:30,]
treatments[["pop_struct"]] <- inference_data[31:35,]

treatments <- c("finite_sites", "infinite_sites", "geno_error", "pol_error", "low_mu_to_rho", "gene_conversion", "pop_struct")

rows <- list(1:5, 6:10, 11:15, 16:20, 21:25, 26:30, 31:35)
names(rows) <- treatments

r = rows[[treatment]]

for (treatment in treatments){
    
    svg(paste0("arg_inference.", treatment, ".svg"), width=8, height=4)
    
    par(mfrow = c(1, 3), mar = c(4,6,1,1))
    
    r = rows[[treatment]]
    
    plot(0, cex = 0, xlim = c(4, 64), ylim = c(0, 1), xlab = "Number of haplotypes", ylab = "Accuracy (1 - quartet distance)", main = "Accuracy", log="x", xaxt="n", bty="l")
    abline(h=seq(0,1,0.1), col="gray90", lty=1)
    abline(v=params[r,"n"], col="gray90", lty=1)
    for (method in methods){
        lines(params[r,"n"], 1-inference_data[r, paste0(method, "_dist")], type="l", col=cols[method])
        points(params[r,"n"], 1-inference_data[r, paste0(method, "_dist")], type="p",, pch=pch[method], col=cols[method], bg=cols[method])
        }
    axis(1, at = params[r,"n"])
    legend("bottom", lty=1, pch=pch, col=cols[methods], pt.bg=cols[methods] , legend=method_names)
    
    
    plot(1, cex = 1, xlim = c(4, 64), ylim = c(2, 10), xlab = "Number of haplotypes", ylab = "Children per node (mean)", main = "Polytomies", xaxt="n", log="x", bty="l")
    abline(h=seq(2,10,1), col="gray90", lty=1)
    abline(v=params[r,"n"], col="gray90", lty=1)
    for (method in methods){
        lines(params[r,"n"], inference_data[r, paste0(method, "_mean_children")], type="l", col=cols[method])
        points(params[r,"n"], inference_data[r, paste0(method, "_mean_children")], type="p",, pch=pch[method], col=cols[method], bg=cols[method])
        }
    axis(1, at = params[r,"n"])
    
    plot(0, cex = 0, xlim = c(4, 64), ylim = c(1, 10000), xlab = "Number of haplotypes", ylab = "Numer of trees", main = "Number of trees", xaxt="n", log="xy", bty="l")
    abline(h=10**(0:3), col="gray90", lty=1)
    abline(v=params[r,"n"], col="gray90", lty=1)
    lines(params[r,"n"], inference_data[r, "sim_trees"], type="l", col="black")
    points(params[r,"n"], inference_data[r, "sim_trees"], type="p",, pch=8, col="black")
    for (method in methods){
        lines(params[r,"n"], inference_data[r, paste0(method, "_trees")], type="l", col=cols[method])
        points(params[r,"n"], inference_data[r, paste0(method, "_trees")], type="p",, pch=pch[method], col=cols[method], bg=cols[method])
        }
    axis(1, at = params[r,"n"])
    
    dev.off()
    }


