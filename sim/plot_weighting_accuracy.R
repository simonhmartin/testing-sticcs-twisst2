

accuracy <- read.csv("admix_n20_weights_accuracy_all.csv", as.is=T)

n = nrow(accuracy)

titles = c("rho/theta = 1", "rho/theta = 10", "rho/theta = 0.1")

################## rooted analysis - excludes NJ50 ####################

svg("admix_model_rooted_weights_accuracy.svg", width=8, height=6)

methods = c("sticcstack", "tsinfer", "relate", "singer", "argweaver")
method_names = c("sticcstack", "tsinfer", "Relate", "Singer", "Argweaver")

scales = c(0,1000,2000,5000,10000,20000)

plot_columns = lapply(methods, function(method) paste0(method, "_rooted_scale", scales)) 

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#4363d8', #azure blue
'#3cb44b', #grass green
'#f58231', #chalk orange
'#e6194b') #pomergranate red

names(cols) <- methods

pch = 21:25
names(pch)<-methods



par(mfrow = c(2, 3), mar = c(4,4,2,1))

for (i in 1:n){
    plot(0, cex = 0, xlim = c(1,length(scales)), ylim = c(0.75, 1), ylab = "", xlab = "", main = "", xaxt="n", bty="l")
    mtext(3, text=titles[i], line=0)
    grid(col="gray90", lty=1)
    for (j in 1:length(plot_columns)){
        lines(1:length(scales), 1-accuracy[i, plot_columns[[j]]], type="l", col=cols[j], lty = 1)
        points(1:length(scales), 1-accuracy[i, plot_columns[[j]]], pch = pch[j], col=cols[j], bg=cols[j], lty = 1)
        }
    axis(1, at = 1:length(scales), labels=scales/1000)
    if (i == 1) mtext(2, text="Weights accuracy (1-error rate)", cex=0.8, line=3)
    }

legend("bottom", legend=method_names, lty = 1, pch = pch, col=cols, pt.bg=cols, cex=0.9)

# dev.off()




################## rooted analysis unphased ####################

methods = c("sticcstack_dip", "sticcstack_tet", "argweaver_dip")
method_names = c("sticcstack diploid unphased", "sticcstack tetraploid unphased", "argweaver diploid unphased")

scales = c(0,1000,2000,5000,10000,20000)

plot_columns = lapply(methods, function(method) paste0(method, "_rooted_scale", scales)) 

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#911eb4', #violet
'#e6194b') #pomergranate red

names(cols) <- methods

pch = c(21,21,25)
names(pch)<-methods

lty = c(2,3,2)

# svg("admix_model_rooted_weights_accuracy.svg", width=8, height=6)

# par(mfrow = c(2, 3), mar = c(4,4,3,1))

for (i in 1:n){
    plot(0, cex = 0, xlim = c(1,length(scales)), ylim = c(0.75, 1), xlab = "", ylab = "", main = "", xaxt="n", bty="l")
    grid(col="gray90", lty=1)
    for (j in 1:length(plot_columns)){
        lines(1:length(scales), 1-accuracy[i, plot_columns[[j]]], type="l", col=cols[j], lty = lty[j])
        points(1:length(scales), 1-accuracy[i, plot_columns[[j]]], pch = pch[j], col=cols[j], bg=cols[j])
        }
    axis(1, at = 1:length(scales), labels=scales/1000)
    if (i == 1) mtext(2, text="Weights accuracy (1-error rate)", cex=0.8, line=3)
    mtext(1, text="Scale of comparison (kb)", cex=0.8, line=3)
    }

legend("bottom", legend=method_names, lty = lty, pch = pch, pt.bg=cols, col=cols, cex=0.9)

dev.off()


######### single plot for phased and unphased

svg("admix_model_weights_accuracy_all.svg", width=8, height=5)

methods = c("sticcstack", "sticcstack_dip", "sticcstack_tet",
            "tsinfer", "relate", "singer", "argweaver", "argweaver_dip")
method_names = c("sticcstack", "sticcstack diploid unphased", "sticcstack tetraploid unphased",
                 "tsinfer", "Relate", "Singer", "Argweaver", "Argweaver diploid unphased")

scales = c(0,1000,2000,5000,10000,20000)

plot_columns = lapply(methods, function(method) paste0(method, "_rooted_scale", scales)) 

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#911eb4', #violet
'#911eb4', #violet
'#4363d8', #azure blue
'#3cb44b', #grass green
'#f58231', #chalk orange
'#e6194b', #pomergranate red
'#e6194b') #pomergranate red

names(cols) <- methods

pch = c(21,21,21,22,23,24,25,25)
names(pch)<-methods

lty = c(1,2,3,1,1,1,1,2)

par(mfrow = c(1, 3), mar = c(4,4,2,1))

for (i in 1:n){
    plot(0, cex = 0, xlim = c(1,length(scales)), ylim = c(0.75, 1), xlab = "", ylab = "", main = "", xaxt="n", bty="l")
    grid(col="gray90", lty=1)
    for (j in 1:length(plot_columns)){
        lines(1:length(scales), 1-accuracy[i, plot_columns[[j]]], type="l", col=cols[j], lty = lty[j])
        points(1:length(scales), 1-accuracy[i, plot_columns[[j]]], pch = pch[j], col=cols[j], bg=cols[j])
        }
    axis(1, at = 1:length(scales), labels=scales/1000)
    if (i == 1) mtext(2, text="Weights accuracy (1-error rate)", cex=0.8, line=3)
    mtext(1, text="Scale of comparison (kb)", cex=0.8, line=3)
    }

legend("bottom", legend=method_names, lty = lty, pch = pch, pt.bg=cols, col=cols, cex=0.9)

dev.off()


######### single plot for phased and unphased, unrooted

svg("admix_model_weights_unrooted_accuracy_all.svg", width=8, height=5)

methods = c("sticcstack", "sticcstack_dip", "sticcstack_tet",
            "tsinfer", "relate", "singer", "argweaver", "argweaver_dip", "NJ50")
method_names = c("sticcstack", "sticcstack diploid unphased", "sticcstack tetraploid unphased",
                 "tsinfer", "Relate", "Singer", "Argweaver", "Argweaver diploid unphased", "NJ 50 SNP windows")

scales = c(0,1000,2000,5000,10000,20000)

plot_columns = lapply(methods, function(method) paste0(method, "_unrooted_scale", scales)) 

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#911eb4', #violet
'#911eb4', #violet
'#4363d8', #azure blue
'#3cb44b', #grass green
'#f58231', #chalk orange
'#e6194b', #pomergranate red
'#e6194b', #pomergranate red
'#9a6324') #coconut bown


names(cols) <- methods

pch = c(21,21,21,22,23,24,25,25,8)
names(pch)<-methods

lty = c(1,2,3,1,1,1,1,2,1)

par(mfrow = c(1, 3), mar = c(4,4,2,1))

for (i in 1:n){
    plot(0, cex = 0, xlim = c(1,length(scales)), ylim = c(0.5, 1), xlab = "", ylab = "", main = "", xaxt="n", bty="l")
    grid(col="gray90", lty=1)
    for (j in 1:length(plot_columns)){
        lines(1:length(scales), 1-accuracy[i, plot_columns[[j]]], type="l", col=cols[j], lty = lty[j])
        points(1:length(scales), 1-accuracy[i, plot_columns[[j]]], pch = pch[j], col=cols[j], bg=cols[j])
        }
    axis(1, at = 1:length(scales), labels=scales/1000)
    if (i == 1) mtext(2, text="Weights accuracy (1-error rate)", cex=0.8, line=3)
    mtext(1, text="Scale of comparison (kb)", cex=0.8, line=3)
    }

legend("bottom", legend=method_names, lty = lty, pch = pch, pt.bg=cols, col=cols, cex=0.9)

dev.off()


################################ runtime #################################

svg("admix_model_rooted_weights_speed.svg", width=3, height=3, bg="white")

par(mfrow = c(1, 1), mar = c(1,4,1,1))

methods = c("sticcstack", "tsinfer", "relate", "singer", "argweaver")

plot_columns = paste0(methods, "_runtime")

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#4363d8', #azure blue
'#3cb44b', #grass green
'#f58231', #chalk orange
'#e6194b') #pomergranate red

names(cols) <- methods

pch = 21:25
names(pch)<-methods


plot(1, cex=0, xlim=c(1, length(methods)), ylim =c(1, 10000), log="y", ylab = "runtime (seconds)", xlab="", xaxt="n")
grid(col="gray90", lty=1)
points(as.numeric(accuracy[1, plot_columns]), type="p", col = cols, bg=cols, pch = pch, cex=2)

dev.off()




svg("admix_model_rooted_weights_unphased_speed.svg", width=3, height=3)

par(mfrow = c(1, 1), mar = c(1,4,1,1))

methods = c("sticcstack_dip", "sticcstack_tet", "argweaver_dip")

plot_columns = paste0(methods, "_runtime")

cols = c( #https://sashamaps.net/docs/resources/20-colors/
'#911eb4', #violet
'#911eb4', #violet
'#e6194b') #pomergranate red

names(cols) <- methods

pch = c(21,21,25)
names(pch)<-methods


plot(1, cex=0, xlim=c(1, length(methods)), ylim =c(1, 10000), log="y", ylab = "runtime (seconds)", xlab="", xaxt="n")
grid(col="gray90", lty=1)
points(as.numeric(accuracy[1, plot_columns]), type="p", col = cols, bg=cols, pch = pch, cex=2)

dev.off()

