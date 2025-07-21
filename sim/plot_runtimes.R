


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

params <- read.csv("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/timetest_params.csv", as.is=T)

runtime_data <- read.csv("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/times_data.csv", as.is=T)


treatments <- list(rho_theta_1=c(mu=1e-8, rec=1e-8), c(mu=1e-8, rec=1e-9), rho_theta_10=c(mu=1e-9, rec=1e-8))

names(treatments) = c("rho=0.01 theta=0.01", "rho=0.001 theta=0.01", "rho=0.01 theta=0.001")


lens = c(1e5, 1e6, 1e7)
names(lens) = c("len=0.1Mb", "len=1Mb", "len=10Mb")


svg(paste0("runtimes.svg"), width=10, height=10)

par(mfrow = c(3, 3), mar = c(6,6,1,1))

    
for (treatment in names(treatments)){
    
    for (len in names(lens)){

        rows = which(params$mutation_rate == treatments[[treatment]]["mu"] & params$recombination_rate == treatments[[treatment]]["rec"] & params$seq_len == lens[len]) 
        
        plot(1, cex = 1, xlim = c(4, 64), ylim = c(0.01, 10000), xlab = "", ylab = "", main = paste(treatment, len), xaxt="n", log="xy", las=2)
        abline(h=10**(-2:4), col="gray90", lty=1)
        abline(v=params[rows,"n"], col="gray90", lty=1)
        
        for (method in methods){
            lines(params[rows,"n"], runtime_data[rows, method], type="l", col=cols[method])
            points(params[rows,"n"], runtime_data[rows, method], type="p",, pch=pch[method], col=cols[method], bg=cols[method])
            }
        axis(1, at = params[rows,"n"])
        mtext(side=2, text="Runtime (seconds)", line=4, cex=0.8)
        mtext(side=1, text="Number of haplotypes", line=2, cex=0.8)
        
        if (treatment == names(treatments)[1] & len == names(lens)[1]) legend("bottomright", lty=1, pch=pch, col=cols[methods], pt.bg=cols[methods], legend=method_names)
        
        }
    }

dev.off()


##################################################### complexity #######################################################################



estimate_complexity <- function(x, runtime) {
  # Create log-log model
  model <- lm(log(runtime) ~ log(x))
  
  # Extract the slope (complexity exponent)
  complexity <- coef(model)[2]
  
  # Get confidence interval for the slope
  ci <- confint(model)[2, ]
  
  # Print results
  cat("Estimated Complexity: O(n^", round(complexity, 3), ")\n", sep = "")
  cat("95% Confidence Interval: [", round(ci[1], 3), ", ", round(ci[2], 3), "]\n", sep = "")
  cat("R-squared: ", round(summary(model)$r.squared, 4), "\n")
  
  return(list(
    complexity = complexity,
    confidence_interval = ci,
    model = model,
    r_squared = summary(model)$r.squared
  ))
}


analyze_complexity <- function(x, runtime) {
  # Estimate complexity
  result <- estimate_complexity(x, runtime)
  
  # Create log-log plot
  plot(log(x), log(runtime), 
       xlab = "log(Sample Size)", 
       ylab = "log(Runtime)",
       main = paste("Complexity Analysis: O(n^", round(result$complexity, 2), ")", sep = ""),
       pch = 19, col = "blue")
  
  # Add regression line
  abline(result$model, col = "red", lwd = 2)
  
  # Add grid
  grid()
  
  # Add legend with R-squared
  legend("topleft", 
         legend = paste("R² =", round(result$r_squared, 4)),
         bty = "n")
  
  return(result)
}


# Common complexity interpretations:
interpret_complexity <- function(exponent) {
  if (exponent < 1.1) {
    return("Linear or sub-linear: O(n) or O(log n)")
  } else if (exponent < 1.5) {
    return("Between linear and quadratic: O(n log n) likely")
  } else if (exponent < 2.1) {
    return("Quadratic: O(n²)")
  } else if (exponent < 3.1) {
    return("Cubic: O(n³)")
  } else {
    return("Higher order polynomial or exponential")
  }
}


treatment = "rho=0.01 theta=0.01"

len = "len=0.1Mb"

rows = which(params$mutation_rate == treatments[[treatment]]["mu"] & params$recombination_rate == treatments[[treatment]]["rec"] & params$seq_len == lens[len]) 


x <- params[rows,"n"]
runtime <- runtime_data[rows,"sticcs_sc"]

analyze_complexity(x, runtime)

interpret_complexity(1.726)


