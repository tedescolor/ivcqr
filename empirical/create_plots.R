path = getwd() #set path to the folder empirical/
setwd(path)
test = FALSE
fullestimation = read.csv("bootstrap_results.csv")
covNames = c("intercept","treatment","age")
fullestimation$cov = rep(covNames,nrow(fullestimation)/3)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
save = TRUE
plot_conf = cbind(covNames[2:3], c(-2,-0.15), c(0.5,0.02))
for(i in 1:nrow(plot_conf)){
  cov = plot_conf[i,1]
  png(filename=paste(cov,".png",sep = ""),width = 800, height = 800)
  layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2), nrow = 5, ncol = 3, byrow = TRUE))
  y <- cbind(fullestimation[fullestimation$cov == cov,"Estimation"])
  # z = cbind(
  #   fullestimation[, c(cov)] + 1.96 * fullestimation[, c(paste("sd_", cov, sep = ""))],
  #   fullestimation[, c(cov)] - 1.96 * fullestimation[, c(paste("sd_", cov, sep = ""))]
  # )
  z = cbind(
    fullestimation[fullestimation$cov == cov,"CI025"],
    fullestimation[fullestimation$cov == cov,"CI975"]
  )
  ylims = c(as.numeric(plot_conf[i,2]),as.numeric(plot_conf[i,3]))
  
  par(mar=c(6,6, 6,6)+2, mgp=c(6.5,3,0))
  plot(1, type="n",ylab="",  xlab="quantile",xlim=c(0.1,0.6), ylim=ylims, 
       main = cov,cex.lab=4, cex.axis = 4, cex.main=4,yaxt="n", xaxs="i")
  
  
  polygon(c(fullestimation[fullestimation$cov == cov,"u"], 
            rev(fullestimation[fullestimation$cov == cov,"u"])), c(z[,1], rev(z[,2])),
          col = "#DADADA",border = NA)
  lines(fullestimation[fullestimation$cov == cov,"u"],y[,1],col="red", lty=1, lwd=4)
  #lines(fullestimation[fullestimation$cov == cov,"u"],rep(0,length(fullestimation[fullestimation$cov == cov,"u"])),col="blue", lty=1, lwd=4)
  #lines(fullestimation$u,y[,2],col="blue", lty=3, lwd=4)
  axis(2,cex.axis = 4)
  plot.new()
  if(cov == "treatment"){
    legend(x = "bottom",
           inset = c(0,0), legend=c("Estimate","Bootstrap CI"),xjust = 1,
           col=  c("red","#DADADA"),
           fill = c(NA,"#DADADA"), bty = "n",
           lty=c(1,NA), border = c(NA,NA),cex=4,x.intersp = c(0.25,-1.4),
           trace = TRUE ,xpd =TRUE, horiz=TRUE, lwd=4, text.width = c(0.25, 0.35))
    
  }# 
  if(cov == "age"){
    legend(x = "bottom",
           inset = c(0, 0), legend=c("Estimate","Bootstrap CI"),xjust = 1,
           col=  c("red","#DADADA"),
           fill = c(NA,"#DADADA"), bty = "n",
           lty=c(1,NA), border = c(NA,NA),cex=4,x.intersp = c(0.25,-1.4), trace = TRUE ,xpd =TRUE, horiz=TRUE, lwd=4, text.width = c(0.25, 0.35))
    
  }
  #   
  # readline(prompt="Press [enter] to continue")
  dev.off()
}


library(dplyr)

# Create a custom function to calculate the color of each bar
get_bar_colors <- function(times, deltas, breaks) {
  # Create a data frame with the times and deltas
  aux <- data.frame(times = times, deltas = deltas)
  
  # Calculate the range for each bin
  bin_ranges <- cut(aux$times, breaks = breaks, right = FALSE, include.lowest = TRUE)
  
  # Calculate the proportion of censored data in each bin
  bin_censored <- aux %>%
    mutate(bin = bin_ranges) %>%
    group_by(bin) %>%
    summarize(censored_proportion = mean(deltas)) %>%
    pull(censored_proportion)
  
  # Calculate the color for each bin based on the proportion of censored data
  bar_colors = sapply(bin_censored, function(value)gray.colors(1, start = value))
  #bar_colors <- gray.colors(breaks - 1, start = 1, end = 0)[match(bin_ranges, unique(bin_ranges))]
  
  return(bar_colors)
}

# Define the number of breaks (bins) for the histogram
breaks <-15

data = read.csv("data.csv")
# Calculate the bar colors
bar_colors <- get_bar_colors(data$days, data$delta, breaks)
png(filename="histogram.png",width = 800, height = 600)
# Create the histogram with the custom bar colors
hist(data$days, breaks = breaks, col = bar_colors, xlab = "Days", main = "",cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
dev.off()
