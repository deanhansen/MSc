# Author(s): Dean Hansen

source("Main/Scripts/Indicies - Sample.R")
colors <- c("#7A003C","#FDBF57")
library('scales')
B <- 10000

get_summary <- function(data) {
  sample_mean <- apply(data, MARGIN = 2, mean)
  sample_var <- apply(data, MARGIN = 2, var)
  sample_di <- sample_var / sample_mean
  indices <- sample_indices(data)
  output <- list(sample_mean=sample_mean,sample_var=sample_var,sample_di=sample_di,indices=indices)
  return(output)
}

# Load the FBReference Datasets
madrid <- read.csv("Main/Data/FBReference/Real Madrid.csv")
madrid <- madrid[,2:4]
madrid_summary <- get_summary(madrid)
madrid_index <- madrid_summary$indices
madrid_summary

barcelona <- read.csv("Main/Data/FBReference/Barcelona.csv")
barcelona <- barcelona[,2:4]
barcelona_summary <- get_summary(barcelona)
barcelona_index <- barcelona_summary$indices
barcelona_summary


# Bootstrap Sample
madrid_boot <- matrix(NA,nrow=B,ncol=3)
madrid_rows <- nrow(madrid)
set.seed(1)

for (i in 1:B) {
  rows <- sample(madrid_rows,50,replace=TRUE)
  indices <- sample_indices(madrid[rows,])
  madrid_boot[i,] <- c(indices$fi,indices$gdi,indices$mdi)
}

madrid_boot_mean <- apply(madrid_boot,2,mean)
round(madrid_boot_mean,4)
madrid_boot_sd <- apply(madrid_boot,2,sd)
round(madrid_boot_sd,4)
apply(madrid_boot,MARGIN=2,function(x) {quantile(x,c(0.025,0.975))}) |> round(4)

madrid_bias <- c(
  as.numeric(madrid_index$fi - madrid_boot_mean[1]), 
  as.numeric(madrid_index$gdi - madrid_boot_mean[2]),
  as.numeric(madrid_index$mdi - madrid_boot_mean[3])
  )
round(madrid_bias,4)

madrid_mse <- c(
  mean((madrid_index$fi - madrid_boot[,1])^2), 
  mean((madrid_index$gdi - madrid_boot[,2])^2), 
  mean((madrid_index$mdi - madrid_boot[,3])^2)
  )
round(madrid_mse,4)


barcelona_boot <- matrix(NA,ncol=3,nrow=B)
barcelona_rows <- nrow(barcelona)
set.seed(1)

for (i in 1:B) {
  rows <- sample(barcelona_rows,50,replace=TRUE)
  indices <- sample_indices(barcelona[rows,])
  barcelona_boot[i,] <- c(indices$fi,indices$gdi,indices$mdi)
}

barcelona_boot_mean <- apply(barcelona_boot,2,mean)
round(barcelona_boot_mean,4)
barcelona_boot_sd <- apply(barcelona_boot,2,sd)
round(barcelona_boot_sd,4)
apply(barcelona_boot,MARGIN=2,function(x) {quantile(x,c(0.025,0.975))}) |> round(4)

barcelona_bias <- c(
  barcelona_index$fi - barcelona_boot_mean[1], 
  barcelona_index$gdi - barcelona_boot_mean[2],
  barcelona_index$mdi - barcelona_boot_mean[3]
  )
round(barcelona_bias,4)

barcelona_mse <- c(
  mean((barcelona_index$fi - barcelona_boot[,1])^2), 
  mean((barcelona_index$gdi - barcelona_boot[2])^2), 
  mean((barcelona_index$mdi - barcelona_boot[3])^2)
  )
round(barcelona_mse,4)


# Plot Raw Data
break_points <- seq(-0.5, 6.5, by = 1)
xlim <- c(0, 5)
ylim <- c(0,1)

barcelona_1_plot <- hist(barcelona[,1], breaks = break_points, plot = FALSE)
barcelona_2_plot <- hist(barcelona[,2], breaks = break_points, plot = FALSE)
barcelona_3_plot <- hist(barcelona[,3], breaks = break_points, plot = FALSE)
real_madrid_1_plot <- hist(madrid[,1], breaks = break_points, plot = FALSE)
real_madrid_2_plot <- hist(madrid[,2], breaks = break_points, plot = FALSE)
real_madrid_3_plot <- hist(madrid[,3], breaks = break_points, plot = FALSE)

png("main/bbc_msn_goal_histograms.png", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 3), mar=c(4, 4, 3, 3), oma = c(1,1,1,1))
plot(real_madrid_1_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Cristiano Ronaldo",ylab="", main="",freq=F)
plot(real_madrid_2_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Gareth Bale",ylab="", main="",freq=F)
plot(real_madrid_3_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Karim Benzema",ylab="", main="",freq=F)
plot(barcelona_1_plot, col=c2,xlim=xlim, ylim=ylim, xlab="Lionel Messi",ylab="", main="",freq=F)
plot(barcelona_2_plot, col=c2,xlim=xlim, ylim=ylim,xlab="Luis Suarez",ylab="", main="",freq=F)
plot(barcelona_3_plot, col=c2,xlim=xlim, ylim=ylim,xlab="Neymar",ylab="", main="",freq=F)
dev.off()

jpeg("main/bbc_msn_goal_histograms.jpeg", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 3), mar=c(4, 4, 3, 3), oma = c(1,1,1,1))
plot(real_madrid_1_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Cristiano Ronaldo",ylab="", main="",freq=F)
plot(real_madrid_2_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Gareth Bale",ylab="", main="",freq=F)
plot(real_madrid_3_plot, col=c1,xlim=xlim, ylim=ylim,xlab="Karim Benzema",ylab="", main="",freq=F)
plot(barcelona_1_plot, col=c2,xlim=xlim, ylim=ylim, xlab="Lionel Messi",ylab="", main="",freq=F)
plot(barcelona_2_plot, col=c2,xlim=xlim, ylim=ylim,xlab="Luis Suarez",ylab="", main="",freq=F)
plot(barcelona_3_plot, col=c2,xlim=xlim, ylim=ylim,xlab="Neymar",ylab="", main="",freq=F)
dev.off()


# Plot Bootstrap Results
xlim <- c(0.5, 2.5)
ylim <- c(0, 2500)

real_madrid_fi_plot <- hist(real_madrid_boot[,1], plot = FALSE)
real_madrid_gdi_plot <- hist(real_madrid_boot[,2], plot = FALSE)
real_madrid_mdi_plot <- hist(real_madrid_boot[,3], plot = FALSE)
barcelona_fi_plot <- hist(barcelona_boot[,1], plot = FALSE)
barcelona_gdi_plot <- hist(barcelona_boot[,2], plot = FALSE)
barcelona_mdi_plot <- hist(barcelona_boot[,3], plot = FALSE)

png("main/bbc_msn_bootstrap_distribution.png", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(3, 1), mar=c(4, 4, 2, 2), oma = c(1,1,1,1))
plot(real_madrid_fi_plot, col = c1,xlim=xlim, ylim=ylim, xlab=expression(FI[3]), ylab="", main="")
plot(barcelona_fi_plot, col = c2, add = TRUE)
plot(real_madrid_gdi_plot, col=c1, xlim=xlim,ylim=ylim, xlab="GDI", ylab="", main="")
plot(barcelona_gdi_plot, col = c2, add = T)
plot(real_madrid_mdi_plot, col=c1, xlim=xlim,ylim=ylim, xlab="MDI", ylab="", main="")
plot(barcelona_mdi_plot, col = c2, add = T)
dev.off()

jpeg("main/bbc_msn_bootstrap_distribution.jpeg", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(3, 1), mar=c(4, 4, 2, 2), oma = c(1,1,1,1))
plot(real_madrid_fi_plot, col = c1,xlim=xlim, ylim=ylim, xlab=expression(FI[3]), ylab="", main="")
plot(barcelona_fi_plot, col =c2, add = TRUE)
plot(real_madrid_gdi_plot, col=c1, xlim=xlim,ylim=ylim, xlab="GDI", ylab="", main="")
plot(barcelona_gdi_plot, col = c2, add = T)
plot(real_madrid_mdi_plot, col=c1, xlim=xlim,ylim=ylim, xlab="MDI", ylab="", main="")
plot(barcelona_mdi_plot, col = c2, add = T)
dev.off()
