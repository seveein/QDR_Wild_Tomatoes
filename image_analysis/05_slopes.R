# Load libraries
rm(list=ls())
library(ggplot2)
library(tidyverse)
library(tidyr)
library(magrittr)
library(segmented)
library(cowplot)
library(progress)

# Define function to compute slope
compute_slope <- function(d, N, myfile) {
  rout <- NULL
  deja <- NULL
  
  # Check if output file exists
  if (file.exists(myfile)) {
    deja <- read.csv(myfile, sep='\t', header=F)
  }
  
  for (uy in levels(as.factor(d$UId))) {
    mpass <- FALSE
    # Check if UId is already processed
    if (length(deja$V1) != 0) {
      if (uy %in% deja$V1) {
        mpass <- TRUE
      }
    }
    # If UId is not processed
    if (mpass == FALSE) {
      sd <- subset(d, UId == uy)
      min <- 0
      mm <- max(sd$S) * 0.5
      mm <- (max(sd$S) - min(sd$S)) * 0.5
      ts <- sd$t[which.max(abs(sd$S))]
      ssd <- subset(sd, t < ts)
      max <- ssd$t[which.min(abs(ssd$S - mm))]
      span <- 50
      mfit <- with(sd, ksmooth(t, S, kernel = "normal", bandwidth = span))
      
      # Plot data and fitted curve
      p <- ggplot() + geom_point(data = sd, aes(x = t, y = S), alpha = 0.3) +
        geom_line(aes(x = sd$t, y = mfit$y), color = 'black') +
        theme_bw() + ggtitle(uy) +
        geom_vline(xintercept = max, color = 'red')
      print(paste0("red line is at ", max))
      print(p)
      sd$Sb <- as.numeric(mfit$y)
      
      # Get user input for limit
      minmax <- readline(paste0("Enter limit [range,[s]=skip,[a]=automatic,[f]=fail] "))
      if (minmax != "f") {
        if (minmax != "s") {
          if (minmax != "a") {
            min <- as.numeric(unlist(strsplit(minmax, ","))[[1]])
            max <- as.numeric(unlist(strsplit(minmax, ","))[[2]])
          }
          # Subset data based on limit
          sd <- subset(sd, UId == uy & t > min & t < max)
          print(N)
          N <- N + 1
          x <- sd$t
          x <- x * 10 # Optional
          y <- sd$Sb
          ylog <- log(sd$Sb + 1)
          ylog[which(ylog == -Inf)] <- NA
          
          # Fit segmented regression
          out.lm <- lm(ylog ~ x)
          o <- segmented(out.lm, seg.Z = ~ x)
          fit <- numeric(length(x)) * NA
          fit[complete.cases(rowSums(cbind(ylog, x)))] <- broken.line(o)$fit
          val <- sd(o$residuals)
          ylog[which(o$residuals > val | o$residuals < -val)] <- NA
          o$coefficients
          data1 <- data.frame(x = x, y = y, ylog = ylog, fit = fit)
          res <- slope(o)
          print(paste0("breakpoint is at ", o$psi[[2]], " seconds"))
          
          # Plot segmented regression
          p1 <- ggplot() + geom_point(data = sd, aes(x = t * 10, y = log(Sb + 1)), alpha = 0.3) +
            geom_line(data = data1, aes(x = x, y = ylog), color = "black") +
            geom_line(aes(x = x, y = fit), color = "steelblue") + theme_bw() +
            geom_vline(xintercept = o$psi[[2]], color = "red") +
            ggtitle(paste(uy, res$x[2, 1], sep = " "))
          p2 <- ggplot() + geom_line(data = data1, aes(x = x, y = y, alpha = 0.3)) +
            geom_point(data = sd, aes(x = t * 10, y = Sb), alpha = 0.3) +
            theme_bw() + geom_vline(xintercept = o$psi[[2]], color = "red")
          
          # Combine plots
          p2
          p <- plot_grid(p1, p2, ncol = 2, nrow = 1, align = "hv")
          print(p)
          
          # Get user input for response
          answer <- readline(paste0("Is it ok: \n \t please >> enter O,N : \n Enter Your Response HERE:   "))
          rout <- rbind(rout, c(uy, o$psi[[2]], res$x[2, 1], sd$mod[1], answer))
          line <- paste(uy, o$psi[[2]], res$x[2, 1], answer, log(((log(800) - summary(o)$coefficients[1, 1]) / slope(o)$x[2, 1] - (log(400) - summary(o)$coefficients[1, 1]) / slope(o)$x[2, 1])), sep = "\t")
          write(line, file = myfile, append = TRUE)
        }
      }
      if (minmax == "s") {
        line <- paste(uy, NA, NA, minmax, sep = "\t")
        write(line, file = myfile, append = TRUE)
      }
      if (minmax == "f") {
        line <- paste(uy, NA, NA, minmax, sep = "\t")
        write(line, file = myfile, append = TRUE)
      }
    }
  }
  return(rout)
}

# Define datasets
d <- read.csv("counts1.txt", header = TRUE, sep = "\t")

# Change column names
colnames(d) <- c('UId', 't', 'S', '')

# Subset data
d <- subset(d, t > 30)

# Define output file
myfile <- 'slope.txt'

compute_slope(d, 1, myfile)
