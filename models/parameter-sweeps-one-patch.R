source("./models/PPdorm-model-one-patch.R")

min.q <- 0.000; max.q <- .5; by.q <- 0.001
param.list <- seq(min.q, max.q, by.q)
param.sweep <- matrix(data = 0, nrow = length(param.list), ncol = 9)
colnames(param.sweep) <- c("Q", "mean.D.A", "sd.D.A", "mean.D.P", "sd.D.P",
                           "mean.NoD.A", "sd.NoD.A", "mean.NoD.P", "sd.NoD.P")
for(i in 1:length(param.list)){
  time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 5)
  colnames(time.dynamics) <- c("t", "A", "D", "P", "R")
  time.dynamics[1, ] <- c(1, A0, D0, P0, R0)
  Q <- param.list[i]
  # Run the model with and without dormancy
  out.dynamics.D <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T)
  out.dynamics.NoD <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F)
  
  param.sweep[i,1] <- Q
  param.sweep[i,2] <- mean(out.dynamics.D[500:timesteps,2])
  param.sweep[i,3] <- sd(out.dynamics.D[500:timesteps,2])
  param.sweep[i,4] <- mean(out.dynamics.D[500:timesteps,4])
  param.sweep[i,5] <- sd(out.dynamics.D[500:timesteps,4])
  param.sweep[i,6] <- mean(out.dynamics.NoD[500:timesteps,2])
  param.sweep[i,7] <- sd(out.dynamics.NoD[500:timesteps,2])
  param.sweep[i,8] <- mean(out.dynamics.NoD[500:timesteps,4])
  param.sweep[i,9] <- sd(out.dynamics.NoD[500:timesteps,4])
  
}

# png("./figures/StabilityPrey.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,2], type = "l", col = "black", lwd = 2,
#      ylim = c(0, max(na.omit(param.sweep[,c(2,4)]))), yaxt = "n",
#      ylab = "Stability, 1/CV", xlab = "Resource Inputs")
# points(param.sweep[,1], param.sweep[,4], type = "l", col = "black", lty = "dashed", lwd = 2)
# axis(side = 1, lwd.ticks = 2)
# axis(side = 2, lwd.ticks = 2, las = 1)
# box(lwd = 2)
# legend(x = "topright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabilityPrey.png"))
# 
# png("./figures/StabilityPred.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,3], type = "l", col = "red", lwd = 2,
#      ylim = c(0, max(na.omit(param.sweep[,c(3,5)]))), yaxt = "n",
#      ylab = "Stability, 1/CV", xlab = "Resource Inputs")
# points(param.sweep[,1], param.sweep[,5], type = "l", col = "red", lty = "dashed", lwd = 2)
# axis(side = 1, lwd.ticks = 2)
# axis(side = 2, lwd.ticks = 2, las = 1)
# box(lwd = 2)
# legend(x = "topright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        col = "red", bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabilityPred.png"))
