source("./models/PPdorm-model-one-patch.R")

param <- "Q"
min.par <- 0.001; max.par <- 1; by.par <- 0.01
param.list <- seq(min.par, max.par, by.par)
param.sweep <- matrix(data = 0, nrow = length(param.list), ncol = 11)
colnames(param.sweep) <- c(param, "mean.D.A", "sd.D.A", "mean.D.D", "sd.D.D", "mean.D.P", "sd.D.P",
                           "mean.NoD.A", "sd.NoD.A", "mean.NoD.P", "sd.NoD.P")

for(i in 1:length(param.list)){
  time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 5)
  colnames(time.dynamics) <- c("t", "A", "D", "P", "R")
  time.dynamics[1, ] <- c(1, A0, D0, P0, R0)
  assign(param, param.list[i])  # change selected parameter value
  
  # Run the model with and without dormancy
  out.dynamics.D <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T, stochastic = F)
  out.dynamics.NoD <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F, stochastic = F)
  
  param.sweep[i,1] <- eval(as.name(param))  # write the value of the parameter 
  param.sweep[i,2] <- mean(out.dynamics.D[500:timesteps,2])  # Dorm: Active
  param.sweep[i,3] <- sd(out.dynamics.D[500:timesteps,2])
  param.sweep[i,4] <- mean(out.dynamics.D[500:timesteps,3])  # Dorm: Dorm
  param.sweep[i,5] <- sd(out.dynamics.D[500:timesteps,3])
  param.sweep[i,6] <- mean(out.dynamics.D[500:timesteps,4])  # Dorm: Pred
  param.sweep[i,7] <- sd(out.dynamics.D[500:timesteps,4])
  param.sweep[i,8] <- mean(out.dynamics.NoD[500:timesteps,2])  # No Dorm: Active
  param.sweep[i,9] <- sd(out.dynamics.NoD[500:timesteps,2])
  param.sweep[i,10] <- mean(out.dynamics.NoD[500:timesteps,4]) # No Dorm: Pred
  param.sweep[i,11] <- sd(out.dynamics.NoD[500:timesteps,4])
  
  print(paste("completed ",i," of ",nrow(param.sweep)," iterations.", sep = ""))
}

### Equilibrium Densities
png(paste("./figures/OnePatchEquilDensities_",param,"fd",f.d,
          "ed",e.d,".png",sep = ""), width = 1200, height = 1000, res = 2*96)
par(mfrow = c(2,1))
par(mar = c(2,5,3,3))
plot(param.sweep[,1], param.sweep[,6], type = "l", col = "red", lwd = 2,
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     ylim = c(0, max(param.sweep[,c(6,4,2)])))
points(param.sweep[,1], param.sweep[,4], type = "l", col = "green", lwd = 2)
points(param.sweep[,1], param.sweep[,2], type = "l", col = "blue", lwd = 2)
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Equilibrium Density\n(With Dormancy)", line = 2.75, cex = 1.2)

par(mar = c(5,5,0,3))
plot(param.sweep[,1], param.sweep[,10], type = "l", col = "red", lwd = 2,
     yaxt = "n", xaxt = "n",
     ylab = "", xlab = "", cex.lab = 1.5,
     ylim = c(0, max(param.sweep[,c(10,8)])))
points(param.sweep[,1], param.sweep[,8], type = "l", col = "blue", lwd = 2)
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Equilibrium Density\n(No Dormancy)", line = 2.75, cex = 1.2)
mtext(side = 1, param, line = 3, cex = 1.5)
legend("topright", c("Active", "Dormant", "Predators"),
       lty = c("solid", "solid", "solid"),
       col = c("blue", "green", "red"), cex = 1, bty = "n")

dev.off()
grid::grid.raster(png::readPNG(paste("./figures/OnePatchEquilDensities_",param,"fd",f.d,".png",sep = "")))
