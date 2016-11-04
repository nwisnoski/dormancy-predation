source("./models/PPdorm-model-one-patch.R")

param <- "Q"
min.par <- 0.001; max.par <- 1; by.par <- 0.01
param.list <- seq(min.par, max.par, by.par)
param2 <- "f.d"
param2.list <- seq(min.par, max.par, by.par)
param.sweep <- matrix(data = 0, nrow = length(param.list) * length(param2.list), 
                      ncol = 12)
colnames(param.sweep) <- c(param, "mean.D.A", "sd.D.A", "mean.D.D", "sd.D.D", "mean.D.P", "sd.D.P",
                           "mean.NoD.A", "sd.NoD.A", "mean.NoD.P", "sd.NoD.P", param2)
k <- 1
for(i in 1:length(param.list)){
  
  for(j in 1:length(param2.list)){
    

    time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 5)
    colnames(time.dynamics) <- c("t", "A", "D", "P", "R")
    time.dynamics[1, ] <- c(1, A0, D0, P0, R0)
 
    assign(param, param.list[i], envir = .GlobalEnv)  # change selected parameter value
    assign(param2, param2.list[j], envir = .GlobalEnv)


    # Run the model with and without dormancy
    out.dynamics.D <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T, stochastic = F)
    out.dynamics.NoD <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F, stochastic = F)

    param.sweep[k,1] <- eval(as.name(param))  # write the value of the parameter
    param.sweep[k,2] <- mean(out.dynamics.D[500:timesteps,2])  # Dorm: Active
    param.sweep[k,3] <- sd(out.dynamics.D[500:timesteps,2])
    param.sweep[k,4] <- mean(out.dynamics.D[500:timesteps,3])  # Dorm: Dorm
    param.sweep[k,5] <- sd(out.dynamics.D[500:timesteps,3])
    param.sweep[k,6] <- mean(out.dynamics.D[500:timesteps,4])  # Dorm: Pred
    param.sweep[k,7] <- sd(out.dynamics.D[500:timesteps,4])
    param.sweep[k,8] <- mean(out.dynamics.NoD[500:timesteps,2])  # No Dorm: Active
    param.sweep[k,9] <- sd(out.dynamics.NoD[500:timesteps,2])
    param.sweep[k,10] <- mean(out.dynamics.NoD[500:timesteps,4]) # No Dorm: Pred
    param.sweep[k,11] <- sd(out.dynamics.NoD[500:timesteps,4])
    param.sweep[k,12] <- eval(as.name(param2))
    
    k <- k + 1

  }
  
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
grid::grid.raster(png::readPNG(paste("./figures/OnePatchEquilDensities_",param,".png",sep = "")))


# png("./figures/StabilityPrey.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,2]/param.sweep[,3], type = "l", col = "black", lwd = 2,
#      ylim = c(0, 15), yaxt = "n", xaxt = "n",
#      ylab = "Stability, 1/CV", xlab = "Resource Inputs", cex.lab = 1.5)
# points(param.sweep[,1], param.sweep[,8]/param.sweep[,9], type = "l", col = "black", lty = "dashed", lwd = 2)
# axis(side = 1, lwd.ticks = 2, cex = 1.2)
# axis(side = 2, lwd.ticks = 2, las = 1, cex = 1.2)
# box(lwd = 2)
# legend(x = "topright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabilityPrey.png"))
# 
# png("./figures/StabilityPred.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,6]/param.sweep[,7], type = "l", col = "red", lwd = 2,
#      ylim = c(0, 20), yaxt = "n", xaxt = "n",
#      ylab = "Stability, 1/CV", xlab = "Resource Inputs", cex.lab = 1.5)
# points(param.sweep[,1], param.sweep[,10]/param.sweep[,11], type = "l", col = "red", lty = "dashed", lwd = 2)
# axis(side = 1, lwd.ticks = 2, cex.lab = 1.2)
# axis(side = 2, lwd.ticks = 2, las = 1, cex.lab = 1.2)
# box(lwd = 2)
# legend(x = "topright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        col = "red", bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabilityPred.png"))
# 
# 
# png("./figures/StabDiffPrey.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,2]/param.sweep[,3] - param.sweep[,8]/param.sweep[,9], type = "l", col = "black", lwd = 2,
#      yaxt = "n", xaxt = "n", main = "Prey Stability",
#      ylab = "Effect on Stability", xlab = "Resource Inputs", cex.lab = 1.5)
# abline(h = 0, lwd = 2, lty = "dotted")
# axis(side = 1, lwd.ticks = 2, cex = 1.2)
# axis(side = 2, lwd.ticks = 2, las = 1, cex = 1.2)
# box(lwd = 2)
# legend(x = "bottomright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabDiffPrey.png"))
# 
# png("./figures/StabDiffPred.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,6]/param.sweep[,7] - param.sweep[,10]/param.sweep[,11], type = "l", col = "red", lwd = 2,
#      yaxt = "n", xaxt = "n", main = "Predator Stability",
#      ylab = "Effect on Stability", xlab = "Resource Inputs", cex.lab = 1.5)
# abline(h = 0, lwd = 2, lty = "dotted")
# axis(side = 1, lwd.ticks = 2, cex.lab = 1.2)
# axis(side = 2, lwd.ticks = 2, las = 1, cex.lab = 1.2)
# box(lwd = 2)
# legend(x = "bottomright", legend = c("With dormancy", "Without dormancy"), lwd = 2, lty = c("solid", "dashed"),
#        col = "red", bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabDiffPred.png"))
# 
# png("./figures/StabDiff.png", height = 1200, width = 1200, res = 192)
# plot(param.sweep[,1], param.sweep[,2]/param.sweep[,3] - param.sweep[,8]/param.sweep[,9], type = "l", col = "black", lwd = 2,
#      yaxt = "n", xaxt = "n", 
#      ylim = c(
#        min(
#         min(na.omit(param.sweep[,2]/param.sweep[,3] - param.sweep[,8]/param.sweep[,9])),
#         min(na.omit(param.sweep[,6]/param.sweep[,7] - param.sweep[,10]/param.sweep[,11]))
#         ),
#        max(
#          max(na.omit(param.sweep[,2]/param.sweep[,3] - param.sweep[,8]/param.sweep[,9])),
#          max(na.omit(param.sweep[,6]/param.sweep[,7] - param.sweep[,10]/param.sweep[,11]))
#          )
#        ),
#      ylab = "Effect on Stability", xlab = "Resource Inputs", cex.lab = 1.5)
# abline(h = 0, lwd = 2, lty = "dotted")
# points(param.sweep[,1], param.sweep[,6]/param.sweep[,7] - param.sweep[,10]/param.sweep[,11], type = "l", col = "red", lwd = 2)
# axis(side = 1, lwd.ticks = 2, cex = 1.2)
# axis(side = 2, lwd.ticks = 2, las = 1, cex = 1.2)
# box(lwd = 2)
# legend(x = "bottomright", legend = c("Prey", "Predator"), lwd = 2, lty = c("solid", "solid"),
#        col = c("black", "red"), bty = "n")
# dev.off()
# grid::grid.raster(png::readPNG("./figures/StabDiff.png"))




### Stable/unstable parameters (with dormancy)
stab.dorm <- param.sweep[which(param.sweep[,2] > 0 & param.sweep[,4] > 0 & param.sweep[,6] > 0),c(1,12)]
unstab.dorm <- param.sweep[which((param.sweep[,2] == 0 & param.sweep[,4] == 0) | param.sweep[,6] == 0),c(1,12)]
# stab.mat.dorm <- rbind(
#   cbind(stab.dorm, rep(1, nrow(stab.dorm))),
#   cbind(unstab.dorm, rep(0, nrow(unstab.dorm))))
stab.dorm <- cbind(stab.dorm, rep(1, nrow(stab.dorm)))
unstab.dorm <- cbind(unstab.dorm, rep(0, nrow(unstab.dorm)))

### Stable/unstable parameters (without dormancy)
stab.nodorm <- param.sweep[which(param.sweep[,8] > 0 & param.sweep[,10] > 0),c(1,12)]
unstab.nodorm <- param.sweep[which(param.sweep[,8] == 0 | param.sweep[,10] == 0),c(1,12)]
# stab.mat.nodorm <- rbind(
#   cbind(stab.nodorm, rep(1, nrow(stab.nodorm))),
#   cbind(unstab.nodorm, rep(0, nrow(unstab.nodorm))))
stab.nodorm <- cbind(stab.nodorm, rep(1, nrow(stab.nodorm)))
unstab.nodorm <- cbind(unstab.nodorm, rep(0, nrow(unstab.nodorm)))

plot.new()
points(stab.nodorm[,1], stab.nodorm[,2], pch = 21, bg = "black")
points(unstab.nodorm[,1], unstab.nodorm[,2], pch = 21, bg = "white")
box(lwd = 2)
axis(side = 1, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 2, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 3, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 4, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
mtext(paste(param), side = 1, line = 3, cex = 1.2)
mtext(paste(param2), side = 2, line = 3, cex = 1.2)
mtext("No dormancy", side = 3, line = 1.5, cex = 1.5)

plot.new()
points(stab.dorm[,1], stab.dorm[,2], pch = 21, bg = "black")
points(unstab.dorm[,1], unstab.dorm[,2], pch = 21, bg = "white")
box(lwd = 2)
axis(side = 1, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 2, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 3, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 4, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
mtext(paste(param), side = 1, line = 3, cex = 1.2)
mtext(paste(param2), side = 2, line = 3, cex = 1.2)
mtext("Dormancy", side = 3, line = 1.5, cex = 1.5)
