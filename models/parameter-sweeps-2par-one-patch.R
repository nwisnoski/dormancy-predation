source("./models/PPdorm-model-one-patch.R")

param <- "dorm.max"
min.par <- 0.001; max.par <- 1; by.par <- 0.01
param.list <- seq(min.par, max.par, by.par)
param2 <- "react.max"
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
    
    print(paste("completed ",k," of ",nrow(param.sweep)," iterations.", sep = ""))
    k <- k + 1

  }
  
}



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


png(paste("./figures/StabPlot_nodorm_",param,eval(as.name(param)),
          param2,eval(as.name(param2)),".png",sep = ""), 
    width = 1200, height = 1000, res = 192)
plot.new()
points(stab.nodorm[,1], stab.nodorm[,2], pch = 21, bg = "black", cex = 0.5)
#points(unstab.nodorm[,1], unstab.nodorm[,2], pch = 21, bg = "white", cex = 0.5)
box(lwd = 2)
axis(side = 1, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 2, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 3, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 4, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
mtext(paste(param), side = 1, line = 3, cex = 1.2)
mtext(paste(param2), side = 2, line = 3, cex = 1.2)
mtext("No dormancy", side = 3, line = 1.5, cex = 1.5)
dev.off()
graphics.off()
grid::grid.raster(
  png::readPNG(paste("./figures/StabPlot_nodorm_",param,eval(as.name(param)),
                     param2,eval(as.name(param2)),".png",sep = "")))

png(paste("./figures/StabPlot_dorm_",param,eval(as.name(param)),
          param2,eval(as.name(param2)),".png",sep = ""), 
    width = 1200, height = 1000, res = 192)
plot.new()
points(stab.dorm[,1], stab.dorm[,2], pch = 21, bg = "black", cex = 0.5)
#points(unstab.dorm[,1], unstab.dorm[,2], pch = 21, bg = "white", cex = 0.5)
box(lwd = 2)
axis(side = 1, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 2, labels = T, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 3, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
axis(side = 4, labels = F, las = 1, lwd.ticks = 2, cex.axis = 1)
mtext(paste(param), side = 1, line = 3, cex = 1.2)
mtext(paste(param2), side = 2, line = 3, cex = 1.2)
mtext("Dormancy", side = 3, line = 1.5, cex = 1.5)
dev.off()
graphics.off()
grid::grid.raster(
  png::readPNG(paste("./figures/StabPlot_dorm_",param,eval(as.name(param)),
                     param2,eval(as.name(param2)),".png",sep = "")))

