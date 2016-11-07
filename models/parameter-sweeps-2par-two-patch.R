source("./models/PPdorm-model-two-patch.R")

param <- "d.a"
min.par <- 0.001; max.par <- 1; by.par <- 0.01
param.list <- seq(min.par, max.par, by.par)
param2 <- "d.d"
param2.list <- seq(min.par, max.par, by.par)
param.sweep <- matrix(data = 0, nrow = length(param.list) * length(param2.list), 
                      ncol = 22)
colnames(param.sweep) <- c(param, param2, "p1.mean.D.A", "p1.sd.D.A", "p1.mean.D.D", "p1.sd.D.D", "p1.mean.D.P", "p1.sd.D.P",
                           "p1.mean.NoD.A", "p1.sd.NoD.A", "p1.mean.NoD.P", "p1.sd.NoD.P",
                           "p2.mean.D.A", "p2.sd.D.A", "p2.mean.D.D", "p2.sd.D.D", "p2.mean.D.P", "p2.sd.D.P",
                           "p2.mean.NoD.A", "p2.sd.NoD.A", "p2.mean.NoD.P", "p2.sd.NoD.P")
k <- 1
for(i in 1:length(param.list)){
  
  for(j in 1:length(param2.list)){
    
    
    time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 9)
    colnames(time.dynamics) <- c("t", "A1", "D1", "P1", "R1",
                                 "A2", "D2", "P2", "R2")
    time.dynamics[1, ] <- c(1, A10, D10, P10, R10, A20, D20, P20, R20)
    
    
    assign(param, param.list[i], envir = .GlobalEnv)  # change selected parameter value
    assign(param2, param2.list[j], envir = .GlobalEnv)
    
    
    # Run the model with and without dormancy
    out.dynamics.D <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T, stochastic = F)
    out.dynamics.NoD <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F, stochastic = F)
    
    param.sweep[k,1] <- eval(as.name(param))  # write the value of the parameter
    param.sweep[k,2] <- eval(as.name(param2))
    param.sweep[k,3] <- mean(out.dynamics.D[500:timesteps,2])  # Dorm: Active 1
    param.sweep[k,4] <- sd(out.dynamics.D[500:timesteps,2])
    param.sweep[k,5] <- mean(out.dynamics.D[500:timesteps,3])  # Dorm: Dorm 1
    param.sweep[k,6] <- sd(out.dynamics.D[500:timesteps,3])
    param.sweep[k,7] <- mean(out.dynamics.D[500:timesteps,4])  # Dorm: Pred 1
    param.sweep[k,8] <- sd(out.dynamics.D[500:timesteps,4])
    param.sweep[k,9] <- mean(out.dynamics.NoD[500:timesteps,2])  # No Dorm: Active 1
    param.sweep[k,10] <- sd(out.dynamics.NoD[500:timesteps,2])
    param.sweep[k,11] <- mean(out.dynamics.NoD[500:timesteps,4]) # No Dorm: Pred 1 
    param.sweep[k,12] <- sd(out.dynamics.NoD[500:timesteps,4])
    
    param.sweep[k,13] <- mean(out.dynamics.D[500:timesteps,6])  # Dorm: Active 2
    param.sweep[k,14] <- sd(out.dynamics.D[500:timesteps,6])
    param.sweep[k,15] <- mean(out.dynamics.D[500:timesteps,7])  # Dorm: Dorm 2
    param.sweep[k,16] <- sd(out.dynamics.D[500:timesteps,7])
    param.sweep[k,17] <- mean(out.dynamics.D[500:timesteps,8])  # Dorm: Pred 2
    param.sweep[k,18] <- sd(out.dynamics.D[500:timesteps,8])
    param.sweep[k,19] <- mean(out.dynamics.NoD[500:timesteps,6])  # No Dorm: Active 2
    param.sweep[k,20] <- sd(out.dynamics.NoD[500:timesteps,6])
    param.sweep[k,21] <- mean(out.dynamics.NoD[500:timesteps,8]) # No Dorm: Pred 2
    param.sweep[k,22] <- sd(out.dynamics.NoD[500:timesteps,8])
    
    
    print(paste("completed ",k," of ",nrow(param.sweep)," iterations.", sep = ""))
    k <- k + 1
    
  }
  
}


### Stable/unstable parameters (with dormancy)
stab.dorm <- param.sweep[which((param.sweep[,3] > 0 | param.sweep[,5] > 0 | param.sweep[,13] > 0 | param.sweep[,15] > 0) & 
                                 (param.sweep[,7] > 0 | param.sweep[,17] > 0)),c(1,2)]
unstab.dorm <- param.sweep[which((param.sweep[,3] > 0 & param.sweep[,5] > 0 & param.sweep[,13] > 0 & param.sweep[,15] > 0) | 
                                   (param.sweep[,7] > 0 & param.sweep[,17] > 0)),c(1,2)]
# stab.mat.dorm <- rbind(
#   cbind(stab.dorm, rep(1, nrow(stab.dorm))),
#   cbind(unstab.dorm, rep(0, nrow(unstab.dorm))))
stab.dorm <- cbind(stab.dorm, rep(1, nrow(stab.dorm)))
unstab.dorm <- cbind(unstab.dorm, rep(0, nrow(unstab.dorm)))

### Stable/unstable parameters (without dormancy)
stab.nodorm <- param.sweep[which((param.sweep[,9] > 0 | param.sweep[,19] > 0) & 
                                   (param.sweep[,11] > 0 | param.sweep[,21] > 0)),c(1,2)]
unstab.nodorm <- param.sweep[which((param.sweep[,9] > 0 & param.sweep[,19] > 0) | 
                                     (param.sweep[,11] > 0 & param.sweep[,21] > 0)),c(1,2)]
# stab.mat.nodorm <- rbind(
#   cbind(stab.nodorm, rep(1, nrow(stab.nodorm))),
#   cbind(unstab.nodorm, rep(0, nrow(unstab.nodorm))))
stab.nodorm <- cbind(stab.nodorm, rep(1, nrow(stab.nodorm)))
unstab.nodorm <- cbind(unstab.nodorm, rep(0, nrow(unstab.nodorm)))


png(paste("./figures/TwoPatchStabPlot_nodorm_",param,eval(as.name(param)),
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
  png::readPNG(paste("./figures/TwoPatchStabPlot_nodorm_",param,eval(as.name(param)),
                     param2,eval(as.name(param2)),".png",sep = "")))

png(paste("./figures/TwoPatchStabPlot_dorm_",param,eval(as.name(param)),
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
  png::readPNG(paste("./figures/TwoPatchStabPlot_dorm_",param,eval(as.name(param)),
                     param2,eval(as.name(param2)),".png",sep = "")))



