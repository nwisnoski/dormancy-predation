# Predator-prey dynamics when the prey can go dormant
# Organisms are modeled as a fraction of system resources
# Actively growing organisms use growth plus maintenance energy
# Dormant organisms only require maintenance energy
# Predators have higher energy requirements
# Author: Nathan Wisnoski, Indiana University

# initialize environment
rm(list = ls())
require(png)
require(grid)

# initialize model parameters:
timesteps = 1000
A10 = 1
D10 = 0
P10 = 1
A20 = 1
D20 = 0
P20 = 1

gm = .1             # maintenance energy req./unit biomass
gg = .05            # growth energy req. / unit biomass
gp = .7             # additional predator energy required / unit biomass
r = 2               # intrinsic growth rate of the prey
delta.max = 0.9     # maximum dormancy rate /step
alpha.max = 0.9     # maximum reactivation rate /step
ea =.7              # conversion efficiency on active prey
ed =.0              # conversion efficiency on dormant prey (0 = inedible)
fa =.4              # feeding rate of predators on active prey
fd = .4             # feeding rate of predators on dormant prey
m.d = 0.001         # death rate of dormant prey
m.a = 0.01           # death rate of active prey
m.p = 0.3           # predator death rate
d.a = 0.0           # active prey dispersal
d.d = 0.0           # dormant prey dispersal
d.p = 0.000           # predator dispersal

# define model
PPdorm.energetic.patch <- function(in.matrix = "", timesteps = "", dormancy = ""){
  if(dormancy == F){
    delta.max = 0
    alpha.max = 0
  }
  
  A1 <- in.matrix[1,2]
  D1 <- in.matrix[1,3]
  P1 <- in.matrix[1,4]
  A2 <- in.matrix[1,5]
  D2 <- in.matrix[1,6]
  P2 <- in.matrix[1,7]
  
  for(i in 1:(timesteps-1)){
    A1 <- as.numeric(in.matrix[i,2])
    D1 <- as.numeric(in.matrix[i,3])
    P1 <- as.numeric(in.matrix[i,4])
    A2 <- as.numeric(in.matrix[i,5])
    D2 <- as.numeric(in.matrix[i,6])
    P2 <- as.numeric(in.matrix[i,7])
    R1 <- 1 - (A1*(gg+gm) + D1*(gm) + P1*(gg+gm+gp))
    R2 <- 1 - (A2*(gg+gm) + D2*(gm) + P2*(gg+gm+gp))
    
    delta1 <- runif(1, max = delta.max)
    alpha1 <- runif(1, max = alpha.max)
    delta2 <- runif(1, max = delta.max)
    alpha2 <- runif(1, max = alpha.max)
    disp.A <- runif(1, min = -d.a, max = d.a)
    disp.D <- runif(1, min = -d.d, max = d.d)
    disp.P <- runif(1, min = -d.p, max = d.p)
    
    dA1     <- r*R1*A1 - delta1*(1-R1)*A1* + alpha1*D1 - fa*A1*P1 - runif(1,max=m.a)*A1 + disp.A*(A1+A2)
    dD1     <- delta1*(1-R1)*A1 - alpha1*D1 - fd*D1*P1 - runif(1,max=m.d)*D1 + disp.D*(D1+D2)
    dP1     <- ea*fa*A1*P1 + ed*fd*D1*P1 - runif(1,max=m.p)*P1 + disp.P*(P1+P2)
    dA2     <- r*R2*A2 - delta2*(2-R2)*A2* + alpha2*D2 - fa*A2*P2 - runif(2,max=m.a)*A2 - disp.A*(A1+A2)
    dD2     <- delta2*(2-R2)*A2 - alpha2*D2 - fd*D2*P2 - runif(2,max=m.d)*D2 - disp.D*(D1+D2)
    dP2     <- ea*fa*A2*P2 + ed*fd*D2*P2 - runif(2,max=m.p)*P2 - disp.P*(P1+P2)
    
    in.matrix[i+1, 1] <- i
    in.matrix[i+1, 2] <- max((A1 + dA1), 0)
    in.matrix[i+1, 3] <- max((D1 + dD1), 0)
    in.matrix[i+1, 4] <- max((P1 + dP1), 0)
    in.matrix[i+1, 5] <- max((A2 + dA2), 0)
    in.matrix[i+1, 6] <- max((D2 + dD2), 0)
    in.matrix[i+1, 7] <- max((P2 + dP2), 0)
  }
  
  return(in.matrix)
}

# Initialize time dynamics
time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 7)
colnames(time.dynamics) <- c("t", "A1", "D1", "P1",
                             "A2", "D2", "P2")
time.dynamics[1, ] <- c(1, A10, D10, P10, A20, D20, P20)

# Run the model with and without dormancy
out.dynamics.D <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T)
out.dynamics.NoD <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F)


# Plot the temporal dynamics
png("./figures/EnergyMod_two-patch_Dynamics.png", width = 1200, height = 1000, res = 2*96)
par(mfcol = c(2,2))
par(mar = c(2,5,3,1))

dorm.plot <- plot(out.dynamics.D[,1], out.dynamics.D[,2], 
                  ylim = c(0,max(out.dynamics.NoD[,2:4], out.dynamics.D[,2:4])),
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,3],
       type = "l", lty = "dashed", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,4],
       type = "l", lty = "longdash", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(With Dormancy)", line = 2.5, cex = 1.2)
mtext(side = 3, "Patch 1", line = 1, cex = 1.5)

par(mar = c(5,5,0,1))
no.dorm.plot <- plot(out.dynamics.NoD[,1], out.dynamics.NoD[,2], 
                     ylim = c(0,max(out.dynamics.NoD[,2:7], out.dynamics.D[,2:7])),
                     xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,3],
       type = "l", lty = "dashed", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,4],
       type = "l", lty = "longdash", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(No Dormancy)", line = 2.5, cex = 1.2)
mtext(side = 1, "Timestep", line = 3, cex = 1.5)

par(mar = c(2,1,3,5))
dorm.plot <- plot(out.dynamics.D[,1], out.dynamics.D[,5], 
                  ylim = c(0,max(out.dynamics.NoD[,2:7], out.dynamics.D[,2:7])),
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,6],
       type = "l", lty = "dashed", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,7],
       type = "l", lty = "longdash", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 3, "Patch 2", line = 1, cex = 1.5)

par(mar = c(5,1,0,5))
no.dorm.plot <- plot(out.dynamics.NoD[,1], out.dynamics.NoD[,5], 
                     ylim = c(0,max(out.dynamics.NoD[,2:7], out.dynamics.D[,2:7])),
                     xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,6],
       type = "l", lty = "dashed", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,7],
       type = "l", lty = "longdash", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 1, "Timestep", line = 3, cex = 1.5)

dev.off()
graphics.off()
img <- png::readPNG("./figures/EnergyMod_two-patch_Dynamics.png")
grid.raster(img)
