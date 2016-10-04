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
timesteps = 100000
s = 100

A0 = 1
D0 = 0
P0 = 1
R0 = 1

#gm = .000001             # maintenance energy req./unit biomass
#gg = .000005            # growth energy req. / unit biomass
#gp = .00007             # additional predator energy required / unit biomass
Q = 2 /s                 # resource renewal rate
r = 1.2 /s           # intrinsic growth rate of the prey
delta.max = 1 /s    # maximum dormancy rate /step
alpha.max = 1 /s    # maximum reactivation rate /step
ea =.5 /s             # conversion efficiency on active prey
ed =.0 /s             # conversion efficiency on dormant prey (0 = inedible)
fa =.5 /s             # feeding rate of predators on active prey
fd = .0  /s           # feeding rate of predators on dormant prey
m.d = 0.0001  /s       # death rate of dormant prey
m.a = 0.05   /s        # death rate of active prey
m.p = 0.05   /s        # predator death rate

# define model
PPdorm.energetic <- function(in.matrix = "", timesteps = "", dormancy = ""){
  if(dormancy == F){
    delta.max = 0
    alpha.max = 0
  }
  
  A <- in.matrix[1,2]
  D <- in.matrix[1,3]
  P <- in.matrix[1,4]
  R <- in.matrix[1,5]
  
  for(i in 1:(timesteps-1)){
    A <- in.matrix[i,2]
    D <- in.matrix[i,3]
    P <- in.matrix[i,4]
    R <- in.matrix[i,5]
    
    delta <- runif(1, max = delta.max)
    alpha <- runif(1, max = alpha.max)
    
    dR     <- Q*R - r*R*A
    dA     <- r*R*A - delta*exp(-R)*A* + alpha*D*R - fa*A*P - runif(1,max=m.a)*A
    dD     <- delta*exp(-R)*A - alpha*D*R - fd*D*P - runif(1,max=m.d)*D
    dP     <- ea*fa*A*P + ed*fd*D*P - runif(1,max=m.p)*P
    
    in.matrix[i+1, 1] <- i
    in.matrix[i+1, 2] <- max((A + dA), 0)
    in.matrix[i+1, 3] <- max((D + dD), 0)
    in.matrix[i+1, 4] <- max((P + dP), 0)
    in.matrix[i+1, 5] <- max((R + dR), 0)
  }
  
  return(in.matrix)
}

# Initialize time dynamics
time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 5)
colnames(time.dynamics) <- c("t", "A", "D", "P", "R")
time.dynamics[1, ] <- c(1, A0, D0, P0, R0)

# Run the model with and without dormancy
out.dynamics.D <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T)
out.dynamics.NoD <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F)


# Plot the temporal dynamics
png("./figures/EnergyMod_one-patch_Dynamics.png", width = 1200, height = 1000, res = 2*96)
par(mfrow = c(2,1))
par(mar = c(2,5,3,3))
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

legend("topright", c("Active", "Dormant", "Predators"),
       lty = c("solid", "dashed", "longdash"),
       col = c("black", "black", "red"), cex = 1, bty = "n")

# Plot Without Dormancy
par(mar = c(5,5,0,3))
no.dorm.plot <- plot(out.dynamics.NoD[,1], out.dynamics.NoD[,2], 
                     ylim = c(0,max(out.dynamics.NoD[,2:4], out.dynamics.D[,2:4])),
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


dev.off()
graphics.off()
img <- png::readPNG("./figures/EnergyMod_one-patch_Dynamics.png")
grid.raster(img)
