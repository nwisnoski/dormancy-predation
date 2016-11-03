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
timesteps = 10000


# set initial densities
A0 = .1
D0 = 0
P0 = .1
R0 = 1

# set habitat quality
Q = .07              # mean resource renewal rate
a = 0*Q                # amplitude of fluctuation in resource inflow
w = .001             # periodicity of flucutations
l = 0.00         # loss rate of resource in patch


# set vital rates
c = .5         # resource consumption rate of the prey
dorm.max = .5    # maximum dormancy rate /step
react.max = .4    # maximum reactivation rate /step
e.r = .8          # conversion rate on resource
e.a = .4            # conversion rate on active prey
e.d = .0            # conversion rate on dormant prey (0 = inedible)
f.a = .6           # feeding rate of predators on active prey
f.d = .0        # feeding rate of predators on dormant prey
m.d = 0.001        # death rate of dormant prey
m.a = 0.01           # death rate of active prey
m.p = 0.05          # predator death rate
a.ii = .05       # strength of intraspecific competition



extinct.thresh = 0.00001

# define model
PPdorm.energetic <- function(in.matrix = "", timesteps = "", dormancy = "", stochastic = T){
  if(dormancy == F){
    dorm.max = 0
    react.max = 0
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
    
    dorm <- if(stochastic) runif(1, max = dorm.max*exp(-R)) else dorm.max*exp(-R)
    react <- if(stochastic) runif(1, max = react.max*(1-exp(-R))) else react.max*(1-exp(-R))
    
    dR     <- Q + a*sin(w*i) - c*R*A - l*R
    dA     <- e.r*c*R*A*(1-A*a.ii) - dorm*A + react*D - f.a*A*P - if(stochastic) runif(1,max=m.a)*A else m.a*A
    dD     <- dorm*A - react*D - f.d*D*P - if(stochastic) runif(1,max=m.d)*D else m.d*D
    dP     <- e.a*f.a*A*P + e.d*f.d*D*P - if(stochastic) runif(1,max=m.p)*P else m.p*P
    
    in.matrix[i+1, 1] <- i
    in.matrix[i+1, 2] <- max(((A + dA) > extinct.thresh)*(A+dA), 0) # return 0 if pop below thresh.
    in.matrix[i+1, 3] <- max((D + dD), 0)
    in.matrix[i+1, 4] <- max(((P + dP) > extinct.thresh)*(P+dP), 0)
    in.matrix[i+1, 5] <- max((R + dR), 0)
  }
  
  return(in.matrix)
}

# Initialize time dynamics
time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 5)
colnames(time.dynamics) <- c("t", "A", "D", "P", "R")
time.dynamics[1, ] <- c(1, A0, D0, P0, R0)

# Run the model with and without dormancy
out.dynamics.D <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T, stochastic = F)
out.dynamics.NoD <- PPdorm.energetic(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F, stochastic = F)


# Plot the temporal dynamics
png("./figures/PPdorm_one-patch_Dynamics.png", width = 1200, height = 1000, res = 2*96)
par(mfrow = c(2,1))
par(mar = c(2,5,3,3))
dorm.plot <- plot(out.dynamics.D[,1], out.dynamics.D[,2], 
                  ylim = c(0,max(out.dynamics.NoD[,2:4], out.dynamics.D[,2:4])),
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)

points(out.dynamics.D[,1], out.dynamics.D[,3],
       type = "l", lty = "solid", cex = 0.5, col = "green")

points(out.dynamics.D[,1], out.dynamics.D[,4],
       type = "l", lty = "solid", cex = 0.5, col = "red")

axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(With Dormancy)", line = 2.5, cex = 1.2)

legend("topright", c("Active", "Dormant", "Predators"),
       lty = c("solid", "solid", "solid"),
       col = c("black", "green", "red"), cex = 1, bty = "n")

# Plot Without Dormancy
par(mar = c(5,5,0,3))
no.dorm.plot <- plot(out.dynamics.NoD[,1], out.dynamics.NoD[,2], 
                     ylim = c(0,max(out.dynamics.NoD[,2:4], out.dynamics.D[,2:4])),
                     xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", lty = "solid", cex = 0.5)

points(out.dynamics.NoD[,1], out.dynamics.NoD[,3],
       type = "l", lty = "solid", cex = 0.5, col = "green")

points(out.dynamics.NoD[,1], out.dynamics.NoD[,4],
       type = "l", lty = "solid", cex = 0.5, col = "red")

axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(No Dormancy)", line = 2.5, cex = 1.2)
mtext(side = 1, "Timestep", line = 3, cex = 1.5)


dev.off()
graphics.off()
img <- png::readPNG("./figures/PPdorm_one-patch_Dynamics.png")
grid.raster(img)

mean(out.dynamics.D[500:timesteps,2])/sd(out.dynamics.D[500:timesteps,2])
mean(out.dynamics.D[500:timesteps,4])/sd(out.dynamics.D[500:timesteps,4])
mean(out.dynamics.NoD[500:timesteps,2])/sd(out.dynamics.NoD[500:timesteps,2])
mean(out.dynamics.NoD[500:timesteps,4])/sd(out.dynamics.NoD[500:timesteps,4])


