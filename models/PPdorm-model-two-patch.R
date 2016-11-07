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


# Starting populations
A10 = .1
D10 = 0
P10 = .1
R10 = 1
A20 = .1
D20 = 0
P20 = .1
R20 = 1

# Set patch quality
Q1 = .2              # mean resource renewal rate patch 1
Q2 = .2              # "" patch 2
a1 = 0.*Q1                # amplitude of fluctuation in resource inflow patch 1
a2 = 0.*Q2             # amplitude patch 2
w1 = 0.001             # periodicity of flucutations in patch 1
w2 = 0.001           # periodicity in patch 2
l1 = 0.00         # loss rate of resource in patch 1
l2 = 0.00         # loss rate of resource in patch 2

# set prey and predator vital rates
c = .5         # resource consumption rate of the prey
dorm.max = .5    # maximum dormancy rate /step
react.max = .4    # maximum reactivation rate /step
e.r = .8          # conversion rate on resource
e.a = .4            # conversion rate on active prey
e.d = .01            # conversion rate on dormant prey (0 = inedible)
f.a = .6           # feeding rate of predators on active prey
f.d = .6        # feeding rate of predators on dormant prey
m.d = 0.0001        # death rate of dormant prey
m.a = 0.01           # death rate of active prey
m.p = 0.05          # predator death rate
a.ii = 0.05       # strength of intraspecific competition
d.a = 0.     # active dispersal
d.d = 0.     # dormant dispersal
d.p = 0.     # predator dispersal


extinct.thresh = 0.00001

# define model
PPdorm.energetic.patch <- function(in.matrix = "", timesteps = "", dormancy = "", stochastic = T){
  if(dormancy == F){
    dorm.max = 0
    react.max = 0
  }
  
  A1 <- in.matrix[1,2]
  D1 <- in.matrix[1,3]
  P1 <- in.matrix[1,4]
  R1 <- in.matrix[1,5]
  A2 <- in.matrix[1,6]
  D2 <- in.matrix[1,7]
  P2 <- in.matrix[1,8]
  R2 <- in.matrix[1,9]
  
  for(i in 1:(timesteps-1)){
    A1 <- as.numeric(in.matrix[i,2])
    D1 <- as.numeric(in.matrix[i,3])
    P1 <- as.numeric(in.matrix[i,4])
    R1 <- as.numeric(in.matrix[i,5])
    A2 <- as.numeric(in.matrix[i,6])
    D2 <- as.numeric(in.matrix[i,7])
    P2 <- as.numeric(in.matrix[i,8])
    R2 <- as.numeric(in.matrix[i,9])
   
    dorm1 <- if(stochastic) runif(1, max = dorm.max*exp(-R1)) else dorm.max*exp(-R1)
    react1 <- if(stochastic) runif(1, max = react.max*(1-exp(-R1))) else react.max*(1-exp(-R1))
    dorm2 <- if(stochastic) runif(1, max = dorm.max*exp(-R2)) else dorm.max*exp(-R2)
    react2 <- if(stochastic) runif(1, max = react.max*(1-exp(-R2))) else react.max*(1-exp(-R2))
    disp.A <- if(stochastic) runif(1, min = -d.a, max = d.a) else d.a
    disp.D <- if(stochastic) runif(1, min = -d.d, max = d.d) else d.d
    disp.P <- if(stochastic) runif(1, min = -d.p, max = d.p) else d.p
    
    dR1     <- Q1 + a1*sin(w1*i) - c*R1*A1 - l1*R1
    dA1     <- e.r*c*R1*A1*(1-A1*a.ii) - dorm1*A1 + react1*D1 - f.a*A1*P1 + 
                  (if(stochastic) disp.A*(A1+A2) else disp.A*(A1-A2)) -
                  (if(stochastic) runif(1,max=m.a)*A1 else m.a*A1)
    dD1     <- dorm1*A1 - react1*D1 - f.d*D1*P1 +
                  (if(stochastic) disp.D*(D1+D2) else disp.D*(D1-D2)) -
                  (if(stochastic) runif(1,max=m.d)*D1 else m.d*D1)
    dP1     <- e.a*f.a*A1*P1 + e.d*f.d*D1*P1 +
                  (if(stochastic) disp.P*(P1+P2) else disp.P*(P1-P2)) -
                  (if(stochastic) runif(1,max=m.p)*P1 else m.p*P1)
    dR2     <- Q2 + a2*sin(w2*i) - c*R2*A2 - l2*R2
    dA2     <- e.r*c*R2*A2*(1-A2*a.ii) - dorm2*A2 + react2*D2 - f.a*A2*P2 -
                  (if(stochastic) disp.A*(A1+A2) else disp.A*(A1-A2)) -
                  (if(stochastic) runif(1,max=m.a)*A2 else m.a*A2)
    dD2     <- dorm2*A2 - react2*D2 - f.d*D2*P2 -
                  (if(stochastic) disp.D*(D1+D2) else disp.D*(D1-D2)) -
                  (if(stochastic) runif(1,max=m.d)*D2 else m.d*D2)
    dP2     <- e.a*f.a*A2*P2 + e.d*f.d*D2*P2 -
                  (if(stochastic) disp.P*(P1+P2) else disp.P*(P1-P2)) -
                  (if(stochastic) runif(1,max=m.p)*P2 else m.p*P2)
    
    in.matrix[i+1, 1] <- i
    in.matrix[i+1, 2] <- max(((A1 + dA1) > extinct.thresh)*(A1+dA1), 0)
    in.matrix[i+1, 3] <- max((D1 + dD1), 0)
    in.matrix[i+1, 4] <- max(((P1 + dP1) > extinct.thresh)*(P1+dP1), 0)
    in.matrix[i+1, 5] <- max((R1 + dR1), 0)
    in.matrix[i+1, 6] <- max(((A2 + dA2) > extinct.thresh)*(A2+dA2), 0)
    in.matrix[i+1, 7] <- max((D2 + dD2), 0)
    in.matrix[i+1, 8] <- max(((P2 + dP2) > extinct.thresh)*(P2+dP2), 0)
    in.matrix[i+1, 9] <- max((R2 + dR2), 0)
  }
  
  return(in.matrix)
}

# Initialize time dynamics
time.dynamics <- matrix(data = NA, nrow = timesteps, ncol = 9)
colnames(time.dynamics) <- c("t", "A1", "D1", "P1", "R1",
                             "A2", "D2", "P2", "R2")
time.dynamics[1, ] <- c(1, A10, D10, P10, R10, A20, D20, P20, R20)

# Run the model with and without dormancy
out.dynamics.D <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = T, stochastic = F)
out.dynamics.NoD <- PPdorm.energetic.patch(in.matrix = time.dynamics, timesteps = timesteps, dormancy = F, stochastic = F)


# Plot the temporal dynamics
png("./figures/PPdorm_two-patch_Dynamics.png", width = 1200, height = 1000, res = 2*96)
par(mfcol = c(2,2))
par(mar = c(2,5,3,1))

plot(out.dynamics.D[,1], out.dynamics.D[,2], 
                  ylim = c(0,max(out.dynamics.NoD[,c(2:4,6:8)], out.dynamics.D[,c(2:4,6:8)])),
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", 
                  col = "blue", lty = "solid", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,3],
       type = "l", col = "green", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,4],
       type = "l", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(With Dormancy)", line = 2.5, cex = 1.2)
mtext(side = 3, "Patch 1", line = 1, cex = 1.5)

par(mar = c(5,5,0,1))
plot(out.dynamics.NoD[,1], out.dynamics.NoD[,2], 
                     ylim = c(0,max(out.dynamics.NoD[,c(2:4,6:8)], out.dynamics.D[,c(2:4,6:8)])),
                     xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l",
                     col = "blue", lty = "solid", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,3],
       type = "l", col = "green", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,4],
       type = "l", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 2, "Density\n(No Dormancy)", line = 2.5, cex = 1.2)
mtext(side = 1, "Timestep", line = 3, cex = 1.5)

par(mar = c(2,1,3,5))
dorm.plot <- plot(out.dynamics.D[,1], out.dynamics.D[,6], 
                  ylim = c(0,max(out.dynamics.NoD[,c(2:4,6:8)], out.dynamics.D[,c(2:4,6:8)])),
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l", 
                  col = "blue", lty = "solid", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,7],
       type = "l", col = "green", cex = 0.5)
points(out.dynamics.D[,1], out.dynamics.D[,8],
       type = "l", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 3, "Patch 2", line = 1, cex = 1.5)

par(mar = c(5,1,0,5))
no.dorm.plot <- plot(out.dynamics.NoD[,1], out.dynamics.NoD[,6], 
                     ylim = c(0,max(out.dynamics.NoD[,c(2:4,6:8)], out.dynamics.D[,c(2:4,6:8)])),
                     xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "l",
                     col = "blue", lty = "solid", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,7],
       type = "l", col = "green", cex = 0.5)
points(out.dynamics.NoD[,1], out.dynamics.NoD[,8],
       type = "l", cex = 0.5, col = "red")
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, labels = F)
box(lwd = 2)
mtext(side = 1, "Timestep", line = 3, cex = 1.5)

dev.off()
graphics.off()
img <- png::readPNG("./figures/PPdorm_two-patch_Dynamics.png")
grid.raster(img)
