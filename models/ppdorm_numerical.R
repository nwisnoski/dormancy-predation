rm(list = ls())
library(deSolve)

parameters <- c(Q = .2,
                c = .5,
                d.max = .5,
                r.max = .4,
                e.r = .8,
                e.a = .4,
                e.d = .01,
                f.a = .6,
                f.d = .6,
                m.d = 0.0001,
                m.a = 0.01,
                m.p = 0.05)

state <- c(R = 1,
           A = 1,
           D = 1,
           P = 1)

PPdorm <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)),{
    
    # rate of change
    
    dR <- Q - c*R*A
    dA <- e.r*c*R*A - d.max*A + r.max*D - f.a*A*P - m.a*A
    dD <- d.max*A - r.max*D - f.d*D*P - m.d*D
    dP <- e.a*f.a*A*P + e.d*f.d*D*P - m.p*P
    
    # return rate of change
    list(c(dR, dA, dD, dP))
  })
}

times <- seq(0, 1000, by = 0.001)

out <- ode(y = state, times = times, func = PPdorm, parms = parameters)
head(out)

png(filename = "")
par(oma = c(0,0,3,0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Y"], pch = ".")
mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)
