# Dylan Sosa
# January 17 2020
# Theoretical Ecology, Winter 2020 

library(tidyverse)
library(gridExtra)

# Problem 1
exponential <-
  ggplot(data.frame(x = seq(0, 15)), aes(x)) +
  stat_function(fun = exp) +
  labs(title = 'Density Independent Pop Growth')

logistic <-
  ggplot(data.frame(x = seq(-15, 15)), aes(x)) +
  stat_function(fun = plogis) +
  labs(title = 'Density Dependent Pop  Growth')

grid.arrange(exponential, logistic, nrow = 1)

# Exercise 2
logGrowth <- function(x) {
  # r = 2 = u, c = 0.1, v = 5
  # r*N*(1-(N/K)) # logistic growth equation
  2 * x * (1 - (x / 15))
}

logGrowthPerCapita <- function(x) {
  (2 * x * (1 - (x / 15))) / x
}

alleeEffect <- function(x) {
  # r = 2 = u, c = 0.1, v = 5
  # ((u*N/v+N)-c*N)*N # Allee effect equation
  (((2 * x) / (5 + x)) - (0.1 * x)) * x
}

alleeEffectPerCapita <- function(x) {
  (((2 * x) / (5 + x)) - (0.1 * x))
}

plotPop <-
  ggplot(data.frame(x = seq(0, 15, 0.001)), aes(x)) +
  stat_function(fun = logGrowth, aes(color = 'log')) +
  stat_function(fun = alleeEffect, aes(color = 'Allee')) +
  labs(x = 'N', y = 'Population dN/dt', title = 'Logistic & Allee Effect Population Growth Rate') +
  scale_colour_manual("Legend", values = c("forestgreen", "brown")) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

plotPerCap <-
  ggplot(data.frame(x = seq(0, 15, 0.001)), aes(x)) +
  stat_function(fun = logGrowthPerCapita, aes(color = 'log')) +
  stat_function(fun = alleeEffectPerCapita, aes(color = 'Allee'), show.legend = TRUE) +
  labs(x = 'N', y = 'Per Capita dN/dt', title = 'Logistic & Allee Effect Per Capita Growth Rate') +
  scale_colour_manual("Legend", values = c("forestgreen", "brown")) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

grid.arrange(plotPop, plotPerCap, nrow = 1)

# Exercise 6
################
## Bifurcation.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"

## This function returns the values of the min and max 
peaks <- function(x) {
  if (min(x)==max(x)) return(min(x)) ## Does not oscillate 
  l <- length(x) 
  xm1 <- c(x[-1], x[l])
  xp1 <- c(x[1], x[-l]) 
  z<-x[x > xm1 & x > xp1 | x < xm1 & x < xp1] 
  if (length(z)==0) return(min(x)) ## It has not converged yet 
  return (z)
} 

## This function creates a simulation of the logistic map 
LogisticMap<-function(N0,r,TimeSteps){
  Results<-rep(0,TimeSteps) 
  Results[1]<-N0 
  for (j in 2:TimeSteps){
    Results[j]<-r*Results[j-1]*(1-Results[j-1])
  } 
  return(Results)
} 

# mention derivation here, show in markdown
rickerMap<-function(N0,r,TimeSteps){
  Results<-rep(0,TimeSteps) 
  Results[1]<-N0 
  for (j in 2:TimeSteps){
    #Results[j]<-log(1/r)/-0.01
    #Results[j] <- r*Results[j-1]*exp(-0.01*Results[j-1])
    Results[j]<-exp(r*(1-(Results[j-1]/2)))*Results[j-1]
  } 
  return(Results)
}

## Plot the Diagram 
plot(0,0, xlim=c(0,4), ylim=c(-0.05,10),type="n", xlab="r", ylab="X") 
for (r in seq(0.001,100,0.005)) { # These are the initial and final values for r
  # out <- LogisticMap(0.5,r,2500) # Initial conditions 
  out <- rickerMap(0.5,r,2500) # Initial conditions 
  l <- length(out) %/% 10 # use only the last 250 steps 
  out <- out[(9*l):(10*l)] 
  p <- peaks(out) 
  l <- length(out) 
  points(rep(r, length(p)), p, pch=".")
}

# Exercise 7
###############
# From Euler.R
ExponentialGrowth<-function(x,lambda){
  return (lambda*x)
}

Euler<-function(x, D, FUN,lambda){
  return(x+FUN(x,lambda)*D)
}

ExponentialSolution<-function(N0,t,lambda){
  return (N0*exp(t*lambda))
}

N0<-0.25
lambda<-0.3
errorsDF <- data.frame(dt = c(0.05,0.1,0.2), Error = NA) # init df
plot(c(0,0),col="0",xlim=c(0,10),ylim=c(0,(N0*exp(lambda*10)+1)))
for (Deltat in c(0.05,0.1,0.2)){
  PointsToEstimate<-seq(0,10,by=Deltat)
  Iterations<-length(PointsToEstimate)-1
  Approx<-rep(0,Iterations+1)
  RealSol<-rep(0,Iterations+1)
  Approx[1]<-RealSol[1]<-N0
  for (i in 1:Iterations){
    RealSol[i+1]<-ExponentialSolution(N0,Deltat*i,lambda)		
    Approx[i+1]<-Euler(Approx[i],Deltat,ExponentialGrowth,lambda)		
  }
  if (Deltat==0.05) {
    mycol="red"
    errorsDF[1,2] <- abs(RealSol[length(RealSol)] - Approx[length(Approx)])
    # last element in approx/real vectors is at time = 10
    # update data frame with error at t = 10 for each Deltat
  }
  if (Deltat==0.1) {
    mycol="blue"
    errorsDF[2,2] <- abs(RealSol[length(RealSol)] - Approx[length(Approx)])
  }
  if (Deltat==0.2)  {
    mycol="green"
    errorsDF[3,2] <- abs(RealSol[length(RealSol)] - Approx[length(Approx)])
  }
  
  points(RealSol~PointsToEstimate,col="black",type="l")
  points(Approx~PointsToEstimate,col=mycol,type="l")
 
}

ggplot(errorsDF, aes(errorsDF$dt, errorsDF$Error)) +
  geom_line() +
  labs(x = 'Step Size', y = 'Error at t= 10', title = 'Error as fxn of step size') +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

# Exercise 8
## ExpGrowth.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"

require(deSolve) 
## package for integrating numerically ODEs see 
##http://cran.r-project.org/web/packages/deSolve/

## This function takes a step in time
ExponentialGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<-r*X 
  return(list(dX)) ## for some reason, you have to return a list
}

## This function runs the model and produces the trajectory
RunExponentialGrowth<-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = ExponentialGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

logisticGrowth <-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<-r*X*(1-(X/r))
  return(list(dX)) ## for some reason, you have to return a list
}

runLogisticGrowth<-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = logisticGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l", xlab="time",ylab="Density of X"))
  return(out)
}

alleeGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<- (((r * X) / (5 + X)) - (0.1 * X)) * X # parameters hard-coded from previous Allee exercise
  return(list(dX)) ## for some reason, you have to return a list
}

runAlleeGrowth<-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = alleeGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l", xlab="time",ylab="Density of X"))
  return(out)
}

levinsGrowth<-function(t, state, parameters) { 
  X <- state[1] 
  r <- parameters[1] 
  dX<- 0.1 * X * (1-X) - exp(X) # c paraeter taken from previous Allee exercise
  return(list(dX)) ## for some reason, you have to return a list
}

runLevinsGrowth <-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = ExponentialGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

thetaLogisticGrowth<-function(t, state, parameters) { 
  X <- state[1] 
  r <- parameters[1] 
  dX<- 2 * X * (1 - (X/15)**1.5) # using theta given on pg. 14, k from exercise 2
  return(list(dX)) ## for some reason, you have to return a list
}

runThetaLogisticGrowth <-function(MaxTime=10,GrowthRate=0.1,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(X=InitialX)
  out <- ode(y = state, times = times, func = thetaLogisticGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

RunExponentialGrowth(MaxTime=10,GrowthRate=0.1,InitialX=1.0)
runLogisticGrowth(MaxTime=10,GrowthRate=0.1,InitialX=1.0)
runAlleeGrowth(MaxTime=10,GrowthRate=0.1,InitialX=1.0)
runLevinsGrowth(MaxTime=10,GrowthRate=0.1,InitialX=1.0)
runThetaLogisticGrowth(MaxTime=10,GrowthRate=0.1,InitialX=1.0)
