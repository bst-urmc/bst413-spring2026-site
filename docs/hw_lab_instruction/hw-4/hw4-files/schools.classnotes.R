# file schools.classnotes.R - same as in classnotes
mu.hat <- function (tau, y, sigma.y){
sum(y/(sigma.y^2 + tau^2))/sum(1/(sigma.y^2 + tau^2))
}

V.mu <- function (tau, y, sigma.y){
1/sum(1/(tau^2 + sigma.y^2))
}

# Read data
schools <- read.table ("schools.txt", header=T)
J <- nrow (schools)
y <- schools$estimate
sigma.y <- schools$sd

n.grid <- 2000
tau.grid <- seq (.01, 40, length=n.grid)
log.p.tau <- rep (NA, n.grid)

# Step 1: Calculate log(p(log(tau) | y)).  Do this at mu=mu.hat, V=v.mu
# Note that mu.hat and V.mu are functions of tau.  
# This posterior uses a flat prior on tau, p(tau) \propto 1.

for (i in 1:n.grid){
  mu <- mu.hat (tau.grid[i], y, sigma.y)
  V <- V.mu(tau.grid[i], y, sigma.y)
  log.p.tau[i] <- .5*log(V) - .5*sum(log(sigma.y^2 + tau.grid[i]^2)) 
- .5*sum((y-mu)^2/(sigma.y^2 + tau.grid[i]^2))
}

# Compute the posterior density for tau on the log scale and rescale it to
# eliminate the possibility of computational overflow or underflow 

log.p.tau <- log.p.tau - max(log.p.tau)
p.tau <- exp(log.p.tau)
# the next line just scales p.tau so that the draws sum to 1
p.tau <- p.tau/sum(p.tau)

# Step 2: now that we have p(tau | y), sample from this.
n.sims <- 1000
# In the sample command, the first argument is a vector of possible
# values to sample, and the last argument are the corresponding probabilities
tau <- sample (tau.grid, n.sims, replace=T, prob=p.tau)

# Step 3: draw from p(mu | tau, y) (mean mu.hat, variance V.u)
# Step 4: draw from p(theta | mu, tau, y)
mu <- rep (NA, n.sims)
theta <- array (NA, c(n.sims,J))
theta.mean <- array (NA, c(n.sims,J))
theta.sd <- array (NA, c(n.sims,J))
for (i in 1:n.sims){
  mu[i] <- rnorm (1, mu.hat(tau[i],y,sigma.y), sqrt(V.mu(tau[i],y,sigma.y)))
  theta.mean[i,] <- (mu[i]/tau[i]^2 + y/sigma.y^2)/ (1/tau[i]^2 + 1/sigma.y^2)
  theta.sd[i,] <- sqrt(1/(1/tau[i]^2 + 1/sigma.y^2))
  theta[i,] <- rnorm (J, theta.mean[i,], theta.sd[i,])
}

# We now have created 1000 draws from the joint posterior distribution of tau,
# mu, theta.

#################################
### Part II: make plots

### Plot 1: a density estimate of the posterior of tau
par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(4.1,4.1,2.1,2.1))
plot(x=tau.grid,y=p.tau,type="l",main="Posterior for tau | y")

###  Plot 2: point estimates and 95% posterior intervals for theta_j
###  First get the empirical quantiles from the draws
theta.medians <- rep(NA,8)
theta.lower <- rep(NA,8)
theta.upper <- rep(NA,8)
cat("Quantiles: 2.5% 25% median 75% 97.5% \n")
for (j in 1:8) {
  cat("School ", schools[j,"school"], round(quantile(theta[,j],c(.025, .25, .5, .75, .975)),2),"\n")
  theta.medians[j] <- quantile(theta[,j],.5)
  theta.lower[j] <- quantile(theta[,j],.025)
  theta.upper[j] <- quantile(theta[,j],.975)
}

plot(x=schools[,"estimate"],y=theta.medians,xlab="Data mean",ylab="Posterior for theta",xlim=c(-3,28),ylim=c(-10,32),main="95% posterior intervals for theta")
abline(0,1,lty=2)
abline(h=7.9,lty=2) # overall mean
for (j in 1:8) {
  temp.x <- schools[j,"estimate"]
  segments(x0=temp.x,y0=theta.lower[j], x1=temp.x,y1=theta.upper[j])
}

# Plot 3: mean theta_j vs tau (averaged over mu)
for (j in 1:8) { 
smooth.estimate <- lowess(x=tau,y=theta.mean[,j])
if (j==1) {
  plot(smooth.estimate,type="l",ylim=c(-5,30),xlab="tau",ylab="Posterior mean of theta",main="Posterior for theta | tau")
} else {
  lines(smooth.estimate)
}
}

# Plot 4: mean SD(theta_j) vs tau (averaged over mu)
### (Using lowess here smoothed PAST the range of observed values)
tau.order <- tau[order(tau)]
theta.sd.order <- theta.sd[order(tau),]

plot(x=tau.order,y=theta.sd.order[,1],type="l",ylim=c(0,17),xlab="tau",ylab="Posterior for SD(theta)",main="Posterior for SD(theta) | tau")
for (j in 2:8) {
  lines(x=tau.order,y=theta.sd.order[,j])
}

mtext("School data (Gelman et al)",side=3,outer=T,line=0,cex=2)
