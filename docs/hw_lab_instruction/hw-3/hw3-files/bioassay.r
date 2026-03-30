# bioassay.r 
# Example:analysis of biossay expriment
# Gelman and Rubin book, page 82
	
input <- function(LL=100){
  biossay <- data.frame(doses=c(-.863,-.296,-.053,.727),
		       rats=c(5,5,5,5),deaths=c(0.5,1,3,4.5),freq=c(0.5,1,3,4.5)/5)
  DD	  <- data.frame(y=log((biossay$freq)/(1-biossay$freq)),x=biossay$doses)
  estimates <- lm(y~x,data=DD)
  alpha.hat <- summary(estimates)$coeff[1,1]
  std.alpha <- summary(estimates)$coeff[1,2]
  beta.hat  <- summary(estimates)$coeff[2,1]
  std.beta  <- summary(estimates)$coeff[2,2]
  alpha       <- seq(-2,2,length=LL)
  beta       <-  seq(-5,10,length=LL)
  return(biossay,estimates,alpha,beta)
}

log.post<-function(alpha,beta,data=DD$biossay){
  ldens <- 0
  for (i in 1:length(data$doses)){
    theta <- 1/(1+exp(-alpha-beta*data$doses[i]))
    ldens  <- ldens + data$deaths[i]*log(theta)+(data$rats[i]-data$deaths[i])*log(1-theta)
  }
  ldens
}

plot.joint.post<-function(data,draws){
  contours <- seq(.05,.95,.1)
#  logdens <-outer(DD$alpha,DD$beta,log.post,data)
   logdens <-outer(DD$alpha,DD$beta,log.post)
  dens<-exp(logdens-max(logdens))
  contour(DD$alpha,DD$beta,dens,levels=contours,xlab="alpha",ylab="beta")
  points(draws$post.alpha,draws$post.beta)
  mtext("Posterior density",3,line=1,cex=1.2)
}

grid.value<-function(data=DD$biossay,alpha,beta){
  ll    <- length(alpha)
  PP    <- matrix(NA,ll,ll)
  for(i in 1:ll){ 
    for(j in 1:ll){
      PP[i,j]<-exp(log.post(alpha[i],beta[j],data))	
    }}
  return(PP/sum(PP))
}



sampling<-function(M=100,PP,alpha,beta){
  alpha.mar<-apply(PP,1,sum)
  alpha.cdf  <-  cumsum(alpha.mar)

  post.alpha  <-  rep(0,M)
  post.beta  <-  rep(0,M)
  for( m in 1:M){
    uuu<-runif(1,0,1)
    Fhat.alpha  <-  max( alpha.cdf[ alpha.cdf <= uuu])
#    post.alpha[m]  <- alpha[(1:length(alpha.cdf))[alpha.cdf == Fhat.alpha]]
    post.alpha[m]  <- alpha[alpha.cdf == Fhat.alpha]
    junk  <- length(alpha[alpha <= post.alpha[m]])
    PP[junk, ]  <-  PP[junk,]/sum(PP[junk,])                                   
    beta.cond.cdf  <-   cumsum(PP[junk,])
    uuu<-runif(1,0,1)
    Fhat.beta  <-  max( beta.cond.cdf[ beta.cond.cdf < uuu])
#    post.beta[m]  <-  beta[(1:length(beta.cond.cdf))[beta.cond.cdf == Fhat.beta]]    
     post.beta[m]  <-  beta[beta.cond.cdf == Fhat.beta]    
  }
  return(post.alpha,post.beta)
}



DD <- input(LL=200)
PP <- grid.value(DD$biossay,DD$alpha,DD$beta)
draws <- sampling(M=1000,PP,DD$alpha,DD$beta)
plot.joint.post(DD$bioassay,draws)




