# file bicycles-snippet.R

y <- c(16,9,10,13,19,20,18,17,35,55)
n <- c(74,99,58,70,122,77,104,129,308,119)

calc.logpost <- function(u, v){
    ### Calculate the log(posterior) for a given u and v, where
    ### u = log(alpha/beta) and v=log(alpha+beta)
    ### Since the log posterior uses alpha and beta, have to
    ### convert u and v to alpha (e.g. a) and beta (e.g. b)
    a <- (exp(u)*exp(v))/(exp(u)+1)
    b <- exp(v)/(exp(u)+1)
    J<- length(y)
    log.prior <- log(a)+log(b)-(5/2)*log(a+b)
    log.lik <- 0
    for (j in 1:J) {
        ### The lgamma function is the log-gamma function
        log.lik <-  log.lik+lgamma(a+b) - lgamma(a) - lgamma(b) + lgamma(y[j]+a) + lgamma(n[j]-y[j]+b) -lgamma(a+b+n[j])
    }
    return(log.prior + log.lik)
}

### Set up grids for u and v - each will be a sequence of numbers (can use trial and error to find reasonable bounds).

u.grid <-  # fill in
v.grid <- # fill in

### Note the outer function!
### From documention, outer(vec1, vec2, fun) will evalute function="fun" over
### all values of vec1 and vec2, returning a matrix

log.post <- outer(u.grid, v.grid, calc.logpost)
### Recall that this is the LOG density.
### Now need to exponentiate and subtract the maximum value to prevent numeric overflow
post <- exp(log.post - max(log.post))

### Finally, make a contour plot - the x-axis is u.grid, y-axis is v.grid
contour(u.grid,v.grid,post,levels=seq(.05,.95,by=.1),xlab="log(alpha/beta)",ylab="log(alpha+beta)")
