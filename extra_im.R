
  ##
  ## ZIP series generator
  ##

zip <- function(n=1, lambda, rho){
  r <- rbinom(n,1,rho)
  (1-r)*rpois(n, lambda = lambda)
}


  ##
  ## GOOD DISTRIBUTION pmf
  ##

require(VGAM)

dlerch <- function(x, q, nu, a)  (q^x)*(x+a)^(-nu)/lerch(q, nu, a)

plerch <- function(q, q1, nu, a, lower.tail=TRUE)
{
  res <- ifelse(q>=0, vapply(q, function(x) sum(dlerch(0:x, q1, nu, a)), 1), 0)
  if (lower.tail == FALSE) res <- 1 - res
  return(res)
}

qlerch <- function(p, q, nu, a, lower.tail=TRUE)
{
  if (q >= 1 | q <= 0 | a <= 0)
  {
    warning("improper parameter specification")
    return(NaN)
  }
  q1 <- vector()
  j <- 1
  if (lower.tail == FALSE) p <- 1-p
  while(!is.na(p[j]))
  {
    f <- 0
    q1[j] <- -1
    if ((p[j] > 1 | p[j] < 0))
    {
      warning("NaNs produced")
      q1[j] <- NaN
    }
    if (p[j] == 1)
    {
      q1[j] <- Inf
    } else {
      if(p[j] < 1)
      {
        while(f<p[j])
        {
          q1[j] <- q1[j] + 1
          f <- dlerch(q1[j], q, nu, a) + f
        }
      }
    }
    j <- j + 1
  }
  q1 <- ifelse(is.nan(q1) | q1>=0, q1, 0)
  return(q1)
}

rlerch <- function(n, q1, nu, a)
{
  inverse <- function(y){
    uniroot(function(q){plerch(q, q1, nu, a) - y},
            lower = -100, upper = 100,
            tol=1e-3, extendInt="yes")[1]
  }
  u <- runif(n)
  res <- vector()
  for (i in 1:n)
  {
    res[i] <- round(as.numeric(inverse(u[i])))
  }
  return(res)
}


  ##
  ## GOOD DISTRIBUTION series generator
  ##

rgood <- function(n, q1, nu) rlerch(n, q1, nu, a=1)


  ##
  ## GENERAL POISSON series generator
  ##

top_gp <- function(x,l,t,e=10^-8) dgenp(x,l,t)-e

rgenp <- function(n, lambda, theta, upp=10^8){
  top <- ceiling(uniroot(top_gp, l=lambda, t=theta, lower = 1, upper = upp)$root)
  replicate(n, which.max(cumsum(dgenp(0:top, lambda, theta))>runif(1)))-1
}


  ##
  ## COM-POISSON series generator
  ##

# require(COMPoissonReg)


##
## HERMITE series generator
##

require(hermite)

  ##
  ## INAR(p) series generator
  ##

binthin <- function(x, alpha){
  k <- length(x)
  p <- length(alpha)
  if(!k>p) return( sum(rbinom(n=k, size=x, prob=alpha[1:k])) )
  return( sum(rbinom(n=p, size=x[k:(k-p)], prob=alpha)) )
}

rinn <- function(par, dist){
  if(dist=="pois"){ rpois(1, par)
  }else if(dist=="zip"){ zip(n=1, par[1], par[2])
  }else if(dist=="nb"){ rnbinom(n=1, par[1], par[2])
  }else if(dist=="geom"){ rgeom(n=1, par)
  }else if(dist=="good"){ rgood(n=1, par[1], par[2])
  }else if(dist=="genp"){ rgenp(n=1, par[1], par[2])
  }else if(dist=="herm"){ rhermite(n=1, a=par[1], b=par[2], m=par[3])
  }
}

rip <- function(n, alpha, par, dist, thin="bin", gamma = NULL){
  Y <- NULL
  Y[1] <- rinn(par, dist)
  if(thin=="bin") for (i in 2:n) Y[i] <- binthin(Y,alpha) + rinn(par, dist)
  if(thin=="I2") for (i in 2:n) Y[i] <- sum(rI2(n=Y[i-1], alp = alpha, gam = gamma)) + rinn(par, dist)
  if(thin=="I3") for (i in 2:n) Y[i] <- sum(rI3(n=Y[i-1], alp = alpha, gam = gamma)) + rinn(par, dist)
  return(Y)
}

  ##
  ## Residuals visualization
  ##

resplot <- function(res){
  par(fig=c(0,0.9,0,1))
  plot(res, pch=20, cex=0.35)
  legend("topleft",legend=c(paste("mu =", round(mean(res),2)),
                            paste("s2 =", round(var(res),2))),
         cex=0.6, box.lty=0, inset = -0.01, bg="transparent", col = 3)
  legend("topright",legend=c(paste("p.value =", round(shapiro.test(res)$p,3))),
         cex=0.6, box.lty=0, inset = 0.01, bg="transparent", col = 3)
  par(fig=c(0.70 ,1,0,1),new=TRUE)
  boxplot(res, axes=FALSE)
  par(fig=c(0,1,0,1))
}

  ##
  ## Series visualization
  ##

splot <- function(x){
  plot(x, type = "l", ylim = c(0, max(x)+2))
  abline(h=0, lty=2)
}


# Auxiliar function to estimate INAR(p) by least squares method

cls <- function(x, param){
  n <- length(x)
  p <- length(param)-1
  
  l <- param[p+1]
  a <- param[1:p]
  q <- NULL
  for (i in (p+1):n) {
    q[i-p]<-x[i] - l - sum(a*x[(i-1):(i-p)])
  }
  sum(q^2)
}

# I2 distribution random series generator

pI2 <- function(x, alph, gam){
  sapply(x, function(y)  sum(dI2(0:y, alph, gam)))
}

rI2 <- function(n, alp, gam, e=10^-5){
  top <- ceiling(uniroot(function(x) 1-pI2(x, alp, gam)-e, lower = 1, upper = 200)$root)
  sapply(1:n, function(y) which(pI2(0:top, alp, gam) > runif(1))[1] -1)
}

# alpha <- 0.7
# gamma <- 0.6
# k <- 6
# nsim <- 10^k
# {
#   tm <- Sys.time()
#   d <- rI2(nsim, alpha, gamma)
#   tab <- table(d)
#   Sys.time()-tm
# }
# data.frame(empiric =tab/nsim, real = round(dI2(as.numeric(names(tab)),alpha, gamma),k))

# I3 distribution random series generator

pI3 <- function(x, alph, gam){
  sapply(x, function(y)  sum(dI3(0:y, alph, gam)))
}

rI3 <- function(n, alp, gam, e=10^-5){
  top <- ceiling(uniroot(function(x) 1-pI3(x, alp, gam)-e, lower = 1, upper = 200)$root)
  sapply(1:n, function(y) which(pI3(0:top, alp, gam) > runif(1))[1] -1)
}

# alpha <- 0.7
# gamma <- 0.6
# k <- 5
# nsim <- 10^k
# {
#   tm <- Sys.time()
#   d <- rI3(nsim, alpha, gamma)
#   tab <- table(d)
#   Sys.time()-tm
# }
# data.frame(empiric =tab/nsim, real = round(dI3(as.numeric(names(tab)),alpha, gamma),k))
