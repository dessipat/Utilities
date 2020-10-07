
require(UReported)

##
## Under reported INAR with poisson inn.
## data generator
##

ruri <- function(n, param){
  alpha <- param[1]
  lambda <- param[2]
  omega <- param[3]
  q <- param[4]
  if(length(param)>4) p01 <- param[5]

  X <- Y <- I <- NULL
  X[1] <- Y[1] <- rpois(1, lambda)
  if(length(param)>4){
    I[1] <- rbinom(n=1, size=1, prob=omega)
    p11 <- 1 - p01*(1-omega)/omega
    for (i in 2:n) {
      ifelse(I[i-1]<1,
             I[i] <- rbinom(n=1, size=1, prob=p01),
             I[i] <- rbinom(n=1, size=1, prob=p11))
    }
  }else I <- rbinom(n=n, size=1, prob=omega)

  for (i in 2:n) {
    X[i] <- rbinom(n=1, size=X[i-1], prob=alpha) + rpois(1, lambda)
    if(I[i]<1){ Y[i] <- X[i]
    }else{ Y[i] <- rbinom(n=1, size=X[i], prob=q)}
  }
  data.frame(X=X, I=I, Y=Y)
}


##
## UR series visualization
##

urplot <- function(df){
  plot(df$X, type = "l", ylim=c(0,max(df$X)+3))
  lines(df$Y, col = 2)
  abline(h=0, lty=2)
  legend("topleft",
         legend=c("hidden", "observed"),
         col=c(1,2), lty=1, cex=0.6, box.lty=0, bg="transparent")
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
