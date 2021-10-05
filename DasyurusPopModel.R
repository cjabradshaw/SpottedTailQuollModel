##########################################################################################################################################
## spotted-tail quoll (Dasyurus maculatus) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


##############################
## DASYURUS (maculatus) (DM)

# mass
DM.mass <- 1.68 ## 1.7 kg (Glen 2008-AJZ)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
DM.rm.pred <- 10^(0.6914 - (0.2622*log10(DM.mass*1000)))
DM.lm.pred <- exp(DM.rm.pred)

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
DM.D.pred <- (10^(1.77 + (-1.02*log10(DM.mass))))/2 # divided by 2 for females only
DM.D.pred # animals/km2

# (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227); they said log10, but only makes sense at ln
lmu <- c(-2.4098957, -0.32328582, 2.5504842)
lDu <- c(4.3807473,2.2986877,-0.70201576)
lml <- c(-2.4042623,0.26574665,2.3397434)
lDl <- c(0.80104977, -1.3227539, -2.954604)
plot(lmu, lDu, pch=19, xlim=c(min(lml),max(lmu)), ylim=c(min(lDl),max(lDu)))
lmlDu.fit <- lm(lDu~lmu)
abline(lmlDu.fit, lty=2, col="red")
points(lml, lDl, pch=19)
lmlDl.fit <- lm(lDl~lml)
abline(lmlDl.fit, lty=2, col="red")
DM.D.pred.l <- exp(as.numeric(coef(lmlDl.fit)[1] + coef(lmlDl.fit)[2]*log(DM.mass)))/2
DM.D.pred.u <- exp(as.numeric(coef(lmlDu.fit)[1] + coef(lmlDu.fit)[2]*log(DM.mass)))/2
DM.D.pred <- DM.D.pred.u # better match to Johnson pers. obs.

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
DM.age.max <- round(10^(0.89 + (0.13*log10(DM.mass*1000))), 0)
DM.age.max <- 4 # Glen & Dickman 2013; Moro et al. 2019 & Cremona et al. 2104 (D. hallucatus); changed to 4 to account for a few females persisting into fourth year, but very low survival

## age vector
DM.age.vec <- 0:DM.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## % females breeding = .643 Glen & Dickman (or 0.91 for D maculatus - Cremona et al. 214)
5*.643/2
DM.F.pred <- exp(2.719 - (0.211*log(DM.mass*1000)))/2 # divided by 2 for females

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
DM.alpha <- ceiling(exp(-1.34 + (0.214*log(DM.mass*1000))))
DM.alpha <- 1 # Glen & Dickman 2013; Moro et al. 2019 (D. hallucatus)

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
DM.s.tran <- ln.a.s + b.s*log(DM.mass*1000) + log(1)
DM.s.ad.yr <- exp(-exp(DM.s.tran))

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.99*DM.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.8 # rate of mortality decline (also known as bt)
a2 <- 1 - (0.99*DM.s.ad.yr) # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.6e-04 # initial adult mortality rate (also known as βt)
b3 <- 2.5 # rate of mortality increase
longev <- DM.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DM.Sx <- c(0.95*DM.s.ad.yr, 1 - qx)
plot(x, DM.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DM.s.sd.vec <- 0.05*DM.Sx

## pre-breeding design with 0-1 survival in first row
DM.m.vec <- c(0.4*DM.F.pred, rep(DM.F.pred,3), 0.95*DM.F.pred) # allowed some breeding in < 1 year because sexual maturing < 12 months in D maculatus (Moro et al. 2019)
DM.m.sd.vec <- 0.05*DM.m.vec
plot(DM.age.vec, DM.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

## create matrix
DM.popmat <- matrix(data = 0, nrow=DM.age.max+1, ncol=DM.age.max+1)
diag(DM.popmat[2:(DM.age.max+1),]) <- DM.Sx[-1]
DM.popmat[DM.age.max+1,DM.age.max+1] <- 0 # catastrophic mortality at longevity DM.Sx[DM.age.max+1]
DM.popmat[1,] <- DM.m.vec
colnames(DM.popmat) <- c(0:DM.age.max)
rownames(DM.popmat) <- c(0:DM.age.max)
DM.popmat.orig <- DM.popmat ## save original matrix

## matrix properties
max.lambda(DM.popmat.orig) ## 1-yr lambda
DM.lm.pred
max.r(DM.popmat.orig) # rate of population change, 1-yr
DM.ssd <- stable.stage.dist(DM.popmat.orig) ## stable stage distribution
plot(DM.age.vec, DM.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(DM.popmat.orig, DM.age.max) # reproductive value
DM.gen.l <- G.val(DM.popmat.orig, DM.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
DM.pop.found <- round(area*DM.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
DM.init.vec <- DM.ssd * DM.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*DM.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

DM.tot.F <- sum(DM.popmat.orig[1,])
DM.popmat <- DM.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
DM.n.mat[,1] <- DM.init.vec

## set up projection loop
for (i in 1:t) {
  DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
}

DM.n.pred <- colSums(DM.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(DM.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
DM.K.max <- 1*DM.pop.found
DM.K.vec <- c(1, 0.25*DM.K.max, DM.K.max/2, 0.75*DM.K.max, DM.K.max) 
DM.red.vec <- c(1,0.96,0.76,0.45,0.20)
plot(DM.K.vec, DM.red.vec,pch=19,type="b")
DM.Kred.dat <- data.frame(DM.K.vec, DM.red.vec)

# logistic power function a/(1+(x/b)^c)
DM.param.init <- c(1, DM.K.max, 2)
DM.fit.lp <- nls(DM.red.vec ~ a/(1+(DM.K.vec/b)^c), 
                 data = DM.Kred.dat,
                 algorithm = "port",
                 start = c(a = DM.param.init[1], b = DM.param.init[2], c = DM.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DM.fit.lp.summ <- summary(DM.fit.lp)
plot(DM.K.vec, DM.red.vec, pch=19,xlab="N",ylab="reduction factor")
DM.K.vec.cont <- seq(1,2*DM.pop.found,1)
DM.pred.lp.fx <- coef(DM.fit.lp)[1]/(1+(DM.K.vec.cont/coef(DM.fit.lp)[2])^coef(DM.fit.lp)[3])
lines(DM.K.vec.cont, DM.pred.lp.fx, lty=3,lwd=3,col="red")

DM.a.lp <- coef(DM.fit.lp)[1]
DM.b.lp <- coef(DM.fit.lp)[2]
DM.c.lp <- coef(DM.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
DM.n.mat <- matrix(0, nrow=DM.age.max+1, ncol=(t+1))
DM.n.mat[,1] <- DM.init.vec
DM.popmat <- DM.popmat.orig

## set up projection loop
for (i in 1:t) {
  DM.totN.i <- sum(DM.n.mat[,i])
  DM.pred.red <- as.numeric(DM.a.lp/(1+(DM.totN.i/DM.b.lp)^DM.c.lp))
  diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.Sx[-1])*DM.pred.red
  DM.popmat[DM.age.max+1,DM.age.max+1] <- 0
  DM.popmat[1,] <- DM.m.vec
  DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
}

DM.n.pred <- colSums(DM.n.mat)
plot(yrs, DM.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=DM.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 500
itdiv <- iter/10

DM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DM.s.arr <- DM.m.arr <- array(data=NA, dim=c(t+1, DM.age.max+1, iter))

for (e in 1:iter) {
  DM.popmat <- DM.popmat.orig
  
  DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
  DM.n.mat[,1] <- DM.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    DM.s.alpha <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$alpha
    DM.s.beta <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$beta
    DM.s.stoch <- rbeta(length(DM.s.alpha), DM.s.alpha, DM.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    DM.fert.stch <- rnorm(length(DM.popmat[,1]), DM.m.vec, DM.m.sd.vec)
    DM.m.arr[i,,e] <- ifelse(DM.fert.stch < 0, 0, DM.fert.stch)
    
    DM.totN.i <- sum(DM.n.mat[,i], na.rm=T)
    DM.pred.red <- DM.a.lp/(1+(DM.totN.i/DM.b.lp)^DM.c.lp)
    
    diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.s.stoch[-1])*DM.pred.red
    DM.popmat[DM.age.max+1,DM.age.max+1] <- 0
    DM.popmat[1,] <- DM.m.arr[i,,e]
    DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
    
    DM.s.arr[i,,e] <- DM.s.stoch * DM.pred.red
    
  } # end i loop
  
  DM.n.sums.mat[e,] <- ((as.vector(colSums(DM.n.mat))/DM.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

DM.n.md <- apply(DM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
DM.n.up <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DM.n.lo <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,DM.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(DM.n.lo),1.05*max(DM.n.up)))
lines(yrs,DM.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,DM.n.up,lty=2,col="red",lwd=1.5)

DM.s.add <- DM.m.add  <- rep(0, DM.age.max+1)
for (m in 1:iter) {
  DM.s.add <- rbind(DM.s.add, DM.s.arr[ceiling(DM.gen.l):(t+1),,m])
  DM.m.add <- rbind(DM.m.add, DM.m.arr[ceiling(DM.gen.l):(t+1),,m])
}
DM.s.add <- DM.s.add[-1,]
DM.m.add <- DM.m.add[-1,]

DM.s.md <- apply(DM.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
DM.s.up <- apply(DM.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DM.s.lo <- apply(DM.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(DM.age.vec,DM.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(DM.s.lo),1.05*(max(DM.s.up))))
lines(DM.age.vec,DM.s.lo,lty=2,col="red",lwd=1.5)
lines(DM.age.vec,DM.s.up,lty=2,col="red",lwd=1.5)

DM.m.md <- apply(DM.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
DM.m.up <- apply(DM.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DM.m.lo <- apply(DM.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(DM.age.vec,DM.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(DM.m.lo),1.05*max(DM.m.up)))
lines(DM.age.vec,DM.m.lo,lty=2,col="red",lwd=1.5)
lines(DM.age.vec,DM.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

