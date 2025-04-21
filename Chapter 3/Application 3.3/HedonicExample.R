## Code to produce housing market example for 
## Chapter 3 of 2Tier book with Alecos. 

## Started January 30th, 2024. 
## This version, February 8th, 2024. 

## Clean memory
rm(list=ls())

## Load libraries
library(AER)
library(minqa)
library(latex2exp)
library(maxLik)
library(numDeriv)
library(Hmisc)
library(VGAM)

## Load in data
data("HousePrices")

attach(HousePrices)

n <- dim(HousePrices)[1]

## Set up variables as needed
## Here we go with log-log specification
## dependent variable
y <- log(price)

## independent variables
llot <- log(lotsize)

#lbd <- log(bedrooms)
#lsty <- log(stories)
#lfb  <- log(bathrooms)

lbd  <- bedrooms
lsty <- stories
lfb  <- bathrooms

## Should go back and not log bedrooms, stories or bathrooms
## These variables are all counts.

## Now convert factors to {0,1}
drv  <- ifelse(driveway=="yes", 1, 0)
rec  <- ifelse(recreation=="yes", 1, 0)
ffin <- ifelse(fullbase=="yes", 1, 0)
ghw  <- ifelse(gasheat=="yes", 1, 0)
ca   <- ifelse(aircon=="yes", 1, 0)
reg  <- ifelse(prefer=="yes", 1, 0)

## Garage places is a pure count
gar <- garage

## Now put everything together
xx <- cbind(1, drv, rec, ffin, ghw, ca, 
            gar, reg, llot, lbd, lfb, lsty)

zu <- as.matrix(xx[, 1])
zw <- as.matrix(xx[, 1])

ols <- summary(lm(y~xx-1))

## Starting parameter vector for optimization
p.start  <- c(8.0385, 0.1141, 0.0595, 0.0941, 0.1599, 0.1587, 
				    0.0487, 0.1217, 0.3044, 0.0823, 0.2650, 0.1753, 
					-2.0949, -1.9635, -2.3089)

twotier.ne <- function(p, y, xx, zu, zw){

	nr     <- ncol(xx) ##Calculate number of regressors in regression
	nzu  <- ncol(zu) ##Calculate number of determinants for u component
	nzw <- ncol(zw) ##Calculate number of determinants for w component

	sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  	sigu <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
  	sigw <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

  	if (sigv<= 1e-6){stop("Variance too small")}

  	e <- y - xx%*%p[1:nr]
  	a <- (sigv^2)/(2*sigw^2) - e/sigw
  	b <- e/sigv - sigv/sigw

  	alpha <- e/sigu + (sigv^2)/(2*sigu^2)
  	beta  <- -e/sigv - sigv/sigu

  	denom <- sigu+sigw

  	term1 <- exp(alpha)
  	term2 <- exp(a)

  	##return will send the summation of the log of the 
  	##density of the composed error

  	ll <- -log(denom)+log((pnorm(beta)*term1)+(pnorm(b)*term2))

  	return(ll)  ## The negative here is if I optimize using bobyqa, need to turn off if using maxLik

}

twotier.hn <- function(p, y, xx, zu, zw){

	  nr  <- ncol(xx) ##Calculate number of regressors in regression
	  nzu <- ncol(zu) ##Calculate number of determinants for u component
	  nzw <- ncol(zw) ##Calculate number of determinants for w component

	  sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
    sigu <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
    sigw <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

    if (sigv<= 1e-6){stop("Variance too small")}
    
    ## Calculate reparameterizations
    theta1  <- sigw/sigv
    theta2  <- sigu/sigv
    s       <- sqrt(sigv^2+sigw^2+sigu^2)
    omega1  <- s*sqrt(1+theta2^2)/theta1
    omega2  <- s*sqrt(1+theta1^2)/theta2
    lambda1 <- (theta2/theta1)*sqrt(1+theta1^2+theta2^2)
    lambda2 <- (theta1/theta2)*sqrt(1+theta1^2+theta2^2)

    e <- y - xx%*%p[1:nr]

	  D <- pbinorm(e/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))-
			    pbinorm(e/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2))

	  ll <- log(2*sqrt(2)/sqrt(pi))-log(s)-e^2/(2*s^2)+log(D)

    ##return will send the summation of the log of the 
    ##density of the composed error

    return(ll)  ## The negative here is if I optimize using bobyqa, need to turn off if using maxLik
    	
}

twotier.nce <- function(p, y, xx){
  
  nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  sigv    <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  a.prime <- exp(p[nr+2])
  b.prime <- exp(p[nr+3])
  m       <- pnorm(p[nr+4])
  
  if (sigv<= 1e-6){stop("Variance too small")}
  
  e <- y - xx%*%p[1:nr]
  
  ## Define omega parameters as before equation (12.8)
  omega2 <- (e/sigv)+b.prime*sigv
  omega3 <- (e/sigv)-a.prime*sigv
  
  term1 <- m*b.prime*exp(0.5*omega2^2)*pnorm(-omega2)
  term2 <- (1-m)*a.prime*exp(0.5*omega3^2)*pnorm(omega3)

  f.e   <- sqrt(2*pi)*dnorm(e/sigv)*(term1+term2)
  
  ##return will send the summation of the log of the 
  ##density of the composed error
  
  ll <- log(f.e)
  
  return(ll)  ## The negative here is if I optimize using bobyqa, need to turn off if using maxLik
  
}

twotier.nctn <- function(p, y, xx){
  
  nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  sigu <- exp(p[nr+2])
  sigw <- exp(p[nr+3])
  rho  <- 2*pnorm(p[nr+4])-1
  
  if (sigv<= 1e-6){stop("Variance too small")}
  
  e <- y - xx%*%p[1:nr]
  
  ## Define shorthands
  sige <- sqrt(sigv^2+sigw^2+sigu^2-2*rho*sigu*sigw)
  psiu <- sqrt(sigv^2+(1-rho^2)*sigu^2)
  psiw <- sqrt(sigv^2+(1-rho^2)*sigw^2)
  r    <- (rho*sigv^2+(1-rho^2)*sigu*sigw)/(psiu*psiw)
  aw   <- (rho*sigw-sigu)/(sige*psiw)
  au   <- (rho*sigu-sigw)/(sige*psiu)

  ## Density as defined in Equation 12.20
  f.e   <- (2*pi/acos(-rho))*dnorm(e/sige)*pbinorm(aw*e,-au*e,cov12=r)/sige
  
  ##return will send the summation of the log of the 
  ##density of the composed error
  
  ll <- log(f.e)
  
  return(ll)  ## The negative here is if I optimize using bobyqa, need to turn off if using maxLik
  
}

## Call optmizer, this minimizes so make sure negative of sum(ll) in likelihood function call
opt <- bobyqa(p.start, twotier.ne, lower=-Inf, upper=Inf, 
			 			  control=list(iprint=1), y, xx, zu, zw)

## Invert Hessian to obtain standard errors
H <- hessian(twotier.ne, x=opt$par, y=y, xx=xx, zu=zu, zw=zw)

## Invert to obtain SEs. 
ses <- sqrt(diag(solve(H)))

opt.ml <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
							   method="nm", control=list(iterlim=20000), 
							   start=opt$par, y=y, xx=xx, zu=zu, zw=zw)

## Make table of estimates
ests <- matrix(0, 30, 2)

ests[seq(from=1, to=24, by=2), 1] <- ols$coefficients[, 1]
ests[seq(from=2, to=24, by=2), 1] <- ols$coefficients[, 2]

ests[seq(from=1, to=30, by=2), 2] <- summary(opt.ml)$estimate[, 1]
ests[seq(from=2, to=30, by=2), 2] <- summary(opt.ml)$estimate[, 2]

## Now correct variance estimates for 2Tier
ests <- ests[c(1:25, 27, 29), ]
ests[25:27, 2] <- exp(ests[25:27, 2])
ests[25, 1] <- sqrt(var(residuals(ols)))

ests <- round(ests, 3)

latex(ests)

## Alecos wants us to use E[e^(w)|epsilon] and E[e^(u)|epsilon] to assess
## information deficiency in this example. 

## So strip off parameter estimates and create epsilon series.
par.hat <- coefficients(opt.ml)   ##opt$par
ep.hat  <- y-xx%*%par.hat[1:12]

sig.v <- exp(par.hat[13])
sig.u <- exp(par.hat[14])
sig.w <- exp(par.hat[15])

## Use 8.23 and 8.26 to construct metrics
## Setup necessary parameters needed.
lambda <- 1/sig.w+1/sig.u
a1     <- sig.v^2/(2*sig.u^2)+ep.hat/sig.u
b1     <- -(ep.hat/sig.v+sig.v/sig.u)
a2     <- sig.v^2/(2*sig.w^2)-ep.hat/sig.w
b2     <-  ep.hat/sig.v-sig.v/sig.w
chi1   <- pnorm(b2)+exp(a1-a2)*pnorm(b1)
chi2   <- exp(a2-a1)*chi1

Eew.cond  <- (lambda/(chi2*(lambda-1)))*(pnorm(b1)+exp(0.5*((b2+sig.v)^2-b1^2))*pnorm(b2+sig.v))

Eemw.cond <- (lambda/(chi2*(1+lambda)))*(pnorm(b1)+exp(a2-a1-b2*sig.v+0.5*sig.v^2)*pnorm(b2-sig.v))

Eeu.cond  <- (lambda/(chi1*(lambda-1)))*(pnorm(b2)+exp(0.5*((b1+sig.v)^2-b2^2))*pnorm(b1+sig.v))

Eemu.cond <- (lambda/(chi1*(1+lambda)))*(pnorm(b2)+exp(a1-a2-b1*sig.v+0.5*sig.v)*pnorm(b1-sig.v))

## Now calculate the M1 and M2 metrics (these are information deficiency 
## relative to the actual price) and M5 and M6 metrics (these are information 
## deficiency relative to balanced price).
M1 <- 1-Eemw.cond
M2 <- Eeu.cond-1

M5 <- Eew.cond-1
M6 <- 1-Eemu.cond

## Alecos wanted M7 as well
M7 <- Eew.cond*Eemu.cond-1

quantile(M1, c(0.25, 0.5, 0.75))
quantile(M2, c(0.25, 0.5, 0.75))
quantile(M5, c(0.25, 0.5, 0.75))
quantile(M6, c(0.25, 0.5, 0.75)) 
quantile(M7, c(0.25, 0.5, 0.75))
mean(M7)

## Now lets plot out these two densities
# pdf("EwEuCond.pdf")
# plot(density(Eew.cond), type="l", lwd=2, lty=1, 
		# xlim= c(1, 1.35), main="",
		# xlab="Information Deficiency", ylab="Density")
# lines(density(Eeu.cond), lwd=2, lty=2)
# legend("topright", legend=c(expression(E~"["~e^w~"|"~epsilon~"]"), 
											 # expression(E~"["~e^u~"|"~epsilon~"]")),
			# lwd=rep(2,2), lty=1:2)
# dev.off()

pdf("M1M2_hedonic.pdf")
plot(density(M1), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.35), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M2), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M5M6_hedonic.pdf")
plot(density(M5), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.25), ylim=c(0, 25), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M6), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M7_hedonic.pdf")
plot(density(M7), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density")
dev.off()

## Calculate Vuong (1989) test between NE and NHN 2TSF models
## Estimate both NE and HN 2TSF models
opt.ml.ne <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
                    method="nm", control=list(iterlim=20000), 
                    start=p.start, y=y, xx=xx, zu=zu, zw=zw)

opt.ml.hn <- maxLik(twotier.hn, grad=NULL, hess=NULL, 
                    method="nm", control=list(iterlim=20000), 
                    start=p.start, y=y, xx=xx, zu=zu, zw=zw)

li.NE <- twotier.ne(coefficients(opt.ml.ne), y, xx, zu, zw)
li.HN <- twotier.hn(coefficients(opt.ml.hn), y, xx, zu, zw)

li.diff <- li.NE-li.HN

LRn        <- sum(li.diff)
omegahat.n <- mean(li.diff^2)-mean(li.diff)^2
Tn         <- sqrt(1/n)*LRn/sqrt(omegahat.n)

## Calculate c and -c from standard normal for given significance level
alpha <- 0.05
c.max <- qnorm(1-alpha/2)
c.min <- qnorm(alpha/2)

p.start.nce <- c(7.427694, 0.1349851, 0.05568961, 0.104027, 
                 		 0.2188659, 0.2531409, -0.01036444, 0.1745869, 
                 		 0.307643, -0.009254036, 0.2543505, 0.0582116, 
                 		-3.516158, 0.9112417, -1.065818, -2.572593)

p.start.nctn <- c(8.0385, 0.1141, 0.0595, 0.0941, 0.1599, 0.1587, 
                  		  0.0487, 0.1217, 0.3044, 0.0823, 0.2650, 0.1753, 
                  		 -3.516158, 0.9112417, -1.065818, 1)

opt.ml.nce <- maxLik(twotier.nce, grad=NULL, hess=NULL, 
                     method="nm", control=list(iterlim=20000), 
                     start=p.start.nce, y=y, xx=xx)

opt.ml.nctn <- maxLik(twotier.nctn, grad=NULL, hess=NULL, 
                      method="nm", control=list(iterlim=20000), 
                      start=p.start.nctn, y=y, xx=xx)

twotier.nce <- function(p, y, xx){
  
  nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  sigv       <- p[nr+1]	##Assume homoscedastic two sided component
  a.prime <- p[nr+2]
  b.prime <- p[nr+3]
  m          <- p[nr+4]
  
  if (sigv<= 1e-6){stop("Variance too small")}
  
  e <- y - xx%*%p[1:nr]
  
  ## Define omega parameters as before equation (12.8)
  omega2 <- (e/sigv)+b.prime*sigv
  omega3 <- (e/sigv)-a.prime*sigv
  
  term1 <- m*b.prime*exp(0.5*omega2^2)*pnorm(-omega2)
  term2 <- (1-m)*a.prime*exp(0.5*omega3^2)*pnorm(omega3)

  f.e   <- sqrt(2*pi)*dnorm(e/sigv)*(term1+term2)
  
  ##return will send the summation of the log of the 
  ##density of the composed error
  
  ll <- log(f.e)
  
  return(ll) 
  
}

twotier.nctn <- function(p, y, xx){
  
  nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  sigv  <- p[nr+1] 	##Assume homoscedastic two sided component
  sigu  <- p[nr+2]
  sigw <- p[nr+3]
  rho   <- p[nr+4]
  
  if (sigv<= 1e-6){stop("Variance too small")}
  
  e <- y - xx%*%p[1:nr]
  
  ## Define shorthands
  sige <- sqrt(sigv^2+sigw^2+sigu^2-2*rho*sigu*sigw)
  psiu <- sqrt(sigv^2+(1-rho^2)*sigu^2)
  psiw <- sqrt(sigv^2+(1-rho^2)*sigw^2)
  r       <- (rho*sigv^2+(1-rho^2)*sigu*sigw)/(psiu*psiw)
  aw   <- (rho*sigw-sigu)/(sige*psiw)
  au    <- (rho*sigu-sigw)/(sige*psiu)

  ## Density as defined in Equation 12.20
  f.e   <- (2*pi/acos(-rho))*dnorm(e/sige)*pbinorm(aw*e,-au*e,cov12=r)/sige

  ##return will send the summation of the log of the 
  ##density of the composed error

  ll <- log(f.e)

  return(ll)

}

A <- matrix(0, 16, 5)
B <- matrix(0, 5, 1)

A[13, 1] <- 1		## sigma_v >0
A[14, 2] <- 1		## a.prime>0
A[15, 3] <- 1		## b.prime>0
A[16, 4] <- 1		## m >0
A[16, 5] <- -1	## m<1
A <- t(A)
B[4] <- 1
B[5] <- 1

p.start.nce <- c(7.82477633, 0.11087207, 0.05757073, 0.10282750, 
						 0.15951663, 0.15413679, 0.04375504, 0.12354117, 
						 0.30143483, 0.03092422, 0.16379650, 0.09652393, 
						 0.11538181, 9.07547700, 7.16551247, 0.64971669)

opt.ml.nce <- maxLik(twotier.nce, grad=NULL, hess=NULL, 
                     			   method="nm", control=list(iterlim=20000), 
                     			   constraints=list(ineqA=A, ineqB=B),
                     			   start=p.start.nce, y=y, xx=xx)

A <- matrix(0, 16, 5)
B <- matrix(0, 5, 1)

A[13, 1] <- 1		## sigma_v >0
A[14, 2] <- 1		## sigma_u>0
A[15, 3] <- 1		## sigma_w>0
A[16, 4] <- 1		## rho >-1
A[16, 5] <- -1	## rho<1
A <- t(A)
B[5] <- 1


p.start.nctn <- c(7.96949933, 0.12758384, 0.05329819, 0.11283846, 
						  0.18143522, 0.16622376, 0.04746548, 0.12184049, 
						  0.28011132, 0.04260589, 0.15475334, 0.09307778, 
						  0.04183599, 0.41250827, 0.34719496, 0.81821724)

opt.ml.nctn <- maxLik(twotier.nctn, grad=NULL, hess=NULL, 
                     			   method="nm", control=list(iterlim=20000), 
                     			   constraints=list(ineqA=A, ineqB=B),
                     			   start=p.start.nctn, y=y, xx=xx)

library(randtoolbox)

twotier.copula <- function(p, y, xx, S=200){
	
	nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  	sigv  <- p[nr+1] ##Assume homoscedastic two sided component
  	sigu  <- p[nr+2] ##parameter for u marginal
  	sigw <- p[nr+3] ##parameter for w marginal
  	rho   <- p[nr+4] ##correlation for Normal copula
  
  	## Calculate residuals
  	e <- y - xx%*%p[1:nr]
  	
  	n <- length(e)

	## Construct Halton Draws
	## init=TRUE here guarantees that everytime
	## this function is called the same Halton draws 
	## are produced
	toss <- halton(1, dim=2, init=TRUE)
	
	f.e <- as.numeric()
	
	for (i in 1:n){

		## Construct Halton Draws
		H.draw <- halton(S, dim=2, init=FALSE)

		w.si <- qexp(H.draw[ ,1], rate=1/sigw)
		u.si <- qexp(H.draw[, 2], rate=1/sigu)

		## Construct summands
		## noise density
		fv <- dnorm(e[i]-w.si+u.si, sd=sigv)
		
		## copula density
		Fw.si <- pexp(w.si, rate=1/sigw)
		Fu.si <- pexp(u.si, rate=1/sigu)
		
		w.star <- qnorm(Fw.si)
		u.star <- qnorm(Fu.si)
		
		cdens <- (1/sqrt(1-rho^2))*exp((-(rho^2)*(w.star^2+u.star^2)+2*rho*w.star*u.star)/(2*(1-rho^2)))
		
		f.e[i] <- mean(fv*cdens)

	}

	ll <- log(f.e)
  
	return(ll) 
  
}

A <- matrix(0, 5, 16)
B <- matrix(0, 5, 1)

A[1, 13] <- 1		## sigma_v >0
A[2, 14] <- 1		## sigma_u>0
A[3, 15] <- 1		## sigma_w>0
A[4, 16] <- 1		## rho >-1
A[5, 16] <- -1	## rho<1

B[4] <- 1
B[5] <- 1

p.start.ncop <- c(8.21783174, 0.11862419, 0.06168671, 0.09860413, 
						   0.16665764, 0.16192581, 0.04905254, 0.13257851, 
						   0.25046346, 0.03208156, 0.16703915, 0.09453291, 
						   0.10502202, 0.15003797, 0.11842377, 0.04435116)

p.start.ncop <- c(8.18289252, 0.12085086, 0.05572379, 0.09911819, 
						   0.16866722, 0.15836194, 0.05031340, 0.13189161,
						   0.25498796, 0.03073642, 0.17055509, 0.09301182, 
						   0.12104155, 0.15725649, 0.12347888, 0.27435604)

opt.ml.ncop <- maxLik(twotier.copula, grad=NULL, hess=NULL, 
                     			   method="nm", control=list(iterlim=20000), 
                     			   constraints=list(ineqA=A, ineqB=B),
                     			   start=p.start.ncop, y=y, xx=xx)

twotier.Icopula <- function(p, y, xx, S=200){
	
	nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  	sigv  <- p[nr+1]  ##Assume homoscedastic two sided component
  	sigu  <- p[nr+2]  ##parameter for u marginal
  	sigw <- p[nr+3]  ##parameter for w marginal
  
  	## Calculate residuals
  	e <- y - xx%*%p[1:nr]
  	
  	n <- length(e)

	## Construct Halton Draws
	## init=TRUE here guarantees that everytime
	## this function is called the same Halton draws 
	## are produced
	toss <- halton(1, dim=2, init=TRUE)
	
	f.e <- as.numeric()

	for (i in 1:n){

		## Construct Halton Draws
		H.draw <- halton(S, dim=2, init=FALSE)

		w.si <- qexp(H.draw[ ,1], rate=1/sigw)
		u.si <- qexp(H.draw[, 2], rate=1/sigu)

		## Construct summands
		## noise density
		fv <- dnorm(e[i]-w.si+u.si, sd=sigv)
		
		## copula density
		cdens <- 1
		
		f.e[i] <- mean(fv*cdens)

	}

	ll <- log(f.e)
  
	return(ll) 
  
}

A <- matrix(0, 3, 15)
B <- matrix(0, 3, 1)

A[1, 13] <- 1		## sigma_v >0
A[2, 14] <- 1		## sigma_u>0
A[3, 15] <- 1		## sigma_w>0

p.trial <- c(8.02418886,  0.09831496,  0.05675755,  0.10951549,
				 0.17710960,  0.17526184,  0.04632807,  0.13223598,
				 0.27427445,  0.01512169,  0.17793611,  0.09707987,
				 exp(-2.10234631), exp(-2.03439196), exp(-1.98687040))

sum(twotier.Icopula(p.trial, y, xx, S=200))

## Frank copula density
num     <- theta*(1-exp(-theta))*exp(-theta*(Fw.si+Fu.si))
denom <- exp(-theta)-1+(exp(-theta*Fw.si)-1)*(exp(-theta*Fu.si)-1)
cdens  <- num/(denom^2)

## Clayton copula density
num     <- (Fw.si^(-theta)+Fu.si^(-theta)-1)^(-2-1/theta)
denom <- (Fw.si*Fu.si)^(1+theta)
cdens  <- (1+theta)*num/denom
		