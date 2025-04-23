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
library(randtoolbox)

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
zz  <- as.matrix(xx[, 1])

ols <- summary(lm(y~xx-1))

twotier.copula.g <- function(p, y, xx, S=200){
	
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

p.start.ncop.g <- c(8.21783174, 0.11862419, 0.06168671, 0.09860413, 
						      0.16665764, 0.16192581, 0.04905254, 0.13257851, 
						      0.25046346, 0.03208156, 0.16703915, 0.09453291, 
						      0.10502202, 0.15003797, 0.11842377, 0.04435116)

## p.start.ncop.g <- c(8.18289252, 0.12085086, 0.05572379, 0.09911819, 
##						   0.16866722, 0.15836194, 0.05031340, 0.13189161,
##						   0.25498796, 0.03073642, 0.17055509, 0.09301182, 
##						   0.12104155, 0.15725649, 0.12347888, 0.27435604)

opt.ml.ncop.g <- maxLik(twotier.copula.g, grad=NULL, hess=NULL, 
                     			   		method="nm", control=list(iterlim=20000), 
                     			   		constraints=list(ineqA=A, ineqB=B),
                     			   		start=p.start.ncop.g, y=y, xx=xx)

twotier.copula.f <- function(p, y, xx, S=200){
	
	nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  	sigv   <- p[nr+1] ##Assume homoscedastic two sided component
  	sigu   <- p[nr+2] ##parameter for u marginal
  	sigw  <- p[nr+3] ##parameter for w marginal
  	theta <- p[nr+4] ##Frank copula parameter
  
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
		
		Fw.si <- pexp(w.si, rate=1/sigw)
		Fu.si <- pexp(u.si, rate=1/sigu)

		## Construct summands
		## noise density
		fv <- dnorm(e[i]-w.si+u.si, sd=sigv)
		
		## Frank copula density
		num     <- theta*(1-exp(-theta))*exp(-theta*(Fw.si+Fu.si))
		denom <- exp(-theta)-1+(exp(-theta*Fw.si)-1)*(exp(-theta*Fu.si)-1)
		cdens  <- num/(denom^2)
		
		f.e[i] <- mean(fv*cdens)

	}

	ll <- log(f.e)
  
	return(ll) 
  
}

A <- matrix(0, 3, 16)
B <- matrix(0, 3, 1)

A[1, 13] <- 1		## sigma_v >0
A[2, 14] <- 1		## sigma_u>0
A[3, 15] <- 1		## sigma_w>0

p.start.ncop.f <- c(8.21783174, 0.11862419, 0.06168671, 0.09860413, 
						   0.16665764, 0.16192581, 0.04905254, 0.13257851, 
						   0.25046346, 0.03208156, 0.16703915, 0.09453291, 
						   0.10502202, 0.15003797, 0.11842377, 0.04435116)

## p.start.ncop.f <- c(8.21383298, 0.10626455,  0.06130027,  0.10986204,
##								  0.22004949, 0.17109216,  0.04117120,  0.12136838,
##								  0.25540681,  0.01587652, 0.17429455,  0.09317039,  
##								  0.01728862, 0.15901851,  0.13872828, -1.82474071)

opt.ml.ncop.f <- maxLik(twotier.copula.f, grad=NULL, hess=NULL, 
                     			   		method="nm", control=list(iterlim=20000), 
                     			   		constraints=list(ineqA=A, ineqB=B),
                     			   		start=p.start.ncop.f, y=y, xx=xx)

twotier.copula.c <- function(p, y, xx, S=200){
	
	nr  <- ncol(xx) ##Calculate number of regressors in regression
  
  	sigv   <- p[nr+1] ##Assume homoscedastic two sided component
  	sigu   <- p[nr+2] ##parameter for u marginal
  	sigw  <- p[nr+3] ##parameter for w marginal
  	theta <- p[nr+4] ##Clayton copula parameter
  
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
		
		Fw.si <- pexp(w.si, rate=1/sigw)
		Fu.si <- pexp(u.si, rate=1/sigu)

		## Construct summands
		## noise density
		fv <- dnorm(e[i]-w.si+u.si, sd=sigv)
		
		## Clayton copula density
		num     <- (Fw.si^(-theta)+Fu.si^(-theta)-1)^(-2-1/theta)
		denom <- (Fw.si*Fu.si)^(1+theta)
		cdens  <- (1+theta)*num/denom
		## New code
		## cdens <- ifelse(cdens<0, 0, cdens)
		
		f.e[i] <- mean(fv*cdens)

	}

	ll <- log(f.e)
  
	return(ll) 
  
}

A <- matrix(0, 4, 16)
B <- matrix(0, 4, 1)

A[1, 13] <- 1		## sigma_v >0
A[2, 14] <- 1		## sigma_u>0
A[3, 15] <- 1		## sigma_w>0
A[4, 16] <- 1		## theta >-1

B[4] <- 1

p.start.ncop.c <- c(8.21783174, 0.11862419, 0.06168671, 0.09860413, 
						      0.16665764, 0.16192581, 0.04905254, 0.13257851, 
						      0.25046346, 0.03208156, 0.16703915, 0.09453291, 
						      0.10502202, 0.15003797, 0.11842377, 0.04435116)

## p.start.ncop.c <- c(8.21991101, 0.11993180, 0.06242044, 0.09831884,  0.16715915, 0.16164817, 0.04920436, 0.13186435, 0.24994784, 0.03341981, 0.16500425, 0.09472129, 0.10362734, 0.14808171, 0.11566232, -0.04183639)

opt.ml.ncop.c <- maxLik(twotier.copula.c, grad=NULL, hess=NULL, 
                     			   		method="nm", control=list(iterlim=20000), 
                     			   		constraints=list(ineqA=A, ineqB=B),
                     			   		start=p.start.ncop.c, y=y, xx=xx)
zu <- as.matrix(xx[, 1])zw <- as.matrix(xx[, 1])zz <- as.matrix(xx[,1])ols    <- lm(y~xx+zz+zu+zw+zu:zw-1)ols.g <- olsols.c <- ols
ols.f  <- ols## Tabulate standard errorsse.g  <- sqrt(diag(vcov(opt.ml.ncop.g)))se.f  <- sqrt(diag(vcov(opt.ml.ncop.f)))se.c <- sqrt(diag(vcov(opt.ml.ncop.c)))stargazer(list(ols.g, ols.f, ols.c),          coef=list(coefficients(opt.ml.ncop.g),                     	coefficients(opt.ml.ncop.f),
                    	coefficients(opt.ml.ncop.c)),           se=list(se.g, se.f, se.c),          label="tab:2TSF_COP_ex",          covariate.labels = c("(Intercept)", "drv", "rec", "ffin",                               "ghw", "ca", "gar", "reg", "$\\ln$ lot",                               "bd", "fb", "sty", "$\\sigma_v$",                                "$\\sigma_u$", "$\\sigma_w$", "$\\rho$"),          star.cutoffs = c(0),          out = "2TSF_COP_ex.tex",          no.space=T, digits=3, single.row=F,          omit.stat = c("rsq", "adj.rsq", "f", "n", "ser"),          intercept.bottom=FALSE,          dep.var.labels.include = FALSE,          title="2 Tier Stochastic Frontier Estimates Allowing Dependence Between $u$ and $w$.",          column.labels=c("Gaussian", "Frank", "Clayton")          )

## Now calculate expectation of e^w or e^u
par.hat <- coefficients(opt.ml.ncop.g)

nr  <- ncol(xx)ep.hat  <- y-xx%*%par.hat[1:12]
n <- length(ep.hat)

sigv  <- par.hat[nr+1] ## Assume homoscedastic two sided component
sigu  <- par.hat[nr+2] ## parameter for u marginal
sigw <- par.hat[nr+3] ## parameter for w marginal
rho   <- par.hat[nr+4] ## correlation for Normal copula

## Halton Draws	to construct density
f.e <- as.numeric()
	
for (i in 1:n){

	## Construct Halton Draws
	H.draw <- halton(200, dim=2, init=FALSE)

	w.si <- qexp(H.draw[ ,1], rate=1/sigw)
	u.si <- qexp(H.draw[, 2], rate=1/sigu)

	## Construct summands
	## noise density
	fv <- dnorm(ep.hat[i]-w.si+u.si, sd=sigv)
		
	## copula density
	Fw.si <- pexp(w.si, rate=1/sigw)
	Fu.si <- pexp(u.si, rate=1/sigu)
		
	w.star <- qnorm(Fw.si)
	u.star <- qnorm(Fu.si)
		
	cdens <- (1/sqrt(1-rho^2))*exp((-(rho^2)*(w.star^2+u.star^2)+
														2*rho*w.star*u.star)/(2*(1-rho^2)))
		
	f.e[i] <- mean(fv*cdens)

}

## Now calculate expectation
Ew.sim <- as.numeric()
Eu.sim <- as.numeric()

for (i in 1:n){

	## Construct Halton Draws
	H.draw <- halton(200, dim=2, init=FALSE)

	w.si <- qexp(H.draw[ ,1], rate=1/sigw)
	u.si <- qexp(H.draw[, 2], rate=1/sigu)

	## Construct summands
	## noise density
	fv <- dnorm(ep.hat[i]-w.si+u.si, sd=sigv)
		
	## copula density
	Fw.si <- pexp(w.si, rate=1/sigw)
	Fu.si <- pexp(u.si, rate=1/sigu)
		
	w.star <- qnorm(Fw.si)
	u.star <- qnorm(Fu.si)
		
	cdens <- (1/sqrt(1-rho^2))*exp((-(rho^2)*(w.star^2+u.star^2)+
														2*rho*w.star*u.star)/(2*(1-rho^2)))
		
	Ew.sim[i] <- mean(exp(w.si)*fv*cdens)/f.e[i]
	Eu.sim[i] <- mean(exp(u.si)*fv*cdens)/f.e[i]

}

pdf("Kdens_Simulation_EworEu.pdf")
plot(density(Ew.sim), type="l", lwd=2, 
	   xlim=c(0.75, 2),
	   main="", sub="", ylab="density", xlab="")
lines(density(Eu.sim), lwd=2, lty=2)
legend("topright", legend=c("E[exp(w)|e]", "E[exp(u)|e]"), lty=1:2, lwd=2)
dev.off()

