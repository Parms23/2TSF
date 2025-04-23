## Code to produce wage example for 
## Chapter 10 of 2Tier book with Alecos. 

## Started July 10th, 2024. 
## This version July 11th, 2024.

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
library(stargazer)
library(plm)
library(mvtnorm)

## Set seed for replications purposes
set.seed(123)

## Load in data
data <- read.table("psidq.dat", h=F)

## Organize data into a panel
names(data) <- c("YEAR", "WAGE", "EDUC", "EXP", "TEN", "SEQN68")

## Now make a panel
pdata <- pdata.frame(data, index=c("SEQN68", "YEAR"))

pdata$TEN <- pdata$TEN/12

## Estimate OLS to replicate column 1 -- Table 1 of PY 1996
py.col1 <- lm(WAGE~EDUC+EXP+I(EXP^2)+TEN+I(TEN^2), 
					  data=pdata)

## dependent variable
y <- pdata$WAGE

xx <- as.matrix(cbind(1, pdata$EDUC, pdata$EXP, (pdata$EXP)^2,
									pdata$TEN, (pdata$TEN)^2))

zu <- as.matrix(rep(1, dim(pdata)[1]))
zw <- as.matrix(rep(1, dim(pdata)[1]))

nmulti <- 5

## Starting parameter vector for optimization
p.start.ne  <- c(coefficients(py.col1), log(0.25), log(0.25), log(0.25))

twotier.ne <- function(p, y, xx, zu, zw){

	nr     <- ncol(xx) ##Calculate number of regressors in regression
	nzu  <- ncol(zu) ##Calculate number of determinants for u component
	nzw <- ncol(zw) ##Calculate number of determinants for w component

	sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  	sigu <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
  	sigw <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

  	#if (sigv<= 1e-6){stop("Variance too small")}

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
  	
  	if(any(is.na(ll))){return(-.Machine$double.xmax)}
	if(is.null(ll)){return(-.Machine$double.xmax)}

  	return(ll)  

}

twotier.ne.py <- function(p, y, xx, zu, zw, ybar, xbar){

	nr     <- ncol(xx) ##Calculate number of regressors in regression
	nzu  <- ncol(zu) ##Calculate number of determinants for u component
	nzw <- ncol(zw) ##Calculate number of determinants for w component

	sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  	sigu <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
  	sigw <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

  	#if (sigv<= 1e-6){stop("Variance too small")}

	a.hat <- ybar-xbar%*%t(t(p[1:nr]))

  	e <- y - xx%*%p[1:nr]-a.hat
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
  	
  	if(any(is.na(ll))){return(-.Machine$double.xmax)}
	if(is.null(ll)){return(-.Machine$double.xmax)}

  	return(ll)  

}

twotier.hn <- function(p, y, xx, zu, zw){

	nr    <- ncol(xx)  ##Calculate number of regressors in regression
	nzu <- ncol(zu)  ##Calculate number of determinants for u component
	nzw <- ncol(zw) ##Calculate number of determinants for w component

	sigv <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
    sigu <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
    sigw <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

    #if (sigv<= 1e-6){stop("Variance too small")}
    
    ## Calculate reparameterizations
    theta1      <- sigw/sigv
    theta2      <- sigu/sigv
    s              <- sqrt(sigv^2+sigw^2+sigu^2)
    omega1   <- s*sqrt(1+theta2^2)/theta1
    omega2   <- s*sqrt(1+theta1^2)/theta2
    lambda1  <- (theta2/theta1)*sqrt(1+theta1^2+theta2^2)
    lambda2  <- (theta1/theta2)*sqrt(1+theta1^2+theta2^2)

    e <- y - xx%*%p[1:nr]

	tryCatch(D <- pbinorm(e/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))-
			    pbinorm(e/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2)), error=function(e) NULL)

	ll <- NULL

	tryCatch(ll <- log(2*sqrt(2)/sqrt(pi))-log(s)-e^2/(2*s^2)+
								log(D), error=function(e) NULL)

	if(any(is.na(ll))){return(-.Machine$double.xmax)}
	if(is.null(ll)){return(-.Machine$double.xmax)}

    ##return will send the log of the density of the composed error
    return(ll)  
    	
}

twotier.ne.plugin <- function(p, ehat){

	ehat <- as.numeric(ehat)

	sigv <- exp(p[1])
  	sigu <- exp(p[2])
  	sigw <- exp(p[3])
  	
  	psumsw <- sigu-sigw

  	a1 <- (sigv^2)/(2*sigu^2) + (ehat+psumsw)/sigu
  	b1 <- -((ehat+psumsw)/sigv + sigv/sigu)

  	a2 <- (sigv^2)/(2*sigw^2) - (ehat+psumsw)/sigw
  	b2 <- (ehat+psumsw)/sigv - sigv/sigw

  	denom <- sigu+sigw

  	##return will send the summation of the log of the 
  	##density of the composed error

  	ll <- -log(denom)+log(pnorm(b1)*exp(a1)+pnorm(b2)*exp(a2))
  	
  	if(any(is.na(ll))){return(-.Machine$double.xmax)}
	if(is.null(ll)){return(-.Machine$double.xmax)}

  	return(ll)  

}

max.j <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
						   method="nm", control=list(iterlim=80000), 
					       start=p.start.ne, y=y, xx=xx, zu=zu, zw=zw)

# is.max         <- max.j$maximum
# coef.start1  <- max.j$estimate
# which.multi <- 1

# mu <- as.numeric(coefficients(max.j))
# vcov <- diag(sqrt(diag(vcov(max.j))))

# for (jjj in 1:nmulti){

	# start.p1 <- p.start.ne+rmvnorm(1, mean=mu, sigma=vcov)

	# ## Estimate likelihood
	# max.j1 <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
							     # method="bfgs", control=list(iterlim=80000), 
							     # start=start.p1, y=y, xx=xx, zu=zu, zw=zw)

	# if(is.null(max.j1$maximum)){max.j1$maximum <- is.max-1}

	# if(max.j1$maximum>is.max){
					# max.j          <- max.j1
					# is.max        <- max.j$maximum
					# which.multi <- jjj+1
	# }
# }

py.col2 <- max.j

## Estimate OLS to replicate column 1 -- Table 1 of PY 1996
py.fe <- plm(WAGE~EDUC+EXP+I(EXP^2)+TEN+I(TEN^2), 
				    effect="individual", model="within", data=pdata)

YMEAN <- Between(pdata$WAGE)

XMEANS <- cbind(Between(pdata$EXP), Between((pdata$EXP)^2),
							   Between(pdata$TEN), Between((pdata$TEN)^2))

beta.hat <- as.matrix(coefficients(py.fe))
alpha.hat <- YMEAN-XMEANS%*%beta.hat

## Starting parameter vector for optimization
p.start.py  <- c(coefficients(py.fe), log(0.03), log(0.1939), log(0.17779))

py.col3 <- maxLik(twotier.ne.py, grad=NULL, hess=NULL, 
						     method="nm", control=list(iterlim=80000), 
					         start=p.start.py, y=y, xx=xx[, 3:6], zu=zu, zw=zw,
					         ybar=YMEAN, xbar=XMEANS)

## Now we assume RE
## Estimate OLS to replicate column 1 -- Table 1 of PY 1996
py.re <- plm(WAGE~EDUC+EXP+I(EXP^2)+TEN+I(TEN^2), 
				    effect="individual", model="random", data=pdata)

ehat.re        <- residuals(py.re)
alphahat.re <- ranef(py.re)

## Now separate plug-in MLE for ehat and alphahat
py.re.e <- maxLik(twotier.ne.plugin, grad=NULL, hess=NULL, 
						     method="nm", control=list(iterlim=80000), 
					         start=c(-2, -2, -2), ehat=ehat.re)

py.re.a <- maxLik(twotier.ne.plugin, grad=NULL, hess=NULL, 
						     method="nm", control=list(iterlim=80000), 
					         start=c(-2, -2, -2), ehat=alphahat.re)

## For posterity also do COLS
ehat <- alphahat.re

k2 <- mean(ehat^2)
k3 <- mean(ehat^3)
k4 <- mean(ehat^4)-3*k2^2

fn   <- function(x, k3, k4){return(0.5*k3-((k4-6*x^4)/6)^0.75+x^3)}
#fn1 <- function(x, k3, k4){return(k4-6*(x^3+0.5*k3)^(4/3)+x^4)}

lower <- -(k4/6)^0.25+0.0000000001
upper <- (k4/6)^0.25-0.0000000001

## Solve for the 0 -- Need tryCatch here
sigu.hat <- uniroot(fn, k3=k3, k4=k4, lower=lower, upper=upper)$root 

sigw.hat <- ((k4-6*sigu.hat^4)/6)^0.25

sigv.hat <- sqrt(k2-sigw.hat^2-sigu.hat^2) ## Need tryCatch here

## Now print out table of estimates
## Have to convert estimates from likelihood function
## Ignorning fixed effects
g.d <- diag(9)
diag(g.d)[7:9] <-  exp(coefficients(py.col2)[7:9])
VCV.py.col2 <- t(g.d)%*%vcov(py.col2)%*%g.d

## Fixed effects
g.d <- diag(7)
diag(g.d)[5:7] <-  exp(coefficients(py.col3)[5:7])
VCV.py.col3 <- t(g.d)%*%vcov(py.col3)%*%g.d

## Random effects
g.d <- diag(3)
diag(g.d) <-  exp(coefficients(py.re.e))
VCV.py.re.e <- t(g.d)%*%vcov(py.re.e)%*%g.d

g.d <- diag(3)
diag(g.d) <-  exp(coefficients(py.re.a))
VCV.py.re.a <- t(g.d)%*%vcov(py.re.a)%*%g.d

ests <- matrix(NA, 24, 4)

ests[seq(from=1, to=12, by=2), 1] <- coefficients(py.col1)
ests[seq(from=1, to=18, by=2), 2] <- c(coefficients(py.col2)[1:6], exp(coefficients(py.col2)[7:9]))
ests[seq(from=5, to=18, by=2), 3] <- c(coefficients(py.col3)[1:4], exp(coefficients(py.col3)[5:7]))
ests[seq(from=1, to=12, by=2), 4] <- coefficients(py.re)
ests[seq(from=13, to=18, by=2), 4] <- exp(coefficients(py.re.e))
ests[seq(from=19, to=24, by=2), 4] <- exp(coefficients(py.re.a))

ests[seq(from=2, to=12, by=2), 1] <- summary(py.col1)$coefficients[, 2]
ests[seq(from=2, to=18, by=2), 2] <- sqrt(diag(VCV.py.col2))
ests[seq(from=6, to=18, by=2), 3] <- sqrt(diag(VCV.py.col3))
ests[seq(from=2, to=12, by=2), 4] <- summary(py.re)$coefficients[, 2]
ests[seq(from=14, to=18, by=2), 4] <- sqrt(diag(VCV.py.re.e))
ests[seq(from=20, to=24, by=2), 4] <- sqrt(diag(VCV.py.re.a))

ests <- round(ests, 4)

rownames(ests) <- c("Constant", "", "Educ", "", "Exper", "", "Exper$^{2}$", "", 
								 "Tenure", "", "Tenure${}^2}$", "", "$\\sigma_v$", "", 
								 "$\\sigma_u$", "", "$\\sigma_w$", "", "$\\sigma_{\\alpha}$", "", 
								 "$\\sigma_{\\upsilon}$", "", "$\\sigma_{\\omega}$", "")

latex(ests, file="Estimates_PY.tex")

## Alecos wants us to use E[e^(w)|epsilon] and E[e^(u)|epsilon] to assess
## information deficiency in this example. 

metrics <- function(p, e=NULL, y=NULL, xx=NULL, zu=NULL, zw=NULL, alphahat=NULL){

	if(is.null(e)){

		nr    <- ncol(xx)  ##Calculate number of regressors in regression
		nzu <- ncol(zu)  ##Calculate number of determinants for u component
		nzw <- ncol(zw) ##Calculate number of determinants for w component
	
		ep.hat  <- y-xx%*%p[1:nr]
		
		if(!is.null(alphahat)){

			ep.hat  <- y-xx%*%p[1:nr]-alpha.hat

		}


		sig.v <- exp(p[nr+1]) 	##Assume homoscedastic two sided component
  		sig.u <- exp((zu%*%p[(nr+2):(nr+nzu+1)]))
  		sig.w <- exp((zw%*%p[(nr+nzu+2):(nr+nzu+nzw+1)]))

	}else{
		
		ep.hat  <- e

		sig.v  <- exp(p[1]) 	
  		sig.u  <- exp(p[2])
  		sig.w <- exp(p[3])

		## Have to correct for the shift in our residuals to begin with
		ep.hat <- e-sig.u+sig.w

	}
	## Use 8.23 and 8.26 to construct metrics
	## Setup necessary parameters needed.
	lambda <- 1/sig.w+1/sig.u
	a1     <- sig.v^2/(2*sig.u^2)+ep.hat/sig.u
	b1     <- -(ep.hat/sig.v+sig.v/sig.u)
	a2     <- sig.v^2/(2*sig.w^2)-ep.hat/sig.w
	b2     <-  ep.hat/sig.v-sig.v/sig.w
	chi1   <- pnorm(b2)+exp(a1-a2)*pnorm(b1)
	chi2   <- exp(a2-a1)*chi1

	Eew.cond  <- (lambda/(chi2*(lambda-1)))*(pnorm(b1)+
									exp(0.5*((b2+sig.v)^2-b1^2))*pnorm(b2+sig.v))

	Eemw.cond <- (lambda/(chi2*(1+lambda)))*(pnorm(b1)+
									exp(a2-a1-b2*sig.v+0.5*sig.v^2)*pnorm(b2-sig.v))

	Eeu.cond  <- (lambda/(chi1*(lambda-1)))*(pnorm(b2)+
									exp(0.5*((b1+sig.v)^2-b2^2))*pnorm(b1+sig.v))

	Eemu.cond <- (lambda/(chi1*(1+lambda)))*(pnorm(b2)+
									exp(a1-a2-b1*sig.v+0.5*sig.v)*pnorm(b1-sig.v))

	Eewmu.cond <- (exp((1+sig.u)*(a1+sig.v^2/2/sig.u))*pnorm(b1-sig.v)+
										exp((1-sig.w)*(a2-sig.v^2/2/sig.w))*
										pnorm(b2+sig.v))/(exp(a1)*pnorm(b1)+
												exp(a2)*pnorm(b2))

	Eeumw.cond <- (exp((1-sig.u)*(a1-sig.v^2/2/sig.u))*pnorm(b1+sig.v)+
										exp((1+sig.w)*(a2+sig.v^2/2/sig.w))*
										pnorm(b2-sig.v))/(exp(a1)*pnorm(b1)+
													exp(a2)*pnorm(b2))

	## Now calculate the M1 and M2 metrics (these are information deficiency 
	## relative to the actual price) and M5 and M6 metrics (these are information 
	## deficiency relative to balanced price).
	M1.ne <- 1-Eemw.cond
	M2.ne <- Eeu.cond-1

	M5.ne <- Eew.cond-1
	M6.ne <- 1-Eemu.cond

	## Alecos wanted M7 as well
	M7.ne   <- Eewmu.cond-1
	M10.ne <- 1-Eeumw.cond

	return(list(M1=M1.ne, M2=M2.ne, M5=M5.ne, 
					M6=M6.ne, M7=M7.ne, M10=M10.ne))

}

M.col2 <- metrics(p=coefficients(py.col2), y=y, xx=xx, zu=zu, zw=zw)

alpha.hat <- YMEAN-XMEANS%*%coefficients(py.col3)[1:5]
M.col3 <- metrics(p=coefficients(py.col3), y=y, xx=xx[, 3:6], zu=zu, zw=zw, alphahat=alpha.hat)
M.col4e <- metrics(p=coefficients(py.re.e), e=as.numeric(ehat.re))
M.col4a <- metrics(p=coefficients(py.re.a), e=as.numeric(alphahat.re))

pdf("M1M2_PY.pdf")
par(mfrow=c(2, 2))

plot(density(M.col2$M1, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 1.05), main="a) Pooled",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col2$M2, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col3$M1, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 1.05), ylim=c(0, 13.2), main="b) Fixed Effects",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col3$M2, na.rm=TRUE, bw=0.1), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col4e$M1, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.57), main="c) Random Effects -- Time Varying",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col4e$M2, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col4a$M1, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.57), main="d) Random Effects -- Time Constant",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col4a$M2, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)

dev.off()

pdf("M5M6_PY.pdf")
par(mfrow=c(2,2))

plot(density(M.col2$M5, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0.1, 0.50), ylim=c(0, 18), main="a) Pooled",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col2$M6, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col3$M5, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0.0, 0.705), ylim=c(0, 8.2), main="b) Fixed Effects",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col3$M6, na.rm=TRUE, bw=0.1), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col4e$M5, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0.1, 0.5), ylim=c(0, 25), main="c) Random Effects -- Time Varying",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col4e$M6, na.rm=TRUE, bw=0.01), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)

plot(density(M.col4a$M5, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0.1, 0.5), ylim=c(0, 18), main="d) Random Effects -- Time Constant",
		xlab="Information Deficiency", ylab="Density")
lines(density(M.col4a$M6, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M7_PY.pdf")
plot(density(M.col2$M7, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-1, 1.5), ylim=c(0, 2.5))
lines(density(M.col3$M7, na.rm=TRUE), lwd=2, lty=2)
lines(density(M.col4e$M7, na.rm=TRUE), lwd=2, lty=3)
legend("topright", legend=c("Pooled", "FE", "RE"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

pdf("M10_PY.pdf")
plot(density(M.col2$M10, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-1.5, 1), ylim=c(0, 2.75))
lines(density(M.col3$M10, na.rm=TRUE), lwd=2, lty=2)
lines(density(M.col4e$M10, na.rm=TRUE), lwd=2, lty=3)
legend("topleft", legend=c("Pooled", "FE", "RE"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

## Make table that contains mean and quartiles of the various measures 
## for each of the three models
tab.Ms <- matrix(0, 12, 5)

colnames(tab.Ms) <- c("Mean", "Q$_1$", "Median", "Q$_3$", "Q$_3$-Q$_1$")

rownames(tab.Ms) <- rep(c("$\\widehat E((e^{w-u})|\\varepsilon)-1$",
									  		"1-$\\widehat E((e^{-w})|\\varepsilon))$",
									  		"1-$\\widehat E((e^{-u})|\\varepsilon))$"), 3)

## Place things in correctly
tab.Ms[1, 1]    <- mean(M.col2$M7, na.rm=TRUE)
tab.Ms[1, 2:4] <- quantile(M.col2$M7, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[1, 5]    <- IQR(M.col2$M7, na.rm=TRUE)

tab.Ms[2, 1]    <- mean(M.col2$M1, na.rm=TRUE)
tab.Ms[2, 2:4] <- quantile(M.col2$M1, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[2, 5]    <- IQR(M.col2$M1, na.rm=TRUE)

tab.Ms[3, 1]    <- mean(M.col2$M6, na.rm=TRUE)
tab.Ms[3, 2:4] <- quantile(M.col2$M6, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[3, 5]    <- IQR(M.col2$M6, na.rm=TRUE)

tab.Ms[4, 1]    <- mean(M.col3$M7, na.rm=TRUE)
tab.Ms[4, 2:4] <- quantile(M.col3$M7, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[4, 5]    <- IQR(M.col3$M7, na.rm=TRUE)

tab.Ms[5, 1]    <- mean(M.col3$M1, na.rm=TRUE)
tab.Ms[5, 2:4] <- quantile(M.col3$M1, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[5, 5]    <- IQR(M.col3$M1, na.rm=TRUE)

tab.Ms[6, 1]    <- mean(M.col3$M6, na.rm=TRUE)
tab.Ms[6, 2:4] <- quantile(M.col3$M6, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[6, 5]    <- IQR(M.col3$M6, na.rm=TRUE)

tab.Ms[7, 1]    <- mean(M.col4e$M7, na.rm=TRUE)
tab.Ms[7, 2:4] <- quantile(M.col4e$M7, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[7, 5]    <- IQR(M.col4e$M7, na.rm=TRUE)

tab.Ms[8, 1]    <- mean(M.col4e$M1, na.rm=TRUE)
tab.Ms[8, 2:4] <- quantile(M.col4e$M1, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[8, 5]    <- IQR(M.col4e$M1, na.rm=TRUE)

tab.Ms[9, 1]    <- mean(M.col4e$M6, na.rm=TRUE)
tab.Ms[9, 2:4] <- quantile(M.col4e$M6, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[9, 5]    <- IQR(M.col4e$M6, na.rm=TRUE)

tab.Ms[10, 1]    <- mean(M.col4a$M7, na.rm=TRUE)
tab.Ms[10, 2:4] <- quantile(M.col4a$M7, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[10, 5]    <- IQR(M.col4a$M7, na.rm=TRUE)

tab.Ms[11, 1]    <- mean(M.col4a$M1, na.rm=TRUE)
tab.Ms[11, 2:4] <- quantile(M.col4a$M1, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[11, 5]    <- IQR(M.col4a$M1, na.rm=TRUE)

tab.Ms[12, 1]    <- mean(M.col4a$M6, na.rm=TRUE)
tab.Ms[12, 2:4] <- quantile(M.col4a$M6, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[12, 5]    <- IQR(M.col4a$M6, na.rm=TRUE)

tab.Ms <- round(tab.Ms, 3)

latex(tab.Ms, file="Metrics_PY.tex")
