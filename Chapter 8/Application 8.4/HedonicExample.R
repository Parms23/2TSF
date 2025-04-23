## Code to produce housing market example for 
## Chapter 8 of 2Tier book with Alecos. 

## Started February 15th, 2024. 
## This version September 18th, 2024.

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
library(minpack.lm)
library(nlstools)
library(sandwich)
library(lmtest)

## Set seed for replications purposes
set.seed(123)

## Load in data
data <- read.csv("hrs.csv", h=T)

## What if we were to drop low income buyers?
##id.lowbuy <- which(data$incbuy<1000)
##id.lowsell <- which(data$incsell<1000)

##id <- c(id.lowbuy, id.lowsell)

## data <- data[-id, ]

## scale income and age by mean so that things do not
## go haywire in optimization
data$incbuy  <- as.numeric(scale(data$incbuy))
data$incsell  <- as.numeric(scale(data$incsell))
data$ageb    <- data$ageb/mean(data$ageb)
data$agen    <- data$agen/mean(data$agen)

## dependent variable
y <- data[, 1]

xx <- data[, c(2:28)]
xx <- as.matrix(cbind(1, xx))

zw <- as.matrix(cbind(1, data[, c(29, 30, 31, 33, 35, 37, 39, 41, 43, 45)]))
zu <- as.matrix(cbind(1, data[, c(32, 34, 36, 38, 40, 42, 44, 46)]))

ols <- summary(lm(y~xx-1))

nmulti <- 5

## Starting parameter vector for optimization
p.start  <- c(coefficients(ols)[, 1], 0.3, rnorm(11, sd=0.6), rnorm(9, sd=0.5))

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

twotier.nls <- function(p, y, xx, zu, zw){

	nr    <- ncol(xx)  ##Calculate number of regressors in regression
	nzu  <- ncol(zu) ##Calculate number of determinants for u component
	nzw <- ncol(zw) ##Calculate number of determinants for w component

  	sigu <- exp((zu%*%p[(nr+1):(nr+nzu)]))
  	sigw <- exp((zw%*%p[(nr+nzu+1):(nr+nzu+nzw)]))

  	#if (sigv<= 1e-6){stop("Variance too small")}

  	e <- y - xx%*%p[1:nr]+sigu-sigw

  	##return will send the summation of the log of the 
  	##density of the composed error

  	ll <- e^2

  	return(sum(ll))

}

# ## Call optmizer, this minimizes so make sure negative of 
# ## sum(ll) in likelihood function call
# opt.ne.bobyqa <- bobyqa(p.start, twotier.ne, lower=-Inf, upper=Inf, 
			 			  				 # control=list(iprint=1), y, xx, zu, zw)

# opt.hn.bobyqa <- bobyqa(p.start, twotier.hn, lower=-Inf, upper=Inf, 
			 			  				  # control=list(iprint=1), y, xx, zu, zw)

# p.start.ne <- opt.ne.bobyqa$par
# p.start.hn <- opt.hn.bobyqa$par

p.start.ne <- c(7.691500673, 0.327032453, -0.068484497, 
						0.268946072, -0.015147192, 0.555115741, 
						0.633091590, -0.056427877, 0.042110987, 
						-0.310954222, -0.270540512, 0.202080962, 
  						0.140814725, 0.055432451, -0.008153885, 
  						0.133513415, -0.009516707, 0.105384825, 
  						0.129204421, 0.211550881, 0.266841649, 
  						0.206746471, 0.247245440, 0.233805395,
  						0.616006235, 0.217288823, 0.232090989, 
  						-0.031707802, -1.354005610, -0.029753573, 
  						-0.271683377, -0.229420291, -0.578819354, 
  						0.687324666, -0.140438062, 0.083650276, 
  						-0.438388489, -0.171737240, -2.573261609, 
  						0.039926582, 0.203928115, 0.297875721, 
  						0.245013815, 0.945060586, -0.197824696, 
  						-0.072940758, 0.006549681, 0.416512464, 
  						0.084998366)

max.j <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
						   method="nm", control=list(iterlim=80000), 
					       start=p.start.ne, y=y, xx=xx, zu=zu, zw=zw)

is.max         <- max.j$maximum
coef.start1  <- max.j$estimate
which.multi <- 1

for (jjj in 1:nmulti){

	start.p1 <- p.start.ne+rnorm(length(p.start.ne), sd=0.2)

	## Estimate likelihood
	max.j1 <- maxLik(twotier.ne, grad=NULL, hess=NULL, 
							     method="bfgs", control=list(iterlim=80000), 
							     start=start.p1, y=y, xx=xx, zu=zu, zw=zw)

	if(is.null(max.j1$maximum)){max.j1$maximum <- is.max-1}

	if(max.j1$maximum>is.max){
					max.j          <- max.j1
					is.max        <- max.j$maximum
					which.multi <- jjj+1
	}
}

opt.ml.ne <- max.j

p.start.hn <- c(6.495266492, 0.463166253, -0.164418033, 
					   0.157120700, 0.004469161, 0.630115937, 
					   0.700133684, 0.049847560, 0.099392668, 
					   -0.260606663, -0.137753073, 0.301028777, 
					   0.170054183, 0.113880194, 0.043906415, 
					   0.206740764, -0.001766426, 0.335337277, 
					   0.328506182, 0.365030271, 0.460726279, 
					   0.408882475, 0.449509857, 0.429562458, 
					   0.598028422, 0.255489389, 0.227915878, 
					   -0.022496323, -1.082922364, -0.128509756, 
					   -0.246355987, -0.123166534, -0.181254582, 
					   0.284980405, -0.035690364, -0.014310092, 
					   -0.294604371, 0.013895466, -1.899516015, 
					   0.045132713, 0.038359871, 0.210145627, 
					   0.122275945, 0.581113173, -1.236936239, 
					   0.065543596, 0.180512515, 0.436741303,
					   -0.124784028)

max.j <- maxLik(twotier.hn, grad=NULL, hess=NULL, 
						   method="nm", control=list(iterlim=80000), 
						   start=p.start.hn, y=y, xx=xx, zu=zu, zw=zw)

is.max         <- max.j$maximum
coef.start1  <- max.j$estimate
which.multi <- 1

for (jjj in 1:nmulti){

	start.p1 <- p.start.hn+rnorm(length(p.start.hn), sd=0.2)

	## Estimate likelihood
	max.j1 <- maxLik(twotier.hn, grad=NULL, hess=NULL, 
							     method="bfgs", control=list(iterlim=80000), 
							     start=start.p1, y=y, xx=xx, zu=zu, zw=zw)

	if(is.null(max.j1$maximum)){max.j1$maximum <- is.max-1}

	if(max.j1$maximum>is.max){
					max.j          <- max.j1
					is.max        <- max.j$maximum
					which.multi <- jjj+1
	}
}

opt.ml.hn <- max.j

ests <- matrix(NA, 22, 6)

ests[c(1, seq(from=7, to=22, by=2)), 1] <- coefficients(opt.ml.ne)[30:38]
ests[c(1, seq(from=7, to=22, by=2)), 2] <- coefficients(opt.ml.hn)[30:38]

ests[seq(from=1, to=22, by=2), 4] <- coefficients(opt.ml.ne)[39:49]
ests[seq(from=1, to=22, by=2), 5] <- coefficients(opt.ml.hn)[39:49]

ests[c(2, seq(from=8, to=22, by=2)), 1] <- summary(opt.ml.ne)$estimate[30:38, 2]
ests[c(2, seq(from=8, to=22, by=2)), 2] <- summary(opt.ml.hn)$estimate[30:38, 2]

ests[seq(from=2, to=22, by=2), 4] <- summary(opt.ml.ne)$estimate[39:49, 2]
ests[seq(from=2, to=22, by=2), 5] <- summary(opt.ml.hn)$estimate[39:49, 2]

## Now use NLS
opt.nls <- optim(coefficients(opt.ml.ne)[-29], fn=twotier.nls, 
						  y=y, xx=xx, zu=zu, zw=zw, 
						  hessian=TRUE, 
						  control=list(maxit=100000))

# opt.ne.bobyqa <- bobyqa(opt.nls$par, twotier.nls, lower=-Inf, upper=Inf, 
			 			  				 # control=list(iprint=1,maxfun=1000000), y, xx, zu, zw)

# opt.nls <- optim(opt.ne.bobyqa$par, fn=twotier.nls, 
						  # y=y, xx=xx, zu=zu, zw=zw, 
						  # hessian=TRUE, 
						  # control=list(maxit=100000))
					  
# lambda <- 0.000001
# step <- lambda
# warn.loop <- TRUE

# while(warn.loop){
	
	# result <- tryCatch(sqrt(diag(solve(opt.nls$hessian+
															# lambda*diag(nrow(opt.nls$hessian))))),
								 # warning = function(e){"Warning"})
	
	# if(result=="Warning"){
		# lambda <- lambda+step
		# warn.loop <- TRUE}
	
	# if(result!="Warning"){warn.loop <- FALSE}
	
# }

nls.ests <- cbind(opt.nls$par, sqrt(diag(solve(opt.nls$hessian))))

gradient.nls <- numDeriv::jacobian(twotier.nls, opt.nls$par, y=y, xx=xx, zu=zu, zw=zw)

opg.nls <- t(gradient.nls)%*%gradient.nls

mod.ols <- lm(lprn~lsf+unitsftc+bathstot+roomsn+sfan+
							   sfdn+cencityn+urbsubn+urbann+
							   riuraln+agelt5+age510+age1015+
                 			   agegte30+inadeq+degreen+
             				   s87+s88+s89+s90+s91+
             				   s92+s93+verylg+large+siz1to3+small, 
             			data=data)

# mod.nls <- nls(lprn~b.0+b.1*lsf+b.2*unitsftc+b.3*bathstot+b.4*roomsn+b.5*sfan+
							   # b.6*sfdn+b.7*cencityn+b.8*urbsubn+b.9*urbann+
							   # b.10*riuraln+b.11*agelt5+b.12*age510+b.13*age1015+
                 			   # b.14*agegte30+b.15*inadeq+b.16*degreen+
             				   # b.17*s87+b.18*s88+b.19*s89+b.20*s90+b.21*s91+
             				   # b.22*s92+b.23*s93+b.24*verylg+b.25*large+
             				   # b.26*siz1to3+b.27*small, data=data, 
             				   # start=list(b.0=0, b.1=0, b.2=0 , b.3=0, b.4=0, b.5=0, 
             				   				  # b.6=0, b.7=0, b.8=0 , b.9=0, b.10=0, b.11=0, 
             				   				  # b.12=0, b.13=0, b.14=0, b.15=0, b.16=0, 
             				   				  # b.17=0, b.18=0, b.19=0, b.20=0, b.21=0, 
             				   				  # b.22=0, b.23=0, b.24=0, b.25=0, b.26=0, 
             				   				  # b.27=0))
                 
# ## Now nonlinear least squares
# start.nls=list(b.0=0, b.1=0, b.2=0, b.3=0, b.4=0, b.5=0, 
				    # b.6=0, b.7=0, b.8=0, b.9=0, b.10=0, b.11=0, 
				    # b.12=0, b.13=0, b.14=0, b.15=0, b.16=0, 
				    # b.17=0, b.18=0, b.19=0, b.20=0, b.21=0, 
				    # b.22=0, b.23=0, b.24=0, b.25=0, b.26=0, 
				    # b.27=0, d.0=log(0.75), d.1=0, d.2=0, d.3=0, 
				    # d.4=0, d.5=0, d.6=0, d.7=0, d.8=0, d.9=0, 
				    # d.10=0, e.0=log(0.1), e.1=0, e.2=0, e.3=0, 
				    # e.4=0, e.5=0, e.6=0, e.7=0, e.8=0)

# min.j <- nlsLM(lprn~b.0+b.1*lsf+b.2*unitsftc+b.3*bathstot+b.4*roomsn+b.5*sfan+
							   # b.6*sfdn+b.7*cencityn+b.8*urbsubn+b.9*urbann+
							   # b.10*riuraln+b.11*agelt5+b.12*age510+b.13*age1015+
							   # b.14*agegte30+b.15*inadeq+b.16*degreen+b.17*s87+
							   # b.18*s88+b.19*s89+b.20*s90+b.21*s91+b.22*s92+
							   # b.23*s93+b.24*verylg+b.25*large+b.26*siz1to3+
							   # b.27*small-exp(d.0+d.1*prvlocn+d.2*firbuy+d.3*incbuy+
							   						    # d.4*qbusb+d.5*ageb+d.6*blkbuy+
							   						    # d.7*marbuy+d.8*sfbuy+d.9*edubuy+
							   						    # d.10*kidbuy)+
							   # exp(e.0+e.1*incsell+e.2*qbusn+e.3*agen+e.4*blksell+
							  		 		# e.5*marsell+e.6*sfsell+e.7*edusell+e.8*kidsell),
					 # data=data, start=start.nls, subset=(cenind==0),
					 # control=list(maxiter=1000, maxfev=10000))

# is.min         <- sum((summary(min.j)$residuals)^2)
# coef.start1  <- min.j$estimate
# which.multi <- 1
# nls1 <- min.j

# for (jjj in 1:nmulti){

	# new.nls.p <- start.nls

	# for(k in 1:length(start.nls)){
			# new.nls.p[[k]] <- start.nls[[k]]+rnorm(1, sd=0.2)
	# }

	# ## Estimate NLS
	# min.j1 <- nlsLM(lprn~b.0+b.1*lsf+b.2*unitsftc+b.3*bathstot+b.4*roomsn+
									  # b.5*sfan+b.6*sfdn+b.7*cencityn+b.8*urbsubn+
									  # b.9*urbann+b.10*riuraln+b.11*agelt5+b.12*age510+
									  # b.13*age1015+b.14*agegte30+b.15*inadeq+
									  # b.16*degreen+b.17*s87+b.18*s88+b.19*s89+
									  # b.20*s90+b.21*s91+b.22*s92+b.23*s93+
									  # b.24*verylg+b.25*large+b.26*siz1to3+
							          # b.27*small-exp(d.0+d.1*prvlocn+d.2*firbuy+d.3*incbuy+
							   						    	   # d.4*qbusb+d.5*ageb+d.6*blkbuy+
							   						    	   # d.7*marbuy+d.8*sfbuy+d.9*edubuy+
							   						    	   # d.10*kidbuy)+
							   		  # exp(e.0+e.1*incsell+e.2*qbusn+e.3*agen+e.4*blksell+
							  		 		 # e.5*marsell+e.6*sfsell+e.7*edusell+e.8*kidsell),
					 			# data=data, start=new.nls.p, subset=(cenind==0),
					 			# control=list(maxiter=1000, maxfev=10000))

	# is.min.j1 <- sum((summary(min.j1)$residuals)^2)

	# if(is.min.j1<is.min){
		# is.min         <- is.min.j1
		# nls1            <- min.j1
		# which.multi <- jjj+1
	# }
# }

# nls.use <- c(b.0 = -8.974583e+01, b.1 = 3.149106e-01, 
				    # b.2 = -1.752722e-01, b.3 = 1.763843e-01,
				    # b.4 = 4.552949e-04, b.5 = 7.617212e-01, 
				    # b.6 = 8.322029e-01, b.7 = -7.588971e-02, 
				    # b.8 = -1.037466e-03, b.9 = -2.214382e-01, 
				    # b.10 = -2.645359e-01, b.11 = 2.592629e-01, 
				    # b.12 = 1.256389e-01, b.13 = 5.440331e-02, 
				    # b.14 = -3.117796e-02, b.15 = 4.067217e-03, 
				    # b.16 = -1.419214e-02, b.17 = 7.036677e-03, 
				    # b.18 = 2.165316e-02, b.19 = 6.272671e-02, 
				    # b.20 = 1.448009e-01, b.21 = 9.625930e-02, 
				    # b.22 = 1.087437e-01, b.23 = 1.126926e-01, 
				    # b.24 = 5.513818e-01, b.25 = 1.282240e-01, 
				    # b.26 = 1.768664e-01, b.27 = -1.028090e-01, 
				    # d.0 = -2.129302e-01, d.1 = -1.459361e-01, 
				    # d.2 = 3.820380e-02, d.3 = -1.441480e-01,
				    # d.4 = 3.854125e-02, d.5 = -3.627852e-01, 
				    # d.6 = -3.423160e-01, d.7 = -1.865434e-01, 
				    # d.8 = -2.135687e-01, d.9 = 1.229494e-01, 
				    # d.10 = -1.034234e-01, e.0 = 4.581809e+00, 
				    # e.1 = 1.908382e-04, e.2 = -3.166285e-04, 
				    # e.3 = 9.221284e-04, e.4 = 1.081247e-03, 
				    # e.5 = 1.683435e-03, e.6 = 2.824576e-03, 
				    # e.7 = -1.440893e-03, e.8 = -1.173438e-04)

# nls.use.ml.ne <- coefficients(opt.ml.ne)[-29]
# names(nls.use.ml.ne) <- names(nls.use)

# min.j <- nlsLM(lprn~b.0+b.1*lsf+b.2*unitsftc+b.3*bathstot+b.4*roomsn+
							    # b.5*sfan+b.6*sfdn+b.7*cencityn+b.8*urbsubn+
								# b.9*urbann+b.10*riuraln+b.11*agelt5+b.12*age510+
								# b.13*age1015+b.14*agegte30+b.15*inadeq+
								# b.16*degreen+b.17*s87+b.18*s88+b.19*s89+
								# b.20*s90+b.21*s91+b.22*s92+b.23*s93+
								# b.24*verylg+b.25*large+b.26*siz1to3+
							    # b.27*small-exp(d.0+d.1*prvlocn+d.2*firbuy+d.3*incbuy+
							   						     # d.4*qbusb+d.5*ageb+d.6*blkbuy+
							   						     # d.7*marbuy+d.8*sfbuy+d.9*edubuy+
							   						     # d.10*kidbuy)+
							    # exp(e.0+e.1*incsell+e.2*qbusn+e.3*agen+e.4*blksell+
							    	   # e.5*marsell+e.6*sfsell+e.7*edusell+e.8*kidsell),
					 	# data=data, start=nls.use.ml.ne, subset=(cenind==0),
					 	# control=list(maxiter=1000, maxfev=10000))

# nls1 <- min.j

# stt <- coeftest(nls1, vcov=sandwich)

ests[c(1, seq(from=7, to=22, by=2)), 3] <- nls.ests[29:37, 1]
ests[seq(from=1, to=22, by=2), 6] <- nls.ests[38:48, 1]

ests[c(2, seq(from=8, to=22, by=2)), 3] <- nls.ests[29:37, 2]
ests[seq(from=2, to=22, by=2), 6] <- nls.ests[38:48, 2]

rownames(ests) <- c("Constant", "", "Out of Town", "", 
								  "First Time", "", "Income", "", 
								  "Own Business",  "", "Age", "",
								  "Black", "", "Married", "", 
								  "Single Female", "",
								  "Education", "", "Kids", "")

# Function to add parentheses to specified rows of a matrix
# Updated function to add parentheses to specified rows of a matrix
add_parentheses <- function(mat, rows_to_parenthesize) {
  mat <- as.matrix(mat)  # Ensure the input is a matrix
  
  # Loop through specified rows and add parentheses only to numeric values
  for (row in rows_to_parenthesize) {
    for (col in seq_len(ncol(mat))) {
      # Only add parentheses if the value is not NA, NaN, or blank
      if (!is.na(mat[row, col]) && !is.nan(mat[row, col]) && mat[row, col] != "") {
        mat[row, col] <- paste0("(", mat[row, col], ")")
      }
    }
  }
  
  return(noquote(mat))
}

# Apply rounding while preserving NA and NaN
ests.rounded <- ifelse(is.na(ests) | is.nan(ests), ests, round(ests, 3))

ests.rounded <- add_parentheses(ests.rounded, seq(from=2, to=22, by=2))

latex(ests.rounded, file="Errors_Hedonic_KP.tex")

## Also make table of coefficient estimates to send to latex as well
ols.ne <- mod.ols
ols.hn <- mod.ols
ols.nls <- mod.ols

se.ols <- sqrt(diag(vcov(mod.ols)))
se.ne  <- sqrt(diag(vcov(opt.ml.ne)))
se.hn  <- sqrt(diag(vcov(opt.ml.hn)))
se.nls <- nls.ests[,2]

stargazer(list(mod.ols, ols.ne, ols.hn, ols.nls),
          		coef=list(coefficients(mod.ols), 
                    		      coefficients(opt.ml.ne)[1:28], 
                    		      coefficients(opt.ml.hn)[1:28], 
                    		      nls.ests[1:28,1]),
          		se=list(se.ols, se.ne[1:28], 
          				   se.hn[1:28], se.nls[1:28]),
          		label="tab:betas_Hedonic_KP",
          		star.cutoffs = c(0),
          		out="betas_Hedonic_KP.tex",
          		no.space=T, digits=3, single.row=F,
          		omit.stat = c("rsq","adj.rsq", "f", "n", "ser"),
          		intercept.bottom=FALSE,
          		dep.var.labels.include = FALSE,
          		title="Two-Tier Stochastic Frontier Estimates",
          		column.labels=c("OLS", "Normal-Exponential", 
          							 		"Normal-Half Normal", "NLS")
)

## Alecos wants us to use E[e^(w)|epsilon] and E[e^(u)|epsilon] to assess
## information deficiency in this example. 

k.x <- dim(xx)[2]
k.u <- dim(zu)[2]
k.w <- dim(zw)[2]

## So strip off parameter estimates and create epsilon series.
## Normal-Exponential First
par.hat <- coefficients(opt.ml.ne)   ##opt$par

ep.hat  <- y-xx%*%par.hat[1:k.x]

sig.v  <- exp(par.hat[k.x+1])
sig.u  <- exp(zu%*%par.hat[(k.x+2):(k.x+k.u+1)])
sig.w <- exp(zw%*%par.hat[(k.x+k.u+2):(k.x+k.u+k.w+1)])

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
					       exp((1-sig.w)*(a2-sig.v^2/2/sig.w))*pnorm(b2+sig.v))/
					       (exp(a1)*pnorm(b1)+exp(a2)*pnorm(b2))

Eeumw.cond <- (exp((1-sig.u)*(a1-sig.v^2/2/sig.u))*pnorm(b1+sig.v)+
						  exp((1+sig.w)*(a2+sig.v^2/2/sig.w))*pnorm(b2-sig.v))/
						  (exp(a1)*pnorm(b1)+exp(a2)*pnorm(b2))

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

quantile(M1.ne, c(0.25, 0.5, 0.75), na.rm=TRUE)
quantile(M2.ne, c(0.25, 0.5, 0.75), na.rm=TRUE)
quantile(M5.ne, c(0.25, 0.5, 0.75), na.rm=TRUE)
quantile(M6.ne, c(0.25, 0.5, 0.75), na.rm=TRUE) 
quantile(M7.ne, c(0.25, 0.5, 0.75), na.rm=TRUE)
quantile(M10.ne, c(0.25, 0.5, 0.75), na.rm=TRUE)
mean(M7.ne, na.rm=TRUE)
mean(M10.ne, na.rm=TRUE)

## Now Normal-Half Normal
## So strip off parameter estimates and create epsilon series.
par.hat <- coefficients(opt.ml.hn)   ##opt$par
ep.hat  <- y-xx%*%par.hat[1:k.x]

sig.v  <- exp(par.hat[k.x+1])
sig.u  <- exp(zu%*%par.hat[(k.x+2):(k.x+k.u+1)])
sig.w <- exp(zw%*%par.hat[(k.x+k.u+2):(k.x+k.u+k.w+1)])

## Use 7.42 and 7.43 to construct metrics
## Setup necessary parameters needed.
theta1 <- sig.w/sig.v
theta2 <- sig.u/sig.v
s <- sqrt(sig.v^2+sig.u^2+sig.w^2)
omega1 <- s*sqrt(1+theta2^2)/theta1
omega2 <- s*sqrt(1+theta1^2)/theta2
lambda1 <- (theta2/theta1)*sqrt(1+theta1^2+theta2^2)
lambda2 <- (theta1/theta2)*sqrt(1+theta1^2+theta2^2)

Di <- pbinorm(ep.hat/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))-
			pbinorm(ep.hat/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2))

f1i <- (2/omega1)*dnorm(ep.hat/omega1)*pnorm(-lambda1*ep.hat/omega1)
f2i <- (2/omega2)*dnorm(ep.hat/omega2)*pnorm(lambda2*ep.hat/omega2)

F1i <- 2*pbinorm(ep.hat/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))
F2i <- 2*pbinorm(ep.hat/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2))

psi1 <- f1i/(F1i-F2i)
psi2 <- f2i/(F1i-F2i)
psii <- (f1i-f2i)/(F1i-F2i)

s1 <- sqrt(sig.v^2+sig.w^2)
s2 <- sqrt(sig.v^2+sig.u^2)

omega.w <- sig.w*s2/s
omega.u <- sig.u*s1/s

Eew.cond  <- 2*(F1i-F2i)^(-1)*exp(0.5*omega.w^2+(omega.w/omega1)*ep.hat)*(pnorm(-(ep.hat-sig.w^2)/omega2)-pbinorm(-(ep.hat-sig.w^2)/omega2, -(omega.w+(ep.hat/omega1)), cov12=-sig.w*sig.u/s1/s2))

Eemw.cond <- 2*(F1i-F2i)^(-1)*exp(0.5*omega.w^2-(omega.w/omega1)*ep.hat)*(pnorm(-(ep.hat+sig.w^2)/omega2)-pbinorm(-(ep.hat+sig.w^2)/omega2, (omega.w-(ep.hat/omega1)), cov12=-sig.w*sig.u/s1/s2))

Eeu.cond  <- 2*(F1i-F2i)^(-1)*exp(0.5*omega.u^2-(omega.u/omega2)*ep.hat)*(pnorm((ep.hat+sig.u^2)/omega1)-pbinorm((ep.hat+sig.u^2)/omega1, ((ep.hat/omega2)-omega.u), cov12=-sig.w*sig.u/s1/s2))

Eemu.cond <- 2*(F1i-F2i)^(-1)*exp(0.5*omega.u^2+(omega.u/omega2)*ep.hat)*(pnorm((ep.hat-sig.u^2)/omega1)-pbinorm((ep.hat-sig.u^2)/omega1, ((ep.hat/omega2)+omega.u), cov12=-sig.w*sig.u/s1/s2))

Eewmu.cond <- exp(((sig.w^2+sig.u^2)/s^2)*(ep.hat+0.5*sig.v^2))*(pbinorm((ep.hat+sig.v^2)/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))-pbinorm((ep.hat+sig.v^2)/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2)))/Di

Eeumw.cond <- exp(((sig.w^2+sig.u^2)/s^2)*(0.5*sig.v^2-ep.hat))*(pbinorm((ep.hat-sig.v^2)/omega1, 0, cov12=lambda1/sqrt(1+lambda1^2))-pbinorm((ep.hat-sig.v^2)/omega2, 0, cov12=-lambda2/sqrt(1+lambda2^2)))/Di

## Now calculate the M1 and M2 metrics (these are information deficiency 
## relative to the actual price) and M5 and M6 metrics (these are information 
## deficiency relative to balanced price).
M1.hn <- 1-Eemw.cond
M2.hn <- Eeu.cond-1

M5.hn <- Eew.cond-1
M6.hn <- 1-Eemu.cond

## Alecos wanted M7 as well
M7.hn   <- Eewmu.cond-1
M10.hn <- 1-Eeumw.cond

quantile(M1.hn, c(0.25, 0.5, 0.75))
quantile(M2.hn, c(0.25, 0.5, 0.75))
quantile(M5.hn, c(0.25, 0.5, 0.75))
quantile(M6.hn, c(0.25, 0.5, 0.75)) 
quantile(M7.hn, c(0.25, 0.5, 0.75))
quantile(M10.hn, c(0.25, 0.5, 0.75))
mean(M7.hn)
mean(M10.hn)

## Now NLS estimates
par.nls.u <- nls.ests[(k.x+1):(k.x+k.u), 1]
par.nls.w <- nls.ests[(k.x+k.u+1):(k.x+k.u+k.w), 1]

w.hat <- exp(zw%*%par.nls.w)
u.hat <- exp(zu%*%par.nls.u)

## Now calculate the M1 and M2 metrics (these are information deficiency 
## relative to the actual price) and M5 and M6 metrics (these are information 
## deficiency relative to balanced price).
M1.nls <- 1-exp(-w.hat)
M2.nls <- exp(u.hat)-1

M5.nls <- exp(w.hat)-1
M6.nls <- 1-exp(-u.hat)

## Alecos wanted M7 as well
M7.nls   <- exp(w.hat-u.hat)-1
M10.nls <- 1-exp(u.hat-w.hat)

quantile(M1.nls, c(0.25, 0.5, 0.75))
quantile(M2.nls, c(0.25, 0.5, 0.75))
quantile(M5.nls, c(0.25, 0.5, 0.75))
quantile(M6.nls, c(0.25, 0.5, 0.75)) 
quantile(M7.nls, c(0.25, 0.5, 0.75))
quantile(M10.nls, c(0.25, 0.5, 0.75))
mean(M7.nls)
mean(M10.nls)

pdf("M1M2_hedonic_NE.pdf")
plot(density(M1.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 1.5), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M2.ne, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M1M2_hedonic_HN.pdf")
plot(density(M1.hn, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 1.35), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M2.hn, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M1M2_hedonic_NLS.pdf")
plot(density(M1.nls, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0.0, 1.75), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M2.nls, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($1-E\[e^{-w}|\epsilon\]$)'), 
 										     TeX(r'($E\[e^u|\epsilon\]-1$)')),
  			lwd=rep(2, 2), lty=1:2)
dev.off()

pdf("M5M6_hedonic_NE.pdf")
plot(density(M5.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.85), ylim=c(0, 6.5), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M6.ne, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M5M6_hedonic_HN.pdf")
plot(density(M5.hn, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 1.05),  main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M6.hn, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2,2), lty=1:2)
dev.off()

pdf("M5M6_hedonic_NLS.pdf")
plot(density(M5.nls, na.rm=TRUE), type="l", lwd=2, lty=1, 
		xlim= c(0, 0.85),  ylim=c(0, 5.7), main="",
		xlab="Information Deficiency", ylab="Density")
lines(density(M6.nls, na.rm=TRUE), lwd=2, lty=2)
legend("topright", legend=c(TeX(r'($E\[e^w|\epsilon\]-1$)'), 
 										     TeX(r'($1-E\[e^{-u}|\epsilon\]$)')),
  			lwd=rep(2, 2), lty=1:2)
dev.off()

pdf("M7_hedonic.pdf")
plot(density(M7.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-1, 3.75), ylim=c(0, 2.85))
lines(density(M7.hn, na.rm=TRUE), lwd=2, lty=2)
lines(density(M7.nls, na.rm=TRUE), lwd=2, lty=3)
legend("topright", legend=c("Normal-Exponential",
											 "Normal-Half Normal",
											 "NLS"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

pdf("M10_hedonic.pdf")
plot(density(M10.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-2, 1), ylim=c(0, 1.25))
lines(density(M10.hn, na.rm=TRUE), lwd=2, lty=2)
lines(density(M10.nls, na.rm=TRUE), lwd=2, lty=3)
legend("topleft", legend=c("Normal-Exponential",
										  "Normal-Half Normal",
										  "NLS"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

## Calculate Vuong (1989) test between NE and NHN 2TSF models
## Estimate both NE and HN 2TSF models
n <- dim(xx)[1]
li.NE <- twotier.ne(coefficients(opt.ml.ne), y, xx, zu, zw)
li.HN <- twotier.hn(coefficients(opt.ml.hn), y, xx, zu, zw)

li.diff <- li.NE-li.HN

LRn              <- sum(li.diff)
omegahat.n <- mean(li.diff^2)-mean(li.diff)^2
Tn                <- sqrt(1/n)*LRn/sqrt(omegahat.n)

## Calculate c and -c from standard normal for given significance level
alpha <- 0.05
c.max <- qnorm(1-alpha/2)
c.min <- qnorm(alpha/2)

Tn
c.max
c.min

## Make table that contains mean and quartiles of the various measures 
## for each of the three models
tab.Ms <- matrix(0, 9, 5)

colnames(tab.Ms) <- c("Mean", "Q$_1$", "Median", "Q$_3$", "Q$_3$-Q$_1$")

rownames(tab.Ms) <- c(rep("$\\widehat E(e^{w-u}|\\varepsilon)-1$", 3),
									  rep("1-$\\widehat E(e^{-w}|\\varepsilon)$",3),
									  rep("1-$\\widehat E(e^{-u}|\\varepsilon)$", 3))

## Place things in correctly
tab.Ms[1, 1]    <- mean(M7.ne, na.rm=TRUE)
tab.Ms[1, 2:4] <- quantile(M7.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[1, 5]    <- IQR(M7.ne, na.rm=TRUE)

tab.Ms[2, 1]    <- mean(M7.hn, na.rm=TRUE)
tab.Ms[2, 2:4] <- quantile(M7.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[2, 5]    <- IQR(M7.hn, na.rm=TRUE)

tab.Ms[3, 1]    <- mean(M7.nls, na.rm=TRUE)
tab.Ms[3, 2:4] <- quantile(M7.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[3, 5]    <- IQR(M7.nls, na.rm=TRUE)

tab.Ms[4, 1]    <- mean(M1.ne, na.rm=TRUE)
tab.Ms[4, 2:4] <- quantile(M1.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[4, 5]    <- IQR(M1.ne, na.rm=TRUE)

tab.Ms[5, 1]    <- mean(M1.hn, na.rm=TRUE)
tab.Ms[5, 2:4] <- quantile(M1.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[5, 5]    <- IQR(M1.hn, na.rm=TRUE)

tab.Ms[6, 1]    <- mean(M1.nls, na.rm=TRUE)
tab.Ms[6, 2:4] <- quantile(M1.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[6, 5]    <- IQR(M1.nls, na.rm=TRUE)

tab.Ms[7, 1]    <- mean(M6.ne, na.rm=TRUE)
tab.Ms[7, 2:4] <- quantile(M6.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[7, 5]    <- IQR(M6.ne, na.rm=TRUE)

tab.Ms[8, 1]    <- mean(M6.hn, na.rm=TRUE)
tab.Ms[8, 2:4] <- quantile(M6.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[8, 5]    <- IQR(M6.hn, na.rm=TRUE)

tab.Ms[9, 1]    <- mean(M6.nls, na.rm=TRUE)
tab.Ms[9, 2:4] <- quantile(M6.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[9, 5]    <- IQR(M6.nls, na.rm=TRUE)

tab.Ms <- round(tab.Ms, 3)

latex(tab.Ms, file="Metrics_Hedonic_KP.tex")

## Now adding in code for three things request by Alecos
## 1. Load in NLS3stages and plot M7 and M10 on the graphs 
## 2. Spit out M2/M5 metrics for NE and NHN to .csv file
## 3. Convert NH and NE results to conditional on z instead of e.

## Recreate Series
## So strip off parameter estimates and create epsilon series.
## Normal-Exponential First
par.hat.ne <- coefficients(opt.ml.ne)   ##opt$par

sig.u.ne  <- exp(zu%*%par.hat.ne[(k.x+2):(k.x+k.u+1)])
sig.w.ne <- exp(zw%*%par.hat.ne[(k.x+k.u+2):(k.x+k.u+k.w+1)])

par.hat.hn <- coefficients(opt.ml.hn)   ##opt$par

sig.u.hn  <- exp(zu%*%par.hat.hn[(k.x+2):(k.x+k.u+1)])
sig.w.hn <- exp(zw%*%par.hat.hn[(k.x+k.u+2):(k.x+k.u+k.w+1)])

## Now calculate new metrics
## Normal Exponential
M1.ne   <- 1-1/(1+sig.w.ne)
M2.ne   <- 1/(1-sig.u.ne)-1
M5.ne   <- 1/(1-sig.w.ne)-1
M6.ne   <- 1-1/(1+sig.u.ne)
M7.ne   <- 1/((1-sig.w.ne)*(1+sig.u.ne))-1
M10.ne <- 1-1/((1+sig.w.ne)*(1-sig.u.ne))

## Normal Half Normal
scale <- sqrt(2/pi)
M1.hn   <- 1-scale*pnorm(-sig.w.hn)/dnorm(sig.w.hn)
M2.hn   <- scale*pnorm(sig.u.hn)/dnorm(sig.u.hn)-1
M5.hn   <- scale*pnorm(sig.w.hn)/dnorm(sig.w.hn)-1
M6.hn   <- 1-scale*pnorm(-sig.u.hn)/dnorm(sig.u.hn)
M7.hn   <- (scale^2)*(pnorm(sig.w.hn)/dnorm(sig.w.hn))*(pnorm(-sig.u.hn)/dnorm(sig.u.hn))-1
M10.hn <- 1-(scale^2)*(pnorm(-sig.w.hn)/dnorm(sig.w.hn))*(pnorm(sig.u.hn)/dnorm(sig.u.hn))

## NLS
par.nls.u <- nls.ests[(k.x+1):(k.x+k.u), 1]
par.nls.w <- nls.ests[(k.x+k.u+1):(k.x+k.u+k.w), 1]

w.hat <- exp(zw%*%par.nls.w)
u.hat <- exp(zu%*%par.nls.u)

## Now calculate the M1 and M2 metrics (these are information deficiency 
## relative to the actual price) and M5 and M6 metrics (these are information 
## deficiency relative to balanced price).
M1.nls <- 1-exp(-w.hat)
M2.nls <- exp(u.hat)-1

M5.nls <- exp(w.hat)-1
M6.nls <- 1-exp(-u.hat)

## Alecos wanted M7 as well
M7.nls   <- exp(w.hat-u.hat)-1
M10.nls <- 1-exp(u.hat-w.hat)

## We do #2 first since it is the simplest. 
data.4ap <- data.frame(M1.ne, M1.hn, M1.nls, M2.ne, M2.hn, M2.nls, 
									 M5.ne, M5.hn, M5.nls, M6.ne, M6.hn, M6.nls, 
									 M7.ne, M7.hn, M7.nls, M10.ne, M10.hn, M10.nls)
write.csv(data.4ap, file="Mseries.csv")

## Now we redo Figures for M7 and M10 with Alecos' 3 stage NLS
library(readxl)
ap.data <- read_excel("NLS3stages.xlsx", sheet="M7-M10")

## Now redo figures 
pdf("M7_hedonic_wAP.pdf")
plot(density(M7.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-1, 1.75), ylim=c(0, 3.85))
lines(density(M7.hn, na.rm=TRUE), lwd=2, lty=2)
lines(density(M7.nls, na.rm=TRUE), lwd=2, lty=3, col="blue")
lines(density(ap.data$M7_cnls, na.rm=TRUE), lwd=2, lty=4, col="red")
legend("topright", legend=c("Normal-Exponential",
											 "Normal-Half Normal",
											 "NLS", "CNLS"),
  			lwd=rep(2, 4), lty=1:4, col=c("black", "black", "blue", "red"))
dev.off()

pdf("M10_hedonic_wAP.pdf")
plot(density(M10.ne, na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-2, 1), ylim=c(0, 1.75))
lines(density(M10.hn, na.rm=TRUE), lwd=2, lty=2)
lines(density(M10.nls, na.rm=TRUE), lwd=2, lty=3, col="blue")
lines(density(ap.data$M10_cnls, na.rm=TRUE), lwd=2, lty=4, col="red")
legend("topleft", legend=c("Normal-Exponential",
										  "Normal-Half Normal",
										  "NLS", "CNLS"),
  			lwd=rep(2, 4), lty=1:4, col=c("black", "black", "blue", "red"))
dev.off()

data <- read.csv(file="Mseries.csv", h=T)
attach(data)
all.M7 <- c(M7.ne, M7.hn, M7.nls)
all.M10 <- c(M10.ne, M10.hn, M10.nls)

id.ne <- which(M7.ne<5 & M7.ne>-2)
id.hn <- which(M7.hn<5 & M7.hn>-2)
id.nls <- which(M7.nls<5 & M7.nls>-2)

id.M7 <- which(M5.ne<0 | M5.hn<0 | M5.nls<0)

pdf("M7_hedonic.pdf")
plot(density(M7.ne[-id.M7], na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-0.75, 2.05), ylim=c(0, 2.75))
lines(density(M7.hn[-id.M7], na.rm=TRUE), lwd=2, lty=2)
lines(density(M7.nls[-id.M7], na.rm=TRUE), lwd=2, lty=3)
legend("topright", legend=c("Normal-Exponential",
											 "Normal-Half Normal",
											 "NLS"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

id.ne <- which(M10.ne<1.4 & M10.ne>-5)
id.hn <- which(M10.hn<1.4 & M10.hn>-5)
id.nls <- which(M10.nls<1.4 & M10.nls>-5)

id.M10 <- which(M2.ne<0 | M2.hn<0 | M2.nls<0)

pdf("M10_hedonic.pdf")
plot(density(M10.ne[-id.M10], na.rm=TRUE), type="l", lwd=2, lty=1, 
	    main="", xlab="Information Deficiency", 
	    ylab="Density", xlim=c(-2, 0.75), ylim=c(0, 1.35))
lines(density(M10.hn[-id.M10], na.rm=TRUE), lwd=2, lty=2)
lines(density(M10.nls[-id.M10], na.rm=TRUE), lwd=2, lty=3)
legend("topleft", legend=c("Normal-Exponential",
										  "Normal-Half Normal",
										  "NLS"),
  			lwd=rep(2, 3), lty=1:3)
dev.off()

## Make table that contains mean and quartiles of the various measures 
## for each of the three models
tab.Ms <- matrix(0, 9, 5)

colnames(tab.Ms) <- c("Mean", "Q$_1$", "Median", "Q$_3$", "Q$_3$-Q$_1$")

rownames(tab.Ms) <- c(rep("$\\widehat E(e^{w-u}|\\varepsilon)-1$", 3),
									  rep("1-$\\widehat E(e^{-w}|\\varepsilon)$",3),
									  rep("1-$\\widehat E(e^{-u}|\\varepsilon)$", 3))

## Place things in correctly
tab.Ms[1, 1]    <- mean(M7.ne, na.rm=TRUE)
tab.Ms[1, 2:4] <- quantile(M7.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[1, 5]    <- IQR(M7.ne, na.rm=TRUE)

tab.Ms[2, 1]    <- mean(M7.hn, na.rm=TRUE)
tab.Ms[2, 2:4] <- quantile(M7.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[2, 5]    <- IQR(M7.hn, na.rm=TRUE)

tab.Ms[3, 1]    <- mean(M7.nls, na.rm=TRUE)
tab.Ms[3, 2:4] <- quantile(M7.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[3, 5]    <- IQR(M7.nls, na.rm=TRUE)

tab.Ms[4, 1]    <- mean(M1.ne, na.rm=TRUE)
tab.Ms[4, 2:4] <- quantile(M1.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[4, 5]    <- IQR(M1.ne, na.rm=TRUE)

tab.Ms[5, 1]    <- mean(M1.hn, na.rm=TRUE)
tab.Ms[5, 2:4] <- quantile(M1.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[5, 5]    <- IQR(M1.hn, na.rm=TRUE)

tab.Ms[6, 1]    <- mean(M1.nls, na.rm=TRUE)
tab.Ms[6, 2:4] <- quantile(M1.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[6, 5]    <- IQR(M1.nls, na.rm=TRUE)

tab.Ms[7, 1]    <- mean(M6.ne, na.rm=TRUE)
tab.Ms[7, 2:4] <- quantile(M6.ne, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[7, 5]    <- IQR(M6.ne, na.rm=TRUE)

tab.Ms[8, 1]    <- mean(M6.hn, na.rm=TRUE)
tab.Ms[8, 2:4] <- quantile(M6.hn, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[8, 5]    <- IQR(M6.hn, na.rm=TRUE)

tab.Ms[9, 1]    <- mean(M6.nls, na.rm=TRUE)
tab.Ms[9, 2:4] <- quantile(M6.nls, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
tab.Ms[9, 5]    <- IQR(M6.nls, na.rm=TRUE)

tab.Ms <- round(tab.Ms, 3)

latex(tab.Ms, file="Metrics_Hedonic_KP.tex")
