##################################
# File: Albert and Chib (Probit).r
# Generate Binary Data and run Gibbs Sampler
# Assumes conjugate normal prior on Beta
# Albert and Chib (1993) section 3
# Nov 15, 2007

# Reproducing code with some modification and cleaning from 
# http://people.musc.edu/~brn200/r/
# Abhijit Oka 
# August 15,2016 
##################################
library(mvtnorm)
library(msm)          			# For rtnorm function

# Simulate data
set.seed(250)
n<-1000					                                 # Sample Size
x<-sample(0:4,n,replace=T)		                   # Data Definition
eps <- sample(0:1,n,replace=T,prob=c(0.15,0.85)) # Disturbance
y<-1*(x>3)   	                               # Observations - simulated

X<-matrix(c(rep(1,n),x), ncol=2)   	# Design matrix first column for constant
k<-2 					                	    # Number of coefficients ( constant + linear factor)

# Initial Values
variable.W<-rep(0,k)                     # Linear VAR
variable.Z<-rep(0,n)				             # Latent Normal Variable VAR

# Priors on Weights
mean_w0<-c(0,0)				                # Prior Mean for Beta
cov_w0<-diag(10,2)		                # Prior Cov of Beta (vague)

# Create vectors to store results
nsim <-1500  			                  	# Number of Iterations of Gibbs Sampler

variable.W.collect <-matrix(0,nrow=nsim,ncol=k) 	      # Store Results

###################
# GIBBS SAMPLER	#
###################
# Posterior Variance of W
prec0 <- solve(cov_w0)                  # Dont know why this is needed
variable.W.sd <-solve(prec0+t(X) %*% X) # MATRIX inversion needed only once

for (i in 2:nsim) {
  
  # Compute mean of Z 
  variable.Z.mean <-X%*%variable.W			# Update Mean of Z
  # Draw Latent Variable, z, from its full conditional, given y
  variable.Z[y==0]<-qnorm(runif(n,0,pnorm(0,variable.Z.mean,1)),variable.Z.mean,1)[y==0]  # Bounded to the right 
  variable.Z[y==1]<-qnorm(runif(n,pnorm(0,variable.Z.mean,1),1),variable.Z.mean,1)[y==1]  # Bounded to the left 
  
  # Compute mean of W given Z ( W is independent of Y )
  variable.W.mean <- variable.W.sd%*%(crossprod(X,variable.Z))  # Compute mean  
  # Draw W 
  variable.W<-c(rmvnorm(1,variable.W.mean,variable.W.sd))    # SAMPLE Weights 
  
  variable.W.collect[i,]<-variable.W
  
  if (i%%100==0) print(i)
}

# Get Summaries
burnout<-500
apply(variable.W.collect[burnout:nsim,],2,mean)

fit<-glm(y~x,family=binomial(link=probit))

plot(1:nsim,variable.W.collect[,2], type="l", col="lightgreen")
abline(h=mean(variable.W.collect[250:nsim,2] ),col="blue4")