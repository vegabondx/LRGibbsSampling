library(faraway)
library(R2OpenBUGS)
library(modeest)

# Traditional Logistic Regression 
rm(list=ls())
data(orings)
x<-orings$temp
y<-orings$damage
orings[1,2]=1             # Convert to Bernouli  
model.lr <- glm(damage~temp,data=orings,family=binomial("logit"))
summary(model.lr)

# Gibbs Sampling Model for Bernouli Random Variable  
N=100000
bin=1000
n=23                                  
#y<-c(1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0) #TRUE
#T<-c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
#y<-c(1,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,1,0,0,1) # TEST 
y<-c(orings$damage)
T<-c(orings$temp)
#T <- x
data=list("y","T","n")
parameters<-c("beta")
inits=function(){ list(beta=c(0,0)) }
bugfile<-gsub("[[:space:]]","",paste(getwd(),"/src/bugs/modelbernouli.txt"))
# library to interface winbugs and R 
modelbernouli.sim=bugs(data, inits, parameters,model.file=bugfile,n.chains=1,n.iter=N,n.thin=1,n.burnin=bin)
output<-modelbernouli.sim$sims.matrix
beta1<-output[,1]
beta2<-output[,2]
deviance<-output[,3]
# Plotting the densities
plot(density(beta1))
plot(density(beta2))
modeest::mlv(beta1,method = "naive") 
modeest::mlv(beta2,method ="naive")

"Results"
results<-rbind(c(mean(beta1),mean(beta2)),c(median(beta1),median(beta2)),c(model.lr$coefficients))
resultdf<-as.data.frame(results,row.names=c("Mean GS Estimate","Median GS Estimate","Standrd LR Estimate"))
resultdf  