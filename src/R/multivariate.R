# Car data 
# logistic regression with r 
library(R2OpenBUGS)
library(modeest)
rm(list=ls())
cardata <- read.delim("data/cardata.csv") 
s<-glm(Acc~buying+maint+doors.+person+lug_boot+safety,data=cardata,family=binomial("logit"))
summary(s)

# Doing regression in Logistic counterpart 
n=1728
data=list("y","X1","X2","X3","X4","X5","X6","n")
parameters=c("beta")
inits=function(){ list(beta=c(0,0,0,0,0,0,0)) }
y=as.numeric(cardata$Acc)
X1=as.numeric(cardata$buying)
X2=as.numeric(cardata$maint)
X3=as.numeric(cardata$doors.)
X4=as.numeric(cardata$person)
X5=as.numeric(cardata$lug_boot)
X6=as.numeric(cardata$safety)
bugfile<-gsub("[[:space:]]","",paste(getwd(),"/src/bugs/multivariate1.txt"))
modelmultivariate.sim=bugs(data, inits, parameters,model.file=bugfile,n.chains=1,n.iter=100,n.thin=1)
A<-modelmultivariate.sim$sims.matrix
estimate_GS_median<-apply(A,2,median)[1:7]
estimate_GS_mean<-apply(A,2,mean)[1:7]
estimate_LR<-s$coefficients

X<-apply(A,2,modeest::mlv,method="naive")
result<-list()
for(I in names(X)) {result<-c(result,eval(substitute(X$i$M,list(i=I))))}
estimate_MLV_modal=as.vector(result)[1:7]
rbind(estimate_GS_mean,estimate_GS_median,estimate_MLV_modal,estimate_LR)


