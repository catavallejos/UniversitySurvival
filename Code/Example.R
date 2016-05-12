#################################################################################
# EXAMPLE FILE VALLEJOS AND STEEL (2016)
#################################################################################

## SIMULATING COVARIATE VALUES FOR 200 STUDENTS 
n.students=400
n.var=2; n.effects=4 # k*=3 VARIABLES, k=4 EFFECTS (but last effect is equal to 0)
x1<-rbinom(n=n.students,size=1,prob=0.5) # DISCRETE PREDICTOR: 2 LEVELS
x2<-rbinom(n=n.students,size=2,prob=0.5) # DISCRETE PREDICTOR: 3 LEVELS
x2.1<-as.numeric(I(x2==1)) # DUMMY INDICATOR RELATED TO x2=1
x2.2<-as.numeric(I(x2==2)) # DUMMY INDICATOR RELATED TO x2=2
x3 <- rbinom(n=n.students,size=1,prob=0.5) # EXTRA IRRELEVANT COVARIATE
X.students=cbind(x1,x2.1,x2.2, x3)

## BUILDS THE DESIGN MATRIX BASED ON COVARIATE INCLUSION INDICATORS
## THIS FUNCTION IS SPECIFIC TO THIS COVARIATE STRUCTURE
## HAS TO BE ADJUSTED FOR THE ANALYSIS OTHER DATASETS
X.design<-function(X.Period,inc)
{
  X=X.Period
  if(inc[1]==1){X=cbind(X,x1)}
  if(inc[2]==1){X=cbind(X,x2.1,x2.2)}
  if(inc[3]==1){X=cbind(X,x3)}
  return(X)
}
## INDICATES WHICH COVARIATE EFFECTS ARE ACTIVE
## THIS IS ALSO SPECIFIC TO THE SIMULATED COVARIATE STRUCTURE
ind.var<-function(inc)
{
  aux=1:16
  if(inc[1]==1){aux=c(aux,17)}
  if(inc[2]==1){aux=c(aux,18,19)}
  if(inc[3]==1){aux=c(aux,20)}
  return(aux)
}

## NUMBER OF POSSIBLE OUTCOMES (EXCLUDING CENSORING)
n.outcomes=3

## NUMBER OF PERIOD INDICATORS 
n.deltas=16

## PARAMETERS VALUES (FIRST delta REPRESENTS AN INTERCEPT) 
delta=matrix(0,ncol=3,nrow=n.deltas)
delta[,1]<-c(-8.23,-5.62,-4.82,-4.88,-5.96,-5.76,-4.89, 4.37, 4.06, 5.98, 5.95, 6.58, 6.93, 7.34, 5.47, 6.40)
delta[,2]<-c(-8.51, 6.07, 5.42, 4.89, 4.81, 4.30, 4.46, 4.96, 4.37, 3.54,-4.52,-4.92,-3.90,-4.35,-3.70,-4.48)
delta[,3]<-c(-3.14, 1.31, 0.28, 0.77,-0.80, 0.45,-2.22, 0.64, 0.24,-0.86,-0.06, 0.76, 0.34,-0.14, 1.85,-0.57)
beta=matrix(0,ncol=3,nrow=n.effects)
beta[,1]<-c(0,1,-1,0)
beta[,2]<-c(1,0,0,0)
beta[,3]<-c(0.5,-0.5,0,0)

## CAUSE-SPECIFIC LOG-HAZARDS RATIO (W.R.T. NO EVENT) AT EACH TIME POINT
hr1<-matrix(0,ncol=n.deltas,nrow=n.students)
hr2<-matrix(0,ncol=n.deltas,nrow=n.students)
hr3<-matrix(0,ncol=n.deltas,nrow=n.students)
aux1<-delta[,1]+c(0,rep(delta[1,1],times=n.deltas-1))
aux2<-delta[,2]+c(0,rep(delta[1,2],times=n.deltas-1))
aux3<-delta[,3]+c(0,rep(delta[1,3],times=n.deltas-1))
for(i in 1:n.students)
{
	hr1[i,]<-exp(aux1+rep(as.numeric(X.students[i,]%*%beta[,1]),n.deltas))
	hr2[i,]<-exp(aux2+rep(as.numeric(X.students[i,]%*%beta[,2]),n.deltas))
	hr3[i,]<-exp(aux3+rep(as.numeric(X.students[i,]%*%beta[,3]),n.deltas))
}

## CAUSE-SPECIFIC HAZARDS AT EACH TIME POINT
h0<-1/(1+hr1+hr2+hr3)
h1<-hr1/(1+hr1+hr2+hr3)
h2<-hr2/(1+hr1+hr2+hr3)
h3<-hr3/(1+hr1+hr2+hr3)

## AUXILIARY MATRIX FOR CREATING DUMMY VARIABLES FOR PERIOD INDICATORS
delta.aux<-matrix(0,ncol=n.deltas,nrow=n.deltas)
delta.aux[,1]<-rep(1,times=n.deltas)
for(i in 2:n.deltas){delta.aux[i,i]=1}

## SIMULATED OUTCOMES FOR EACH PERIOD/STUDENT
## IF STUDENT DID NOT EXPERIENCED THE EVENT BEFORE 16 SEMESTERS, IT IS MARKED AS CENSORED
## NO OTHER CENSORING MECHANISM WAS SIMULATED
data.long<-NULL; period=0
for(i in 1:n.students)
{
	outcome=0; period=0
	while(outcome==0 & period<n.deltas)
	{
		period=period+1
		outcome<-sample(x=c(0,1,2,3),size=1,replace=T,prob=c(h0[i,period],h1[i,period],h2[i,period],h3[i,period]))	
		data.long<-rbind(data.long,c(i,outcome,delta.aux[period,],X.students[i,]))
	}
}

## SURVIVAL/CENSORING TIMES 
table(data.long[,1])

## OUTCOMES FREQUENCIES
table(data.long[,2])

## SOURCING CODES
source("Internal_Codes.R")
source("User_Codes_new.R")

## THE DATA
Y=data.long[,2]
X=as.matrix(data.long[,-c(1,2)])

## PARAMETERS REQUIRED FOR MCMC.MLOG
k=dim(X)[2]
beta0=array(0,dim=c(1,k,n.outcomes))
mean.beta=array(0,dim=c(1,k,n.outcomes))

## MCMC CHAIN 
## A LONGER RUN IS REQUIRED FOR MORE ACCURATE RESULTS 
chain=MCMC.MLOG(N=10000,thin=5,Y=Y,X=X,t0=n.deltas,beta0=beta0,mean.beta=mean.beta,prec.delta=1/100,df.delta=1,
                logg0=c(5,5,5),ls.g0=c(3,3,3),prior="Benchmark-Beta",ar=0.44,fix.g=FALSE, ncov = 3)


# Removing burn-in period (change according to the number of iterations and thining period used above)
Burn = 1:1001
chain.b = list()
chain.b$beta = chain$beta[-Burn,,,drop=FALSE]
chain.b$gamma = chain$gamma[-Burn,]

# Visual convergence diagnostics for MPPIs ('true' covariate inclusion shown by horizontal line)
plot(cumsum(chain.b$gamma[,1])/(1:length(chain.b$gamma[,1])), type = "l", 
     xlab = "Iteration", ylab = "MPPI", ylim = c(0,1))
abline(h = 1, col = "red", lwd = 2)
plot(cumsum(chain.b$gamma[,2])/(1:length(chain.b$gamma[,2])), type = "l", 
     xlab = "Iteration", ylab = "MPPI", ylim = c(0,1))
abline(h = 1, col = "red", lwd = 2)
plot(cumsum(chain.b$gamma[,3])/(1:length(chain.b$gamma[,3])), type = "l", 
     xlab = "Iteration", ylab = "MPPI", ylim = c(0,1))
abline(h = 0, col = "red", lwd = 2)

# Traceplots for regression coefficients (e.g. 1st covariate, 1st outcome)
plot(chain$beta[,1,1], type = "l")
# Traceplots for regression coefficients (e.g. 1st covariate, 2nd outcome)
plot(chain$beta[,1,2], type = "l")
# Traceplots for regression coefficients (e.g. 1st covariate, 3rd outcome)
plot(chain$beta[,1,3], type = "l")

# Posterior densities, vertical line located at real simulated effect value
plot(density(chain.b$beta[,17,1])); abline(v=0,col="red",lwd=3)
plot(density(chain.b$beta[,18,1])); abline(v=1,col="red",lwd=3)
plot(density(chain.b$beta[,19,1])); abline(v=-1,col="red",lwd=3)
plot(density(chain.b$beta[,20,1])); abline(v=0,col="red",lwd=3)

plot(density(chain.b$beta[,17,2])); abline(v=1,col="red",lwd=3)
plot(density(chain.b$beta[,18,2])); abline(v=0,col="red",lwd=3)
plot(density(chain.b$beta[,19,2])); abline(v=0,col="red",lwd=3)
plot(density(chain.b$beta[,20,2])); abline(v=0,col="red",lwd=3)

plot(density(chain.b$beta[,17,3])); abline(v=0.5,col="red",lwd=3)
plot(density(chain.b$beta[,18,3])); abline(v=-0.5,col="red",lwd=3)
plot(density(chain.b$beta[,19,3])); abline(v=0,col="red",lwd=3)
plot(density(chain.b$beta[,20,3])); abline(v=0,col="red",lwd=3)


## Marginal probabilities of inclusion
round(colMeans(chain.b$gamma),2)

## Model's posterior probabilties (ordered from highest to lowest)
ModelProbs <- function(Chain, ncov)
{
  if(ncov < 1) {"No covariates are provided"}
  Models = NULL
  if(ncov > 1)
  {
    for(index in 1:(ncov-1))
    {
      Models <- paste0(Models, Chain$gamma[,index],"-")
    }
    Models <- paste0(Models, Chain$gamma[,ncov])
  }
  else{ Models = Chain$gamma[,ncov]}
  
  Probs = table(Models)/length(Models)
  
  Table = cbind.data.frame(Probs)
  
  return(Table[order(Table[,2], decreasing = TRUE),])
}

ModelProbs(chain.b, ncov = 3)
