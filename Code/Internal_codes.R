################################################################################
################################################################################
# IMPLEMENTATION OF MULTINOMIAL LOGISTIC MODEL: INTERNAL CODES 
# BASED ON POLSON ET AL (2013) AND HOLMES AND HELD (2006)
################################################################################
################################################################################

library(BayesLogit)
library(MASS)
library(mvtnorm)
library(Matrix)
library(compiler); enableJIT(3)
library(fastGHQuad)

## HYPER PRIORS FOR g[r]
## TWO DIFFERENT OPTIONS ARE ADMITTED 
prior.g<-function(g,prior,nt,k,Log=TRUE)
{
	if(prior=="Benchmark-Beta")
	{
		aux1=g/(1+g); aux2=-2*log(1+g); aux3=1*(1+g)^(-2)
		hyp1=0.01*max(nt,k^2); hyp2=0.01
	}
	if(prior %in% c("Hyper-g/n", "Bove-Held"))
	{
		aux1=g/(nt+g); aux2=log(nt)-2*log(nt+g); aux3=nt*(nt+g)^(-2)
		hyp1=1; hyp2=1/2
	}
	if(Log==TRUE)
	{
		aux=ifelse(aux1==1,-Inf,dbeta(aux1,shape1=hyp1,shape2=hyp2,log=TRUE)+aux2)
	}
	if(Log==FALSE)
	{
		aux=ifelse(aux1==1,0,dbeta(aux1,shape1=hyp1,shape2=hyp2,log=FALSE)*aux3)
	}
	return(aux)
}

## METROPOLIS HASTINGS SAMPLER FOR log(g[r])
GRWMH.MLOG.logg<-function(N=1,omega2,logg0,beta,mean.beta,X.beta,prior,nt)
{
	logg<-rep(0,times=N+1); logg[1]<-logg0; ind<-rep(0,N+1); k<-length(mean.beta); 
  
  if(prior %in% c("Benchmark-Beta","Hyper-g/n")) {prec.beta=t(X.beta)%*%X.beta}
	if(prior == "Bove-Held") {prec.beta = 4 *t(X.beta)%*%X.beta}

	for(l in 1:N)
	{
    y<-rnorm(n=1,mean=logg[l],sd=sqrt(omega2))
		u.aux<-runif(1,min=0,max=1)

		aux1=(k/2)*(logg[l]-y)-0.5*(exp(-y)-exp(-logg[l]))*t(beta-mean.beta)%*%prec.beta%*%(beta-mean.beta)
		aux2=prior.g(exp(y),prior,nt,k,Log=TRUE)-prior.g(exp(logg[l]),prior,nt,k,Log=TRUE)+(y-logg[l])
		log.aux<-aux1+aux2
		if(is.na(log.aux)==TRUE){print(logg[l]); print(y)}		
		if(is.na(log.aux)==FALSE & log(u.aux)<log.aux) {logg[l+1]<-y; ind[l+1]<-1}
		else {logg[l+1]<-logg[l]; ind[l+1]<-0}
	}
	logg<-logg[-1]; ind=ind[-1]
	list("logg"=logg,"ind"=ind)
}

## METROPOLIS HASTINGS SAMPLER FOR gamma
MH.MLOG.gamma<-function(gamma0, ncov, Y, X, beta.aux, logg.aux, t0, k, prec.delta, df.delta, prior, factor.int)
{
    # Proposal step
    gamma.select = sample(1:ncov, 1)
    gamma1 = gamma0
    gamma1[gamma.select] = 1 - gamma0[gamma.select]
    u.aux = runif(1,min=0,max=1)
  
    # Acceptance step
    inc.gamma0 = 1:k %in% ind.var(gamma0)
    inc.gamma1 = 1:k %in% ind.var(gamma1)
 
    aux1_1 = log.lik.MLOG(Y,X[,inc.gamma1],beta.aux[,inc.gamma1, , drop = F]) 
    aux1_0 = log.lik.MLOG(Y,X[,inc.gamma0],beta.aux[,inc.gamma0, , drop = F])
    
    aux2_1 = log.prior.MLOG(beta.aux[,inc.gamma1, , drop = F], #logg.aux,
                            X[,inc.gamma1],t0,prec.delta,df.delta,prior,factor.int)
    aux2_0 = log.prior.MLOG(beta.aux[,inc.gamma0, , drop = F], #logg.aux,
                            X[,inc.gamma0],t0,prec.delta,df.delta,prior,factor.int)
  
    log.aux = (aux1_1 - aux1_0) + (aux2_1 - aux2_0)

    if(is.na(log.aux)==TRUE){print("Problem"); print(gamma0); print(gamma1)}		
    if(is.na(log.aux)==FALSE & log(u.aux)<log.aux) {gamma<-gamma1; ind<-1}
    else {gamma<-gamma0; ind<-0}
  
    list("gamma"=gamma,"ind"=ind)
}

## METROPOLIS HASTINGS SAMPLER FOR gamma
MH.MLOG.gamma.index<-function(index, gamma0, ncov, Y, X, beta.aux, logg.aux, t0, k, prec.delta, df.delta, prior, factor.int)
{
  # Proposal step
#  gamma.select = sample(1:ncov, 1)
  gamma1 = gamma0
  gamma1[index] = rbinom(1,1,0.5)
  u.aux = runif(1,min=0,max=1)
  
  # Acceptance step
  inc.gamma0 = 1:k %in% ind.var(gamma0)
  inc.gamma1 = 1:k %in% ind.var(gamma1)
  
  aux1_1 = log.lik.MLOG(Y,X[,inc.gamma1],beta.aux[,inc.gamma1, , drop = F]) 
  aux1_0 = log.lik.MLOG(Y,X[,inc.gamma0],beta.aux[,inc.gamma0, , drop = F])
  
  aux2_1 = log.prior.MLOG(beta.aux[,inc.gamma1, , drop = F], #logg.aux,
                          X[,inc.gamma1],t0,prec.delta,df.delta,prior,factor.int)
  aux2_0 = log.prior.MLOG(beta.aux[,inc.gamma0, , drop = F], #logg.aux,
                          X[,inc.gamma0],t0,prec.delta,df.delta,prior,factor.int)
  
  log.aux = (aux1_1 - aux1_0) + (aux2_1 - aux2_0)
  
  if(is.na(log.aux)==TRUE){print("Problem"); print(gamma0); print(gamma1)}  	
  if(is.na(log.aux)==FALSE & log(u.aux)<log.aux) {gamma<-gamma1; ind<-1}
  else {gamma<-gamma0; ind<-0}
  
  list("gamma"=gamma,"ind"=ind)
}



## COMPUTES THE LOG-LIKELIHOOD FOR A GIVEN VALUE OF beta*. 
log.lik.MLOG<-function(Y,X,beta)
{
	CATEGORIES<-dim(beta)[3]; nt=dim(X)[1]
	aux1=matrix(0,ncol=CATEGORIES,nrow=nt); aux2=rep(0,times=nt)
	for(CAT in 1:CATEGORIES)
	{
		aux1[,CAT]=as.vector(X%*%beta[,,CAT])
		aux2=aux2+I(Y==CAT)*aux1[,CAT]
	}  
	aux3=aux2-log(1+apply(exp(aux1),1,sum))
	aux=sum(aux3)
	return(aux)
}

## AUXILIARY FUNCTION USED FOR log.prior.MLOG
auxiliary1<-function(g,hyp1,hyp2,x,k.beta,factor.int)
{
	aux=(g^(hyp1-k.beta/2-1))*((1+g)^(-(hyp1+hyp2)))*exp(-x/g)
	return(aux/factor.int)
}

## AUXILIARY FUNCTION USED FOR log.prior.MLOG
auxiliary1ChangeVar<-function(logg,hyp1,hyp2,x,k.beta,factor.int)
{
  aux=((exp(logg))^(hyp1-k.beta/2-1))*((1+exp(logg))^(-(hyp1+hyp2)))*exp(-x/(exp(logg)))*exp(logg)
  return(aux/factor.int)
}


# COMPUTES THE LOGARIGHTM OF THE PRIOR FOR beta*
log.prior.MLOG<-function(beta,X,t0,prec.delta,df.delta,prior,factor.int)
{
	CATEGORIES<-dim(beta)[3]; aux=rep(0,times=CATEGORIES); k=dim(X)[2]; nt=dim(X)[1]
	for(CAT in 1:CATEGORIES)
	{
		aux1=dmvt(beta[,1:t0,CAT],sigma=(1/prec.delta)*diag(t0),df=1,delta=rep(0,times=t0),log=TRUE)
			#(t0/2)*log(prec.delta)+sum(dt(beta[,1:t0,CAT]*sqrt(prec.delta),df=df.delta,log=TRUE))
		if(t0==k){aux2=0;aux3=0}
		else
		{
			hyp1=0.01*max(nt,(k-t0)^2); hyp2=0.01; logdet=sum(log(eigen(as.matrix(t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]),symmetric=TRUE,only.values=TRUE)$values))
			aux2=lgamma(hyp1+hyp2)-lgamma(hyp1)-lgamma(hyp2)+0.5*logdet-0.5*(k-t0)*log(2*pi)
      
#			auxiliary1<-function(g,hyp1,hyp2,x,k.beta,factor.int)
#			{
#			  aux=(g^(hyp1-k.beta/2-1))*((1+g)^(-(hyp1+hyp2)))*exp(-x/g)
#			  return(aux/factor.int)
#			}
      
			aux3 = log(aghQuad(auxiliary1ChangeVar, muHat = 0, sigmaHat = 1, rule = gaussHermiteData(10),
			        hyp1=hyp1,hyp2=hyp2,
			        x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),
			        k.beta = k-t0,factor.int=factor.int))
      
#			aux3.1=factor.int*integrate(auxiliary1,lower=0,upper=1,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
#			aux3.2=factor.int*integrate(auxiliary1,lower=1,upper=800,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
#			aux3.3=factor.int*integrate(auxiliary1,lower=1000,upper=Inf,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
#      aux3.3=factor.int*integrate(auxiliary1,lower=1000,upper=2000,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value

#      aux3=log(aux3.1+aux3.2+aux3.3)
#      aux3=log(aux3.1+aux3.2)
		}
		aux[CAT]=aux1+aux2+aux3
	}
	return(sum(aux))
}

# COMPUTES THE LOGARIGHTM OF THE PRIOR FOR beta*
log.prior.MLOG1<-function(beta.aux.gamma,logg.aux,X.gamma,t0,prec.delta,df.delta,prior,factor.int)
{
  CATEGORIES<-dim(beta.aux.gamma)[3]; aux=rep(0,times=CATEGORIES); k=dim(X.gamma)[2]; nt=dim(X.gamma)[1]
  for(CAT in 1:CATEGORIES)
  {
    aux1=dmvt(beta.aux[,1:t0,CAT],sigma=(1/prec.delta)*diag(t0),df=1,delta=rep(0,times=t0),log=TRUE)
    #(t0/2)*log(prec.delta)+sum(dt(beta[,1:t0,CAT]*sqrt(prec.delta),df=df.delta,log=TRUE))
    if(t0==k){aux2=0;aux3=0}
    else
    {
      hyp1=0.01*max(nt,(k-t0)^2); hyp2=0.01; logdet=sum(log(eigen(as.matrix(t(X.gamma[,(t0+1):k,drop=FALSE])%*%X.gamma[,(t0+1):k,drop=FALSE]),symmetric=TRUE,only.values=TRUE)$values))
      aux2=lgamma(hyp1+hyp2)-lgamma(hyp1)-lgamma(hyp2)+0.5*logdet-0.5*(k-t0)*log(2*pi)
      
      auxiliary2<-function(logg,hyp1,hyp2,x,k.beta,factor.int)
      {
        aux=((exp(logg))^(hyp1-k.beta/2-1))*((1+(exp(logg)))^(-(hyp1+hyp2)))*exp(-x/(exp(logg))) * (exp(logg))
        return(aux/factor.int)
      }
      
      aux3 = log(aghQuad(auxiliary2, muHat = logg.aux[,,CAT], sigmaHat = 1, rule = gaussHermiteData(10),
                         hyp1=hyp1,hyp2=hyp2,
                         x=0.5*as.numeric(t(beta.aux.gamma[,(t0+1):k,CAT])%*%t(X.gamma[,(t0+1):k,drop=FALSE])%*%X.gamma[,(t0+1):k,drop=FALSE]%*%beta.aux.gamma[,(t0+1):k,CAT]),
                         k.beta = k-t0,factor.int=factor.int))
      
      #aux3.1=factor.int*integrate(auxiliary1,lower=0,upper=1,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
      #aux3.2=factor.int*integrate(auxiliary1,lower=1,upper=1000,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
      #aux3.3=factor.int*integrate(auxiliary1,lower=1000,upper=Inf,hyp1=hyp1,hyp2=hyp2,x=0.5*as.numeric(t(beta[,(t0+1):k,CAT])%*%t(X[,(t0+1):k,drop=FALSE])%*%X[,(t0+1):k,drop=FALSE]%*%beta[,(t0+1):k,CAT]),k.beta=k-t0,factor.int=factor.int)$value
      #			aux3=log(aux3.1+aux3.2+aux3.3)
    }
    aux[CAT]=aux1+aux2+aux3
  }
  return(sum(aux))
}


## AUXILIARY FUNCTIONS USED FOR BRIDGE SAMPLING ESTIMATOR
logq1<-function(beta,Y,X,t0,prec.delta,df.delta,prior,factor.lik,factor.int)
{
	aux<-factor.lik+log.lik.MLOG(Y,X,beta)+log.prior.MLOG(beta,X,t0,prec.delta,df.delta,prior,factor.int)
	return(aux)
}
logq2<-function(chain1.beta,chain0.beta,t0)
{
	CATEGORIES<-dim(chain1.beta)[3]; k<-dim(chain1.beta)[1]
	aux<-rep(0,times=CATEGORIES)
	for(i in 1:CATEGORIES)
	{
		aux[i]<-dmvnorm(chain1.beta[,,i],mean=apply(chain0.beta[,,i],2,median),sigma=0.2*var(chain0.beta[,,i]),log=TRUE)
	}
	return(sum(aux))
}