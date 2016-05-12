################################################################################
################################################################################
# IMPLEMENTATION OF MULTINOMIAL LOGISTIC MODEL: USER CODES
# BASED ON POLSON ET AL (2013) AND HOLMES AND HELD (2006)
################################################################################
################################################################################


## BUILDS THE DESIGN MATRIX BASED ON COVARIATE INCLUSION INDICATORS
## THIS FUNCTION WAS USED FOR THE ANALYSIS OF THE PUC DATASET
## BUT IT HAS TO BE ADJUSTED FOR THE ANALYSIS OTHER DATASETS
#X.design<-function(X.Period,inc)
#{
#  X=X.Period
#  if(inc[1]==1){X=cbind(X,x1)}
#  if(inc[2]==1){X=cbind(X,x2)}
#  if(inc[3]==1){X=cbind(X,x3)}
#  if(inc[4]==1){X=cbind(X,x4.1,x4.2)}
#  if(inc[5]==1){X=cbind(X,x5.1,x5.2,x5.3)}
#  if(inc[6]==1){X=cbind(X,x6)}
#  if(inc[7]==1){X=cbind(X,x7)}
#  if(inc[8]==1){X=cbind(X,x8)}
#  return(X)
#}

## MAKING A MAP BETWEEN COVARIATE INCLUSION INDICATORS AND COVARIATE EFFECTS (FOR CATEGORIC VARIABLES)
## THIS FUNCTION WAS USED FOR THE ANALYSIS OF THE PUC DATASET
## BUT IT HAS TO BE ADJUSTED FOR THE ANALYSIS OTHER DATASETS
#ind.var<-function(inc)
#{
#  aux=1:16
#  if(inc[1]==1){aux=c(aux,17)}
#  if(inc[2]==1){aux=c(aux,18)}
#  if(inc[3]==1){aux=c(aux,19)}
#  if(inc[4]==1){aux=c(aux,20,21)}
#  if(inc[5]==1){aux=c(aux,22,23,24)}
#  if(inc[6]==1){aux=c(aux,25)}
#  if(inc[7]==1){aux=c(aux,26)}
#  if(inc[8]==1){aux=c(aux,27)}
#  return(aux)
#}

## MCMC SAMPLER
## g[r]'s CAN BE EITHER FIXED OR ASSIGNED AN HYPER-PRIOR (TWO DIFFERENT HYPER OPTIONS ARE ADMITTED) 
MCMC.MLOG<-function(N,thin,Y,X,t0,beta0,mean.beta,prec.delta,df.delta=1,logg0,ls.g0,prior,ar=0.44,fix.g=FALSE, ncov = 8, include = NULL)
{
	# FOR CAT+1 CATEGORIES
	# t0: INDICATES THE NUMBER OF DELTA'S 
	# beta0; mean.beta: MUST BE ARRAYS WITH DIMENSION 1 \times k \times CAT
	# prec.beta: MUST BE AN ARRAY WITH DIMENSION k \times k \times CAT
	CATEGORIES<-dim(beta0)[3] 
	k<-dim(X)[2]; nt<-dim(X)[1]; N.aux=round(N/thin,0); X1=Matrix(X) 
  prec.mat.delta = prec.delta*diag(t0)
	
	# ARRAYS WHERE CHAINS ARE STORED 
	beta<-array(0,dim=c(N.aux+1,k,CATEGORIES))
#	omega<-array(1,dim=c(N.aux+1,nt,CATEGORIES))
	logg<-array(1,dim=c(N.aux+1,1,CATEGORIES))
	lambda<-array(1,dim=c(N.aux+1,1,CATEGORIES))
	accept.g=array(0,dim=c(1,1,CATEGORIES))
	pg.aux=array(0,dim=c(1,1,CATEGORIES))
	ls.g=array(0,dim=c(N.aux+1,1,CATEGORIES)) 
  gamma = matrix(1, nrow = N.aux+1, ncol = ncov) 
  accept.gamma = rep(0, ncov)

	for(CAT in 1:CATEGORIES)
	{
		# INITIALIZATION OF PARAMETERS
		beta[1,,CAT]<-beta0[,,CAT]
		logg[1,,CAT]<-logg0[CAT]
		ls.g[1,,CAT]<-ls.g0[CAT]
	}

	beta.aux<-beta0
	omega.aux<-array(1,dim=c(1,nt,CATEGORIES)) 
	lambda.aux<-lambda[1,,,drop=FALSE]; logg.aux<-logg[1,,,drop=FALSE]; ls.g.aux=ls.g[1,,,drop=FALSE]
  if(is.null(include)) {gamma.aux <- gamma[1,]}
  else{ gamma.aux <- include}

	i_batch=rep(0,times=CATEGORIES);

	for(iter in 2:(N+1))
	{
	  # The update of gamma (only if include == NULL)
	  if(is.null(include))
	  {	    
	    for(index in 1:ncov)
	    {
	      MH.gamma = MH.MLOG.gamma.index(index, gamma.aux, ncov, Y, X, beta.aux, logg.aux, t0, k, prec.delta, df.delta, prior, factor.int = 1)
	      gamma.aux[index] = MH.gamma$gamma[index]
	      accept.gamma[index] = accept.gamma[index] + MH.gamma$ind
	    }
	    if(iter %% 20 == 0)
	    {
	      gamma.select = sample(1:ncov, 1)
	      gamma.aux[gamma.select] = rbinom(1,1,0.5)
	    }
	    print(gamma.aux)      
	  }
        
    # To avoid getting stuck, every 20 iterations, make 1 forced move (with prob 1)
      
	  # Updated design matrix
	  active.gamma = ind.var(gamma.aux)
	  inc.gamma = 1:k %in% active.gamma   
    inc.gamma.beta = inc.gamma
    inc.gamma.beta[1:t0] = F
	  X.gamma = X[,inc.gamma]
	  X1.gamma = Matrix(X.gamma) 
	  k.gamma = ncol(X.gamma)

	  # Updated precision matrix for beta  	  
	  prec.beta=array(0,dim=c(k.gamma,k.gamma,CATEGORIES))
	  prec.beta[1:t0,1:t0,1:CATEGORIES] = prec.mat.delta
	  if(t0<k.gamma)
    {
      if(prior %in% c("Benchmark-Beta","Hyper-g/n"))
      {
        prec.beta[(t0+1):k.gamma,
                  (t0+1):k.gamma,1:CATEGORIES]=t(X.gamma[,(t0+1):k.gamma,drop=FALSE])%*%X.gamma[,(t0+1):k.gamma,drop=FALSE]
      }
      if(prior == "Bove-Held")
      {
        prec.beta[(t0+1):k.gamma,
                  (t0+1):k.gamma,1:CATEGORIES]= 4 * t(X.gamma[,(t0+1):k.gamma,drop=FALSE])%*%X.gamma[,(t0+1):k.gamma,drop=FALSE]
      }
	  
	  }
    
		for(CAT in 1:CATEGORIES)
		{      
			i_batch[CAT]=i_batch[CAT]+1;

			Omega=Diagonal(x=omega.aux[,,CAT])

      # This doesn't change 
			if(CATEGORIES>1)
			{
				beta.aux.const=beta.aux[, inc.gamma,-CAT,drop=FALSE]
				AUX.beta=sapply(1:(CATEGORIES-1),function(i,M,b){exp(M%*%b[,,i])},M=X.gamma,b=beta.aux.const) # ALL BUT ONE beta'S ARE FIXED
				Const=log(1+apply(AUX.beta,1,sum))			
			}
			if(CATEGORIES==1){Const=0}
			kappa=I(Y==CAT)-1/2+diag(Omega*Const)

			prec.beta.aux=prec.beta[,,CAT]
			prec.beta.aux[1:t0,1:t0]=lambda.aux[CAT]*prec.beta[1:t0,1:t0,CAT]
			if(t0<k.gamma){prec.beta.aux[(t0+1):k.gamma,
                                   (t0+1):k.gamma]=exp(-logg.aux[CAT])*prec.beta[(t0+1):k.gamma,(t0+1):k.gamma,CAT]}

			#################################
			# UPDATING BETA'S
			V<-solve(prec.beta.aux+as.matrix(t(X1.gamma)%*%Omega%*%X1.gamma))
			B<-V%*%(prec.beta.aux%*%mean.beta[,inc.gamma,CAT]+t(X.gamma)%*%kappa)    
			beta.aux[,inc.gamma,CAT]<- mvrnorm(n=1,mu=B,Sigma=V)

			#################################
			# UPDATING OMEGA'S
			M.AUX=as.vector(X.gamma%*%beta.aux[,inc.gamma,CAT])-Const 
			omega.aux[,,CAT]<-rpg(num=nt,h=rep(1,times=nt),M.AUX)

			if(fix.g==FALSE)
			{
				#################################
				# UPDATING G
				if(t0<k)
				{
					MH.g=GRWMH.MLOG.logg(N=1,omega2=exp(ls.g.aux[,,CAT]),logg0=logg.aux[,,CAT],
                               beta=beta.aux[,inc.gamma.beta,CAT],mean.beta[,inc.gamma.beta,CAT],
                               X.beta=X.gamma[,-(1:t0)],prior=prior,nt)						
					if(MH.g$ind==1) {logg.aux[,,CAT]<-MH.g$logg; accept.g[,,CAT]=accept.g[,,CAT]+1; pg.aux[,,CAT]=pg.aux[,,CAT]+1}
   				if(i_batch[CAT]==50)
					{
						pg.aux[,,CAT]=pg.aux[,,CAT]/50; Pg.aux=as.numeric(pg.aux[,,CAT])<ar
						ls.g.aux[,,CAT]=ls.g.aux[,,CAT]+((-1)^Pg.aux)*min(0.01,1/sqrt(iter))
						i_batch[CAT]=0; pg.aux[,,CAT]=0; 
					}
				}
			}

			#################################
			# UPDATING LAMBDA'S
			z.aux=beta.aux[,1:t0,CAT]-mean.beta[,1:t0,CAT]
			lambda.aux[,,CAT]<-rgamma(1,shape=(df.delta+t0)/2,rate=0.5*(t(z.aux)%*%prec.beta[1:t0,1:t0,CAT]%*%z.aux+df.delta))	
		
			if(iter%%thin==0)
			{	
				beta[iter/thin+1,,CAT]=beta.aux[,,CAT]
				logg[iter/thin+1,,CAT]=logg.aux[,,CAT]
				ls.g[iter/thin+1,,CAT]=ls.g.aux[,,CAT]
#				omega[iter/thin+1,,CAT]=omega.aux[,,CAT]
				lambda[iter/thin+1,,CAT]=lambda.aux[,,CAT]
        gamma[iter/thin+1, ] = gamma.aux
			}

		}
		if((iter-1)%%20==0){print(paste("Iteration :",iter))}
	}

	print(paste("AR g",1:CATEGORIES,":",round(accept.g[,,1:CATEGORIES]/N,2)))
  print(paste("AR gamma:", round(accept.gamma / N, 2)))

	return(list("beta"=beta,"logg"=logg,"ls.g"=ls.g,"lambda"=lambda, "gamma" = gamma))
}



