enetss = function(y,x,c, max.steps) {
  n <- nrow(x)
  p <- ncol(x)
  q1 <- ncol(c)
  
  XtX <- t(x) %*% x	#Time saving
  xy <- t(x) %*% y
  
  betaSamples <- matrix(0, max.steps, p)
  sigma2Samples <- matrix(0, max.steps,1) 
  TauSamples <- matrix(0, max.steps, p)
  piSamples <- matrix(0, max.steps,1)
  SS <- matrix(0, max.steps, p)
  eta1sample = matrix(0,max.steps,1)
  eta2sample = matrix(0,max.steps,1)
  GammaSamples <- matrix(0,max.steps,q1)
  
  #ls = lm(y ~ 0 + x)
  #sigma2hat = var(ls$resid)
  betahat = rep(1, p)
  lambda1 = 1
  lambda2 = 1
  
  beta <- rep(1, p)
  sigma2 <- 1
  Tau <- rep(1.5, p)
  r=1
  u=1
  pi <- 1/2
  eta1=lambda1^2/(4*lambda2)
  eta2 = lambda2
  Gamma = rep(1,q1)
  c1 = 1e-1
  d1 = 1e-1
  c2 = 1e-1
  d2 = 1e-1
  
  
  k <- 0
  while (k < max.steps) 
  {
    k <- k + 1
    
    z <- rep(0,p)
    sg <- rep(0,p)
    
    # sample beta
    for (j in 1:p) 
    {
      inv_Tau = (Tau[j]*eta2)/(Tau[j]-1)
      A <- t(x[,j])%*%x[,j] + inv_Tau
      inv_A <- 1/A
      res <- (y-x[,-j]%*%beta[-j]-c%*%Gamma)
      mean <- inv_A*t(x[,j])%*%res
      var <- sigma2 * inv_A
      L <- inv_A*t(res)%*%x[,j]%*%t(x[,j])%*%res
      #L <- sqrtm(inv_A)%*%t(x[,sub])%*%(y-x%*%Beta+x[,sub]%*%Beta[sub])
      l0 <- pi+(1-pi)*(inv_Tau^(1/2))*sqrt(abs(inv_A))*exp((1/(2*sigma2))*L)
      l <- pi/l0
      u<-stats::runif(1)
      
      if(u<=l) 
      { beta[j] <- 0; sg[j]=0; z[j]=0
      }else {
        beta[j] <- stats::rnorm(1, mean, sqrt(var)); sg[j]=1; z[j]=1}
    
    }
    
    betaSamples[k,] <- beta
    SS[k,] <- sg
    
    # sample Gamma
    invsig1 = solve(diag(1,nrow=q1))
    B = invsig1+t(c)%*%c/sigma2
    varcov1 = solve(B)
    res1 = y-x%*%beta
    mean1 = varcov1%*%t(t(res1)%*%c/sigma2)
    Gamma = MASS::mvrnorm(1,mean1,varcov1)
    
    GammaSamples[k,] <- Gamma
  
    # sample tau
    
    for (j in 1:p) 
    {
      nu.prime = sqrt(eta1 / (eta2*(beta[j]^2)))
      lambda.prime = eta1 / sigma2
      
      if(z[j]==0) 
      {
        flag = 1
        while (flag) {
        temp <-hbmem::rtgamma(1,shape= 1/2,scale = (2*sigma2)/eta1,a=1, b=Inf)
         # rtrunc(1,"gamma",shape= 1/2,scale = (8*lambda2*sigma2)/lambda1^2,a=1, b=Inf)
        #flag = (temp = Inf)
         flag = (temp <=1)|(temp == Inf)
        }
        Tau[j]=temp
      }else {
        flag = 1
        while (flag) {
          temp = VGAM::rinv.gaussian(n=1, mu=nu.prime, lambda=lambda.prime)
          flag = (temp <= 0)
        }
        Tau[j] = 1 + 1 / temp
      }
        
    }
    TauSamples[k, ] <- Tau
    
    #sample pi
    shape1 <- r + p - sum(z)
    shape2 <- u + sum(z)
    pi <- stats::rbeta(1,shape1, shape2)
    piSamples[k] <- pi
    
    
    # sample sigma2
    a.temp = n/2 +p/2+ sum(z)/2
    b.temp = 1/2*sum( (y - x%*%beta-c%*%Gamma)^2 ) + 1/2*eta2*sum((Tau / (Tau - 1)) * (beta^2)*z) + 1/2*eta1*sum(Tau)
    
    flag.temp = 1
    while(flag.temp) {
      z.temp = 1/stats::rgamma(1, shape = a.temp, rate = b.temp)
      u.temp = stats::runif(1)
      if (log(u.temp) <= p*log(base::gamma(0.5))-p*log(gsl::gamma_inc(1/2, eta1/(2*z.temp)))) {
        sigma2 = z.temp
        flag.temp = 0
      }
    }
    sigma2Samples[k] <- sigma2
    
    #sample eta1
    rejections = 0
    temp.shape = p+c1
    temp.rate = sum(Tau-1) + (2*sigma2)*d1
    temp.eta1 = stats::rgamma(1, shape=temp.shape, rate=temp.rate)
    temp =	p*log(stats::pgamma(eta1,shape=1/2,lower=F) / stats::pgamma(temp.eta1,shape=1/2,lower=F)) + 
      p/2 * log(eta1/temp.eta1) + p*(eta1 - temp.eta1)
    temp = min(temp, 0)
    u = log(stats::runif(1))
    if(u <= temp) {eta1 = temp.eta1} else  {rejections = rejections + 1}
    
    eta1sample[k] = eta1*2*sigma2
    
    #sample eta2
    shape2 = c2 + 1/2*sum(z)
    rate2 = sum(Tau/(Tau-1)*beta^2*z)/(2*sigma2) + d2
    eta2 = stats::rgamma(1, shape=shape2, rate=rate2)
    eta2sample[k] = eta2
    
    
  }
    
  t1 = betaSamples
  t2 = sigma2Samples
  t3 = TauSamples
  t4 = piSamples
  t5 = SS
  
  dat = list(GS.beta = betaSamples,sigma2=sigma2Samples,tau=TauSamples,GS.b=GammaSamples,pi=piSamples, GS.SS=SS,eta1=eta1sample, eta2=eta2sample)
  return(dat) 
    
  
}
