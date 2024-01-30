enet = function(y,x,c, max.steps) {
  n <- nrow(x)
  p <- ncol(x)
  q1 <- ncol(c)
  
  XtX <- t(x) %*% x	#Time saving
  xy <- t(x) %*% y
  
  betaSamples <- matrix(0, max.steps, p)
  sigma2Samples <- matrix(0, max.steps,1) 
  TauSamples <- matrix(0, max.steps, p)
  eta1sample = matrix(0,max.steps,1)
  eta2sample = matrix(0,max.steps,1)
  GammaSamples <- matrix(0,max.steps,q1)
  
  #ls = lm(y ~ 0 + x)
  #sigma2hat = var(ls$resid)
  #betahat = na.omit(ls$coef)
  lambda1 = 1
  lambda2 = 1
  
  beta <- rep(1, p)
  sigma2 <- 1
  Tau <- rep(1.5, p)
  eta1=1
  eta2 = 1
  Gamma = rep(1,q1) # coeffecient of clinical factor
  c1 = 1e-1
  d1 = 1e-1
  c2 = 1e-1
  d2 = 1e-1
  
  k <- 0
  while (k < max.steps) 
{
    k <- k + 1
    
    # sample beta
    AA = XtX + eta2 * diag(Tau / (Tau - 1))
    invA = solve(AA)
    mean_beta = invA %*% t(x)%*%(y-c%*%Gamma)
    varcov_beta = sigma2 * invA
    beta <- MASS::mvrnorm(1, mean_beta, varcov_beta)
    betaSamples[k,] <- beta
    
    # sample Gamma
    invsig1 = solve(diag(1,nrow=q1))
    B = invsig1+t(c)%*%c/sigma2
    varcov1 = solve(B)
    res1 = y-x%*%beta
    mean1 = varcov1%*%t(t(res1)%*%c/sigma2)
    Gamma = MASS::mvrnorm(1,mean1,varcov1)
    
    GammaSamples[k,] <- Gamma
    
    # sample tau
    nu.prime = sqrt(eta1 / (eta2*(beta^2)))
    lambda.prime = eta1 / sigma2
    for (j in 1:p) {
      flag = 1
      while (flag) {
        temp = VGAM::rinv.gaussian(n=1, mu=nu.prime[j], lambda=lambda.prime)
        flag = (temp <= 0)
      }
      Tau[j] = 1 + 1 / temp
    }
    
    TauSamples[k, ] <- Tau
    
    # sample sigma2
    a.temp = n/2 + p
    b.temp = 1/2*sum( (y - x%*%beta-c%*%Gamma)^2 ) + 1/2*eta2*sum((Tau / (Tau - 1)) * beta^2) + 1/2*eta1*sum(Tau)
    
    flag.temp = 1
    while(flag.temp) {
      z.temp = MCMCpack::rinvgamma(n=1, shape = a.temp, scale = b.temp)
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
    shape2 = p/2 + c2 
    rate2 = sum(Tau/(Tau-1)*beta^2)/(2*sigma2) + d2
    eta2 = stats::rgamma(1, shape=shape2, rate=rate2)
    eta2sample[k] = eta2
    
    
    
    
  }
  
  
  dat = list(GS.beta=betaSamples,sigma2=sigma2Samples,tau=TauSamples,GS.b=GammaSamples,eta1=eta1sample,eta2=eta2sample)
  return(dat)
  
  
}
