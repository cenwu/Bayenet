Qenetss = function(y,x,c, theta, max.steps)
{
  n = nrow(x)
  p = ncol(x)
  q1 <- ncol(c)
  
  betasample = matrix(0,max.steps,p)
  tsample = matrix(0,max.steps,p)
  tausample = matrix(0,max.steps,1)
  eta1sample = matrix(0,max.steps,1)
  eta2sample = matrix(0,max.steps,1)
  vsample = matrix(0,max.steps,n)
  pisample = matrix(0,max.steps,1)
  SS = matrix(0,max.steps,p)
  Gammasamples <- matrix(0,max.steps,q1)
  
  #theta = 0.5
  xi1 = (1 - 2*theta) / (theta*(1-theta))
  xi2 = sqrt(2 / (theta*(1-theta)))
  beta = rep(1,p)
  t = rep(2,p)
  v = rep(1,n)
  tau = 1
  eta1 = 1
  eta2 = 1
  a = 1e-1
  b = 1e-1
  c1 = 1e-1
  d1 = 1e-1
  c2 = 1e-1
  d2 = 1e-1
  r1=1
  u1=1
  pi <- 1/2
  Gamma = rep(1,q1)
  gamma0 = 1
  
  for(k in 1:max.steps){
    
    #sample Gamma
    for(j in 1:q1)
    {
      A = c[,j]^2/v
      invsigma2 = tau*sum(A)/xi2+1/gamma0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -c[,-j]%*%Gamma[-j]-x%*%beta-xi1*v)
      B = y_j*c[,j]/v
      mu = tau*sum(B)*sigma2/xi2
      Gamma[j] = stats::rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    Gammasamples[k,] <- Gamma
    
    #sample v
    res = y-x%*%beta-c%*%Gamma
    lambda = xi1^2*tau/(xi2^2) + 2*tau
    mu = sqrt(lambda * xi2^2/(tau*res^2))
    index = 1:n
    flag=1
    while(flag){
      inv_v= SuppDists::rinvGauss(length(index),nu=mu[index],lambda=lambda)
      flag = any(inv_v<=0)|any(is.na(inv_v))
      v[index[inv_v>0]] = 1/inv_v[inv_v>0]
      index = base::setdiff(index,index[inv_v>0])
    }
    vsample[k,] = v
    
    #sample beta
    z <- rep(0,p)
    sg <- rep(0,p)
    for(j in 1:p){
      A = x[,j]^2/v
      invsigma2 = tau*sum(A)/xi2^2 + 2*eta2*t[j]/(t[j]-1)
      sigma2 = 1/invsigma2
      y_j = as.vector(y -xi1*v -x[,-j]%*%beta[-j]-c%*%Gamma)
      B = y_j*x[,j]/v
      mu = tau*sum(B)*sigma2/xi2^2
      BB = tau*sum(B)/xi2^2
      d = (2*eta2*t[j]/(t[j]-1))^(1/2)*sqrt(sigma2)*exp((1/2)*(sqrt(sigma2)*BB)^2)
      l = pi/(pi+(1-pi)*d)
      u = stats::runif(1)
      if(u<l){
        beta[j] = 0; z[j]=0; sg[j]=0
      }
      else{
      beta[j] = stats::rnorm(1,mean=mu,sd=sqrt(sigma2)); z[j]=1; sg[j]=1
      }
    }
    betasample[k,]=beta
    SS[k,] <- sg
    
    #sample tau
    res = y-x%*%beta-xi1*v-c%*%Gamma
    vec = res^2/(2*xi2^2*v)+v
    rate = sum(vec)+b
    shape = a+3*n/2
    tau = stats::rgamma(1,shape=shape,rate=rate)
    tausample[k,] = tau
    
    #sample t
    for(j in 1:p){
      if(beta[j]==0){
        flag = 1
        while (flag) {
          temp.t = hbmem::rtgamma(1,shape= 1/2,scale = 1/eta1,a=1, b=Inf)
          flag = (temp.t <=1)|(temp.t == Inf)
        }
        t[j]=temp.t
      }
      else{
    temp.lambda = 2 * eta1
    temp.nu =  sqrt(temp.lambda / (2*eta2*beta[j]^2))
    flag = 1
    while (flag) {
      temp.s = SuppDists::rinvGauss(1, lambda = temp.lambda, nu = temp.nu)
      flag = (temp.s<=0)|(is.na(temp.s))|(temp.s == Inf)
    }
    t[j] = 1 / temp.s + 1
      }
    }
    tsample[k,] = t
    
    #sample eta1
    rejections = 0
    temp.shape = p + c1 
    temp.rate = sum(t-1) + d1
    temp.eta1 = stats::rgamma(1, shape=temp.shape, rate=temp.rate)
    temp =	p*log(stats::pgamma(eta1,shape=1/2,lower=F) / stats::pgamma(temp.eta1,shape=1/2,lower=F)) + 
      p/2 * log(eta1/temp.eta1) + p*(eta1 - temp.eta1)
    temp = base::min(temp, 0)
    u = log(stats::runif(1))
    if(u <= temp) {eta1 = temp.eta1} else  {rejections = rejections + 1}
    
    eta1sample[k] = eta1
    
    #sample eta2
    shape2 = c2 + 1/2*sum(z) 
    rate2 = sum(t/(t-1)*beta^2*z) + d2
    eta2 = stats::rgamma(1, shape=shape2, rate=rate2)
    eta2sample[k] = eta2
    
    #sample pi
    shape1 <- r1 + p - sum(z)
    shape2 <- u1 + sum(z)
    pi <- stats::rbeta(1,shape1, shape2)
    pisample[k] <- pi
    
  }
  
  
  dat = list(GS.beta=betasample,t=tsample,tau=tausample,GS.b=Gammasamples,eta1=eta1sample,eta2=eta2sample,
             v=vsample,pi=pisample, GS.SS=SS)
  return(dat)
  
  
  
}