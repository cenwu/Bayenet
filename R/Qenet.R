Qenet = function(y,x,c, theta,max.steps)
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
  Gamma = rep(1,q1)
  gamma0 = 1
  
  for(k in 1:max.steps){
    
    #sample Gamma
    for(j in 1:q1)
    {
      A = c[,j]^2/v
      invsigma2 = tau*sum(A)/(xi2^2)+1/gamma0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -c[,-j]%*%Gamma[-j]-x%*%beta-xi1*v)
      B = y_j*c[,j]/v
      mu = tau*sum(B)*sigma2/(xi2^2)
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
    for(j in 1:p){
      A = x[,j]^2/v
      invsigma2 = tau*sum(A)/xi2^2 + 2*eta2*t[j]/(t[j]-1)
      sigma2 = 1/invsigma2
      y_j = as.vector(y -xi1*v -x[,-j]%*%beta[-j]-c%*%Gamma)
      B = y_j*x[,j]/v
      mu = tau*sum(B)*sigma2/xi2^2
      beta[j] = stats::rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    betasample[k,]=beta
    
    #sample tau
    res = y-x%*%beta-xi1*v-c%*%Gamma
    vec = res^2/(2*xi2^2*v)+v
    rate = sum(vec)+b
    shape = a+3*n/2
    tau = stats::rgamma(1,shape=shape,rate=rate)
    tausample[k,] = tau
    
    #sample t
    temp.lambda = 2 * eta1
    temp.nu =  sqrt(temp.lambda / (2*eta2*beta^2))
    index = 1:p
    flag = 1
    while (flag) {
      temp.s = SuppDists::rinvGauss(length(index), lambda = temp.lambda, nu = temp.nu[index])
      flag = any(temp.s <= 0) | any(is.na(temp.s))
      t[index[temp.s>0]] = 1 / temp.s[temp.s>0] + 1
      index = base::setdiff(index, index[temp.s>0])
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
    shape2 = p/2 + c2 
    rate2 = sum(t/(t-1)*beta^2) + d2
    eta2 = stats::rgamma(1, shape=shape2, rate=rate2)
    eta2sample[k] = eta2
    
  }
  
  
  dat = list(GS.beta=betasample,t=tsample,tau=tausample,GS.b=Gammasamples,eta1=eta1sample,eta2=eta2sample,v=vsample)
  return(dat)
  

  
}