
nonrobust <- function(X, Y, clin, max.steps, sparse, penalty,debugging=FALSE)
{
  dat = DataMatrix(X, Y, clin, intercept=TRUE, debugging=FALSE)
  c=dat$c; g=dat$g; y=dat$y; beta_true=dat$coef
  n = dat$n; p= dat$p; q=ncol(c)
  
  G.names = dat$G.names
  clin.names = dat$clin.names
  
  hatb = rep(1,q); hatEta= rep(1,p); 
  invSigb0 = diag(rep(1,q))
  hatInvTauSq2 = rep(1,p)
  sg2 = rep(1,p); hatLambdaSqStar2=1; hatSigmaSq=1
  aStar=1; bStar=1; alpha=1; gamma=1
  progress = 0; hatPiEta=1/2; 
  mu0=1; nu0=1
  
  if(sparse){
    fit=switch (penalty,
                "lasso" = BLSS(y,g,c,max.steps,hatEta,hatb,hatInvTauSq2,sg2,hatPiEta,invSigb0,hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma,mu0,nu0, progress),
                "elastic net" = enetss(y,g,c, max.steps)
    )
  }else{
    fit=switch (penalty,
                "lasso" = BL(y,g,c,max.steps,hatEta,hatb,hatInvTauSq2,invSigb0, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma, progress),
                "elastic net" = enet(y,g,c, max.steps)
    )
  }
  
  out = list( GS.alpha = fit$GS.b,
              GS.beta = fit$GS.beta)
  
  if(sparse){
    class(out)=c("Sparse", "BVS")
  }else{
    class(out)=c("NonSparse", "BVS")
  }
  
  out
  
  
}  