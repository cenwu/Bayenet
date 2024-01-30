
robust <- function(X, Y, clin, max.steps, sparse, penalty,debugging=FALSE)
{
  dat = DataMatrix(X, Y, clin, intercept=TRUE, debugging=FALSE)
  c=dat$c; g=dat$g; y=dat$y; beta_true=dat$coef
  n = dat$n; p= dat$p; q=ncol(c)
  
  G.names = dat$G.names
  clin.names = dat$clin.names
  
  theta=0.5
  hatb = rep(1,q); hatEta= rep(1,p); hatTau=1; hatV = rep(1,n)
  invSigb0 = diag(rep(1,q))
  hatSg2 = rep(1,p); hatEtaSq2=1;
  r=1; a=1; b=1; ss2 = rep(1,p)
  hatPiBeta=1/2; hatPiEta=1/2
  sh0=1; sh1=1
  progress = 0
  
  if(sparse){
    fit=switch (penalty,
                "lasso" = QBLSS(y,g,c,max.steps,hatb,hatEta,hatTau,hatV,hatSg2,ss2,invSigb0, hatPiEta,hatEtaSq2, theta, r,a ,b ,sh0,sh1, progress),
                "elastic net" = Qenetss(y,g,c, theta, max.steps)
    )
  }else{
    fit=switch (penalty,
                "lasso" = QBL(y,g,c,max.steps,hatb,hatEta,hatTau,hatV,hatSg2,invSigb0, hatEtaSq2, theta, r,a ,b, progress),
                "elastic net" = Qenet(y,g,c, theta, max.steps)
    )
  }
  
  out = list( GS.alpha = fit$GS.b,
              GS.beta = fit$GS.beta)
  
  if(sparse){
    class(out)=c("Sparse", "RBVS")
  }else{
    class(out)=c("NonSparse", "RBVS")
  }
  
  out
  
  
}  