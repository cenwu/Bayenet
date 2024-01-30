
DataMatrix <- function(X, Y, clin,intercept=TRUE, debugging=FALSE)
{
  g = as.matrix(X);  y = Y
  n = nrow(g); p = ncol(g)
  c=clin; 
  
  clin.names = G.names = NULL
  g = scale(g, center = TRUE, scale=FALSE)
  
  if(!is.null(y)){
    if(length(y) != n)  stop("Length of Y does not match the number of rows of X.");
  }
  
  if(!is.null(clin)){
    clin = as.matrix(clin)
    if(nrow(clin) != n)  stop("clin has a different number of rows than X.");
    if(is.null(colnames(clin))){colnames(clin)=paste("clin.", 1:ncol(clin), sep="")}
    CLC = clin
    noClin = FALSE
    clin.names = colnames(clin)
  }
  
  if(intercept){ # add intercept
    CLC = cbind(matrix(1,n,1,dimnames=list(NULL, "IC")), CLC)
  }
  
  
  if(is.null(colnames(g))){
    G.names = paste("G", 1:p, sep="")
  }else{
    G.names = colnames(g)
  }
 
   
  
  dat = list(y=y, c=CLC, g=g, n=n, p=p, G.names=G.names, clin.names=clin.names)
  return(dat) 
}
