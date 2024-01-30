
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bayenet

> Bayesian Quantile Elastic Net for Genetic Study 

As heavy-tailed error distribution and outliers in the response variable widely exist, models which are robust to data contamination are highly demanded. Here, We develop a novel robust Bayesian variable selection method with elastic net penalty for quantile regression in genetic analysis. In particular, the spike-and-slab priors have been incorporated to impose sparsity. An efficient Gibbs sampler has been developed to facilitate computation. The algorithms of the proposed and alternative methods are efficiently implemented in 'C++'.
## How to install

 - To install from github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("cenwu/Bayenet")

  - Released versions of Bayenet are available on CRAN
    [(link)](https://cran.r-project.org/package=Bayenet), and can be
    installed within R via

<!-- end list -->

    install.packages("Bayenet")




## Examples

#### Example.1 (default method: Bayesian quantile elastic net with spike-and-slab priors)

    library(Bayenet)
    data(dat)
    
    max.steps=10000
    fit= Bayenet(X, Y, clin, max.steps, penalty="elastic net")
    ## coefficients of parameters
    fit$coefficient
    ## Estimated values of main G effects 
    fit$coefficient$G
    ## Estimated values of clincal effects 
     fit$coefficient$clin
        
     ## Prediction   
     test=sample((1:nrow(X)), floor(nrow(X)/5))
     fit=Bayenet(X[-test,], Y[-test], clin[-test,], max.steps=10000, penalty="elastic net")  
     predict.Bayenet(fit, X[test,], clin[test,], Y[test,])

#### Example.2 (alternative: Bayesian quantile lasso)

    fit= Bayenet(X, Y, clin, max.steps, sparse="False", penalty="lasso")
    

## Methods

This package provides implementation for methods proposed in

  - Lu, X., Ren, J., Ma, S and Wu, C. (2023+). Bayesian quantile elastic net with spike-and-slab priors. 
