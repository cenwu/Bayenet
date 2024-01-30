#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BL (arma::vec y, arma:: mat x, arma:: mat c, int maxSteps, arma:: vec hatEta, arma::vec hatb, arma:: vec hatInvTauSq2, arma:: mat invSigb0, double hatLambdaSqStar2, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, int progress)
{
    unsigned int n = x.n_rows, p = x.n_cols, q2 = c.n_cols;
    arma::mat
            gsb(maxSteps, q2),
            gseta(maxSteps,p),
            gsInvTauSq2(maxSteps, p);
        
    arma::vec
            gsLambdaStar2(maxSteps),
            gsSigmaSq(maxSteps);
            

    arma::mat tCC = c.t()*c, varb;
    arma::vec res, meanb, tRsRs2, muInvTauSq2;
    double tempS2, meanRs2, varRs2, lInvTauSq2;
    
    arma:: mat txx = x.t()*x;
    arma::vec tBrBr2Diag = txx.diag();
    
    
    for (int k = 0; k < maxSteps; k++) {
        
        res = y - c*hatb - x*hatEta;
        
        // b|
        varb = arma::inv(tCC/hatSigmaSq + invSigb0);
        res += c*hatb;
        meanb = varb * (c.t() * res/hatSigmaSq);
        hatb = mvrnormCpp(meanb, varb);
        res -= c* hatb;
        gsb.row(k) = hatb.t();
        
        
        // eta|
        
        for(unsigned int j=0;j<p;j++){
            tempS2 = 1/(tBrBr2Diag(j) + hatInvTauSq2(j));
            varRs2 = hatSigmaSq * tempS2;
            res += x.col(j) * hatEta(j);
            meanRs2 = arma::as_scalar(tempS2 * x.col(j).t() * res);
            hatEta(j) = R::rnorm(meanRs2, std::sqrt(varRs2));
            res -= x.col(j) * hatEta(j);
        }

        gseta.row(k) = hatEta.t();
        
        // sigma.sq|
        double shapeSig = alpha + (n+p)/2;
        double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                      arma::accu(square(hatEta) % hatInvTauSq2));
        hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
        gsSigmaSq(k) = hatSigmaSq;
        
        
        
        // invTAUsq.star2|
    
        lInvTauSq2 = hatLambdaSqStar2;
        tRsRs2 = arma::square(hatEta);
        muInvTauSq2 = arma::sqrt(hatLambdaSqStar2 * hatSigmaSq / tRsRs2);
        for(unsigned int j = 0; j<p; j++){
            hatInvTauSq2(j) = rinvgaussian(muInvTauSq2(j), lInvTauSq2);
        }
            
        gsInvTauSq2.row(k) = hatInvTauSq2.t();
        
        
        // lambda.star2|
        double shapeS2 = aStar + p;
        double rateS2 = bStar + arma::accu(1/hatInvTauSq2)/2;
        hatLambdaSqStar2 = R::rgamma(shapeS2, 1/rateS2);
        gsLambdaStar2(k) = hatLambdaSqStar2;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("GS.b") = gsb,
                            
                            //Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.beta") = gseta,
                            Rcpp::Named("GS.invTAUsq2") = gsInvTauSq2,
                            Rcpp::Named("GS.lambda.sq2") = gsLambdaStar2,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
                          
}
