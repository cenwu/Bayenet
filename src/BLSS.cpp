#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BLSS (arma::vec y, arma:: mat x, arma:: mat c, int maxSteps, arma:: vec hatEta, arma::vec hatb, arma:: vec hatInvTauSq2,arma::vec sg2,double hatPiEta, arma:: mat invSigb0, double hatLambdaSqStar2, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, double mu0, double nu0,int progress)
{
    unsigned int n = x.n_rows, p = x.n_cols, q2 = c.n_cols;
    arma::mat
            gsb(maxSteps, q2),
            gseta(maxSteps,p),
            gsInvTauSq2(maxSteps, p),
            gsSS2(maxSteps,p);
        
    arma::vec
            gsLambdaStar2(maxSteps),
            gsSigmaSq(maxSteps),
            gsPiEta(maxSteps);


    arma::mat tCC = c.t()*c, varb;
    arma::vec res, meanb, tRsRs2, muInvTauSq2;
    double tempS2, meanRs2, varRs2, lE, t, lInvTauSq2, WjtRes;
    
    arma:: mat txx = x.t()*x;
    arma::vec tBrBr2Diag = txx.diag();
    
    
    for (int k = 0; k < maxSteps; k++) {
        
        res = y - c*hatb- x*hatEta;
        
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
            WjtRes = arma::as_scalar(x.col(j).t() * res);
            meanRs2 = tempS2*WjtRes;
            
            lE = hatPiEta/(hatPiEta + (1-hatPiEta)*sqrt(hatInvTauSq2(j)*tempS2)*exp(0.5/hatSigmaSq*tempS2*pow(WjtRes,2)));
            
            t = R::runif(0, 1);
            if(t<lE){
                hatEta(j) = 0;sg2(j)=0;
            }else{
                hatEta(j) = R::rnorm(meanRs2, sqrt(varRs2));sg2(j)=1;
                                 
            }
            res -= x.col(j) * hatEta(j);
        }

        gseta.row(k) = hatEta.t();
        gsSS2.row(k) = sg2.t();
        
        // sigma.sq|
        double shapeSig = alpha + n/2 + arma::accu(hatEta!=0)/2;
        double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                      arma::accu(square(hatEta) % hatInvTauSq2));
        hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
        gsSigmaSq(k) = hatSigmaSq;
        
        
        
        // invTAUsq.star2|
    
        lInvTauSq2 = hatLambdaSqStar2;
        tRsRs2 = arma::square(hatEta);
        muInvTauSq2 = arma::sqrt(hatLambdaSqStar2 * hatSigmaSq / tRsRs2);
        for(unsigned int j = 0; j<p; j++){
            if(hatEta(j)==0){
                hatInvTauSq2(j) = 1/R::rgamma(1,2/lInvTauSq2);
            }else{
                hatInvTauSq2(j) = rinvgaussian(muInvTauSq2(j), lInvTauSq2);
            }
            
        }
            
        gsInvTauSq2.row(k) = hatInvTauSq2.t();
        
        
        // lambda.star2|
        double shapeS2 = aStar + p;
        double rateS2 = bStar + arma::accu(1/hatInvTauSq2)/2;
        hatLambdaSqStar2 = R::rgamma(shapeS2, 1/rateS2);
        gsLambdaStar2(k) = hatLambdaSqStar2;
        
    
        // pi.star|
        double shape1_e = mu0 + arma::accu(hatEta == 0);
        double shape2_e = nu0 + arma::accu(hatEta != 0);
        hatPiEta = R::rbeta(shape1_e, shape2_e);
        gsPiEta(k) = hatPiEta;
        
    }
    
    
    return Rcpp::List::create(Rcpp::Named("GS.b") = gsb,
                            Rcpp::Named("GS.beta") = gseta,
                            Rcpp::Named("GS.invTAUsq2") = gsInvTauSq2,
                            Rcpp::Named("GS.lambda.sq2") = gsLambdaStar2,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                            Rcpp::Named("GS.SS") = gsSS2,
                            Rcpp::Named("GS.Pi.beta") = gsPiEta);
                          
}
