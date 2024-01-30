#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List QBL (arma::vec y, arma:: mat x, arma:: mat c, int maxSteps, arma::vec hatb, arma::vec hatEta, double hatTau, arma::vec hatV, arma:: vec hatSg2, arma:: mat invSigb0, double hatEtaSq2, double theta, double r, double a, double b, int progress)
{
    unsigned int n = x.n_rows, p = x.n_cols, q2 = c.n_cols;
    arma::mat
            gsb(maxSteps,q2),
            gsEta(maxSteps,p),
            gsV(maxSteps, n),
            gsSg2(maxSteps, p);
    
    arma::vec
              gsEtaSq2(maxSteps),
              gsTau(maxSteps);
           
    
    arma::mat varb, tCCoV(q2,q2);
    arma::vec res, RCoV(q2), muV, muS2, meanb;

    double lambV, xi1=(1-2*theta)/(theta*(1-theta)),xi2=std::sqrt(2/(theta*(1-theta))),xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), XgXgoV2, RXgoV2, meanG2,varG2, ResSqoV;
    
    
    for (int k = 0; k < maxSteps; k++) {
        
        res = y - c*hatb - xi1*hatV -x*hatEta;
        
        
        
        // Rcpp::Rcout << "b" << std::endl;
        
        tCCoV = (c.each_col()/hatV).t() * c;
        res += c*hatb;
        RCoV = arma::sum(c.each_col()% (res/hatV), 0).t();
        varb = arma::inv_sympd(tCCoV*hatTau/xi2Sq + invSigb0);
        meanb = varb * RCoV * hatTau / xi2Sq;
        hatb = mvrnormCpp(meanb, varb);
        res -= c* hatb;
        gsb.row(k) = hatb.t();
 
        
        // Rcpp::Rcout << "v" << std::endl;
        res += xi1*hatV;
        lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
        muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
        for(unsigned int i = 0; i<n; i++){
            bool flag = true;
            while(flag){
                hatV(i) = 1/rinvGauss(muV(i), lambV);
                if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
                    if(progress != 0) Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl;
                    Rcpp::checkUserInterrupt();
                }else{
                    flag = false;
                }
            }
        }
        res -= xi1*hatV;
        gsV.row(k) = hatV.t();
        
        

    
        // Rcpp::Rcout << "S2" << std::endl;
       muS2 = std::sqrt(hatEtaSq2)/ arma::abs(hatEta);
       for(unsigned int j = 0; j<p; j++){
           bool flag = true;
           while(flag){
               hatSg2(j) = 1/rinvGauss(muS2(j), hatEtaSq2);
               if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
                   if(progress != 0) Rcpp::Rcout << "hatSg2(j): " << hatSg2(j) << std::endl;
                   Rcpp::checkUserInterrupt();
               }else{
                   flag = false;
               }
           }
       }
       gsSg2.row(k) = hatSg2.t();
        
        
        // Rcpp::Rcout << "Eta" << std::endl;
        for(unsigned int j=0; j<p; j++){
            res += x.col(j) * hatEta(j);
            XgXgoV2 = arma::as_scalar((x.col(j)/hatV).t() * x.col(j));
            varG2 = 1/(XgXgoV2*hatTau/xi2Sq + 1/hatSg2(j));
            
            RXgoV2 = arma::sum(x.col(j) % (res/hatV));
            meanG2 = varG2 * RXgoV2 * hatTau / xi2Sq;
            hatEta(j) = R::rnorm(meanG2, sqrt(varG2));
            res -= x.col(j) * hatEta(j);
        }
        gsEta.row(k) = hatEta.t();
        
        
        // Rcpp::Rcout << "tau" << std::endl;
        double shape = a + 3*n/2;
        ResSqoV = arma::accu(arma::square(res)/hatV);
        double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
        hatTau = R::rgamma(shape, 1/rate);
        gsTau(k) = hatTau;
        
        // Rcpp::Rcout << "eta2Sq2" << std::endl;
        double shape2 = p+1;
        double rate2 = arma::accu(hatSg2)/2 + r;
        hatEtaSq2 = R::rgamma(shape2, 1/rate2);
        gsEtaSq2(k) = hatEtaSq2;
        
        
    }
    
    return Rcpp::List::create(
                            Rcpp::Named("GS.b")=gsb,
                            Rcpp::Named("GS.beta") = gsEta,
                            Rcpp::Named("GS.tau") = gsTau,
                            Rcpp::Named("GS.v") = gsV,
                            Rcpp::Named("GS.s2") = gsSg2,
                            Rcpp::Named("GS.eta22.sq") = gsEtaSq2);
}
