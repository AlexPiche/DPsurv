// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "RcppCommon.h"
using namespace Rcpp;


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//' @title
//' Random Number From a Truncated Normal
//' @description
//' Draw from a truncated normal
//' 
//' @param y is the value at which we want to left-truncated our normal distribution
//' @param mu is the mean of the truncated normal
//' @param sigma is the s.d. of the truncated normal
//' @param K is the number of replication of draws
//' @export
// [[Rcpp::export]]
double rtruncnorm(double y, double mu, double sigma, int K=1){
  NumericVector yy(1);
  yy(0) = (y-mu)/sigma;
  NumericVector cdf = pnorm(yy, 0, 1, true, false);
  NumericVector temp = rep(cdf, K)+runif(K)*rep((1-cdf),K);
  NumericVector draws = rep(mu, K)+rep(sigma, K)*qnorm(temp);
  return mean(draws) == arma::datum::inf ? y : mean(draws);
}


//' @title
//' eStep
//' @description
//' Determines to which cluster belongs an obervation, and a group of observations.
//' 
//' @param theta
//' @param phi
//' @param w
//' @param DataStorage
//' @export
// [[Rcpp::export]]
NumericMatrix eStep(NumericVector data, NumericVector censoring, NumericVector theta, NumericVector phi, NumericVector w) {
  //NumericVector censoring = DataStorage.slot("censoring");
  //NumericVector data = DataStorage.slot("computation");
  const int N = data.length();
  const int L = theta.length();
  NumericMatrix myMat(L, N);
  
  NumericVector log_w = log(w);
  NumericVector ones(N, 1.0);
  
  for(int l=0; l < L; l++){
    //function (x, mean = 0, sd = 1, log = FALSE)
    // function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    myMat(l,_) = log_w(l)+(censoring)*(dnorm(data, theta[l], std::sqrt(phi[l]), true))+(ones-censoring)*(pnorm(data, theta[l], std::sqrt(phi[l]), false, true));
  }
  
  return myMat;
}

//' @title
//' mStep
//' @description
//' Return the parameters' posterior.
//' 
//' @param DP
//' @param DataStorage
//' @param xi
//' @param zeta
//' @export
// [[Rcpp::export]]
arma::cube mStep(arma::cube prior, NumericVector data, IntegerVector xi, IntegerVector zeta){
  //NumericVector data = DataStorage.slot("simulation");
  const int N = data.length();
  
  //arma::cube prior = DP.slot("prior");
  
  arma::mat mu = prior.slice(0);
  arma::mat Nmat = prior.slice(1);
  arma::mat v = prior.slice(2);
  arma::mat vs2 = prior.slice(3);
  
  std::fill(v.begin(), v.end(), 3.0);
  std::fill(vs2.begin(), vs2.end(), 3.0);
  std::fill(Nmat.begin(), Nmat.end(), 0.1);
  
  int temp;
  double update_mu;
  double update_vs2;
  double update_n;
  double update_v;
  
  for(int n=0; n < N; n++){
    // if n is not NA
    if(data(n) == data(n)){
      temp = xi[n];
      const int l = temp - 1;
      temp = zeta[n];
      const int k = temp - 1;
      
      update_mu = (Nmat(l,k)*mu(l,k)+data(n))/(1+Nmat(l,k));
      update_vs2 = vs2(l,k) + (Nmat(l,k)*std::pow((mu(l,k)-data(n)),2)/(1+Nmat(l,k)));
      update_n = Nmat(l,k) + 1.0;
      update_v = v(l,k) + 1.0;
      mu(l,k) = update_mu;
      Nmat(l,k) = update_n;
      v(l,k) = update_v;
      vs2(l,k) = update_vs2;
    }
  }
  arma::cube posterior = prior;
  posterior.slice(0) = mu;
  posterior.slice(1) = Nmat;
  posterior.slice(2) = v;
  posterior.slice(3) = vs2;
  
  //DP.slot("prior") = posterior;
  return posterior;
}

//' @title
//' Augment censored data using a Gibbs step
//' @description
//' Augment censored data by drawing them from a truncated normal
//' 
//' @param DP is an S4 object of type DP, HDP, or NDP.
//' @param DataStore is an S4 object of the same name.
//' @param xi is an integer vector that describes to which cluster belong an observation
//' @param zeta is an integer vector that describes in which cluster belong a group of observations
//' @export
// [[Rcpp::export]]
S4 gibbsStep(S4 DP, S4 DataStorage, IntegerVector xi, IntegerVector zeta){
  NumericVector RealData = DataStorage.slot("computation");
  IntegerVector censoring = DataStorage.slot("censoring");
  const int N = RealData.length();
  NumericMatrix Theta = DP.slot("theta");
  NumericMatrix Phi = DP.slot("phi");
  NumericVector toRetData(N, NumericVector::get_na());
  int temp;
  
  for(int n=0; n < N; n++){
    // if n is not NA
    if(RealData(n) == RealData(n)){
        temp = xi[n];
        const int l = temp - 1;
        
        temp = zeta[n];
        const int k = temp - 1;
        
        toRetData(n) = censoring(n) ? RealData(n) : rtruncnorm(RealData(n), Theta(l,k), std::sqrt(Phi(l,k)));
    }
  }
  DataStorage.slot("simulation") = toRetData;
  return DataStorage;
}