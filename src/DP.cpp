// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadilloExtensions/sample.h>
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double rtruncnorm(double y, double mu, double sigma, int K=1){
  NumericVector yy(1);
  yy(0) = (y-mu)/sigma;
  NumericVector cdf = pnorm(yy, 0, 1, true, false);
  NumericVector temp = rep(cdf, K)+runif(K)*rep((1-cdf),K);
  NumericVector draws = rep(mu, K)+rep(sigma, K)*qnorm(temp);
  return mean(draws) == arma::datum::inf ? y : mean(draws);
}

// [[Rcpp::export]]
NumericMatrix vec2mat(NumericVector x, int nb) {
  NumericMatrix toRet(nb, x.size());
  for (int i = 0; i < nb; i++) {
    toRet.row(i) = x;
  }
  return toRet;
}

NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}


NumericVector colSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(ncol);
  
  for (int i = 0; i < ncol; i++) {
    double total = 0;
    for (int j = 0; j < nrow; j++) {
      total += x(j, i);
    }
    out[i] = total;
  }
  return out;
}


NumericVector MaxPerCol(NumericMatrix X){
  NumericVector toRet(X.ncol());
  for(int i = 0; i < X.ncol(); i++){
    NumericVector Xi = X.column(i);
    toRet[i] = Rcpp::max(Xi);
  }
  return toRet;
}

NumericMatrix elementWiseSubstraction(NumericMatrix X, NumericMatrix Y) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericMatrix out(nrow, ncol);
  
  for (int i = 0; i < ncol; i++) {
    for (int j = 0; j < nrow; j++) {
      out(j,i) = X(j,i) - Y(j,i) ;
    }
  }
  return out;
}



NumericMatrix elementExp(NumericMatrix X) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericMatrix out(nrow, ncol);
  
  for (int i = 0; i < ncol; i++) {
    for (int j = 0; j < nrow; j++) {
      double Xji = X(j,i);
      out(j,i) = exp(Xji);
    }
  }
  return out;
}

NumericVector dropNA(NumericVector data, int mask){
  NumericVector out(mask);
  for(int i = 0; i<mask; i++){
    out(i) = data(i);
  }
  return out;
}



// [[Rcpp::export]]
NumericMatrix eStep(NumericVector theta, NumericVector phi, NumericVector w, S4 DataStorage) {
  NumericVector censoring = DataStorage.slot("censoring");
  NumericVector data = DataStorage.slot("computation");
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


// [[Rcpp::export]]
S4 mStep(S4 DP, S4 DataStorage, IntegerVector xi, IntegerVector zeta){
  NumericVector data = DataStorage.slot("simulation");
  const int N = data.length();
  
  arma::cube prior = DP.slot("prior");
  
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
  
  DP.slot("prior") = posterior;
  return DP;
}

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

/*** R
if(F){
  library(DPsurv)
  
  data <- sim.data()
  G1 <- new("DP")
  G1 <- init.DP(G1, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=55)
  G1 <- MCMC.DP(G1, data, 50)
  plot.ICDF(G1@theta, G1@phi, G1@weights, G1@L, grid=0:50,
            distribution=data@presentation, xlim=50)
  
  G2 <- new("NDP")
  G2 <- init.NDP(G2, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=5, L=15)
  
  G3 <- new("HDP")
  G3 <- init.HDP(G3, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=15, J=2)
  G3 <- MCMC.HDP(G3, data, 50)
  plot.ICDF(G3@theta, G3@phi, G3@weights[,1], G3@L, grid=0:50,
            distribution=data@presentation, xlim=50)
  set.seed(123)
  data <- sim.data()
  
  set.seed(123)
  G2 <- new("NDP")
  G2 <- init.NDP(G2, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=3, L=55)
  G2 <- MCMC.NDP(G2, data, 25)
  plot.ICDF(G2@theta[,2], G2@phi[,2], G2@weights[,2], G2@L, grid=0:50,
            distribution=data@presentation, xlim=50)
  G3 <- init.HDP(G3, prior=list(mu=0, n=0.1, v=3, vs2=1*3), 55, 2)
  profvis({MCMC.NDP(G2,data, 5000)})
}
print("dne")
*/
