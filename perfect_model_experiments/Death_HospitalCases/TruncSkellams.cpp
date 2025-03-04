// [[Rcpp::depends(sitmo, RcppArmadillo)]]

#include <iostream>
#include <stdlib.h>
#include <RcppArmadillo.h>
#include <sitmo.h>

// headers for bessel functions
#include "bessel.h"

using namespace Rcpp;

// log-sum-exp function to prevent numerical overflow
double log_sum_exp(arma::vec &x, int mn = 0) {
  double maxx = max(x);
  double y = maxx + log(sum(exp(x - maxx)));
  if(mn == 1) y -= log(x.n_elem);
  return y;
}

// log Skellam CDF (note, no checks on inputs)
// adapted from "pskellam" source code by Patrick Brown (mistakes are mine)
// [[Rcpp::export]]
double lpskellam_cpp(int x, double lambda1, double lambda2) {
  // different formulas for negative & nonnegative x (zero lambda is OK)
  double ldens = 0.0;
  if(x < 0) {
    ldens = R::pnchisq(2.0 * lambda2, -2.0 * x, 2.0 * lambda1, 1, 1);
  } else {
    ldens = R::pnchisq(2.0 * lambda1, 2.0 * (x + 1), 2.0 * lambda2, 0, 1);
  }
  //    if(!arma::is_finite(ldens)) {
  //        Rprintf("Non-finite pdensity in lpskellam = %f x = %d lambda1 = %f lambda2 = %f\n", ldens, x, lambda1, lambda2);
  //    }
  return ldens;
}

// log truncated Skellam density (note, no checks on inputs)
// adapted from "dskellam" source code by Patrick Brown (mistakes are mine)
// [[Rcpp::export]]
double ldtskellam_cpp(int x, double lambda1, double lambda2, int LB = 0, int UB = -1) {
  
  // if UB < LB then returns untruncated density
  
  // check for this instance only
  if(LB > UB) {
    stop("ldt LB > UB\n");
  }
  
  // create array for Bessel (trying to deal with recursive
  // gc error message that may be coming from memory allocation
  // in bessel_i - hence swapping to bessel_i_ex and controlling
  // vector allocation directly here through calloc and free -
  // OK to use std::calloc here because object not going
  // back to R I think)
  double *bi = (double *) calloc (floor(abs(x)) + 1, sizeof(double));
  
  // from "skellam" source code: log(besselI(y, nu)) == y + log(besselI(y, nu, TRUE))
  double ldens = -(lambda1 + lambda2) + (x / 2.0) * (log(lambda1) - log(lambda2)) +
    log(R::bessel_i_ex(2.0 * sqrt(lambda1 * lambda2), abs(x), 2, bi)) + 2.0 * sqrt(lambda1 * lambda2);
  
  // free memory from the heap
  free(bi);
  
  //    // if non-finite density, then try difference of CDFs
  //    if(!arma::is_finite(ldens)) {
  //        Rprintf("Trying difference of CDFs in ldtskellam\n");
  //        double norm0 = 0.0, norm1 = 0.0;
  //        norm0 = lpskellam_cpp(x - 1, lambda1, lambda2);
  //        norm1 = lpskellam_cpp(x, lambda1, lambda2);
  //        if(!arma::is_finite(norm0) && !arma::is_finite(norm1)) {
  //            ldens = norm0;
  //        } else {
  //            ldens = norm1 + log(1.0 - exp(norm0 - norm1));
  //        }
  //    }
  
  // if truncated then adjust log-density
  if(UB >= LB) {
    if(x < LB || x > UB) {
      Rprintf("'x' must be in bounds in dtskellam_cpp: x = %d LB = %d UB = %d\n", x, LB, UB);
      stop("");
    }
    if(arma::is_finite(ldens)) {
      // declare variables
      int normsize = UB - LB + 1;
      if(normsize > 1) {
        double norm0 = 0.0, norm1 = 0.0;
        norm0 = lpskellam_cpp(LB - 1, lambda1, lambda2);
        norm1 = lpskellam_cpp(UB, lambda1, lambda2); 
        if(!arma::is_finite(norm0) && !arma::is_finite(norm1)) {
          stop("Issue with truncation bounds in Skellam\n");
        }
        double lnorm = norm1 + log(1.0 - exp(norm0 - norm1));
        ldens -= lnorm;
      } else {
        ldens = 0.0;
      }
    }
  }
  //    if(!arma::is_finite(ldens)) {
  //        Rprintf("Non-finite density in ldtskellam = %f x = %d l1 = %f l2 = %f LB = %d UB = %d\n", ldens, x, lambda1, lambda2, LB, UB);
  //    }
  return ldens;
}


// truncated Skellam sampler
// [[Rcpp::export]]
int rtskellam_cpp(double lambda1, double lambda2, int LB = 0, int UB = -1) {
  
  // if UB < LB then returns untruncated density
  
  // check in this setting:
  if(LB > UB) {
    stop("rdt LB > UB\n");
  }
  
  // declare variables
  int x = 0;
  //double mx = sitmo::prng::max();
  
  // if bounds are the same, then return the
  // only viable value
  if(LB == UB) {
    x = LB;
  } else {
    // draw untruncated sample
    // x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
    x = R::rpois(lambda1) - R::rpois(lambda2);
    
    // rejection sample if necessary
    if(UB > LB) {
      int k = 0;
      int ntries = 1000;
      while((x < LB || x > UB) && k < ntries) {
        // x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
        x = R::rpois(lambda1) - R::rpois(lambda2);
        k++;
      }
      // if rejection sampling doesn't work
      // then try inverse transform sampling
      if(k == ntries) {
        // Rprintf("LB = %d UB = %d\n", LB, UB);
        //double u = eng() / mx;
        NumericVector u = runif(1);
        k = 0;
        double xdens = exp(ldtskellam_cpp(LB + k, lambda1, lambda2, LB, UB));
        if(!arma::is_finite(xdens)) xdens = 0.0;
        while(u[0] > xdens && k < (UB - LB)) {
          k++;
          xdens += exp(ldtskellam_cpp(LB + k, lambda1, lambda2, LB, UB));
          if(!arma::is_finite(xdens)) xdens += 0.0;
        }
        if(k == (UB - LB + 1) && u[0] > xdens) {
          stop("Something wrong in truncated Skellam sampling\n");
        }
        x = LB + k;
      }
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericVector rtskellam_n_cpp(DoubleVector lambda1, DoubleVector lambda2, IntegerVector LB, IntegerVector UB){
  int n = lambda1.size();
  NumericVector x(n);
  for(int i = 0; i<n; ++i){
    x[i] = rtskellam_cpp(lambda1[i], lambda2[i], LB[i], UB[i]);
  }
  return x;
}

// [[Rcpp::export]]
NumericVector ldtskellam_n_cpp(IntegerVector y, DoubleVector lambda1, DoubleVector lambda2, IntegerVector LB, IntegerVector UB){
  int n = lambda1.size();
  NumericVector x(n);
  for(int i = 0; i<n; ++i){
    x[i] = ldtskellam_cpp(y[i],lambda1[i], lambda2[i], LB[i], UB[i]);
  }
  return x;
}