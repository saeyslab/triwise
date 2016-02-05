#include "RcppArmadillo.h"
#include "RcppArmadilloExtensions/sample.h"
#include <math.h>
#include "prob.hpp"

using namespace Rcpp;

//'
//' @export
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List backgroundModel(NumericVector angles, int nsamples, int samplesize, int nangles, double bw) {
  NumericVector sample(samplesize);
  NumericVector A(nsamples);
  NumericVector Z(nsamples);
  double x;
  double y;

  for (int i=0;i<nsamples; ++i) {
    sample = RcppArmadillo::sample(angles, samplesize, true);
    x = 0;
    y = 0;
    for (int j=0;j<samplesize; ++j) {
      x += cos(sample[j]);
      y += sin(sample[j]);
    }
    A[i] = atan2(y, x);
    Z[i] = sqrt(x*x + y*y);
  }

  NumericMatrix W(nangles, nsamples);
  double angle;
  for (int i=0;i<nangles; ++i) {
    angle = i/nangles * 2 * M_PI;
    for (int j=0;j<nsamples; ++j) {
      W[i, j] = von_mises_pdf(angle, A[j], bw);
    }
  }


  return List::create(Rcpp::Named("angles") = A, Rcpp::Named("z") = Z, Rcpp::Named("weights") = W);
}
