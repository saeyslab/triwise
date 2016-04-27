#include "RcppArmadillo.h"
#include "RcppArmadilloExtensions/sample.h"
#include <math.h>
#include <cmath>
#include "prob.hpp"

using namespace Rcpp;

//' @title
//' backgroundModel2
//' @description Randomly samples a set of `samplesize` genes `nsamples` times, and return the distance to the origin `z` and the mean direction `dir` as a unit vector.
//' @param angles
//' @param nsamples
//' @param samplesize
//' @param nangles
//' @param bw
//'
// [[Rcpp::depends("RcppArmadillo")]]
//' @export
// [[Rcpp::export]]
List backgroundModel2(NumericMatrix P, int nsamples, int samplesize) {
  NumericVector sample(samplesize);
  NumericVector Psample(nsamples);
  NumericVector Z(nsamples);
  IntegerVector samplefrom(P.nrow());
  samplefrom = seq_len(samplefrom.length());

  int dim = P.ncol();

  NumericMatrix D(nsamples,dim);

  NumericVector p(dim);

  for (int i=0;i<nsamples; ++i) {
    std::fill(p.begin(), p.end(), 0);
    sample = RcppArmadillo::sample(samplefrom, samplesize, true);
    for (int j=0;j<samplesize; ++j) {
      for (int k=0;k<dim;++k){
        p(k) += P(sample(j),k);
      }
    }

    Z(i) = sqrt(sum(pow(p, 2)))/(float)samplesize;

    D(i,_) = p/Z(i)/(float)samplesize;
  }

  return List::create(
    Rcpp::Named("z") = Z,
    Rcpp::Named("dir") = D
  );
}
