#include "RcppArmadillo.h"
#include "RcppArmadilloExtensions/sample.h"
#include <math.h>
#include <cmath>
#include "prob.h"

using namespace Rcpp;

//' @title backgroundModel2
//' @description Generates a triwise background model, containing weights, z-values and angle p-values for a set of uniformly distributed reference angles
//' @param angles
//' @param nsamples
//' @param samplesize
//' @param nangles
//' @param bw
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List backgroundModel2(NumericVector angles, NumericVector R, int nsamples, int samplesize, NumericVector anglesoi, double bw) {
  NumericVector sample(samplesize);
  NumericVector A(nsamples);
  NumericVector Z(nsamples);
  IntegerVector samplefrom(angles.size());
  samplefrom = seq_len(angles.size());
  double x;
  double y;
  double z;
  double basepval;
  basepval = 0;
  int nangles = anglesoi.size();

  for (int i=0;i<nsamples; ++i) {
    sample = RcppArmadillo::sample(samplefrom, samplesize, true);
    x = 0;
    y = 0;
    z = 0;
    for (int j=0;j<samplesize; ++j) {
      x += cos(angles(sample(j)-1)) * R(sample(j)-1);
      y += sin(angles(sample(j)-1)) * R(sample(j)-1);
    }
    A(i) = atan2(y, x);
    Z(i) = sqrt(x*x + y*y)/samplesize;

    if (Z(i) > 0) {
      basepval += 1;
    }
  }

  NumericMatrix W(nangles, nsamples);
  NumericVector AD(nangles); // average angle density
  double angle;
  double delta;
  for (int i=0;i<nangles; ++i) {
    angle = anglesoi(i);
    for (int j=0;j<nsamples; ++j) {
      if (Z(j) > 0) { // zero weight if angle could not be determined
        delta = angle - A(j);
        if (delta > M_PI) {
          delta = delta - 2*M_PI;
        }
        // TODO: pre-calculate the expensive bessel function here to speed up everything
        W(i, j) = von_mises_pdf(delta, 0, bw);
        AD(i) += W(i,j);
      }
    }

    for (int j=0;j<nsamples; ++j) {
      W(i,j) = W(i,j)/AD(i); // make sure weights sum to 1
    }
    AD(i) = AD(i)/(float)basepval;
  }


  // now normalize every weight based on the density of its nearest angle of interest
  // this is inactive for now
  double density;
  int id;
  for (int i=0;i<nsamples; ++i) {
    // get the id of the nearest angle and get its density
    id = A(i) / (2 * M_PI) * nangles;
    if (id < 0) {
      id += nangles;
    }
    density = AD(id);

    for (int j=0;j<nangles; ++j) {
      //W(j, i) /= density;
    }
  }

  // normalize the weights of every angle so that they sum to one
  double wsum;
  for (int i=0;i<nangles; ++i) {
    wsum = 0;
    for (int j=0;j<nsamples; ++j) {
      wsum += W(i,j);
    }
    for (int j=0;j<nsamples; ++j) {
      W(i,j) /= wsum;
    }
  }


  // now calculate the p-values of a particular angle using the trapezoid rule (for speed)
  NumericVector AP(nangles);
  float maxdens;
  int k;
  float dens1;
  float dens2;
  float deltangle;
  for (int i=0;i<nangles; ++i) {
    maxdens = AD(i);
    for (int j=0;j<nangles; ++j) {
      if (j == nangles-1) {
        k = 0;
      } else {
        k = j+1;
      }

      // densities
      dens1 = AD(j);
      if (dens1 > maxdens) {
        dens1 = maxdens;
      }
      dens2 = AD(k);
      if (dens2 > maxdens) {
        dens2 = maxdens;
      }

      // difference in angle
      deltangle = (anglesoi(k) - anglesoi(j));
      if (deltangle < 0) {
        deltangle += 2 * M_PI;
      }

      AP(j) += deltangle * (dens2 + dens1)/2;
    }
  }

  return List::create(
    Rcpp::Named("angles") = A,
    Rcpp::Named("z") = Z,
    Rcpp::Named("weights") = W,
    Rcpp::Named("anglesp") = AP,
    Rcpp::Named("anglesd") = AD,
    Rcpp::Named("anglesoi") = anglesoi,
    Rcpp::Named("basepval") = basepval/nsamples
  );
}
