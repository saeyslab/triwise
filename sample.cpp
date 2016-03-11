#include "RcppArmadillo.h"
#include "RcppArmadilloExtensions/sample.h"
#include <math.h>
#include <cmath>
#include "prob.hpp"

using namespace Rcpp;

//' @title
//' backgroundModel
//' @description Generates a triwise background model, containing weights, z-values and angle p-values for a set of uniformly distributed reference angles
//' @param angles
//' @param nsamples
//' @param samplesize
//' @param nangles
//' @param bw
//'
// [[Rcpp::depends("RcppArmadillo")]]
//' @export
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
      x += cos(sample(j));
      y += sin(sample(j));
    }
    A(i) = atan2(y, x);
    Z(i) = sqrt(x*x + y*y)/samplesize;
  }

  NumericVector anglesoi(nangles);
  NumericMatrix W(nangles, nsamples);
  NumericVector AD(nangles); // average angle density
  double angle;
  double delta;
  for (int i=0;i<nangles; ++i) {
    angle = (float)i/(float)nangles * 2 * M_PI;
    anglesoi(i) = angle;
    for (int j=0;j<nsamples; ++j) {
      delta = angle - A(j);
      if (delta > M_PI) {
        delta = delta - 2*M_PI;
      }
      W(i, j) = von_mises_pdf(delta, 0, bw);
      AD(i) += W(i,j);
    }

    for (int j=0;j<nsamples; ++j) {
      //W(i,j) = W(i,j)/AD(i); // make sure weights sum to 1, is not necessary because we will adjust the weights later on anyway
    }
    AD(i) = AD(i)/(float)nsamples;
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
    //Rcout << "The value is " << density << std::endl;
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

      //Rcout << "The value is " << deltangle << std::endl;

      AP(j) += deltangle * (dens2 + dens1)/2;
    }
  }

  return List::create(
    Rcpp::Named("angles") = A,
    Rcpp::Named("z") = Z,
    Rcpp::Named("weights") = W,
    Rcpp::Named("anglep") = AP,
    Rcpp::Named("angled") = AD,
    Rcpp::Named("anglesoi") = anglesoi
  );
}
