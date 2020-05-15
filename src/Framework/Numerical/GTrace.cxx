//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 8, 2016 - CA
   Added Cholesky's method functions
*/
//____________________________________________________________________________

#include <float.h>

#include <TMath.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"

//____________________________________________________________________________
TMatrixD genie::utils::math::CholeskyDecomposition(const TMatrixD& cov_matrix)
{
// Perform a Cholesky decomposition of the input covariance matrix and
// return the lower triangular matrix
//
   const double epsilon = 1E-12;

   int ncols = cov_matrix.GetNcols();
   int nrows = cov_matrix.GetNrows();

   assert(ncols==nrows);

   int n = nrows;

   TMatrixD L(n, n);

   for (int i = 0; i < n; ++i) {

     // calculate the diagonal term first
     L(i,i) = cov_matrix(i,i);
     for (int k = 0; k < i; ++k) {
        double tmp = L(k,i);
        L(i,i) -= tmp*tmp;
     }//k

     if(L(i,i) <= 0) {
       if(fabs(L(i,i)) < epsilon){
         L(i,i)=epsilon;
         LOG("Cholesky", pINFO) 
           << "Changed element (" << i << ", " << i << ") to " << L(i,i);
       }
       else{
         LOG("Cholesky", pERROR) 
            << "Decomposed covariance matrix not positive-definite";
         LOG("Cholesky", pERROR) 
            << "L(" << i << "," << i << ") = " << L(i,i);
         exit(1);
       }
     }
     L(i,i) = TMath::Sqrt(L(i,i));
     // then the off-diagonal terms
     for (int j = i+1; j < n; ++j) {
        L(i,j) = cov_matrix(i,j);
        for (int k = 0; k < i; ++k) {
           L(i,j) -= L(k,i)*L(k,j);
        }
        L(i,j) /= L(i,i);
     }//j
  }//i

  // create the transpose of L
  TMatrixD LT(TMatrixD::kTransposed,L);

  return LT;
}
//____________________________________________________________________________
TVectorD  genie::utils::math::CholeskyGenerateCorrelatedParams (
    const TMatrixD& cholesky_triangular, TVectorD& mean_params) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();

  if(ncols != nrows) {
    LOG("Cholesky", pERROR) 
        << "Mismatch between number of columns (" << ncols 
        << ") & rows (" << nrows << ")";
    exit(1);
  }
  if(npars != nrows) {
    LOG("Cholesky", pERROR) 
        << "Mismatch between number of parameters (" << npars
        << ") & array size (" << nrows << ")";
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = RandomGen::Instance()->RndNum().Gaus();
  }
  g *= cholesky_triangular;

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g(i);
     correlated_params[i] = v;
  }

  return correlated_params;
}
//____________________________________________________________________________
TVectorD  genie::utils::math::CholeskyGenerateCorrelatedParams (
 const TMatrixD& cholesky_triangular, TVectorD& mean_params, TVectorD& g_uncorrelated) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();
  int nunco = g_uncorrelated.GetNrows();

  if(ncols != nrows) {
    LOG("Cholesky", pERROR) 
        << "Mismatch between number of columns (" << ncols 
        << ") & rows (" << nrows << ")";
    exit(1);
  }
  if(npars != nrows) {
    LOG("Cholesky", pERROR) 
        << "Mismatch between number of parameters (" << npars
        << ") & array size (" << nrows << ")";
    exit(1);
  }
  if(nunco != nrows) {
    LOG("Cholesky", pERROR) 
        << "Mismatch between size of uncorrelated parameter vector (" << nunco
        << ") & array size (" << nrows << ")";
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  g_uncorrelated *= cholesky_triangular;

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g_uncorrelated(i);
     correlated_params[i] = v;
  }

  return correlated_params;
}
//____________________________________________________________________________
TVectorD  genie::utils::math::CholeskyGenerateCorrelatedParamVariations (
    const TMatrixD& cholesky_triangular) 
{
  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();

  assert(ncols==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = RandomGen::Instance()->RndNum().Gaus();
  }
  g *= cholesky_triangular;

  return g;
}
//____________________________________________________________________________
TVectorD  genie::utils::math::CholeskyCalculateCorrelatedParamVariations (
    const TMatrixD& cholesky_triangular, TVectorD & g_uncorrelated) 
{
  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = g_uncorrelated.GetNrows();

  assert(ncols==nrows);
  assert(npars==nrows);

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(g_uncorrelated);
  g *= cholesky_triangular;

  return g;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(double x[], unsigned int n)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<n; i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(const vector<double> & x)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<x.size(); i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(double x1, double x2)
{
  double err = 0.001*DBL_EPSILON;
  double dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(float x1, float x2)
{
  float err = FLT_EPSILON;
  float dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(double x, Range1D_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(float x, Range1F_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(int i, Range1I_t range)
{
  return ( i >= range.min && i <= range.max );
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(double x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max(0., x);
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(float x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max( (float)0., x);
}
//____________________________________________________________________________

