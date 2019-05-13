//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Boyd ( s.b.boyd \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <cstdlib>
#include <complex>

#include <TMath.h>

#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/IntegrationTools.h"

namespace genie {
namespace alvarezruso {

typedef std::complex<double> cdouble;

// All numerical values from Abscissas and weights for Gaussian quadratures of high order (1956)
// https://archive.org/details/jresv56n1p35

//
// Routine to work out the evaluation points for a Gaussian integral given the
// end points, and the number of sampling points
//
void integrationtools::SG20R(const double a, const double b, const unsigned int n, const unsigned int nsamp, 
           double* x, unsigned int& np, double* /*w*/)
{
  static const double y[10] = {.9931285991, .9639719272, .9122344282, .8391169718, .7463319064, .6360536807, .5108670019, .3737060887, .2277858511, .0765265211 };
  np = 20 * n;
  double dint = (b - a) / double(n);
  double delt = dint * 0.5;
  double orig = a - delt;
  int i1 = -nsamp;
  int i2, j1, j2;
  double dorig;
  for(unsigned int i = 1; i <= n; i++)
  {
    orig += dint;
    dorig = orig + orig;
    i1 += nsamp;
    i2 = i1 + nsamp+1;
    for(unsigned int j = 1; j <= 10; j++)
    {
      j1 = i1 + j;
      j2 = i2 - j;
      x[j1-1] = orig - delt * y[j-1];
      x[j2-1] = dorig - x[j1-1];
    }
  }
}

//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
cdouble integrationtools::RG201D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const cdouble CF[])
{
  // Gaussian-Legendre integration of the function defined by CF
  const double W[10] = {.0176140071,.0406014298,.0626720483,.0832767415,.1019301198,.1181945319,.1316886384,.1420961093,.1491729864,.1527533871};
  cdouble CR(0.0, 0.0);
  int I1 = -nsamp;
  int I2, J1, J2;
  for(unsigned int i = 1; i <= N; ++i)
  {
    I1 += nsamp;
    I2 = I1 + nsamp+1;
    for(unsigned int j = 1; j <= 10; ++j)
    {
      J1 = I1 + j;
      J2 = I2 - j;
      CR += W[j-1] * (CF[J1-1]+CF[J2-1]);
    }
  } 
  //CRES=CR*0.5*(B-A)/float(N)
  cdouble CRES = CR*0.5*(B-A)/Double_t(N);
  return CRES;
}
//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
void integrationtools::RG202D(const double a, const double b, unsigned int n, unsigned int l, 
              unsigned int m, std::vector< std::vector<cdouble> >& cf, 
              const unsigned int nsamp, std::vector<cdouble>& cres)
{
  // This is a fast integrator based on a Gauss-Legendre method. This only support two-dimensional integration
  n = 2; l = 0; m = 3;

  static const double w[10] = {.0176140071,.0406014298,.0626720483,.0832767415,.1019301198,.1181945319,.1316886384,.1420961093,.1491729864,.1527533871};

  std::vector<cdouble> cr(4, cdouble(0.0,0.0));

  int i1 = -nsamp;
  int i2;
  int j1;
  int j2;

  for(unsigned int i = 0; i != n; ++i)
  {
    i1 += nsamp; 
    i2 = i1 + nsamp-1;
    for(unsigned int j = 0; j != 10; ++j)
    {
      j1 = i1 + j;
      j2 = i2 - j;

      for(unsigned int ll = l; ll <= m; ++ll)
      {
        cr[ll] += w[j] * ( cf[ll][j1] + cf[ll][j2] );
      }
    }
  }

  for(unsigned int i = 0; i != 4; ++i)
  {
    cres[i] = cr[i] * 0.5 * (b-a) / static_cast<double>(n);
  }
}
//-----------------------------------------------------------------------------------------------------------
//
// Routine to work out the evaluation points for a Gaussian integral given the
// end points, and the number of sampling points
//
void integrationtools::SG48R(const double a, const double b, const unsigned int n, const unsigned int nsamp, 
           double* x, unsigned int& np, double* /*w*/)
{
  static const double y[24] = {0.9987710072, 0.9935301722, 0.9841245873, 0.9705915925, 0.9529877031,
   0.9313866907, 0.9058791367, 0.8765720202, 0.8435882616, 0.8070662040, 0.7671590325, 0.7240341309,
   0.6778723796, 0.6288673967, 0.5772247260, 0.5231609747, 0.4669029047, 0.4086864819, 0.3487558862,
   0.2873624873, 0.2247637903, 0.1612223560, 0.0970046992, 0.0323801709};
  np = 48 * n;
  double dint = (b - a) / double(n);
  double delt = dint * 0.5;
  double orig = a - delt;
  int i1 = -nsamp;
  int i2, j1, j2;
  double dorig;
  for(unsigned int i = 1; i <= n; i++)
  {
    orig += dint;
    dorig = orig + orig;
    i1 += nsamp;
    i2 = i1 + nsamp+1;
    for(unsigned int j = 1; j <= 24; j++)
    {
      j1 = i1 + j;
      j2 = i2 - j;
      x[j1-1] = orig - delt * y[j-1];
      x[j2-1] = dorig - x[j1-1];
    }
  }
}
//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
cdouble integrationtools::RG481D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const cdouble CF[])
{
  // Gaussian-Legendre integration of the function defined by CF
  static const double W[24] = {0.0031533460, 0.0073275539, 0.0114772345, 0.0155793157, 0.0196161604,
   0.0235707608, 0.0274265097, 0.0311672278, 0.0347772225, 0.0382413510, 0.0415450829, 0.0446745608,
   0.0476166584, 0.0503590355, 0.0528901894, 0.0551995036, 0.0572772921, 0.0591148396, 0.0607044391,
   0.0620394231, 0.0631141922, 0.0639242385, 0.0644661644, 0.0647376968 };
  cdouble CR(0.0, 0.0);
  int I1 = -nsamp;
  int I2, J1, J2;
  for(unsigned int i = 1; i <= N; ++i)
  {
    I1 += nsamp;
    I2 = I1 + nsamp+1;
    for(unsigned int j = 1; j <= 24; ++j)
    {
      J1 = I1 + j;
      J2 = I2 - j;
      CR += W[j-1] * (CF[J1-1]+CF[J2-1]);
    }
  } 
  //CRES=CR*0.5*(B-A)/float(N)
  cdouble CRES = CR*0.5*(B-A)/Double_t(N);
  return CRES;
}
//-----------------------------------------------------------------------------------------------------------
// Gaussian-Legendre integration of the function defined by CF
void integrationtools::RG482D(const double a, const double b, unsigned int n, unsigned int l, 
              unsigned int m, std::vector< std::vector<cdouble> >& cf, 
              const unsigned int nsamp, std::vector<cdouble>& cres)
{
  // This is a fast integrator based on a Gauss-Legendre method. This only support two-dimensional integration
  n = 2; l = 0; m = 3;

  static const double w[24] = {0.0031533460, 0.0073275539, 0.0114772345, 0.0155793157, 0.0196161604,
   0.0235707608, 0.0274265097, 0.0311672278, 0.0347772225, 0.0382413510, 0.0415450829, 0.0446745608,
   0.0476166584, 0.0503590355, 0.0528901894, 0.0551995036, 0.0572772921, 0.0591148396, 0.0607044391,
   0.0620394231, 0.0631141922, 0.0639242385, 0.0644661644, 0.0647376968 };

  std::vector<cdouble> cr(4, cdouble(0.0,0.0));

  int i1 = -nsamp;
  int i2;
  int j1;
  int j2;

  for(unsigned int i = 0; i != n; ++i)
  {
    i1 += nsamp; 
    i2 = i1 + nsamp-1;
    for(unsigned int j = 0; j != 24; ++j)
    {
      j1 = i1 + j;
      j2 = i2 - j;

      for(unsigned int ll = l; ll <= m; ++ll)
      {
        cr[ll] += w[j] * ( cf[ll][j1] + cf[ll][j2] );
      }
    }
  }

  for(unsigned int i = 0; i != 4; ++i)
  {
    cres[i] = cr[i] * 0.5 * (b-a) / static_cast<double>(n);
  }
}

//______________________________________________________________________
// Calls correct integration tool for current sampling
void integrationtools::SGNR(const double a, const double b, const unsigned int n, const unsigned int nsamp, 
                            double* x, unsigned int& np, double* w)
{
  if (nsamp==20) {
    integrationtools::SG20R(a, b, n, nsamp, x, np, w);
  }
  else if (nsamp==48) {
    integrationtools::SG48R(a, b, n, nsamp, x, np, w);
  }
  else {
    std::cerr<<"IntegrationTools/SGNR >> Unsupported nuclear sampling: "<<nsamp<<std::endl;
    exit(-1);
  }
}

cdouble integrationtools::RGN1D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const cdouble CF[])
{
  if (nsamp==20) {
    return integrationtools::RG201D(A, B, N, nsamp, CF);
  }
  else if (nsamp==48) {
    return integrationtools::RG481D(A, B, N, nsamp, CF);
  }
  else {
    std::cerr<<"IntegrationTools/RGN1D >> Unsupported nuclear sampling: "<<nsamp<<std::endl;
    exit(-1);
  }
}

void integrationtools::RGN2D (const double a, const double b, unsigned int n, unsigned int l, 
              unsigned int m, std::vector< std::vector<cdouble> >& cf, 
              const unsigned int nsamp, std::vector<cdouble>& cres)
{
  if (nsamp==20) {
    integrationtools::RG202D(a, b, n, l, m, cf, nsamp, cres);
  }
  else if (nsamp==48) {
    integrationtools::RG482D(a, b, n, l, m, cf, nsamp, cres);
  }
  else {
    std::cerr<<"IntegrationTools/RGN2D >> Unsupported nuclear sampling: "<<nsamp<<std::endl;
    exit(-1);
  }
}

} // alvarezruso namespace
} // genie namespace
