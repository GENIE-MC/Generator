//____________________________________________________________________________
/*!

\namespace  genie::alvarezruso::integrationTools

\brief      Some fast integration tools

\author     Steve Boyd (s.b.boyd \at warwick.ac.uk)
            University of Warwick

\created    Oct 24th, 2013

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AR_INT_TOOLS_H_
#define _AR_INT_TOOLS_H_

#include <vector>
#include <complex>

namespace genie {
namespace alvarezruso {
namespace integrationtools
{
  void SG20R(const double a, const double b, const unsigned int n, const unsigned int nsamp, double* x, unsigned int& np, double* w);
  std::complex<double>  RG201D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const std::complex<double>  CF[]);
  void RG202D(const double a, const double b, unsigned int n, unsigned int l,
              unsigned int m, std::vector< std::vector<std::complex<double> > >& cf,
              const unsigned int nsamp, std::vector<std::complex<double> >& cres);
              
  void SG48R(const double a, const double b, const unsigned int n, const unsigned int nsamp, double* x, unsigned int& np, double* w);
  std::complex<double>  RG481D(const double A, const double B, const unsigned int N, const unsigned int nsamp, const std::complex<double>  CF[]);
  void RG482D(const double a, const double b, unsigned int n, unsigned int l,
              unsigned int m, std::vector< std::vector<std::complex<double> > >& cf,
              const unsigned int nsamp, std::vector<std::complex<double> >& cres);
              
  void SGNR (const double a, const double b, const unsigned int n, const unsigned int nsamp, double* x, unsigned int& np, double* w);
  std::complex<double>  RGN1D (const double A, const double B, const unsigned int N, const unsigned int nsamp, const std::complex<double>  CF[]);
  void RGN2D (const double a, const double b, unsigned int n, unsigned int l,
              unsigned int m, std::vector< std::vector<std::complex<double> > >& cf,
              const unsigned int nsamp, std::vector<std::complex<double> >& cres);

} // IntegrationTools namespace
} // alvarezruso namespace
} // genie namespace

#endif // AR_INT_TOOLS_H_
