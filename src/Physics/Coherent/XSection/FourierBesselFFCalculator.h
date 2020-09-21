//____________________________________________________________________________
/*!

\class    genie::FourierBesselFFCa

\brief    This class impelments the calculation of De Vries Form factor
          as a function of Q, for a given nucleus radius and a series of
          coeffictients.

\ref      Atom.Data Nucl.Data Tabl. 36 (1987) 495-536
          DOI: 10.1016/0092-640X(87)90013-1


\author  Marco Roda <mroda@liverpool.ac.uk>
         University of Liverpool

\created July 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FOURIER_BESSEL_FF_CALCULATOR_
#define _FOURIER_BESSEL_FF_CALCULATOR_

#include <vector>

namespace genie {

class FourierBesselFFCalculator {

public:

  FourierBesselFFCalculator( const std::vector<double> & coeffs, double radius, double Qmax, double Qmin ) noexcept :
    fFBCs(coeffs),
    fRadius(radius), 
    fQmax(Qmax),
    fQmin(Qmin) {;}

  FourierBesselFFCalculator( const FourierBesselFFCalculator & ) = default ;

  ~FourierBesselFFCalculator() = default ;

  double FormFactor( double Q ) const ;
  // The Q has to be in GeV
  // The returned FF is in fm^3

  const std::vector<double> & Coefficients() const noexcept { return fFBCs; }
  double Radius() const noexcept { return fRadius; }
  double QMax() const noexcept { return fQmax; }
  double QMin() const noexcept { return fQmin; }

protected:
  FourierBesselFFCalculator( ) noexcept : fFBCs(), fRadius(0.), fQmax(0.), fQmin(0.) {;}

private:

  std::vector<double> fFBCs ;  // Fourier-Bessel Coeffictients
  double fRadius ;   // this is the radius of the nucleus in GeV^-1
  double fQmax ;   // this is the max q for the DeVries coeffs in GeV^-2
  double fQmin ;   // this is the min q for the DeVries coeffs in GeV^-2

};

}       // genie namespace
#endif  // _FOURIER_BESSEL_FF_CALCULATOR_
