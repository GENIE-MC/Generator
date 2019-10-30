//____________________________________________________________________________
/*!

\class    genie::DeVriesFormFactor

\brief    De Vries Form factor interfaces for COH Production model
          The class is develope specifically for the NC COH Gamma
          But in principle these Form Factors could be reused.

\ref      Atom.Data Nucl.Data Tabl. 36 (1987) 495-536
          DOI: 10.1016/0092-640X(87)90013-1


\author  Marco Roda <mroda@liverpool.ac.uk>
         University of Liverpool

\created October 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DEVRIES_FORM_FACTOR_H_
#define _DEVRIES_FORM_FACTOR_H_

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

class DeVriesFormFactor : public Algorithm {

public:
  DeVriesFormFactor();
  DeVriesFormFactor(string config);
  virtual ~DeVriesFormFactor();

  int NucleusPDG() const { return fNuclPDG ; }
  double FormFacotor( double Q ) const ;
  // The Q has to be in GeV
  // The returned FF is in fm^3

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);


private:

  void LoadConfig(void);

  std::vector<double> fFBCs ;  // Fourier-Bessel Coeffictients
  double fRadius ;   // this is the radius of the nucleus in GeV^-1
  int fPDG ;


};

}       // genie namespace
#endif  // _COHERENT_ELASTIC_PXSEC_H_
