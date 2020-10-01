//____________________________________________________________________________
/*!

\class    genie::COHGammaIntegrationLimits

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

#ifndef _COH_GAMMA_INTEGRATION_LIMITS_H_
#define _COH_GAMMA_INTEGRATION_LIMITS_H_

#include "Framework/Utils/Range1.h"
#include "Framework/Algorithm/Algorithm.h"

namespace genie {

class COHGammaIntegrationLimits : public Algorithm {

public:

  COHGammaIntegrationLimits();
  COHGammaIntegrationLimits(string config);

  virtual ~COHGammaIntegrationLimits();

  virtual Range1D_t EGamma( const Interaction & ) const ;
  virtual Range1D_t ThetaGamma( const Interaction & ) const ;
  virtual Range1D_t PhiGamma( const Interaction & ) const ;

  virtual Range1D_t ThetaLepton( const Interaction & ) const ;
  virtual Range1D_t t( const Interaction & ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  void LoadConfig(void);

 private:

};

}       // genie namespace
#endif  // _COH_GAMMA_INTEGRATION_LIMITS_H_ 
