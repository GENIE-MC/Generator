//____________________________________________________________________________
/*!

\class    genie::COHGammaIntegrationLimits

\brief    Algorithm to define ina single place the hyper-cube to be used to 
          integrate the complicated phase space of the COH Gamma cross sections

\author  Marco Roda <mroda@liverpool.ac.uk>
         University of Liverpool

\created October 2020

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_GAMMA_INTEGRATION_LIMITS_H_
#define _COH_GAMMA_INTEGRATION_LIMITS_H_

#include "Framework/Utils/Range1.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

#include "Physics/Coherent/XSection/COHFormFactorI.h" 

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

  const COHFormFactorI * fFF ;
  double fMaxEg ;

};

}       // genie namespace
#endif  // _COH_GAMMA_INTEGRATION_LIMITS_H_ 
