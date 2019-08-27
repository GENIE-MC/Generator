//____________________________________________________________________________
/*!

  \class    genie::AlvarezRusoSaulSalaCOHGammaPXsec

  \brief    Implementation of the Alvarez-Ruso & Saul-Sala coherent gamma production model

  Is a concrete implementation of the XSecAlgorithmI interface.

  \ref      

  \author   Marco Roda
  University of Liverpool

  \created  August 23, 2019

  \cpright  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALVAREZ_RUSO_SAUL_SALA_COH_GAMMA_XSEC_H_
#define _ALVAREZ_RUSO_SAUL_SALA_COH_GAMMA_XSEC_H_ 

#include "Framework/EventGen/XSecAlgorithmI.h"


namespace genie {

  class XSecIntegratorI; 
  class Interaction;

  class AlvarezRusoSaulSalaCOHGammaPXSec : public XSecAlgorithmI {

  public:
    AlvarezRusoSaulSalaCOHGammaPXSec() ;
    AlvarezRusoSaulSalaCOHGammaPXSec(string config);
    virtual ~AlvarezRusoSaulSalaCOHGammaPXSec();

    //-- XSecAlgorithmI interface implementation
    double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
    double Integral        (const Interaction * i) const;
    bool   ValidProcess    (const Interaction * i) const;

    //-- overload the Algorithm::Configure() methods to load private data
    //   members from configuration options
    void Configure(const Registry & config);
    void Configure(string config);

  private:
  
    void LoadConfig(void);

    //-- private data members loaded from config Registry or set to defaults

    const XSecIntegratorI * fXSecIntegrator;
 

  };

}      // genie namespace

#endif  // _REIN_SEGHAL_COHPI_PXSEC_H_
