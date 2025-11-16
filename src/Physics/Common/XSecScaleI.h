//____________________________________________________________________________
/*!

\class    genie::XSecScaleI

\brief    This class is responsible to compute a scaling factor for the XSec

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _XSEC_SCALE_I_H_
#define _XSEC_SCALE_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {
  
  class XSecScaleI: public Algorithm {
    
  public:    
    virtual ~XSecScaleI();
    // This function returns the scaling value for a given interaction:
    virtual double GetScaling( const Interaction & ) const = 0 ; 

  protected:
    XSecScaleI( string name, string config = "Default" );

    void Configure(const Registry & config) override ;
    virtual void Configure (string config) override ;
    
    virtual void LoadConfig(void) = 0 ;
    
  };
  
}       // genie namespace
#endif  // _XSEC_SCALE_I_H_
