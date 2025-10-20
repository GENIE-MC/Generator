//____________________________________________________________________________
/*!

\class    genie::XSecScaleMap

\brief    This class is responsible to compute a scaling factor for the XSec

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _XSEC_SCALE_MAP_H_
#define _XSEC_SCALE_MAP_H_

#include "Physics/Common/XSecScaleI.h"
#include <map>

namespace genie {
  
  class XSecScaleMap: public XSecScaleI {
    
  public:
    XSecScaleMap();
    XSecScaleMap(string config);
    virtual ~XSecScaleMap();
    
    // This function returns the scaling value for a given interaction:
    virtual double GetScaling( const Interaction & ) const override ; 
    
  protected:
    
    // Load algorithm configuration
    virtual void LoadConfig (void) override ;
    
 private: 
    const XSecScaleI * fXSecScaleDefault ; 
    std::map<int,const XSecScaleI *> fXSecScaleMap ;

  };
  
}       // genie namespace
#endif  // _XSEC_SCALE_MAP_H_
