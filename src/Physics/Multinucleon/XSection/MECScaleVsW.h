//____________________________________________________________________________
/*!

\class    genie::MECScaleVsW

\brief    This class is responsible to compute the MEC scaling factor given 
          Q0, Q3. The scaling is done as a function of the hadronic invariant
	  mass.

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MEC_SCALE_VS_W_H_
#define _MEC_SCALE_VS_W_H_

#include "Physics/Common/XSecScaleI.h"
#include <TSpline.h>

namespace genie {
  
  class MECScaleVsW: public XSecScaleI {

    using weight_type_map = std::map<double,double> ;
    using weight_type_pair = std::pair<double,double> ;

  public:
    MECScaleVsW();
    MECScaleVsW(string config);
    virtual ~MECScaleVsW();
    
    // This function returns the scaling value at a given Q0 Q3:
    virtual double GetScaling( const Interaction & ) const override ; 
    
  protected:
    
    // Load algorithm configuration
    virtual void LoadConfig (void) override ;

    // This function returns the scaling value at a given Q0 Q3:
    virtual double GetScaling( const double Q0, const double Q3 ) const ; 

    // This function adds the limits of the phase space if they are not set by the user
    weight_type_map GetMapWithLimits( const double Q0, const double Q3 ) const ;

    // Thist function calculates the scale factor value at W as a linear interpolation
    // between two W values (Wmin,Wmax) with weights (scale_min,scale_max).
    virtual double ScaleFunction( const double W, const weight_type_pair min, const weight_type_pair max ) const ;

 private: 
    double fDefaultWeight ; 
    weight_type_map fWeightsMap ;
    // Adding Spline to handle the limits of W1:
    TSpline3 fW1_Q0Q3_limits ; 

    double fLowLimitWeight ; 
    double fUpperLimitWeight ; 

  };
  
}       // genie namespace
#endif  // _MEC_SCALE_VS_W_H_
