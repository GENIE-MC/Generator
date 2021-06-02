//____________________________________________________________________________
/*!

\class    genie::MECScaleVsW

\brief    This class is responsible to compute the MEC scaling factor given 
          q0, q3. The scaling is done as a function of the hadronic invariant
	  mass.

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MEC_SCALE_VS_W_H_
#define _MEC_SCALE_VS_W_H_

#include "Physics/Common/XSecScaleMap.h"
#include <TSpline.h>

namespace genie {
  
  class MECScaleVsW: public XSecScaleMap {
    
  public:
    MECScaleVsW();
    MECScaleVsW(string config);
    virtual ~MECScaleVsW();
    
    // This function returns the scaling value at a given q0 q3:
    virtual double GetScaling( const Interaction & ) const ; 
    
    void Configure (const Registry & config);
    void Configure (string config);
    
  protected:
    
    // Load algorithm configuration
    void LoadConfig (void);

    // This function returns the scaling value at a given q0 q3:
    virtual double GetScaling( const double q0, const double q3 ) const ; 

    // This function adds the limits of the phase space if they are not set by the user
    void GetVectorWithLimits( std::vector<double> & W_limits, std::vector<double> & weights,
			      const double q0, const double q3 , const double weight ) const ;

    // Thist function calculates the scale factor value at W as a linear interpolation
    // between two W values (Wmin,Wmax) with weights (scale_min,scale_max).
    virtual double ScaleFunction( const double W, const double Win, const double Wmax, 
				  const double scale_min, const double scale_max ) const ; 

 private: 
    double fDefaultWeight ; 
    std::vector<double> fWeights ; 
    std::vector<double> fWValues ; 
    TSpline3 * fW1_q0q3_limits ; 

  };
  
}       // genie namespace
#endif  // _MEC_SCALE_VS_W_H_
