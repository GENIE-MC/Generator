//____________________________________________________________________________
/*!

\class    genie::ScaleMECW

\brief    This class is responsible to compute the MEC scaling factor given 
          q0, q3. The scaling is done as a function of the hadronic invariant
	  mass.

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SCALE_MEC_W_H_
#define _SCALE_MEC_W_H_

#include "Framework/Algorithm/Algorithm.h"
#include <TSpline.h>

using std::map; 

namespace genie {
  
  class ScaleMECW: public Algorithm {
    
  public:
    ScaleMECW();
    ScaleMECW(string config);
    virtual ~ScaleMECW();
    
    // This function returns the scaling value at a given q0 q3:
    virtual double GetScaling( const double q0, const double q3 ) const ; 
    
    void Configure (const Registry & config);
    void Configure (string config);
    
  protected:
    
    // Load algorithm configuration
    void LoadConfig (void);
    
    // Thist function calculates the scale factor value at W as a linear interpolation
    // between two W values (Wmin,Wmax) with weights (scale_min,scale_max).
    virtual double ScaleFunction( double W, double Win, double Wmax, 
				  double scale_min, double scale_max ) const ; 

 private: 
    std::vector<double> fWeights ; 
    std::vector<double> fLimits ; 
    TSpline3 * fW1_q0q3_limits ; 

  };
  
}       // genie namespace
#endif  // _SCALE_MEC_W_H_
