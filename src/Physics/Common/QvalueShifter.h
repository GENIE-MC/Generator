//____________________________________________________________________________
/*!

\class    genie::QvalueShfiter

\brief    This class is responsible to compute a relative shift to a Qvalue

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _QVALUE_SHIFTER_H_
#define _QVALUE_SHIFTER_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/Target.h"
#include <map>

using std::map; 

namespace genie {
  
  class QvalueShifter: public Algorithm {
    
  public:
    QvalueShifter();
    QvalueShifter(string config);
    virtual ~QvalueShifter();
    
    virtual double Shift( const Target & target ) const ; 
    virtual double Shift( const Interaction & interaction ) const ; 
    
    void Configure (const Registry & config);
    void Configure (string config);
    
  protected:
    
    // Load algorithm configuration
    void LoadConfig (void);
    
 private: 
    double fRelShiftDefault ; 
    std::map<int,double> fRelShift ; // map from pdf_target to the Qvalue relative shift 
    
  };
  
}       // genie namespace
#endif  // _QVALUE_SHIFTER_H_
