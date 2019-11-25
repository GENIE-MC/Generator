//____________________________________________________________________________
/*!

\class    genie::COHDeltaCurrent

\brief    Class for COH Hadronic Current specification for Delta Resonance

\author   Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  November 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_DELTA_CURRENT_H_
#define _COH_DELTA_CURRENT_H_

#include "Physics/Coherent/XSection/COHHadronicCurrentI.h"

namespace genie {

class COHDeltaCurrent : public COHHadronicCurrentI {

public:

  COHDeltaCurrent() ;
  COHDeltaCurrent( string config );

  virtual ~COHDeltaCurrent() { ; } 

  virtual std::complex<double> R( uint8_t l, uint8_t m, 
				  const Interaction * i, 
				  const COHFormFactorI * ff ) const ;

  virtual Resonance_t Resonance() const { return kP33_1232 ; } 

  // methods to implemented from Algorithm 
  void Configure (const Registry & config);
  void Configure (string param_set);


 protected: 

  void LoadConfig(void);

};

}       // genie namespace
#endif  // #ifndef _COH_DELTA_CURRENT_H_

