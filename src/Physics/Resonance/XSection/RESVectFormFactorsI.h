//____________________________________________________________________________
/*!
\class    genie::RESVectFormFactorsI

\brief    Pure abstract base class. Defines the RESVectFormFactorsI interface.

\author   Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
          Tel Aviv University

\created  January 2023

\cpright  Copyright (c) 2023-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _RES_VECT_FF_AMP_I_H_
#define _RES_VECT_FF_AMP_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Resonance/XSection/RESVectFFAmplitude.h"

namespace genie {

class RESVectFormFactorsI : public Algorithm
{
public:
  virtual ~RESVectFormFactorsI();

  // define the RESVectFormFactorsI interface
  virtual RESVectFFAmplitude Compute( const Interaction interaction ) = 0;

protected:
  RESVectFormFactorsI();
  RESVectFormFactorsI(string name);
  RESVectFormFactorsI(string name, string config);
};

}        // namespace

#endif   // _RES_VECT_FF_AMP_I_H_
