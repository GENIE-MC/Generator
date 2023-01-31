//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmplModelEMp

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free neutrons, as computed in MAID
	  parameterization.

\author   Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
          Tel Aviv Universtiy

\created  January 2023

\cpright  Copyright (c) 2023-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _MAIDHELICITY_AMPL_MODEL_EM_P_H_
#define _MAIDHELICITY_AMPL_MODEL_EM_P_H_

#include "Physics/Resonance/XSection/MAIDHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelEMp : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelEMp();
  MAIDHelicityAmplModelEMp(string config);
  virtual ~MAIDHelicityAmplModelEMp();

  const MAIDHelicityAmpl & Compute( const Interaction * interaction ) const;

private:

mutable MAIDHelicityAmpl fAmpl; 
};

}        // genie namespace
#endif
