//____________________________________________________________________________
/*!

  \class    genie::MAIDRESVectFormFactorsEMn

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

#ifndef _MAID_RESVECTFF_EM_N_H_
#define _MAID_RESVECTFF_EM_N_H_

#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"

namespace genie {

  class MAIDRESVectFormFactorsEMn : public RESVectFormFactorsI {

  public:
    MAIDRESVectFormFactorsEMn();
    MAIDRESVectFormFactorsEMn(string config);
    virtual ~MAIDRESVectFormFactorsEMn();

    RESVectFFAmplitude Compute( const Interaction interaction ) const;

  };

}        // genie namespace
#endif
