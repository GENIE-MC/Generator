//____________________________________________________________________________
/*!

  \class    genie::MAIDHelicityAmplModelEMn

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

#ifndef _MAIDHELICITY_AMPL_MODEL_EM_N_H_
#define _MAIDHELICITY_AMPL_MODEL_EM_N_H_

#include "Physics/Resonance/XSection/MAIDHelicityAmplModelI.h"

namespace genie {

  class MAIDHelicityAmplModelEMn : public MAIDHelicityAmplModelI {

  public:
    MAIDHelicityAmplModelEMn();
    MAIDHelicityAmplModelEMn(string config);
    virtual ~MAIDHelicityAmplModelEMn();

    const MAIDHelicityAmpl & Compute( const Interaction * interaction ) const;

  private:
    mutable MAIDHelicityAmpl fAmpl;

  };

}        // genie namespace
#endif
