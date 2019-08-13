//____________________________________________________________________________
/*!

\class    genie::MKHelicityAmplModelNCp

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free protons, modified for MK-model.

          Concrete implementation of the MKHelicityAmplModelI interface.

\authors  Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research \n
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab


\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_MK_NC_P_H_
#define _HELICITY_AMPL_MODEL_MK_NC_P_H_

#include "Physics/Resonance/XSection/MKHelicityAmplModelI.h"

namespace genie {

class MKHelicityAmplModelNCp : public MKHelicityAmplModelI {

public:
  MKHelicityAmplModelNCp();
  MKHelicityAmplModelNCp(string config);
  virtual ~MKHelicityAmplModelNCp();

  // RSHelicityAmplModelI interface implementation
  const MKHelicityAmpl & Compute(Resonance_t res, const FKR_MK & fkr) const;
 
  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  mutable MKHelicityAmpl fAmpl;
  
  double fSin28w; 
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_MK_NC_P_H_


