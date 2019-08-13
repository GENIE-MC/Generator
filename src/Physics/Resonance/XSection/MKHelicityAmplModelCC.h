//____________________________________________________________________________
/*!

\class    genie::MKHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, modified for MK-model.

          Concrete implementation of the RSHelicityAmplModelI interface.

\authors  Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research \n
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab


\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_MK_CC_H_
#define _HELICITY_AMPL_MODEL_MK_CC_H_

#include "Physics/Resonance/XSection/MKHelicityAmplModelI.h"

namespace genie {

class MKHelicityAmplModelCC : public MKHelicityAmplModelI {

public:
  MKHelicityAmplModelCC();
  MKHelicityAmplModelCC(string config);
  virtual ~MKHelicityAmplModelCC();

  // MKHelicityAmplModelI interface implementation
 const MKHelicityAmpl & Compute(Resonance_t res, const FKR_MK & fkr) const;

private:
  mutable MKHelicityAmpl fAmpl; 
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_MK_CC_H_


