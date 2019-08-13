//____________________________________________________________________________
/*!
\class    genie::MKHelicityAmplModelI

\brief    Pure abstract base class. Defines the MKHelicityAmplModelI interface.

\authors  Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research \n
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MK_HELICITY_AMPL_MODEL_I_H_
#define _MK_HELICITY_AMPL_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/Resonance/XSection/FKR_MK.h"
#include "Physics/Resonance/XSection/MKHelicityAmpl.h"

namespace genie {

class MKHelicityAmplModelI : public Algorithm
{
public:
  virtual ~MKHelicityAmplModelI();

  // define the MKHelicityAmplModelI interface
  virtual const MKHelicityAmpl & Compute(Resonance_t res, const FKR_MK & fkr) const = 0;

protected:
  MKHelicityAmplModelI();
  MKHelicityAmplModelI(string name);
  MKHelicityAmplModelI(string name, string config);
};

}        // namespace

#endif   // _MK_HELICITY_AMPL_MODEL_I_H_



