//____________________________________________________________________________
/*!

\class    genie::DecayModelI

\brief    Pure abstract base class. Defines the DecayModelI interface to be
          implemented by any algorithmic class decaying a particle.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 20, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________


#ifndef _DECAY_MODEL_I_H_
#define _DECAY_MODEL_I_H_

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <TLorentzVector.h>

#include "Algorithm/Algorithm.h"
#include "Decay/DecayerInputs.h"

namespace genie {

class DecayModelI : public Algorithm {

public:

  virtual ~DecayModelI();

  //! define DecayModelI interface

  virtual bool           IsHandled  (int pdgc)                    const = 0; ///< can this particle be decayed?
  virtual void           Initialize (void)                        const = 0; ///< decayer initialization
  virtual TClonesArray * Decay      (const DecayerInputs_t & inp) const = 0; ///< return a TClonesArray of TMCParticle objects (NOTE: all TMCParticle units in GeV^n [hbar=c=1])
  virtual double         Weight     (void)                        const = 0; ///< last decay weight

protected:

  DecayModelI();
  DecayModelI(string name);
  DecayModelI(string name, string config);
};

}         // genie namespace

#endif    // _DECAY_MODEL_I_H_
