//____________________________________________________________________________
/*!

\class    genie::DecayModelI

\brief    Pure abstract base class. Defines the DecayModelI interface to be
          implemented by any algorithmic class decaying a particle.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 20, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________


#ifndef _DECAY_MODEL_I_H_
#define _DECAY_MODEL_I_H_

class TDecayChannel;

#include "Framework/Algorithm/Algorithm.h"
#include "Physics/Decay/DecayerInputs.h"

namespace genie {

class DecayModelI : public Algorithm {

public:

  virtual ~DecayModelI();

  //
  // define the DecayModelI interface
  //

  //! can this particle be decayed?
  virtual bool IsHandled  (int pdgc) const = 0; 

  //! decayer initialization
  virtual void Initialize (void) const = 0; 

  //! return a TClonesArray of TMCParticle objects (NOTE: all TMCParticle units in GeV^n [hbar=c=1])
  virtual TClonesArray * Decay (const DecayerInputs_t & inp) const = 0; 

  //! last decay weight
  virtual double Weight(void) const = 0; 

  //! inhibit input decay channel for the input particle (inhibit all decays if dc is null)
  virtual void InhibitDecay(int pdgc, TDecayChannel * dc = 0) const = 0;

  //! clear inhibit flags & re-enable all decay channel (enable all if dc is null)
  virtual void UnInhibitDecay(int pdgc, TDecayChannel * dc = 0) const = 0;

protected:

  DecayModelI();
  DecayModelI(string name);
  DecayModelI(string name, string config);
};

}         // genie namespace

#endif    // _DECAY_MODEL_I_H_
