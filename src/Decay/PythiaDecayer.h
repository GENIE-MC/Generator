//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayer. \n
          The PythiaDecayer is a concrete implementation of the DecayModelI
          interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 20, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PYTHIA_DECAYER_I_H_
#define _PYTHIA_DECAYER_I_H_

#include <TPythia6.h>

#include "Decay/DecayModelI.h"

namespace genie {

class PythiaDecayer : public DecayModelI {

public:

  PythiaDecayer();
  PythiaDecayer(string config);
  virtual ~PythiaDecayer();

  //! implement the DecayModelI interface
  
  bool           IsHandled  (int pdgc)                    const;
  void           Initialize (void)                        const;
  TClonesArray * Decay      (const DecayerInputs_t & inp) const;
  double         Weight     (void)                        const;

  //! overload the Algorithm::Configure() methods to load private data
  //!  members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig                 (void);
  void SwitchOnAllChannels        (int pdgc) const;
  void SwitchOffInhibitedChannels (int pdgc, const TClonesArray * inhibited) const;
  bool MatchDecayChannel          (int ichannel, TDecayChannel & dc) const;

  TPythia6 * fPythia;
  bool       fForceDecay;
};

}         // genie namespace

#endif    // _PYTHIA_DECAYER_H_
