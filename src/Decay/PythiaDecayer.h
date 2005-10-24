//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayers.

          The PythiaDecayer is a concrete implementation of the DecayModelI
          interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 20, 2004

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

  //-- implement the DecayModelI interface
  
  bool           IsHandled  (int pdgc)                    const;
  void           Initialize (void)                        const;
  TClonesArray * Decay      (const DecayerInputs_t & inp) const;
  
private:

  void SwitchOnAllChannels        (int pdgc) const;
  void SwitchOffInhibitedChannels (int pdgc, const TClonesArray * inhibited) const;
  bool MatchDecayChannel          (int ichannel, TDecayChannel & dc) const;

  TPythia6 * fPythia;
};

}         // genie namespace

#endif    // _PYTHIA_DECAYER_H_
