//____________________________________________________________________________
/*!

\class    genie::DecayerInputs

\brief    A primitive class with public data members to keep a short argument
          list in Particle Decayer algorithms

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  December 02, 2004

*/
//____________________________________________________________________________

#ifndef _DECAYER_INPUTS_H_
#define _DECAYER_INPUTS_H_

#include <TLorentzVector.h>
#include <TClonesArray.h>

namespace genie {

const UInt_t kDChIsInhibited = 1<<17;

class DecayerInputs_t {

public:

  DecayerInputs_t()  { this->Init(); }
  ~DecayerInputs_t() { }

  int                    PdgCode;
  const TLorentzVector * P4;
  const TClonesArray *   InhibitedChannels;

private:

  void Init(void) {
    PdgCode           = 0;
    P4                = 0;
    InhibitedChannels = 0; 
  }

};

}         // genie namespace

#endif    // _DECAYER_INPUTS_H_
