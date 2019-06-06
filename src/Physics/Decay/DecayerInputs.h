//____________________________________________________________________________
/*!

\class    genie::DecayerInputs

\brief    A primitive class with public data members to keep a short argument
          list in Particle Decayer algorithms

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 02, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DECAYER_INPUTS_H_
#define _DECAYER_INPUTS_H_

#include <TLorentzVector.h>
#include <TVector3.h>

namespace genie {

class DecayerInputs_t {

public:

  DecayerInputs_t()  { this->Init(); }
  ~DecayerInputs_t() { }

  int                    PdgCode; ///< pdg code
  const TLorentzVector * P4;      ///< 4-momentum
  const TVector3 *       Polz;    ///< polarization

private:

  void Init(void) {
    PdgCode = 0;
    P4      = 0;
    Polz    = 0; 
  }

};

}         // genie namespace

#endif    // _DECAYER_INPUTS_H_
