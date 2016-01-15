//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayer. \n
          The PythiaDecayer is a concrete implementation of the DecayModelI
          interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 20, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

  // implement the DecayModelI interface  
  bool           IsHandled      (int pdgc)                      const;
  void           Initialize     (void)                          const;
  TClonesArray * Decay          (const DecayerInputs_t & inp)   const;
  double         Weight         (void)                          const;
  void           InhibitDecay   (int pdg, TDecayChannel * dc=0) const;
  void           UnInhibitDecay (int pdg, TDecayChannel * dc=0) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig             (void);
  double SumBR                  (int kc) const;
  int    FindPythiaDecayChannel (int kc, TDecayChannel* dc) const;
  bool   MatchDecayChannels     (int ichannel, TDecayChannel * dc) const;

  mutable TPythia6 * fPythia;  ///< PYTHIA6 wrapper class
  mutable double fWeight;
//bool fForceDecay;
};

}         // genie namespace

#endif    // _PYTHIA_DECAYER_H_
