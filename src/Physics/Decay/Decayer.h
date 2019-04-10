//____________________________________________________________________________
/*!

\class    genie::Decayer

\brief    Base class for decayer classes.
          Implements common configuration, allowing users to toggle on/off
          flags for particles and decay channels.
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 14, 2018

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DECAYER_H_
#define _DECAYER_H_

class TDecayChannel;

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/GHEP/GHepStatus.h"

namespace genie {

class GHepParticle;

class Decayer : public EventRecordVisitorI {

public:
  virtual ~Decayer();

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

protected:
  Decayer();
  Decayer(string name);
  Decayer(string name, string config);

  virtual void LoadConfig    (void);
  virtual bool ToBeDecayed   (int pdgc, GHepStatus_t ist) const;
  virtual bool IsUnstable    (int pdgc) const;
  virtual bool IsHandled     (int pdgc) const = 0;
  virtual void InhibitDecay  (int pdgc, TDecayChannel * dc=0) const = 0;
  virtual void UnInhibitDecay(int pdgc, TDecayChannel * dc=0) const = 0;

  bool        fGenerateWeighted;    ///< generate weighted or unweighted decays?
  bool        fRunBefHadroTransp;   ///< is invoked before or after FSI?
  PDGCodeList fParticlesToDecay;    ///< list of particles to be decayed
  PDGCodeList fParticlesNotToDecay; ///< list of particles for which decay is inhibited
};

}      // genie namespace
#endif // _DECAYER_H_
