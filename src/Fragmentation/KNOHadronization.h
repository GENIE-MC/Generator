//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    The KNO hadronization model.
          This hadronization scheme is similar to the one originally used
          in NeuGEN by G.Barr, G.F.Pearce, H.Gallagher et al. \n
          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_H_
#define _KNO_HADRONIZATION_H_

#include <vector>

#include <TGenPhaseSpace.h>

#include "Fragmentation/HadronizationModelI.h"

using std::vector;

namespace genie {

class MultiplicityProbModelI;
class DecayModelI;

class KNOHadronization : public HadronizationModelI {

public:

  KNOHadronization();
  KNOHadronization(string config);
  virtual ~KNOHadronization();

  //-- implement the HadronizationModelI interface
  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void          LoadConfig            (void);
  vector<int> * GenerateFSHadronCodes (int mult, int maxQ, double W) const;
  int           GenerateBaryonPdgCode (int mult, int maxQ)           const;
  int           HadronShowerCharge    (const Interaction * proc)     const;
  void          HandleDecays          (TClonesArray * particle_list) const;

  const MultiplicityProbModelI * fMultProbModel;
  const DecayModelI *            fDecayer;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;

  double       fPpi0;         ///< pi0 pi0 production probability
  double       fPpic;         ///< pi+ pi- production probability
  double       fPKc;          ///< K+  K- production probability
  double       fPK0;          ///< K0  K0bar production probability
  bool         fForceDecays;  ///< force decays of unstable hadrons produced?
  bool         fForceMinMult; ///< force minimum multiplicity if (at low W) generated less
  unsigned int fMaxMult;      ///< maximum allowed multiplicity
};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_H_

