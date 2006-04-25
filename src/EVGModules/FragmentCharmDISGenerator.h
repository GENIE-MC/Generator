//____________________________________________________________________________
/*!

\class    genie::FragmentCharmDISGenerator

\brief    Generates the final state hadronic system in v DIS charm production
          interactions. The charm hadron is generated according to the input
          charm fractions, its longitudinal momentum is generated from the
          input fragmentation function, its transverse momentum is generated
          from an exponential function for the input pT^2 scale. The hadronic
          remnants are generated so as to conserve the hadronic shower charge
          and are distributed uniformly in the remaining phase space.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 22, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _FRAGMENT_CHARM_DIS_GENERATOR_H_
#define _FRAGMENT_CHARM_DIS_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class Interaction;
class FragmentationFunctionI;

class FragmentCharmDISGenerator : public HadronicSystemGenerator {

public :
  FragmentCharmDISGenerator();
  FragmentCharmDISGenerator(string config);
  ~FragmentCharmDISGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * evr) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   GenerateHadronicSystem   (GHepRecord * evr) const;
  bool   GenerateCharmHadronOnly  (GHepRecord * evr, bool ign) const;
  int    CharmedHadronPdgCode     (double E) const;
  double GeneratePT2              (double pT2max) const;
  void   LoadConfig               (void);

  const FragmentationFunctionI * fFragmFunc;
  double fpT2scale;
  bool   fCharmOnly;
};

}      // genie namespace

#endif // _FRAGMENT_CHARM_DIS_GENERATOR_H_
