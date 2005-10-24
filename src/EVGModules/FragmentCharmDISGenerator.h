//____________________________________________________________________________
/*!

\class   genie::FragmentCharmDISGenerator

\brief   Generates the final state hadronic system in v DIS charm production
         interactions. The charm hadron is generated according to the input
         charm fractions, its longitudinal momentum is generated from the
         input fragmentation function, its transverse momentum is generated
         from an exponential function for the input pT^2 scale. The hadronic
         remnants are generated so as to conserve the hadronic shower charge
         and are distributed uniformly in the remaining phase space.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 22, 2005

*/
//____________________________________________________________________________

#ifndef _FRAGMENT_CHARM_DIS_GENERATOR_H_
#define _FRAGMENT_CHARM_DIS_GENERATOR_H_

#include <TLorentzVector.h>
#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class Interaction;
class FragmentCharmDISGenerator : public HadronicSystemGenerator {

public :

  FragmentCharmDISGenerator();
  FragmentCharmDISGenerator(string config);
  ~FragmentCharmDISGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * evr) const;

private:

  TLorentzVector HadronicSystemP4 (GHepRecord * evr) const;
  void   GenerateHadronicSystem   (GHepRecord * evr) const;
  bool   GenerateCharmHadronOnly  (GHepRecord * evr, bool ign) const;
  int    CharmedHadronPdgCode     (double E) const;
  int    HadronShowerCharge       (const Interaction * i) const;
  double GeneratePT2              (double pT2max) const;
};

}      // genie namespace

#endif // _FRAGMENT_CHARM_DIS_GENERATOR_H_
