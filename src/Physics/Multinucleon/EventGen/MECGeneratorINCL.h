//____________________________________________________________________________
/*!

\class    genie::MECGeneratorINCL

\brief    INCL-aware MEC primary-interaction generator.

          Derives from MECGenerator. Overrides only the two methods that
          touch the nuclear model (GenerateFermiMomentum and
          GenerateNSVInitialHadrons) plus LoadConfig. Everything else --
          ProcessEventRecord dispatch, AddTargetRemnant, all
          Select*LeptonKinematics, AddFinalStateLepton, RecoilNucleonCluster,
          DecayNucleonCluster, NucleonClusterConstituents -- is inherited
          from MECGenerator and stays automatically in sync with future
          changes to the base class.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _MEC_GENERATOR_INCL_H_
#define _MEC_GENERATOR_INCL_H_

#include "Physics/Multinucleon/EventGen/MECGenerator.h"

namespace genie {

class NucleusGenI;

class MECGeneratorINCL : public MECGenerator {

public :
  MECGeneratorINCL();
  MECGeneratorINCL(string config);
  virtual ~MECGeneratorINCL();

  // ProcessEventRecord, Configure(const Registry &), Configure(string)
  // are all inherited unchanged from MECGenerator.

protected:

  // Overrides: the two methods that depend on the nuclear model, plus
  // LoadConfig to additionally wire fNucleusGen.
  void LoadConfig                    (void) override;
  void GenerateFermiMomentum         (GHepRecord * event) const override;
  void GenerateNSVInitialHadrons     (GHepRecord * event) const override;

private:

  const NucleusGenI * fNucleusGen;   ///< INCL-aware nucleus generator
                                     ///  (position + momentum + binding 4-momentum
                                     ///   for the di-nucleon cluster)
};

}      // genie namespace

#endif // _MEC_GENERATOR_INCL_H_
#endif // __GENIE_INCL_ENABLED__
