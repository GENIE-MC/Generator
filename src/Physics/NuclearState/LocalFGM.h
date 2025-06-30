//____________________________________________________________________________
/*!

\class    genie::LocalFGM

\brief    local Fermi gas model. Implements the NuclearModelI
          interface.

\ref

\author   Joe Johnston, Steven Dytman

\created  December 2015

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _LOCAL_FGM_H_
#define _LOCAL_FGM_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class LocalFGM : public NuclearModelI {

public:
  LocalFGM();
  LocalFGM(string config);
  virtual ~LocalFGM();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- allow methods to be called with a radius
  bool   GenerateNucleon (const Target & t, double hitNucleonRadius) const;
  double Prob            (double p, double w, const Target & t,
			  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmLocalFermiGas;
  }

  virtual double LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set) ;

 protected:
  void   LoadConfig (void);


private:
  TH1D * ProbDistro (const Target & t, double r) const;

  /// Throw a value from the Maxwell-Boltzmann distribution with the configured
  /// parameters
  double MaxwellBoltzmannRemovalE(const Target & t, double Ermv_min,
    double Ermv_max) const;

  map<int, double> fNucRmvE;

  double fPMax;
  bool fMomDepErmv;
  bool fForcePositiveErmv;
  bool fUseMBDist;

  // options related to SRC pairs
  double fSRC_Fraction;
  double fPCutOff;

  /// Center of Maxwell-Boltmann distribution used for SRC removal energy
  /// distribution, GeV
  double fSRC_Ermv_C;

  /// Sigma of Maxwell-Boltmann distribution used for SRC removal energy
  /// distribution, GeV
  double fSRC_Ermv_sigma;
};

}         // genie namespace
#endif    // _LOCAL_FGM_H_

