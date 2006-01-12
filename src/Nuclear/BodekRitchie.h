//____________________________________________________________________________
/*!

\class    genie::BodekRitchie

\brief    The Bodek-Ritchie model for the probability distribution of nucleon
          momenta within a nucleus.
          Implements the NuclMomentumModelI interface.

\ref      Bodek and Ritchie, Phys.Rev.D23:1070 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_RITCHIE_H_
#define _BODEK_RITCHIE_H_

#include "Nuclear/NuclMomentumModelI.h"

namespace genie {

class BodekRitchie : public NuclMomentumModelI {

public:

  BodekRitchie();
  BodekRitchie(string config);
  virtual ~BodekRitchie();

  //-- implement the NuclMomentumModelI interface
  TH1D * ProbabilityDistribution(const Target & target) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfigData();

  int    fNPBins;
  double fPMax;
  double fPCutOff;
};

}         // genie namespace

#endif    // _BODEK_RITCHIE_H_

