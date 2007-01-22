//____________________________________________________________________________
/*!

\class    genie::FGMBodekRitchie

\brief    The Bodek Richie Fermi Gass model. Implements the NuclearModelI 
          interface.

\ref      

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _FGM_BODEK_RITCHIE_H_
#define _FGM_BODEK_RITCHIE_H_

#include <map>

#include <TH1D.h>
#include "Nuclear/NuclearModelI.h"

using std::map;

namespace genie {

class FGMBodekRitchie : public NuclearModelI {

public:
  FGMBodekRitchie();
  FGMBodekRitchie(string config);
  virtual ~FGMBodekRitchie();

  //-- implement the NuclearModelI interface
  bool     GenerateNucleon (const Target & t) const;
  double   Prob            (double p, double w, const Target & t) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void   LoadConfig (void);
  TH1D * ProbDistro (const Target & t) const;

  mutable map<string, TH1D *> fProbDistroMap;

  int    fNPBins;
  double fPMax;
  double fPCutOff;
  string fKFTable;
};

}         // genie namespace
#endif    // _FGM_BODEK_RITCHIE_H_

