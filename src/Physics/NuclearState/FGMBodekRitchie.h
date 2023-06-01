//____________________________________________________________________________
/*!

\class    genie::FGMBodekRitchie

\brief    The Bodek Richie Fermi Gass model. Implements the NuclearModelI 
          interface.

\ref      

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  October 09, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _FGM_BODEK_RITCHIE_H_
#define _FGM_BODEK_RITCHIE_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class FGMBodekRitchie : public NuclearModelI {

public:
  FGMBodekRitchie();
  FGMBodekRitchie(string config);
  virtual ~FGMBodekRitchie();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double mom, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmFermiGas; 
  }

  virtual double         FermiMomentum( const Target & t, int nucleon_pdg ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

 protected:
  void   LoadConfig (void);
  
private:

  TH1D * ProbDistro (const Target & t) const;

  mutable map<string, TH1D *> fProbDistroMap;

  map<int, double> fNucRmvE;

  double fPMax;
  double fPCutOff;
  bool fUseParametrization;

};

}         // genie namespace
#endif    // _FGM_BODEK_RITCHIE_H_

