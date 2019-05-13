//____________________________________________________________________________
/*!

\class    genie::TransverseEnhancementFFModel

\brief    Modification of magnetic form factors to match observed enhancement
          in transverse cross section of the quasi-elastic peak.
          Implements ElFormFactorsModelI.  Requires another subclass of
          ElFormFactorsModelI to calculate original form factors, which
          are then enhances.

\ref      http://arxiv.org/pdf/1106.0340
          http://arxiv.org/abs/1405.0583

\author   Brian Coopersmith, University of Rochester

\created  10/22/2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _TRANSVERSE_ENHANCEMENT_FF_MODEL_H_
#define _TRANSVERSE_ENHANCEMENT_FF_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {
class Target;

class TransverseEnhancementFFModel : public ELFormFactorsModelI {

public:
  TransverseEnhancementFFModel();
  TransverseEnhancementFFModel(string config);
  virtual ~TransverseEnhancementFFModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;
  
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);
  
  void SetElFFBaseModel(const ELFormFactorsModelI* ffBase) const {
    fElFormFactorsBase = ffBase;
  }

private:

  void GetTransEnhParams(const Target& target, double* transEnhA,
                         double* transEnhB) const;
  void LoadConfig(void);
  double GetTransEnhMagFF(
      double magFF, const Interaction * interaction) const;
 
  mutable ELFormFactorsModelI const* fElFormFactorsBase;
  map<int, double> fNucMagFF_RT_A;
  map<int, double> fNucMagFF_RT_B;
  
  map<pair<int, int>, double> fRangeMagFF_RT_A;
  map<pair<int, int>, double> fRangeMagFF_RT_B;
};

}         // genie namespace
#endif    // _TRANSVERSE_ENHANCEMENT_FF_MODEL_H_
