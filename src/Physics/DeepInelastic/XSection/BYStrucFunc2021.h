//____________________________________________________________________________
/*!

\class    genie::BYStrucFunc2021

\brief    2021 update of the Bodek and Yang structure function model

\ref      arXiv:2108.09240v2 [hep-ph]

\author   Júlia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
          Tel Aviv University
          
          Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 20, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_2021_H_
#define _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_2021_H_

#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDFModelI.h"

namespace genie {

class BYStrucFunc2021 : public QPMDISStrucFuncBase {

public:
  BYStrucFunc2021();
  BYStrucFunc2021(string config);
  virtual ~BYStrucFunc2021();

  // overload Algorithm::Configure() to read the config. registry
  // at the algorithm initialization and set private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  void Init         (void);
  void ReadBYParams (void);

  // override part of the DISStructureFuncModel implementation
  // to compute all the corrections applied by the Bodek-Yang model.
  double ScalingVar (const Interaction * i) const;
  void   KFactors   (const Interaction * i, double & kuv,
                         double & kdv, double & kus, double & kds) const;

  // Bodek-Yang model-specific parameters

  double fA;     ///< better scaling var parameter A
  double fB;     ///< better scaling var parameter B
  double fCsU;   ///< U-sea K factor parameter
  double fCsD;   ///< D-sea K factor parameter
  double fCv1U;  ///< U-val K factor parameter
  double fCv2U;  ///< U-val K factor parameter
  double fCv1D;  ///< D-val K factor parameter
  double fCv2D;  ///< D-val K factor parameter
};

}         // genie namespace

#endif    // _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_2021_H_
