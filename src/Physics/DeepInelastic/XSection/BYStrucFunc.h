//____________________________________________________________________________
/*!

\class    genie::BYStrucFunc

\brief    Bodek Yang structure function model

\ref      hep-ph/0411202

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  September 28, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
#define _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_

#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDFModelI.h"

namespace genie {

class BYStrucFunc : public QPMDISStrucFuncBase {

public:
  BYStrucFunc();
  BYStrucFunc(string config);
  virtual ~BYStrucFunc();

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

#endif    // _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
