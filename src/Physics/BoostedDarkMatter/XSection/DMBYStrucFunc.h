//____________________________________________________________________________
/*!

\class    genie::BYStrucFunc

\brief    Bodek Yang structure function model

\ref      hep-ph/0411202

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 28, 2004

*/
//____________________________________________________________________________

#ifndef _DM_BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
#define _DM_BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_

#include "Physics/BoostedDarkMatter/XSection/QPMDMDISStrucFuncBase.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDFModelI.h"

namespace genie {

class DMBYStrucFunc : public QPMDMDISStrucFuncBase {

public:
  DMBYStrucFunc();
  DMBYStrucFunc(string config);
  virtual ~DMBYStrucFunc();

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

#endif    // _DM_BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
