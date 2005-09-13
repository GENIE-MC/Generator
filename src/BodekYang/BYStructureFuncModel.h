//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModel

\brief    Abstract class. Provides common implementation for concrete
          DISStructureFuncModelI objects computing the Bodek Yang structure
          functions.

\ref      hep-ph/0411202

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
#define _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_

#include "PartonModel/DISStructureFuncModel.h"
#include "Interaction/Interaction.h"
#include "PDF/PDFModelI.h"

namespace genie {

class BYStructureFuncModel : public DISStructureFuncModel {

public:

  virtual ~BYStructureFuncModel();

  // Overload Algorithm::Configure() to read the config. registry
  // at the algorithm initialization and set private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  // protected constructors - abstract class
  BYStructureFuncModel();
  BYStructureFuncModel(const char * param_set);

  // override part of the DISStructureFuncModel implementation
  // to compute all the corrections applied by the Bodek-Yang model.
  double ScalingVar (const Interaction * interaction) const;
  double KSea       (const Interaction * interaction) const;
  double KVal       (const Interaction * interaction) const;

  // Bodek-Yang model-specific parameters A,B,Csea,Cv1,Cv2
  void   Init(void);
  void   ReadBYParams(void);
  double fA;
  double fB;
  double fCs;
  double fCv1;
  double fCv2;
};

}         // genie namespace

#endif    // _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_H_
