//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Provides common implementation for concrete
          DISStructureFuncModelI objects

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_MODEL_H_
#define _DIS_STRUCTURE_FUNCTIONS_MODEL_H_

#include "Base/DISStructureFuncModelI.h"
#include "Interaction/Interaction.h"
#include "PDF/PDF.h"

namespace genie {

class DISStructureFuncModel : public DISStructureFuncModelI {

public:

  virtual ~DISStructureFuncModel();

  //-- common code for all DISFormFactorsModelI interface implementations
  virtual double F1 (void) const { return fF1; }
  virtual double F2 (void) const { return fF2; }
  virtual double F3 (void) const { return fF3; }
  virtual double F4 (void) const { return fF4; }
  virtual double F5 (void) const { return fF5; }
  virtual double F6 (void) const { return fF6; }

  virtual void Calculate (const Interaction * interaction) const;

  //-- Overload Algorithm's Configure() to set the PDF data member
  //   from the configuration registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

protected:

  DISStructureFuncModel();
  DISStructureFuncModel(const char * param_set);

  //-- commom code for SF calculation for all DISFormFactorsModelI
  //   interface implementations inheriting DISFormFactorsModel
  virtual void   ConfigPDF  (void);
  virtual double Q2         (const Interaction * interaction) const;
  virtual double ScalingVar (const Interaction * interaction) const;
  virtual double KSea       (const Interaction * interaction) const;
  virtual double KVal       (const Interaction * interaction) const;
  virtual void   CalcPDFs   (const Interaction * interaction) const;
  virtual double Q          (const Interaction * interaction) const;
  virtual double QBar       (const Interaction * interaction) const;
  virtual double NuclMod    (const Interaction * interaction) const;

  PDF * fPDF;
  PDF * fPDFc;  // to take into account charm suppression

  mutable double fF1;
  mutable double fF2;
  mutable double fF3;
  mutable double fF4;
  mutable double fF5;
  mutable double fF6;
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNCTIONS_MODEL_H_

