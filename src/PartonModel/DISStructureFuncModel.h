//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Implements the DISFormFactorsModelI interface
          but can not be instantiated. Its mere role of existence is to factor
          out common implementation from concrete implementations like
          the PartonModelNC, PartonModelCCAboveCharmThr and
          PartonModelCCBelowCharmThr.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_H_
#define _DIS_STRUCTURE_FUNCTIONS_H_

#include "Base/DISStructureFuncModelI.h"
#include "Interaction/Interaction.h"
#include "PDF/PDF.h"

namespace genie {

class DISStructureFuncModel : public DISStructureFuncModelI {

public:

  virtual ~DISStructureFuncModel();

  //-- common code for all DISFormFactorsModelI interface implementations
  virtual double xF1 (const Interaction * interaction) const
  {
     return this->F2(interaction)/2;
  }
  virtual double F4  (const Interaction * interaction) const
  {
     return 0.;
  }
  virtual double xF5 (const Interaction * interaction) const
  {
     return this->F2(interaction);
  }
  virtual double F6 (const Interaction * interaction) const
  {
     return 0.;
  }
  virtual bool CalculatePDFs (
                      const Interaction * interaction) const;


  //-- Overload Algorithm's Configure() to set the PDF data member
  //   from the configuration registry

  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

protected:

  DISStructureFuncModel();
  DISStructureFuncModel(const char * param_set);

  //-- commom code for PDF set calculation for all
  //   DISFormFactorsModelI interface implementations
  virtual void   ConfigPDF (void);
  //virtual void   CalcPDFs (const Interaction * interaction) const;
  //virtual void   CalcPDFs (double x, double Q2)             const;
  virtual double CalcQ2   (const Interaction * interaction) const;

  PDF * fPDF; // PDF (values + PDG model)
};

}         // genie namespace

#endif    // _DIS_STRUCTURE_FUNCTIONS_H_

