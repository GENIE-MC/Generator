//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Provides common implementation for concrete
          DISStructureFuncModelI objects

\ref      For a discussion of DIS SF see for eaxample E.A.Paschos and J.Y.Yu, 
          Phys.Rev.D 65.033002 and R.Devenish and A.Cooper-Sarkar, OUP 2004.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
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
  DISStructureFuncModel(string name);
  DISStructureFuncModel(string name, string config);

  //-- commom code for SF calculation for all DISFormFactorsModelI
  //   interface implementations inheriting DISFormFactorsModel
  virtual void   LoadConfig (void);
  virtual void   InitPDF    (void);
  virtual double Q2         (const Interaction * i) const;
  virtual double ScalingVar (const Interaction * i) const;
  virtual void   CalcPDFs   (const Interaction * i) const;
  virtual double NuclMod    (const Interaction * i) const;
  virtual double R          (const Interaction * i) const;
  virtual void   QQBar      (const Interaction * i, double & q, double & qbar) const;
  virtual void   KFactors   (const Interaction * i, double & kuv, 
                                     double & kdv, double & kus, double & kds) const;

  mutable double fF1;
  mutable double fF2;
  mutable double fF3;
  mutable double fF4;
  mutable double fF5;
  mutable double fF6;

  PDF *  fPDF;           ///<
  PDF *  fPDFc;          ///< to take into account charm contribution
  double fQ2min;         ///< min Q^2 allowed for PDFs: PDF(Q2<Q2min):=PDF(Q2min)
  bool   fIncludeFL;     ///<
  bool   fIncludeNuclMod;///<
  bool   fCorrectF3;     ///<
  double fMc;            ///<
  double fVcd;           ///<
  double fVcs;           ///<
  double fVud;           ///<
  double fVus;           ///<
  double fVcd2;          ///<
  double fVcs2;          ///<
  double fVud2;          ///<
  double fVus2;          ///<
};

}         // genie namespace
#endif    // _DIS_STRUCTURE_FUNCTIONS_MODEL_H_

