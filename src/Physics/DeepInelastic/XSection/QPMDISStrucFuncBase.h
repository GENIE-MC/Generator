//____________________________________________________________________________
/*!

\class    genie::QPMDISStrucFuncBase

\brief    Abstract base class. 
          Provides common implementation for concrete objects implementing the
          DISStructureFuncModelI interface.

\ref      For a discussion of DIS SF see for example E.A.Paschos and J.Y.Yu, 
          Phys.Rev.D 65.033002 and R.Devenish and A.Cooper-Sarkar, OUP 2004.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Adapted from neugen 3.
          Primary authors: 
          D.Naples (Pittsburgh U.), H.Gallagher (Tufts U), CA (RAL)

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QPM_DIS_STRUCTURE_FUNCTIONS_BASE_H_
#define _QPM_DIS_STRUCTURE_FUNCTIONS_BASE_H_

#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDF.h"

namespace genie {

class QPMDISStrucFuncBase : public DISStructureFuncModelI {

public:
  virtual ~QPMDISStrucFuncBase();

  // common code for all DISFormFactorsModelI interface implementations
  virtual double F1 (void) const { return fF1; }
  virtual double F2 (void) const { return fF2; }
  virtual double F3 (void) const { return fF3; }
  virtual double F4 (void) const { return fF4; }
  virtual double F5 (void) const { return fF5; }
  virtual double F6 (void) const { return fF6; }

  virtual void Calculate (const Interaction * interaction) const;

  // overload Algorithm's Configure() to set the PDF data member
  // from the configuration registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

protected:
  QPMDISStrucFuncBase();
  QPMDISStrucFuncBase(string name);
  QPMDISStrucFuncBase(string name, string config);

  // commom code for SF calculation for all DISFormFactorsModelI
  // interface implementations inheriting from QPMDISStrucFuncBase
  virtual void   LoadConfig (void);
  virtual void   InitPDF    (void);
  virtual double Q2         (const Interaction * i) const;
  virtual double ScalingVar (const Interaction * i) const;
  virtual void   CalcPDFs   (const Interaction * i) const;
  virtual double NuclMod    (const Interaction * i) const;
  virtual double R          (const Interaction * i) const;
  virtual void   KFactors   (const Interaction * i, double & kuv, 
                                     double & kdv, double & kus, double & kds) const;
  // configuration
  //
  double fQ2min;             ///< min Q^2 allowed for PDFs: PDF(Q2<Q2min):=PDF(Q2min)
  bool   fCharmOff;          ///< turn charm production off?
  bool   fIncludeR;          ///< include R (~FL) in DIS SF calculation?
  bool   fIncludeNuclMod;    ///< include nuclear factor (shadowing, anti-shadowing,...)?
  double fMc;                ///< charm mass used
  double fVcd;               ///< CKM element Vcd used
  double fVcs;               ///< CKM element Vcs used
  double fVud;               ///< CKM element Vud used
  double fVus;               ///< CKM element Vcs used
  double fVcd2;              ///<
  double fVcs2;              ///<
  double fVud2;              ///<
  double fVus2;              ///<
  double fSin2thw;           ///<
  bool   fUse2016Corrections;///< Use 2016 SF relation corrections
  double fLowQ2CutoffF1F2;   ///< Set min for relation between 2xF1 and F2

  mutable double fF1;
  mutable double fF2;
  mutable double fF3;
  mutable double fF4;
  mutable double fF5;
  mutable double fF6;
  PDF *  fPDF;           ///< computed PDFs @ (x,Q2)
  PDF *  fPDFc;          ///< computed PDFs @ (slow-rescaling-var,Q2)
  mutable double fuv;
  mutable double fus; 
  mutable double fdv; 
  mutable double fds; 
  mutable double fs;
  mutable double fc;
  mutable double fuv_c; 
  mutable double fus_c; 
  mutable double fdv_c;
  mutable double fds_c;
  mutable double fs_c; 
  mutable double fc_c; 

};

}         // genie namespace
#endif    // _QPM_DIS_STRUCTURE_FUNCTIONS_BASE_H_

