/*!

\class    genie::BY21StrucFunc

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

class BY21StrucFunc : public QPMDISStrucFuncBase {

public:
  BY21StrucFunc();
  BY21StrucFunc(string config);
  virtual ~BY21StrucFunc();

  // overload Algorithm::Configure() to read the config. registry
  // at the algorithm initialization and set private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  void Init         (void);
  void ReadBYParams (void);

  // override part of the DISStructureFuncModel implementation
  // to compute all the corrections applied by the Bodek-Yang model.
  double ScalingVar (const Interaction * i, double Mf = 0 ) const;
  void   KVectorFactors (const Interaction * i, double & kuv, double & kdv, double & kus, double & kds, double & kss ) const;
  void   KAxialFactors(const Interaction * i, double & kuv, double & kdv, double & kus, double & kds, double & kss ) const ;
  
  double R(const Interaction * interaction) const ; // overrides QPMDISStrucFuncBase implementation
  double H(const Interaction * interaction) const ; // overrides QPMDISStrucFuncBase implementation
  double KCharm(const Interaction * i, double Mf = 0) const; // overrides QPMDISStrucFuncBase implementation
  // Bodek-Yang model-specific parameters

  double fMv;    ///< Vector Mass
  double fMv2;   ///< Vector Mass Squared
  double fA;     ///< better scaling var parameter A
  double fB;     ///< better scaling var parameter B
  double fCvLW;   ///< C low-nu vector paramter
  double fCsU;   ///< U-sea K factor parameter
  double fCsD;   ///< D-sea K factor parameter
  double fCsS;   ///< S-sea K factor parameter
  double fCv1U;  ///< U-val K factor parameter
  double fCv2U;  ///< U-val K factor parameter
  double fCv1D;  ///< D-val K factor parameter
  double fCv2D;  ///< D-val K factor parameter
  double fPsA;   ///< P-axial sea parameter
  double fPvA;   ///< P-axial valance paramter
  double fCsA;   ///< C-axial sea parameter
  double fCaLW_nubar; ///< C-axial neutrino LW parameter
  double fCaLW_nu; ///< C-axial anti-neutrino LW parameter
  double fH0;      /// high order QCD parameter
  double fH1;      /// high order QCD parameter
  double fH2;      /// high order QCD parameter
  double fH3;      /// high order QCD parameter
  double fRQ2min;    /// Q2 below corrections are applied  
  bool   fIncludeH; ///< Include H correction 
  bool   fIncludeAxial; ///< Include difference between A and V;
  bool   fIncludeKCharm; ///< Include KCharm correction

};

}         // genie namespace

#endif    // _BODEK_YANG_STRUCTURE_FUNCTION_MODEL_2021_H_
