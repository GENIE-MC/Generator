//____________________________________________________________________________
/*!

\class    genie::ProcessInfo

\brief    A class encapsulating an enumeration of interaction types (EM,
          Weak-CC, Weak-NC) and scattering types (Elastic, Quasi Elastic,
          Deep Inelastic, Resonant Single Pion Production, Coherent Pion
          Production).

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PROCESS_INFO_H_
#define _PROCESS_INFO_H_

#include <iostream>
#include <string>

#include <TObject.h>

#include "Interaction/InteractionType.h"
#include "Interaction/ScatteringType.h"

using std::ostream;
using std::string;

namespace genie {

class ProcessInfo : public TObject {

public:
  ProcessInfo();
  ProcessInfo(ScatteringType_t sc_type, InteractionType_t  int_type);
  ProcessInfo(const ProcessInfo & proc);
 ~ProcessInfo();

  //-- set process information
  void Set(ScatteringType_t sc_type, InteractionType_t  int_type);

  //-- query for process information
  bool IsQuasiElastic      (void) const;
  bool IsDeepInelastic     (void) const;
  bool IsResonant          (void) const;
  bool IsCoherent          (void) const;
  bool IsCoherentElas      (void) const;
  bool IsSingleKaon        (void) const;
  bool IsElectronScattering(void) const;
  bool IsNuElectronElastic (void) const;
  bool IsInverseMuDecay    (void) const;
  bool IsIMDAnnihilation   (void) const;
  bool IsInverseBetaDecay  (void) const;
  bool IsGlashowResonance  (void) const;
  bool IsAMNuGamma         (void) const;
  bool IsMEC               (void) const;
  bool IsDiffractive       (void) const;
  bool IsEM                (void) const;
  bool IsWeak              (void) const;
  bool IsWeakCC            (void) const;
  bool IsWeakNC            (void) const;
  bool IsWeakMix           (void) const;

  //-- get scattering and interaction type enumerations
  ScatteringType_t  ScatteringTypeId  (void) const;
  InteractionType_t InteractionTypeId (void) const;

  //-- get scattering and interaction types as strings
  string ScatteringTypeAsString  (void) const;
  string InteractionTypeAsString (void) const;

  //-- Copy, reset, compare, print itself and build string code
  void   Reset    (void);
  void   Copy     (const ProcessInfo & proc);
  bool   Compare  (const ProcessInfo & proc) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const ProcessInfo & proc) const;
  ProcessInfo &    operator =  (const ProcessInfo & proc);
  friend ostream & operator << (ostream& stream, const ProcessInfo & proc);

private:

  //-- Private data members
  ScatteringType_t  fScatteringType;  ///< scattering type  (QEL, RES, DIS, ...)
  InteractionType_t fInteractionType; ///< interaction type (Weak CC/NC, E/M, ...)

ClassDef(ProcessInfo,1)
};

}        // genie namespace

#endif   // _PROCESS_INFO_H_
