//____________________________________________________________________________
/*!

\class    genie::ProcessInfo

\brief    A class encapsulating an enumeration of interaction types (EM,
          Weak-CC, Weak-NC) and scattering types (Elastic, Quasi Elastic,
          Deep Inelastic, Resonant Single Pion Production, Coherent Pion
          Production).

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

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

  void Set(ScatteringType_t sc_type, InteractionType_t  int_type);

  bool IsElastic        (void) const;
  bool IsQuasiElastic   (void) const;
  bool IsDeepInelastic  (void) const;
  bool IsResonant       (void) const;
  bool IsCoherent       (void) const;
  bool IsInverseMuDecay (void) const;
  bool IsEM             (void) const;
  bool IsWeak           (void) const;
  bool IsWeakCC         (void) const;
  bool IsWeakNC         (void) const;

  ScatteringType_t  ScatteringTypeId  (void) const;
  InteractionType_t InteractionTypeId (void) const;

  string AsString                (void) const;
  string ScatteringTypeAsString  (void) const;
  string InteractionTypeAsString (void) const;

  bool Compare (const ProcessInfo & proc) const;
  void Copy    (const ProcessInfo & proc);
  void Print   (ostream & stream)         const;

  ProcessInfo &    operator =  (const ProcessInfo & proc);
  friend ostream & operator << (ostream& stream, const ProcessInfo & proc);

private:

  ScatteringType_t  fScatteringType;
  InteractionType_t fInteractionType;

ClassDef(ProcessInfo,1)
};

}        // genie namespace

#endif   // _PROCESS_INFO_H_
