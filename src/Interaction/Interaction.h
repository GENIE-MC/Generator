//____________________________________________________________________________
/*!

\class    genie::Interaction

\brief    Summary information for an interaction.

          It is a container of an InitialState, a ProcessInfo, an XclsTag
          and a ScatteringParams object.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 25, 2004

*/
//____________________________________________________________________________

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <ostream>
#include <string>

#include "Conventions/RefFrame.h"
#include "Interaction/InitialState.h"
#include "Interaction/ProcessInfo.h"
#include "Interaction/ScatteringParams.h"
#include "Interaction/XclsTag.h"

using std::ostream;
using std::string;

namespace genie {

class Interaction {

public:

  Interaction();
  Interaction(const InitialState & state, const ProcessInfo & proc_info);
  Interaction(const Interaction & interaction);
  virtual ~Interaction();

  const InitialState &     GetInitialState     (void) const { return *fInitialState;     }
  const ProcessInfo &      GetProcessInfo      (void) const { return *fProcInfo;         }
  const ScatteringParams & GetScatteringParams (void) const { return *fScatteringParams; }
  const XclsTag &          GetExclusiveTag     (void) const { return *fExclusiveTag;     }

  ScatteringParams *       GetScatParamsPtr    (void) const { return fScatteringParams;  }
  InitialState *           GetInitialStatePtr  (void) const { return fInitialState;      }
  
  void   SetExclusiveTag (const XclsTag & xcls_tag);
  void   ResetExclusive  (void);
  bool   IsExclusive     (void) const { return (fExclusiveTag != 0); }
  double CrossSection    (void) const { return fXSec; }

  TParticlePDG * GetFSPrimaryLepton (void) const; // tmp
  
  string AsString (void)             const;
  void   Print    (ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const Interaction & interaction);
  
private:

  InitialState *     fInitialState;     ///< Initial State info
  ProcessInfo *      fProcInfo;         ///< Process info (scattering, weak current,...)
  ScatteringParams * fScatteringParams; ///< Scattering parameters
  XclsTag *          fExclusiveTag;     ///< Additional info for exclusive channels
  double             fXSec;             ///< XSec for the process as set during event generation
};

}      // genie namespace

#endif // _INTERACTION_H_
