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

  void Copy  (const Interaction & interaction);
  void Reset (void);

  //! get read-only interaction information
  const InitialState &     GetInitialState     (void) const { return *fInitialState;     }
  const ProcessInfo &      GetProcessInfo      (void) const { return *fProcInfo;         }
  const ScatteringParams & GetScatteringParams (void) const { return *fScatteringParams; }
  const XclsTag &          GetExclusiveTag     (void) const { return *fExclusiveTag;     }

  //! get read/write interaction information
  ScatteringParams *       GetScatParamsPtr    (void) const { return fScatteringParams;  }
  InitialState *           GetInitialStatePtr  (void) const { return fInitialState;      }
  
  //! methods to set/reset and examine whether an XclsTag object has been attached
  void   SetExclusiveTag (const XclsTag & xcls_tag);
  void   ResetExclusive  (void);
  bool   IsExclusive     (void) const { return (fExclusiveTag != 0); }

  //! get final state primary lepton / uniquely determined from the inputs
  TParticlePDG * GetFSPrimaryLepton (void) const; // tmp
  
  //! setting/getting cross sections set during event generation
  void   SetXSec     (double xsec) { fXSec  = xsec; } // to be set when selecting interaction
  void   SetDiffXSec (double xsec) { fdXSec = xsec; } // to be set when selecting kinematics
  double XSec        (void) const  { return fXSec;  }
  double DiffXSec    (void) const  { return fdXSec; }
  
  //! printing itself and the interaction string code
  string AsString (void)             const;
  void   Print    (ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const Interaction & interaction);
  
private:

  void Init (void);

  //! initial state, process info, scattering parameters and exclusive information
  InitialState *     fInitialState;     ///< Initial State info
  ProcessInfo *      fProcInfo;         ///< Process info (scattering, weak current,...)
  ScatteringParams * fScatteringParams; ///< Scattering parameters
  XclsTag *          fExclusiveTag;     ///< Additional info for exclusive channels

  //! cross section for this interaction
  double fXSec;    ///< xsec for the given energy
  double fdXSec;   ///< diff. xsec for the given kinematical parameters
};

}      // genie namespace

#endif // _INTERACTION_H_
