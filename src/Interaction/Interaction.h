//____________________________________________________________________________
/*!

\class    genie::Interaction

\brief    Summary information for an interaction.

          It is a container of an InitialState, a ProcessInfo, an XclsTag
          and a Kinematics object.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 25, 2004

*/
//____________________________________________________________________________

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <ostream>
#include <string>

#include <TObject.h>

#include "Conventions/RefFrame.h"
#include "Interaction/InitialState.h"
#include "Interaction/ProcessInfo.h"
#include "Interaction/Kinematics.h"
#include "Interaction/XclsTag.h"

using std::ostream;
using std::string;

namespace genie {

const UInt_t kISkipProcessChk   = 1<<17;
const UInt_t kISkipKinematicChk = 1<<16;

class Interaction : public TObject {

public:

  Interaction();
  Interaction(const InitialState & init, const ProcessInfo & proc);
  Interaction(const Interaction & i);
  ~Interaction();

  //! Get read-only interaction information
  const InitialState & GetInitialState (void) const { return *fInitialState; }
  const ProcessInfo &  GetProcessInfo  (void) const { return *fProcInfo;     }
  const Kinematics &   GetKinematics   (void) const { return *fKinematics;   }
  const XclsTag &      GetExclusiveTag (void) const { return *fExclusiveTag; }

  //! Get read/write interaction information
  InitialState * GetInitialStatePtr (void) const { return fInitialState; }
  ProcessInfo *  GetProcessInfoPtr  (void) const { return fProcInfo;     }
  Kinematics *   GetKinematicsPtr   (void) const { return fKinematics;   }
  XclsTag *      GetExclusiveTagPtr (void) const { return fExclusiveTag; }

  //! Methods to 'block' set interaction's properties
  void SetInitialState (const InitialState & init);
  void SetProcessInfo  (const ProcessInfo &  proc);
  void SetKinematics   (const Kinematics &   kine);
  void SetExclusiveTag (const XclsTag &      xcls);

  //! Get final state primary lepton / uniquely determined from the inputs
  TParticlePDG * GetFSPrimaryLepton (void) const;

  //! Copy, reset, print itself and build string code
  void   Reset    (void);
  void   Copy     (const Interaction & i);
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  Interaction &    operator =  (const Interaction & i);
  friend ostream & operator << (ostream & stream, const Interaction & i);

private:

  //! Methods for Interaction initialization and clean up
  void Init    (void);
  void CleanUp (void);

  //! Private data members
  InitialState * fInitialState;  ///< Initial State info
  ProcessInfo *  fProcInfo;      ///< Process info (scattering, weak current,...)
  Kinematics *   fKinematics;    ///< kinematical variables
  XclsTag *      fExclusiveTag;  ///< Additional info for exclusive channels

ClassDef(Interaction,1)
};

}      // genie namespace

#endif // _INTERACTION_H_
