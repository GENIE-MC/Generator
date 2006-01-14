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
  Interaction(const InitialState & state, const ProcessInfo & proc_info);
  Interaction(const Interaction & interaction);
  virtual ~Interaction();

  void Copy  (const Interaction & interaction);
  void Reset (void);

  //! get read-only interaction information
  const InitialState & GetInitialState (void) const { return *fInitialState; }
  const ProcessInfo &  GetProcessInfo  (void) const { return *fProcInfo;     }
  const Kinematics &   GetKinematics   (void) const { return *fKinematics;   }
  const XclsTag &      GetExclusiveTag (void) const { return *fExclusiveTag; }

  //! get read/write interaction information
  InitialState * GetInitialStatePtr (void) const { return fInitialState; }
  ProcessInfo *  GetProcessInfoPtr  (void) const { return fProcInfo;     }
  Kinematics *   GetKinematicsPtr   (void) const { return fKinematics;   }
  XclsTag *      GetExclusiveTagPtr (void) const { return fExclusiveTag; }

  //! methods to 'block' set interaction's properties
  void SetInitialState (const InitialState & init_state);
  void SetProcessInfo  (const ProcessInfo & proc_info);
  void SetKinematics   (const Kinematics & kinematics);
  void SetExclusiveTag (const XclsTag & xcls_tag);

  //! get final state primary lepton / uniquely determined from the inputs
  TParticlePDG * GetFSPrimaryLepton (void) const;

  //! printing itself and the interaction string code
  string AsString (void)             const;
  void   Print    (ostream & stream) const;

  friend ostream & operator<< (ostream & stream, const Interaction & interaction);

private:

  void Init    (void);
  void CleanUp (void);

  //! initial state, process info, scattering parameters and exclusive information
  InitialState * fInitialState;  ///< Initial State info
  ProcessInfo *  fProcInfo;      ///< Process info (scattering, weak current,...)
  Kinematics *   fKinematics;    ///< kinematical variables describing the scattering
  XclsTag *      fExclusiveTag;  ///< Additional info for exclusive channels

ClassDef(Interaction,1)
};

}      // genie namespace

#endif // _INTERACTION_H_
