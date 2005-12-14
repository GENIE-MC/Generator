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

#include <sstream>

#include "Conventions/Constants.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

using std::endl;
using std::ostringstream;

ClassImp(Interaction)

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const Interaction & interaction)
 {
   interaction.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
Interaction::Interaction()
{
  this->Init();
}
//___________________________________________________________________________
Interaction::Interaction(
              const InitialState & init_state, const ProcessInfo & proc_info)
{
  this->Init();

  fInitialState -> Copy (init_state);
  fProcInfo     -> Copy (proc_info);
}
//___________________________________________________________________________
Interaction::Interaction(const Interaction & interaction)
{
  this->Init();
  this->Copy(interaction);
}
//___________________________________________________________________________
Interaction::~Interaction()
{
  this->Reset();
}
//___________________________________________________________________________
void Interaction::Copy(const Interaction & interaction)
{
  fInitialState = new InitialState ( *interaction.fInitialState  );
  fProcInfo     = new ProcessInfo  ( *interaction.fProcInfo      );
  fKinematics   = new Kinematics   ( *interaction.fKinematics    );
  fExclusiveTag = new XclsTag      ( * interaction.fExclusiveTag );
}
//___________________________________________________________________________
void Interaction::Reset(void)
{
  if ( fInitialState ) delete fInitialState;
  if ( fProcInfo     ) delete fProcInfo;
  if ( fKinematics   ) delete fKinematics;
  if ( fExclusiveTag ) delete fExclusiveTag;

  fInitialState = 0;
  fProcInfo     = 0;
  fKinematics   = 0;
  fExclusiveTag = 0;

  this->Init();
}
//___________________________________________________________________________
void Interaction::Init(void)
{
  fInitialState = new InitialState ();
  fProcInfo     = new ProcessInfo  ();
  fKinematics   = new Kinematics   ();
  fExclusiveTag = new XclsTag      ();
}
//___________________________________________________________________________
TParticlePDG * Interaction::GetFSPrimaryLepton(void) const
{
  const ProcessInfo &  proc_info  = this -> GetProcessInfo();
  const InitialState & init_state = this -> GetInitialState();

  int pdgc = init_state.GetProbePDGCode();

  LOG("Interaction", pDEBUG) << "Probe PDG code: " << pdgc;

  // -- vN (Weak-NC) or eN (EM)
  if ( proc_info.IsWeakNC() || proc_info.IsEM() )
  {
     return PDGLibrary::Instance()->Find(pdgc);
  }

  // -- vN (Weak-CC)
  else if ( proc_info.IsWeakCC() )
  {
     int clpdgc = pdg::Neutrino2ChargedLepton(pdgc);
     return PDGLibrary::Instance()->Find(clpdgc);
  }
  LOG("Interaction", pWARN)
        << "Could not figure out the final state primary lepton pdg code!!";

  return 0;
}
//___________________________________________________________________________
void Interaction::SetInitialState(const InitialState & init_state)
{
  if (!fInitialState) fInitialState = new InitialState();
  fInitialState->Copy(init_state);
}
//___________________________________________________________________________
void Interaction::SetProcessInfo(const ProcessInfo & proc_info)
{
  if (!fProcInfo) fProcInfo = new ProcessInfo();
  fProcInfo->Copy(proc_info);
}
//___________________________________________________________________________
void Interaction::SetKinematics(const Kinematics & kinematics)
{
  if (!fKinematics) fKinematics = new Kinematics();
  fKinematics->Copy(kinematics);
}
//___________________________________________________________________________
void Interaction::SetExclusiveTag(const XclsTag & xcls_tag)
{
  if (!fExclusiveTag) fExclusiveTag = new XclsTag();
  fExclusiveTag->Copy(xcls_tag);
}
//___________________________________________________________________________
void Interaction::Print(ostream & stream) const
{
  const string line(110, '-');

  stream << endl;
  stream << line << endl;

  stream << "GENIE Interaction Summary" << endl;
  stream << line << endl;

  stream << *fInitialState << endl; // print initial state
  stream << *fProcInfo;             // print process info
  stream << *fKinematics;           // print scattering parameters
  stream << *fExclusiveTag;         // print exclusive process tag

  stream << line << endl;
}
//___________________________________________________________________________
string Interaction::AsString(void) const
{
// Code-ify the interaction in a string to be used as (part of a) cache
// branch key.
// Template:
// nu:x;tgt:x;N:x;q:x(s/v);intp:x;sctp:x;xclv:x

  const Target & tgt = fInitialState->GetTarget();

  ostringstream interaction;

  interaction << "nu:"  << fInitialState->GetProbePDGCode() << ";";
  interaction << "tgt:" << tgt.PDGCode() << ";";

  if(tgt.StruckNucleonIsSet()) {
    interaction << "N:" << tgt.StruckNucleonPDGCode() << ";";
  }
  if(tgt.StruckQuarkIsSet()) {
    interaction << "q:" << tgt.StruckQuarkPDGCode() 
                << (tgt.StruckQuarkIsFromSea() ? "(s)" : "(v)") << ";";
  }

  interaction << "intp:" << fProcInfo->InteractionTypeAsString() << ";";
  interaction << "sctp:" << fProcInfo->ScatteringTypeAsString()  << ";";

  interaction << "xclv:" << fExclusiveTag->AsString() << ";";

  return interaction.str();
}
//___________________________________________________________________________

