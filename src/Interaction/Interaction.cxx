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

  fInitialState = new InitialState (init_state);
  fProcInfo     = new ProcessInfo  (proc_info);

  fScatteringParams = new ScatteringParams ();
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
  fInitialState = new InitialState (*interaction.fInitialState);
  fProcInfo     = new ProcessInfo  (*interaction.fProcInfo    );

  try {
     const Registry & reg =
              dynamic_cast<const Registry &>(*interaction.fScatteringParams);

    fScatteringParams = new ScatteringParams(reg);

  } catch( std::bad_cast ) {
     LOG("Interaction", pERROR) << "dynamic_cast to const Registry & failed";
  }

  if( interaction.IsExclusive() )
       fExclusiveTag = new XclsTag(* interaction.fExclusiveTag);
  else fExclusiveTag = 0;

  fXSec  = interaction.fXSec;
  fdXSec = interaction.fdXSec;
}
//___________________________________________________________________________
void Interaction::Reset(void)
{
  if ( fInitialState     ) delete fInitialState;
  if ( fProcInfo         ) delete fProcInfo;
  if ( fScatteringParams ) delete fScatteringParams;
  if ( fExclusiveTag     ) delete fExclusiveTag;

  this->Init();
}
//___________________________________________________________________________
void Interaction::Init(void)
{
  fInitialState     = 0;
  fProcInfo         = 0;
  fScatteringParams = 0;
  fExclusiveTag     = 0;
  fXSec             = 0;
  fdXSec            = 0;
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
void Interaction::SetExclusiveTag(const XclsTag & xcls_tag)
{
  fExclusiveTag = new XclsTag(xcls_tag);
}
//___________________________________________________________________________
void Interaction::ResetExclusive(void)
{
  if(fExclusiveTag) {
     LOG("Interaction", pDEBUG)  << "Reseting exclusive tag";
     delete fExclusiveTag;
     fExclusiveTag = 0;
  }
}
//___________________________________________________________________________
void Interaction::Print(ostream & stream) const
{
  const string line(100, '-');

  stream << endl;
  stream << line << endl;

  stream << *fInitialState << endl; // print initial state
  stream << *fProcInfo;             // print process info
  stream << *fScatteringParams;     // print scattering parameters

  // print exclusive process tag - if exists
  if( this->IsExclusive() ) stream << *fExclusiveTag;

  // print cross section information - if exists
  if(fXSec>0 || fdXSec>0) {
    stream << "[-] [Cross Sections]" << endl;
    if(fXSec>0) {
       stream << " |--> xsec (@ given E) = "
              << fXSec << endl;
    }
    if(fdXSec>0) {
       stream << " |--> diff. xsec (@ given E & kinematical vars) = "
              << fdXSec << endl;
    }
  }
  stream << line << endl;
}
//___________________________________________________________________________
string Interaction::AsString(void) const
{
// Code-ify the interaction in a string to be used as (part of a) cache
// branch key.
//
// Template:
//     nu_pdg:code;tgt-pdg:code;nucl-pdg:code;intype:name;sctype:name

  ostringstream interaction;

  interaction << "nu-pdg:"
              << fInitialState->GetProbePDGCode() << ";";
  interaction << "tgt-pdg:"
              << fInitialState->GetTarget().PDGCode() << ";";
  interaction << "nucl-pdg:"
              << fInitialState->GetTarget().StruckNucleonPDGCode() << ";";
  interaction << "intype:"
              << fProcInfo->InteractionTypeAsString() << ";";
  interaction << "sctype:"
              << fProcInfo->ScatteringTypeAsString() << ";";

  if( this->IsExclusive() )
                interaction << "xclv:" << fExclusiveTag->AsString() << ";";

  return interaction.str();
}
//___________________________________________________________________________

