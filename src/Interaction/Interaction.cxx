//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - April 25, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
Interaction::Interaction(const InitialState & ist, const ProcessInfo & prc)
{
  this->Init();

  fInitialState -> Copy (ist);
  fProcInfo     -> Copy (prc);
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
  this->CleanUp();
}
//___________________________________________________________________________
void Interaction::Reset(void)
{
  this->CleanUp();
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
void Interaction::CleanUp(void)
{
  if ( fInitialState ) delete fInitialState;
  if ( fProcInfo     ) delete fProcInfo;
  if ( fKinematics   ) delete fKinematics;
  if ( fExclusiveTag ) delete fExclusiveTag;

  fInitialState = 0;
  fProcInfo     = 0;
  fKinematics   = 0;
  fExclusiveTag = 0;
}
//___________________________________________________________________________
void Interaction::Copy(const Interaction & interaction)
{
  const InitialState & init = *interaction.fInitialState;
  const ProcessInfo &  proc = *interaction.fProcInfo;
  const Kinematics &   kine = *interaction.fKinematics;
  const XclsTag &      xcls = *interaction.fExclusiveTag;

  fInitialState -> Copy (init);
  fProcInfo     -> Copy (proc);
  fKinematics   -> Copy (kine);
  fExclusiveTag -> Copy (xcls);
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
int Interaction::RecoilNuclPDGCode(void) const
{
// Determine the recoil nucleon PDG code

  const Target & target = fInitialState->GetTarget();

  int recoil_nuc = 0;
  int struck_nuc = target.StruckNucleonPDGCode();

  if(fProcInfo->IsQuasiElastic()) {
    assert(pdg::IsNeutronOrProton(struck_nuc) && fProcInfo->IsWeak());
    if(fProcInfo->IsWeakCC()) 
       recoil_nuc = pdg::SwitchProtonNeutron(struck_nuc); // CC
    else 
       recoil_nuc = struck_nuc; // NC
  }

  LOG("Interaction", pDEBUG) << "Recoil nucleon PDG = " << recoil_nuc;
  return recoil_nuc;
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
string Interaction::AsString(void) const
{
// Code-ify the interaction in a string to be used as (part of a) cache
// branch key.
// Template:
// nu:x;tgt:x;N:x;q:x(s/v);proc:x;xclv_tag

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

  interaction << "proc:" << fProcInfo->InteractionTypeAsString() 
              << "," << fProcInfo->ScatteringTypeAsString()  << ";";

  string xcls = fExclusiveTag->AsString();
  interaction << xcls;
  if(xcls.size()>0) interaction << ";";

  return interaction.str();
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
Interaction & Interaction::operator = (const Interaction & interaction)
{
  this->Copy(interaction);
  return (*this);
}
//___________________________________________________________________________
//
//       **** Methods using the "named constructor" C++ idiom ****
//
//___________________________________________________________________________
Interaction * Interaction::Create(
            int target, int probe, ScatteringType_t st, InteractionType_t it)
{
  InitialState init_state (target, probe);
  ProcessInfo  proc_info  (st, it);

  Interaction * interaction = new Interaction(init_state, proc_info);
  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
       int target, int hitnuc, int hitqrk, bool fromsea, int probe, double E)
{
  Interaction* interaction = Interaction::DISCC(target,hitnuc,probe,E);

  Target * tgt = interaction->GetInitialStatePtr()->GetTargetPtr();
  tgt -> SetStruckQuarkPDGCode (hitqrk);
  tgt -> SetStruckSeaQuark     (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
             int target, int hitnuc, int hitqrk, bool fromsea, 
                                   int probe, const TLorentzVector & p4probe)
{
  Interaction* interaction = Interaction::DISCC(target,hitnuc,probe,p4probe);

  Target * tgt = interaction->GetInitialStatePtr()->GetTargetPtr();
  tgt -> SetStruckQuarkPDGCode (hitqrk);
  tgt -> SetStruckSeaQuark     (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
       int target, int hitnuc, int hitqrk, bool fromsea, int probe, double E)
{
  Interaction* interaction = Interaction::DISNC(target,hitnuc,probe,E);

  Target * tgt = interaction->GetInitialStatePtr()->GetTargetPtr();
  tgt -> SetStruckQuarkPDGCode (hitqrk);
  tgt -> SetStruckSeaQuark     (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
             int target, int hitnuc, int hitqrk, bool fromsea, 
                                   int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = Interaction::DISNC(target,hitnuc,probe,p4probe);

  Target * tgt = interaction->GetInitialStatePtr()->GetTargetPtr();
  tgt -> SetStruckQuarkPDGCode (hitqrk);
  tgt -> SetStruckSeaQuark     (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELCC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
              Interaction::Create(target,probe,kScQuasiElastic, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELCC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
              Interaction::Create(target,probe,kScQuasiElastic, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
              Interaction::Create(target,probe,kScQuasiElastic, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELNC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
              Interaction::Create(target,probe,kScQuasiElastic, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESCC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
                  Interaction::Create(target,probe,kScResonant, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESCC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
              Interaction::Create(target,probe,kScResonant, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
                   Interaction::Create(target,probe,kScResonant, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESNC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
                   Interaction::Create(target,probe,kScResonant, kIntWeakNC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->GetTargetPtr()->SetStruckNucleonPDGCode(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::IMD(int target, double E)
{
  Interaction * interaction = 
          Interaction::Create(target,kPdgNuMu,kScInverseMuDecay, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::IMD(int target, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
          Interaction::Create(target,kPdgNuMu,kScInverseMuDecay, kIntWeakCC);

  InitialState * init_state = interaction->GetInitialStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}
//___________________________________________________________________________




