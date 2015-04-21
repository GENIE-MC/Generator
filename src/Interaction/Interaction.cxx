//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2007 - CA
   Added == operator and Compare() method.
 @ Nov 21, 2007 - CA
   In order to handle the introduction of a new type of coherent interactions 
   (coherent elastic) renamed the COHCC() and COHNC() methods to COHPiCC() 
   and COHPiNC() respectivelly, and added COHEl() methods.
 @ Dec 01, 2007 - CA
   For ve- 'weak mix' interactions (ve+e->ve+e) the neutrino is always set
   as the primary final state lepton
 @ Feb 15, 2008 - CA
   Added named ctors for anomaly mediated nu-gamma interactions.
 @ Sep 24, 2008 - CA
   Added named ctors for MEC interactions.
 @ Feb 09, 2009 - CA
   Added named ctors for diffractive interactions.
 @ Mar 03, 2009 - CA
   Adapted COH name ctors (COHPi -> COH) in anticipation of including coherent
   vector meson production.
 @ Aug 21, 2009 - CR
   Added IBD() named ctors for Inverse Beta Decay interactions.
 @ Sep 18, 2009 - CA
   Added QELEM() named ctors for charged lepton QEL interactions. 
   In RecoilNucleonPdg() allow EM interactions.
 @ Oct 19, 2009 - CA
   Added RESEM() named ctors for charged lepton RES interactions.
 @ Dec 14, 2009 - CA
   Added GLR() named ctors for Glashow resonance interactions.
 @ May 05, 2010 - CR
   Adding special ctor for ROOT I/O purposes so as to avoid memory leak due to
   memory allocated in the default ctor when objects of this class are read by 
   the ROOT Streamer. 
 @ Nov 17, 2011 - CA
   Added NDecay() named ctor. Removed unused Compare() method and operator.
 @ Nov 24, 2011 - CA
   Tweaked RecoilNucleonPdg() so that it works with MEC wgere the hit object
   is a nucleon-cluster and not a single nucleon. The MEC named ctors now 
   have an argument for specifying the hit nucleon cluster PDG.
 @ Feb 29, 2012 - CA
   Added MECNC() and MECEM().
 @ Apr 20, 2012 - CA
   Added DISEM().
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   In elastic neutrino-electron scattering, always set the electron as the 
   final state primary lepton. Handle the IMD annihilation channel.
*/
//____________________________________________________________________________

#include <sstream>

#include <TRootIOCtor.h>

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
Interaction::Interaction() :
TObject()
{
  this->Init();
}
//___________________________________________________________________________
Interaction::Interaction(const InitialState & ist, const ProcessInfo & prc) :
TObject()
{
  this->Init();

  fInitialState -> Copy (ist);
  fProcInfo     -> Copy (prc);
}
//___________________________________________________________________________
Interaction::Interaction(const Interaction & interaction) :
TObject()
{
  this->Init();
  this->Copy(interaction);
}
//___________________________________________________________________________
Interaction::Interaction(TRootIOCtor*) :
TObject(),
fInitialState(0), 
fProcInfo(0),
fKinematics(0), 
fExclusiveTag(0), 
fKinePhSp(0)
{

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
  fKinePhSp     = new KPhaseSpace  (this);
}
//___________________________________________________________________________
void Interaction::CleanUp(void)
{
  if ( fInitialState ) delete fInitialState;
  if ( fProcInfo     ) delete fProcInfo;
  if ( fKinematics   ) delete fKinematics;
  if ( fExclusiveTag ) delete fExclusiveTag;
  if ( fKinePhSp     ) delete fKinePhSp;

  fInitialState = 0;
  fProcInfo     = 0;
  fKinematics   = 0;
  fExclusiveTag = 0;
  fKinePhSp     = 0;
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
TParticlePDG * Interaction::FSPrimLepton(void) const
{
  int pdgc = this->FSPrimLeptonPdg();

  if(pdgc) return PDGLibrary::Instance()->Find(pdgc);
  else     return 0;
}
//___________________________________________________________________________
int Interaction::FSPrimLeptonPdg(void) const
{
  const ProcessInfo &  proc_info  = this -> ProcInfo();
  const InitialState & init_state = this -> InitState();

  int pdgc = init_state.ProbePdg();

  LOG("Interaction", pDEBUG) << "Probe PDG code: " << pdgc;

  if (proc_info.IsNuElectronElastic())
    return kPdgElectron;

  // vN (Weak-NC) or eN (EM)
  if (proc_info.IsWeakNC() || proc_info.IsEM() || proc_info.IsWeakMix()) return pdgc;

  // vN (Weak-CC)
  else if (proc_info.IsWeakCC()) {
    int clpdgc;
    if (proc_info.IsIMDAnnihilation())
      clpdgc = kPdgMuon;
    else
      clpdgc = pdg::Neutrino2ChargedLepton(pdgc);
    return clpdgc;
  }

  LOG("Interaction", pWARN)
        << "Could not figure out the final state primary lepton pdg code!!";

  return 0;
}
//___________________________________________________________________________
TParticlePDG * Interaction::RecoilNucleon(void) const
{
  int rnuc = this->RecoilNucleonPdg();

  if(rnuc) return PDGLibrary::Instance()->Find(rnuc);
  else     return 0;
}
//___________________________________________________________________________
int Interaction::RecoilNucleonPdg(void) const
{
// Determine the recoil nucleon PDG code

  const Target & target = fInitialState->Tgt();

  int recoil_nuc = 0;
  int struck_nuc = target.HitNucPdg();

  if(fProcInfo->IsQuasiElastic() || fProcInfo->IsInverseBetaDecay()) {
    bool struck_is_nuc = pdg::IsNucleon(struck_nuc);
    bool is_weak = fProcInfo->IsWeak();
    bool is_em   = fProcInfo->IsEM();
    assert(struck_is_nuc && (is_weak || is_em));
    if(fProcInfo->IsWeakCC()) {
       recoil_nuc = pdg::SwitchProtonNeutron(struck_nuc); // CC
    } else {
       recoil_nuc = struck_nuc; // NC, EM
    }
  }

  if(fProcInfo->IsMEC()) {
    bool struck_is_2nuc_cluster = pdg::Is2NucleonCluster(struck_nuc);
    bool is_weak = fProcInfo->IsWeak();
    bool is_em   = fProcInfo->IsEM();
    assert(struck_is_2nuc_cluster && (is_weak || is_em));
    if(fProcInfo->IsWeakCC()) {
       bool isnu = pdg::IsNeutrino(fInitialState->ProbePdg());
       // nucleon cluster charge should be incremented by +1 for 
       // neutrino CC and by -1 for antineutrino CC
       int dQ = (isnu) ? +1 : -1;
       recoil_nuc = pdg::ModifyNucleonCluster(struck_nuc,dQ); // CC
    }
    else {
       recoil_nuc = struck_nuc; // NC, EM
    }
  }

  LOG("Interaction", pDEBUG) << "Recoil nucleon PDG = " << recoil_nuc;
  return recoil_nuc;
}
//___________________________________________________________________________
void Interaction::SetInitState(const InitialState & init_state)
{
  if (!fInitialState) fInitialState = new InitialState();
  fInitialState->Copy(init_state);
}
//___________________________________________________________________________
void Interaction::SetProcInfo(const ProcessInfo & proc_info)
{
  if (!fProcInfo) fProcInfo = new ProcessInfo();
  fProcInfo->Copy(proc_info);
}
//___________________________________________________________________________
void Interaction::SetKine(const Kinematics & kinematics)
{
  if (!fKinematics) fKinematics = new Kinematics();
  fKinematics->Copy(kinematics);
}
//___________________________________________________________________________
void Interaction::SetExclTag(const XclsTag & xcls_tag)
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

  const Target & tgt = fInitialState->Tgt();

  ostringstream interaction;

  interaction << "nu:"  << fInitialState->ProbePdg() << ";";
  interaction << "tgt:" << tgt.Pdg() << ";";

  if(tgt.HitNucIsSet()) {
    interaction << "N:" << tgt.HitNucPdg() << ";";
  }
  if(tgt.HitQrkIsSet()) {
    interaction << "q:" << tgt.HitQrkPdg()
                << (tgt.HitSeaQrk() ? "(s)" : "(v)") << ";";
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

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, double E)
{
  Interaction* interaction = Interaction::DISCC(target,hitnuc,probe,E);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScDeepInelastic, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISCC(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, 
   const TLorentzVector & p4probe)
{
  Interaction* interaction = Interaction::DISCC(target,hitnuc,probe,p4probe);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, double E)
{
  Interaction* interaction = Interaction::DISNC(target,hitnuc,probe,E);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScDeepInelastic, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISNC(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, 
   const TLorentzVector & p4probe)
{
  Interaction * interaction = Interaction::DISNC(target,hitnuc,probe,p4probe);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISEM(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
             Interaction::Create(target,probe,kScDeepInelastic, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISEM(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, double E)
{
  Interaction* interaction = Interaction::DISEM(target,hitnuc,probe,E);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISEM(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScDeepInelastic, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DISEM(
   int target, int hitnuc, int hitqrk, bool fromsea, int probe, 
   const TLorentzVector & p4probe)
{
  Interaction * interaction = Interaction::DISEM(target,hitnuc,probe,p4probe);

  Target * tgt = interaction->InitStatePtr()->TgtPtr();
  tgt -> SetHitQrkPdg (hitqrk);
  tgt -> SetHitSeaQrk (fromsea);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELCC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELCC(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELNC(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELEM(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::QELEM(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScQuasiElastic, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::IBD(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction =
     Interaction::Create(target,probe,kScInverseBetaDecay,kIntWeakCC);
   
  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);
   
  return interaction;
} 
//___________________________________________________________________________ 
Interaction * Interaction::IBD(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction =
      Interaction::Create(target,probe,kScInverseBetaDecay,kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);
  
  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESCC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESCC(
           int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESNC(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESNC(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESEM(int target, int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::RESEM(
   int target, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,probe,kScResonant, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DFRCC(int tgt,int hitnuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScDiffractive, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::DFRCC(
   int tgt, int hitnuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScDiffractive, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(hitnuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHCC(int tgt, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherent, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHCC(
    int tgt, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherent, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHNC(int tgt, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherent, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHNC(
   int tgt, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherent, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHEl(int tgt, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherentElas, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::COHEl(
   int tgt, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScCoherentElas, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::IMD(int target, double E)
{
  Interaction * interaction = 
     Interaction::Create(target,kPdgNuMu,kScInverseMuDecay, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::IMD(int target, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(target,kPdgNuMu,kScInverseMuDecay, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::AMNuGamma(int tgt, int nuc, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScAMNuGamma, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(nuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::AMNuGamma(
   int tgt, int nuc, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScAMNuGamma, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(nuc);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECCC(int tgt, int ncluster, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECCC(
   int tgt, int ncluster, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECNC(int tgt, int ncluster, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECNC(
   int tgt, int ncluster, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntWeakNC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECEM(int tgt, int ncluster, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::MECEM(
   int tgt, int ncluster, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt, probe, kScMEC, kIntEM);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(ncluster);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::GLR(int tgt, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt, kPdgAntiNuE, kScGlashowResonance, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);
  init_state->TgtPtr()->SetHitNucPdg(0);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::GLR(int tgt, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt, kPdgAntiNuE, kScGlashowResonance, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);
  init_state->TgtPtr()->SetHitNucPdg(0);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::NDecay(int tgt, int decay_mode)
{
  Interaction * interaction = 
     Interaction::Create(tgt, 0, kScNull, kIntNDecay);
  interaction->ExclTagPtr()->SetDecayMode(decay_mode);
  return interaction;
}
//___________________________________________________________________________
//___________________________________________________________________________
Interaction * Interaction::ASK(int tgt, int probe, double E)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScSingleKaon, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeE(E);

  return interaction;
}
//___________________________________________________________________________
Interaction * Interaction::ASK(
    int tgt, int probe, const TLorentzVector & p4probe)
{
  Interaction * interaction = 
     Interaction::Create(tgt,probe,kScSingleKaon, kIntWeakCC);

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4(p4probe);

  return interaction;
}

