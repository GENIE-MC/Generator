//____________________________________________________________________________
/*!

\namespace  genie::interaction_utils

\brief      Interaction Utilities 

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/ 
//____________________________________________________________________________

#include <cassert>

#include "Interaction/IUtils.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
int genie::interaction_utils::RecoilNucleonPdgCode(
                                              const Interaction * interaction)
{
  if( interaction->GetProcessInfo().IsQuasiElastic() ) 
  {
       return interaction_utils::QELRecoilNucleonPdgCode(interaction);
  } 
  else return 0;
}
//____________________________________________________________________________
int genie::interaction_utils::QELRecoilNucleonPdgCode(
                                              const Interaction * interaction)
{
  const InitialState & init_state = interaction->GetInitialState();

  //-- Determine the pdg code of the recoil nucleon
  int recoil_nuc_pdgc = 0;
  int struck_nuc_pdgc = init_state.GetTarget().StruckNucleonPDGCode();

  assert( pdg::IsProton(struck_nuc_pdgc) || pdg::IsNeutron(struck_nuc_pdgc) );
  assert( interaction->GetProcessInfo().IsWeak() );

  if ( interaction->GetProcessInfo().IsWeakCC() )
  {
    recoil_nuc_pdgc = pdg::SwitchProtonNeutron(struck_nuc_pdgc); // CC  
  }
  else
  {
    recoil_nuc_pdgc = struck_nuc_pdgc; // NC
  }

  LOG("Interaction", pDEBUG) << "Recoil nucleon PDG = " << recoil_nuc_pdgc;
  return recoil_nuc_pdgc;
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDis(int Z, int A, 
          int probe, const TLorentzVector & p4probe, InteractionType_t intype)
{
  Target       target(Z,A);
  InitialState init_state(target, probe);
  ProcessInfo  proc(kScDeepInelastic, intype);

  init_state.SetProbeP4(p4probe);
  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDis(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return interaction_utils::GetDis(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDisCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetDis(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDisCC(int Z, int A, int probe)
{
  return interaction_utils::GetDis(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDisNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetDis(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetDisNC(int Z, int A, int probe)
{
  return interaction_utils::GetDis(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQel(int Z, int A, 
          int probe, const TLorentzVector & p4probe, InteractionType_t intype)
{
  Target       target(Z,A);
  InitialState init_state(target, probe);
  ProcessInfo  proc(kScDeepInelastic, intype);

  init_state.SetProbeP4(p4probe);
  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQel(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return interaction_utils::GetQel(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQelCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetQel(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQelCC(int Z, int A, int probe)
{
  return interaction_utils::GetQel(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQelNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetQel(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetQelNC(int Z, int A, int probe)
{
  return interaction_utils::GetQel(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetRes(int Z, int A, 
          int probe, const TLorentzVector & p4probe, InteractionType_t intype)
{
  Target       target(Z,A);
  InitialState init_state(target, probe);
  ProcessInfo  proc(kScResonant, intype);

  init_state.SetProbeP4(p4probe);
  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetRes(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return interaction_utils::GetQel(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetResCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetQel(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetResCC(int Z, int A, int probe)
{
  return interaction_utils::GetQel(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetResNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return interaction_utils::GetQel(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetResNC(int Z, int A, int probe)
{
  return interaction_utils::GetQel(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetIMD(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  // IMD: Inverse Muon Decay

  Target       target(Z,A,0);
  InitialState init_state(target, probe);
  ProcessInfo  proc(kScInverseMuDecay, kIntWeakCC);

  init_state.SetProbeP4(p4probe);
  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________
Interaction * genie::interaction_utils::GetEl(int Z, int A, 
                                     int probe, const TLorentzVector & p4probe)
{
  Target target(Z,A);

  InitialState init_state(target, probe);

  init_state.SetProbeP4(p4probe);

  ProcessInfo proc(kScElastic, kIntWeakNC);

  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________






