//____________________________________________________________________________
/*!

\namespace  genie::utils::interaction

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
int genie::utils::interaction::RecoilNucleonPdgCode(
                                              const Interaction * interaction)
{
  if( interaction->GetProcessInfo().IsQuasiElastic() )
  {
       return utils::interaction::QELRecoilNucleonPdgCode(interaction);
  }
  else return 0;
}
//____________________________________________________________________________
int genie::utils::interaction::QELRecoilNucleonPdgCode(
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
Interaction * genie::utils::interaction::GetDis(int Z, int A,
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
Interaction * genie::utils::interaction::GetDis(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return utils::interaction::GetDis(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetDisCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetDis(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetDisCC(int Z, int A, int probe)
{
  return utils::interaction::GetDis(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetDisNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetDis(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetDisNC(int Z, int A, int probe)
{
  return utils::interaction::GetDis(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQel(int Z, int A,
          int probe, const TLorentzVector & p4probe, InteractionType_t intype)
{
  Target       target(Z,A);
  InitialState init_state(target, probe);
  ProcessInfo  proc(kScQuasiElastic, intype);

  init_state.SetProbeP4(p4probe);
  Interaction * interaction = new Interaction(init_state, proc);

  return interaction;
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQel(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return utils::interaction::GetQel(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQelCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetQel(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQelCC(int Z, int A, int probe)
{
  return utils::interaction::GetQel(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQelNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetQel(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetQelNC(int Z, int A, int probe)
{
  return utils::interaction::GetQel(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetRes(int Z, int A,
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
Interaction * genie::utils::interaction::GetRes(
                            int Z, int A, int probe, InteractionType_t intype)
{
  TLorentzVector p4(0,0,0,0); // null probe 4-momentum

  return utils::interaction::GetRes(Z,A,probe,p4,intype);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetResCC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetRes(Z,A,probe,p4probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetResCC(int Z, int A, int probe)
{
  return utils::interaction::GetRes(Z,A,probe,kIntWeakCC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetResNC(
                      int Z, int A, int probe, const TLorentzVector & p4probe)
{
  return utils::interaction::GetRes(Z,A,probe,p4probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetResNC(int Z, int A, int probe)
{
  return utils::interaction::GetRes(Z,A,probe,kIntWeakNC);
}
//____________________________________________________________________________
Interaction * genie::utils::interaction::GetIMD(
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
Interaction * genie::utils::interaction::GetEl(int Z, int A,
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






