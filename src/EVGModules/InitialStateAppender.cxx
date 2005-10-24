//____________________________________________________________________________
/*!

\class   genie::InitialStateAppender

\brief   Appends the initial state information to the event record.

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

*/
//____________________________________________________________________________

#include <TLorentzVector.h>

#include "EVGModules/InitialStateAppender.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
InitialStateAppender::InitialStateAppender() :
EventRecordVisitorI("genie::InitialStateAppender")
{

}
//___________________________________________________________________________
InitialStateAppender::InitialStateAppender(string config) :
EventRecordVisitorI("genie::InitialStateAppender", config)
{

}
//___________________________________________________________________________
InitialStateAppender::~InitialStateAppender()
{

}
//___________________________________________________________________________
void InitialStateAppender::ProcessEventRecord(GHepRecord * evrec) const
{
// Adds the initial state particles at the event record (the order is
// significant)

  LOG("ISApp", pINFO) << "Adding the initial state to the event record";

  //-- add the incoming neutrino to the event record
  this->AddNeutrino(evrec);

  //-- add the nuclear target at the event record (if any)
  this->AddNucleus(evrec);

  //-- add the struck nucleon to the event record (if any)
  //   It is added with status-code = 0 (init state) if the target was a
  //   free nucleon, or with a status-code = 11 (nucleon target) if the
  //   target was a nucleus.
  //   If the interaction was an IMD it will add the target e- instead.
  this->AddStruckParticle(evrec);
}
//___________________________________________________________________________
void InitialStateAppender::AddNeutrino(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfLab);
  const TLorentzVector v4(0.,0.,0.,0.);

  int pdgc = init_state.GetProbePDGCode();

  LOG("ISApp", pINFO) << "Adding neutrino [pdgc = " << pdgc << "]";

  evrec->AddParticle(pdgc,kIStInitialState, -1,-1,-1,-1, *p4, v4);

  delete p4;
}
//___________________________________________________________________________
void InitialStateAppender::AddNucleus(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  bool is_nucleus = init_state.GetTarget().IsNucleus();
  if(!is_nucleus) {
    LOG("ISApp", pINFO)
         << "Not an interaction with a nuclear target - no nucleus to add";
    return;
  }
  int    A    = init_state.GetTarget().A();
  int    Z    = init_state.GetTarget().Z();
  int    pdgc = pdg::IonPdgCode(A, Z);
  double M    = PDGLibrary::Instance()->Find(pdgc)->Mass();

  LOG("ISApp", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                              << ", pdg = " << pdgc << "]";

  evrec->AddParticle(pdgc,kIStInitialState,-1,-1,-1,-1, 0,0,0,M, 0,0,0,0);
}
//___________________________________________________________________________
void InitialStateAppender::AddStruckParticle(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo & proc_info   = interaction->GetProcessInfo();

  if(proc_info.IsInverseMuDecay()) {
    int    pdgc = kPdgElectron;
    double mass = PDGLibrary::Instance()->Find(pdgc)->Mass();
    const TLorentzVector p4(0,0,0, mass);
    const TLorentzVector v4(0.,0.,0.,0.);

    LOG("ISApp", pINFO) << "Adding struck electron";
    evrec->AddParticle(pdgc, kIStInitialState, 1, -1, -1, -1, p4, v4);
    return;
  }

  int pdgc = init_state.GetTarget().StruckNucleonPDGCode();

  if(pdgc != 0) {

    bool is_nucleus = init_state.GetTarget().IsNucleus();

    GHepStatus_t ist   = (is_nucleus) ? kIstNucleonTarget : kIStInitialState;
    int          imom1 = (is_nucleus) ? 1 : -1;
    int          imom2 = -1;

    const TLorentzVector p4(*init_state.GetTarget().StruckNucleonP4());
    const TLorentzVector v4(0.,0.,0.,0.);

    LOG("ISApp", pINFO)<< "Adding struck nucleon [pdgc = " << pdgc << "]";

    evrec->AddParticle(pdgc, ist, imom1, imom2, -1, -1, p4, v4);

  }//if struck nucleon was set
}
//___________________________________________________________________________
