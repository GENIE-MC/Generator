//____________________________________________________________________________
/*!

\class   genie::FermiMover

\brief   It visits the event record & computes a Fermi motion momentum for
         initial state nucleons bound in nuclei.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 08, 2004

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Conventions/Constants.h"
#include "EVGModules/FermiMover.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclMomentumGenerator.h"
#include "Nuclear/NuclMomentumModelI.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
FermiMover::FermiMover() :
EventRecordVisitorI("genie::FermiMover")
{

}
//___________________________________________________________________________
FermiMover::FermiMover(string config) :
EventRecordVisitorI("genie::FermiMover", config)
{

}
//___________________________________________________________________________
FermiMover::~FermiMover()
{

}
//___________________________________________________________________________
void FermiMover::ProcessEventRecord(GHepRecord * event_rec) const
{
  Interaction *  interaction = event_rec   -> GetInteraction();
  InitialState * init_state  = interaction -> GetInitialStatePtr();
  Target *       tgt         = init_state  -> GetTargetPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  TLorentzVector * p4 = tgt->StruckNucleonP4();

  // do nothing if the struct nucleon 4-momentum was set (eg as part of the
  // initial state selection)
  if(p4->Px()>0 || p4->Py()>0 || p4->Pz()>0) return;

  // generate a Fermi momentum
  assert(fNuclPModel);
  NuclMomentumGenerator * nucp_gen = NuclMomentumGenerator::Instance();
  nucp_gen->UseProbDistribution(fNuclPModel, *tgt);
  TVector3 p3 = nucp_gen->RandomMomentum3();
  LOG("FermiMover", pINFO) << "Generated nucleon momentum: ("
                  << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << ")";
  
  // access ths hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = event_rec->StruckNucleon();
  GHepParticle * nucleus = event_rec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  //-- compute A,Z for final state nucleus & get its PDG code and its mass
  int nucleon_pdgc = nucleon->PdgCode();
  bool is_p  = pdg::IsProton(nucleon_pdgc);
  int Z = (is_p) ? nucleus->Z()-1 : nucleus->Z();
  int A = nucleus->A() - 1;

  TParticlePDG * pdgp = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  pdgp = PDGLibrary::Instance()->Find(ipdgc);
  if(!pdgp) {
      LOG("FermiMover", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
      exit(1);
  }
  double m_remnant_nucleus = pdgp->Mass();

  //-- compute the mass difference between the two nuclei

  double m_original_nucleus = nucleus->Mass();
  double dm = m_original_nucleus - m_remnant_nucleus;

  //-- update the struck nucleon 4p at the interaction summary and at
  //   the GHEP record

  p4->SetPx( p3.Px() );
  p4->SetPy( p3.Py() );
  p4->SetPz( p3.Pz() );
  p4->SetE ( dm      ); // off the mass-shell

  nucleon->SetMomentum(*p4);

  // sometimes, for interactions near threshold, Fermi momentum might bring
  // the neutrino energy in the nucleon rest frame below threshold (for the
  // selected interaction). In this case mark the event as unphysical and
  // abort the current thread.
  double E    = init_state->GetProbeE(kRfStruckNucAtRest);
  double Ethr = utils::kinematics::EnergyThreshold(interaction);
  if(E<=Ethr) {
     LOG("FermiMover", pNOTICE)
                  << "Event below threshold after generation Fermi momentum";
     event_rec->EventFlags()->SetBitNumber(kBelowThrNucleonRestFrame, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E < Ethr after generating nucleon Fermi momentum");
     exception.SwitchOnFastForward();
     throw exception;
  }
}
//___________________________________________________________________________
void FermiMover::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::LoadConfig(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed

  fNuclPModel = dynamic_cast<const NuclMomentumModelI *>
         (this->SubAlg("nucl-p-distribution-alg","nucl-p-distribution-conf"));
}
//____________________________________________________________________________

