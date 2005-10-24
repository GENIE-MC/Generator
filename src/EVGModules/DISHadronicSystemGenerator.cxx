//____________________________________________________________________________
/*!

\class   genie::DISHadronicSystemGenerator

\brief   Generates the final state hadronic system in v DIS interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "Conventions/Constants.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISHadronicSystemGenerator.h"
#include "Fragmentation/HadronizationModelI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepOrder.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::print;

//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator() :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator")
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::~DISHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void DISHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the recoil nucleon
  //   Its 4-momentum is computed by requiring the energy + momentum to be
  //   conserved.
  this->AddFragmentationProducts(evrec);
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::AddFragmentationProducts(
                                                    GHepRecord * evrec) const
{
  //-- Get the requested hadronization model
  const HadronizationModelI * hadr_model =
            dynamic_cast<const HadronizationModelI *> (this->SubAlg(
                       "hadronization-alg-name", "hadronization-param-set"));

  //-- Compute the hadronic system invariant mass
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E     = p4->Energy();
  double x     = interaction->GetScatteringParams().x();
  double y     = interaction->GetScatteringParams().y();
  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = Mnuc * Mnuc;

  double W2 = TMath::Max(0., Mnuc2 + 2*E*Mnuc*y*(1-x));
  double W  = TMath::Sqrt(W2);

  interaction->GetScatParamsPtr()->Set("W", W);
  delete p4;

  //-- Run the hadronization model and get the fragmentation products:
  //   A collection of ROOT TMCParticles (equivalent to a LUJETS record)
  TClonesArray * particle_list = hadr_model->Hadronize(interaction);

  if(!particle_list) {
     LOG("DISHadronicVtx", pWARN) 
                    << "Got an empty particle list. Hadronizer failed!";
     LOG("DISHadronicVtx", pWARN) 
                      << "Quitting the current event generation thread";

     evrec->SwitchGenericErrFlag(true);

     genie::exceptions::EVGThreadException exception;
     exception.SetReason(
                "Unphysical Event [Not enough phase space for hadronizer]");
     exception.SwitchOnFastForward();
     throw exception;

     return;
  }

  //-- Velocity for the [Hadronic CM] -> [LAB] active Lorentz transform
  TVector3 beta = this->HCM2LAB(evrec);

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.

  int mom = GHepOrder::StruckNucleonPosition(interaction);

  TMCParticle * p = 0;
  TIter particle_iter(particle_list);

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

     // the fragmentation products are generated in the final state
     // hadronic CM Frame - take each particle back to the LAB frame
     TLorentzVector p4(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());
     p4.Boost(beta);

     // copy the particle to the event record
     int          pdgc   = p->GetKF();
     GHepStatus_t status = GHepStatus_t(p->GetKS());
     TLorentzVector v4(0,0,0,0); // dummy position 4-vector

     evrec->AddParticle(pdgc, status, mom,-1,-1,-1, p4, v4);

  } // fragmentation-products-iterator

  particle_list->Delete();

  delete particle_list;
}
//___________________________________________________________________________
TVector3 DISHadronicSystemGenerator::HCM2LAB(GHepRecord * evrec) const
{
// Velocity for the Hadronic CM -> LAB active Lorentz transform

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  //incoming v:
  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfLab);
  assert(nu_p4);

  //struck nucleon:
  TLorentzVector * nucl_p4 = init_state.GetTarget().StruckNucleonP4();
  assert(nucl_p4);

  //final state primary lepton:
  int fsl_pdgc = interaction->GetFSPrimaryLepton()->PdgCode();
  GHepParticle * fsl = evrec->FindParticle(
                                     fsl_pdgc, kIStStableFinalState, 0);
  assert(fsl);

  TLorentzVector * fsl_p4 = fsl->GetP4();

  LOG("DISHadronicVtx", pINFO)
                 << "\n v [LAB]: " << P4AsString( nu_p4   )
                 << "\n N [LAB]: " << P4AsString( nucl_p4 )
                 << "\n l [LAB]: " << P4AsString( fsl_p4  );

  //-- Compute the velocity of the LAB frame in the Final State Hadronic
  //   CM Frame (Pv + Pnuc = Pfsl + Sum{Phad_{i}})

  double PX = nu_p4->Px()     + nucl_p4->Px()     - fsl_p4->Px();
  double PY = nu_p4->Py()     + nucl_p4->Py()     - fsl_p4->Py();
  double PZ = nu_p4->Pz()     + nucl_p4->Pz()     - fsl_p4->Pz();
  double E  = nu_p4->Energy() + nucl_p4->Energy() - fsl_p4->Energy();

  delete fsl_p4;
  delete nu_p4;

  assert(E>0);
  TVector3 beta( PX/E, PY/E, PZ/E );

  LOG("DISHadronicVtx", pINFO)
                   << "\n beta (HCM --> LAB): " << Vec3AsString(&beta);
  return beta;
}
//___________________________________________________________________________
