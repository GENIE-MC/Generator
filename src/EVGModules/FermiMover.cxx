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

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Conventions/Constants.h"
#include "EVGModules/FermiMover.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearPDistribution.h"
#include "Nuclear/NuclearPDistributionModelI.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
FermiMover::FermiMover() :
EventRecordVisitorI()
{
  fName = "genie::FermiMover";
}
//___________________________________________________________________________
FermiMover::FermiMover(const char * param_set) :
EventRecordVisitorI(param_set)
{
  fName = "genie::FermiMover";

  this->FindConfig();
}
//___________________________________________________________________________
FermiMover::~FermiMover()
{

}
//___________________________________________________________________________
void FermiMover::ProcessEventRecord(GHepRecord * event_rec) const
{
  Interaction * interaction = event_rec->GetInteraction();

  const InitialState & init_state = interaction->GetInitialState();

  if(init_state.GetTarget().IsNucleus()) {

     const NuclearPDistributionModelI * nucl_p_model =
       dynamic_cast<const NuclearPDistributionModelI *> (this->SubAlg(
                       "nucl-p-distribution-alg","nucl-p-distribution-conf"));

     NuclearPDistribution nucleon_momentum_generator;

     nucleon_momentum_generator.AttachModel(nucl_p_model);
     nucleon_momentum_generator.BuildProbDistribution(init_state.GetTarget());

     TVector3 p3 = nucleon_momentum_generator.RandomMomentum3();

     LOG("FermiMover", pINFO) << "Generated nucleon momentum: ("
                      << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << ")";

     //-- update the 4-vector in the Interaction object

     TLorentzVector * nucleon_p4 = init_state.GetTarget().StruckNucleonP4();

     double M = init_state.GetTarget().StruckNucleonMass();

     nucleon_p4->SetPx( p3.Px() );
     nucleon_p4->SetPy( p3.Py() );
     nucleon_p4->SetPz( p3.Pz() );
     nucleon_p4->SetE ( sqrt(p3.Mag2() + M*M ) ); // should be kept on-shell?

     //-- update the event record (particle with Ist == kIstNucleonTarget)

     int nucleon_pdgc = init_state.GetTarget().StruckNucleonPDGCode();

     GHepParticle * nucleon =
                event_rec->FindParticle(nucleon_pdgc, kIstNucleonTarget, 0);

     nucleon->SetMomentum(*nucleon_p4);
  }
}
//___________________________________________________________________________
