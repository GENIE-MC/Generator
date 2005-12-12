//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models.

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "Fragmentation/PythiaHadronization.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDF/PDF.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;

//-- the actual PYTHIA call

extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
PythiaHadronization::PythiaHadronization() :
HadronizationModelI("genie::PythiaHadronization")
{
  fPythia = new TPythia6();
}
//____________________________________________________________________________
PythiaHadronization::PythiaHadronization(string config) :
HadronizationModelI("genie::PythiaHadronization", config)
{
  fPythia = new TPythia6();
}
//____________________________________________________________________________
PythiaHadronization::~PythiaHadronization()
{

}
//____________________________________________________________________________
void PythiaHadronization::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * PythiaHadronization::Hadronize(
                                        const Interaction * interaction) const
{
  LOG("PythiaHad", pINFO) << "Running PYTHIA hadronizer";

  const Kinematics & kinematics = interaction->GetKinematics();
  double W = kinematics.W();

  const InitialState & init_state = interaction->GetInitialState();

  int hit_nucleon = init_state.GetTarget().StruckNucleonPDGCode();
  int hit_quark   = init_state.GetTarget().StruckQuarkPDGCode();
  int diquark     = 0;

  LOG("PythiaHad", pINFO)
         << "Struck quark pdgc = " << hit_quark
                   << ", nucleon pdgc = " << hit_nucleon << " at W = " << W;

  bool isp = pdg::IsProton(hit_nucleon);
  bool isn = pdg::IsNeutron(hit_nucleon);
  if( !(isp||isn) ) {
    LOG("PythiaHad", pERROR) << "Can not handle nucleon: " << hit_nucleon;
    exit(1);
  }

  //-- check for hit-quark assignment

  if(hit_quark==0) {

     LOG("PythiaHad", pINFO)
          << "No input hit quark. Selecting one based on specified PDF set";

     const PDFModelI * pdf_model =
                  dynamic_cast<const PDFModelI *>
                             (this->SubAlg("pdf-alg-name", "pdf-param-set"));

     double x  = kinematics.x();
     double Q2 = utils::kinematics::CalcQ2(interaction);

     PDF pdf;
     pdf.SetModel(pdf_model);
     pdf.Calculate(x,Q2);

     double xu  = pdf.UpValence();
     double xd  = pdf.DownValence();
     double xud = xu + xd;

     RandomGen * rnd = RandomGen::Instance();
     double t = xud * (rnd->Random1().Rndm()); // e [0, xu+xd]

     hit_quark = kPdgUQuark;
     if( (isp && t > xu/xud) || (isn && t > xd/xud) ) hit_quark = kPdgDQuark;

     LOG("PythiaHad", pINFO) << "Selected hit quark pdgc = " << hit_quark;
  }

  bool isu = pdg::IsUQuark(hit_quark);
  bool isd = pdg::IsDQuark(hit_quark);
  if( !(isu||isd) ) {
    LOG("PythiaHad", pERROR) << "Can not handle quark: " << hit_quark;
    exit(1);
  }

  //-- select diquark
  if (isp && isu) diquark = kPdgUDDiquarkS1; /* u + ud */
  if (isp && isd) diquark = kPdgUUDiquarkS1; /* d + uu */
  if (isn && isu) diquark = kPdgDDDiquarkS1; /* u + dd */
  if (isn && isd) diquark = kPdgUDDiquarkS1; /* d + ud */

  //-- PYTHIA->HADRONIZE:
  int ip = 0;
  py2ent_(&ip, &hit_quark, &diquark, &W);

  //-- get LUJETS record
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles =
                      (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("TMCParticle", np);
  register unsigned int i = 0;
  TMCParticle * particle = 0;
  TIter particle_iter(pythia_particles);

  while( (particle = (TMCParticle *) particle_iter.Next()) ) {
       LOG("PythiaHad", pINFO)
               << "Adding final state particle pdgc = " << particle->GetKF();
       new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  particle_list->SetOwner(true);
  return particle_list;
}
//____________________________________________________________________________
