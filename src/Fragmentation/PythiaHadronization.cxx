//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "Algorithm/AlgConfigPool.h"
#include "Fragmentation/PythiaHadronization.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
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
  this->Initialize();
}
//____________________________________________________________________________
PythiaHadronization::PythiaHadronization(string config) :
HadronizationModelI("genie::PythiaHadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaHadronization::~PythiaHadronization()
{

}
//____________________________________________________________________________
void PythiaHadronization::Initialize(void) const
{
  fPythia = TPythia6::Instance();
}
//____________________________________________________________________________
TClonesArray * PythiaHadronization::Hadronize(
                                        const Interaction * interaction) const
{
  LOG("PythiaHad", pNOTICE) << "Running PYTHIA hadronizer";

  this->SyncSeeds();

  //-- get kinematics / init-state / process-info

  const Kinematics &   kinematics = interaction->GetKinematics();
  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();
  const Target &       target     = init_state.GetTarget();

  assert(target.StruckQuarkIsSet()); 

  double W = kinematics.W();

  int  probe       = init_state.GetProbePDGCode();
  int  hit_nucleon = target.StruckNucleonPDGCode();
  int  hit_quark   = target.StruckQuarkPDGCode();
  bool from_sea    = target.StruckQuarkIsFromSea();

  LOG("PythiaHad", pNOTICE)
          << "Hit nucleon pdgc = " << hit_nucleon << ", W = " << W;
  LOG("PythiaHad", pNOTICE)
            << "Selected hit quark pdgc = " << hit_quark
                           << ((from_sea) ? "[sea]" : "[valence]");

  //-- check hit-nucleon assignment, input neutrino & weak current
  bool isp  = pdg::IsProton(hit_nucleon);
  bool isn  = pdg::IsNeutron(hit_nucleon);
  bool isv  = pdg::IsNeutrino(probe);
  bool isvb = pdg::IsAntiNeutrino(probe);
  bool iscc = proc_info.IsWeakCC();
  bool isnc = proc_info.IsWeakNC();
  if( !(isp||isn) ) {
    LOG("PythiaHad", pERROR) << "Can not handle nucleon: " << hit_nucleon;
    exit(1);
  }
  if( !(iscc||isnc) ) {
    LOG("PythiaHad", pERROR) << "Can only handle weak interactions";
    exit(1);
  }
  if( !(isv||isvb) ) {
    LOG("PythiaHad", pERROR)
                      << "Can not handle non-neutrino probe: " << probe;
    exit(1);
  }

  //-- assert that the interaction mode is allowed
  bool isu  = pdg::IsUQuark     (hit_quark);
  bool isd  = pdg::IsDQuark     (hit_quark);
  bool isub = pdg::IsUAntiQuark (hit_quark);
  bool isdb = pdg::IsDAntiQuark (hit_quark);

  bool allowed = (iscc && isv  && (isd||isub)) ||
                 (iscc && isvb && (isu||isdb)) ||
                 (isnc && (isv||isvb) && (isu||isd||isub||isdb));
  assert(allowed);

  //-- generate the quark system (q + qq) initiating the hadronization

  int  final_quark = 0; // leading quark (hit quark after the interaction)
  int  diquark     = 0; // remnant diquark (xF<0 at hadronic CMS)

  //-- figure out the what happens to the hit quark after the interaction

  if (proc_info.IsWeakNC()) final_quark = hit_quark;
  else {
    if      (isv  && isd ) final_quark = kPdgUQuark;
    else if (isv  && isub) final_quark = kPdgDQuarkBar;
    else if (isvb && isu ) final_quark = kPdgDQuark;
    else if (isvb && isdb) final_quark = kPdgUQuarkBar;
    else {
      LOG("PythiaHad", pERROR)
        << "Not allowed mode. Refused to make a final quark assignment!";
      exit(1);
    }
  }//CC

  //-- figure out what the remnant diquark is

  //   Note from Hugh, following a conversation with his local HEP theorist 
  //   (Gary Goldstein): "I am told that the probability of finding the diquark 
  //   in the singlet vs. triplet states is 50-50."  

  // hit quark = valence quark
  if(!from_sea) {
    if (isp && isu) diquark = kPdgUDDiquarkS1; /* u(->q) + ud */
    if (isp && isd) diquark = kPdgUUDiquarkS1; /* d(->q) + uu */
    if (isn && isu) diquark = kPdgDDDiquarkS1; /* u(->q) + dd */
    if (isn && isd) diquark = kPdgUDDiquarkS1; /* d(->q) + ud */
  }
  // hit quark = sea quark
  else {

    if(isp && isu) diquark = kPdgUDDiquarkS1; /* u(->q) + bar{u} uud (=ud) */
    if(isp && isd) diquark = kPdgUUDiquarkS1; /* d(->q) + bar{d} uud (=uu) */
    if(isn && isu) diquark = kPdgDDDiquarkS1; /* u(->q) + bar{u} udd (=dd) */
    if(isn && isd) diquark = kPdgUDDiquarkS1; /* d(->q) + bar{d} udd (=ud) */

    // This section needs revisiting.
    // The neutrino is scattered off a sea antiquark, materializing its quark
    // partner and leaving me with a 5q system ( <qbar + q> + qqq(valence) )
    // I will force few qbar+q annhilations below to get my quark/diquark system
    // Probably it is best to leave the qqq system in the final state and then
    // just do the fragmentation of the qbar q system? But how do I figure out
    // how to split the available energy?

    /* bar{u} (-> bar{d}) + u uud => u + uu */
    if(isp && isub && iscc) {final_quark = kPdgUQuark; diquark = kPdgUUDiquarkS1;}
    /* bar{u} (-> bar{u}) + u uud => u + ud */
    if(isp && isub && isnc) {final_quark = kPdgUQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{d} (-> bar{u}) + d uud => d + ud */
    if(isp && isdb && iscc) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{d} (-> bar{d}) + d uud => d + uu */
    if(isp && isdb && isnc) {final_quark = kPdgDQuark; diquark = kPdgUUDiquarkS1;}
    /* bar{u} (-> bar{d}) + u udd => u + ud */
    if(isn && isub && iscc) {final_quark = kPdgUQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{u} (-> bar{u}) + u udd => u + dd */
    if(isn && isub && isnc) {final_quark = kPdgUQuark; diquark = kPdgDDDiquarkS1;}
    /* bar{d} (-> bar{u}) + d udd => d + dd */
    if(isn && isdb && iscc) {final_quark = kPdgDQuark; diquark = kPdgDDDiquarkS1;}
    /* bar{d} (-> bar{d}) + d udd => d + ud */
    if(isn && isdb && isnc) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}

    // if the diquark is a ud, switch it to the singlet state with 50% probability
    if(diquark == kPdgUDDiquarkS1) {
      RandomGen * rnd = RandomGen::Instance();
      double Rqq = rnd->RndHadro().Rndm();
      if(Rqq<0.5) diquark = kPdgUDDiquarkS0;
    }
  }
  assert(diquark!=0);

  //-- PYTHIA->HADRONIZE:

  LOG("PythiaHad", pNOTICE)
        << "Fragmentation / Init System: "
                      << "q = " << final_quark << ", qq = " << diquark;
  int ip = 0;
  py2ent_(&ip, &final_quark, &diquark, &W);

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
double PythiaHadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
void PythiaHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaHadronization::LoadConfig(void)
{
  // the configurable PYTHIA parameters used here are the ones used by NUX 
  // (see A.Rubbia's talk @ NuINT-01)
  // The defaults are the values used by PYTHIA
  // Use the NUX config set to set the tuned values as used in NUX.

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fSSBarSuppression = fConfig->GetDoubleDef(
              "ssbar-suppression", gc->GetDouble("PYTHIA-SSBarSuppression"));
  fGaussianPt2 = fConfig->GetDoubleDef(
                        "gaussian-pt2", gc->GetDouble("PYTHIA-GaussianPt2"));
  fNonGaussianPt2Tail = fConfig->GetDoubleDef(
        "non-gaussian-pt2-tail", gc->GetDouble("PYTHIA-NonGaussianPt2Tail"));
  fRemainingECutoff = fConfig->GetDoubleDef(
   "remaining-energy-cutoff", gc->GetDouble("PYTHIA-RemainingEnergyCutoff"));

  fPythia->SetPARJ(2,  fSSBarSuppression);
  fPythia->SetPARJ(21, fGaussianPt2);
  fPythia->SetPARJ(23, fNonGaussianPt2Tail);
  fPythia->SetPARJ(33, fRemainingECutoff);
}
//____________________________________________________________________________
void PythiaHadronization::SyncSeeds(void) const
{
// Keep PYTHIA6 random number seed in sync with GENIE's random number seed
//
  long int cs = RandomGen::Instance()->GetSeed();
  if(fCurrSeed != cs) {
     fCurrSeed = cs;
     fPythia->SetMRPY(1,fCurrSeed);
  }
}
//____________________________________________________________________________


