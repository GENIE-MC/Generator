//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TH1D.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h" 
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Hadronization/Pythia6Hadronization.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Hadronization/FragmRecUtils.h"


using namespace genie;
using namespace genie::constants;

// the actual PYTHIA call
extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
Pythia6Hadronization::Pythia6Hadronization() :
EventRecordVisitorI("genie::Pythia6Hadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadronization::Pythia6Hadronization(string config) :
EventRecordVisitorI("genie::Pythia6Hadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadronization::~Pythia6Hadronization()
{

}
//____________________________________________________________________________
void Pythia6Hadronization::Initialize(void) const
{
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
}
//____________________________________________________________________________
void Pythia6Hadronization::ProcessEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  TClonesArray * particle_list = this->Hadronize(interaction);

  if(! particle_list ) {
    LOG("Pythia6Hadronization", pWARN) << "Got an empty particle list. Hadronizer failed!";
    LOG("Pythia6Hadronization", pWARN) << "Quitting the current event generation thread";
    
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    
    return;
   }

  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = is_nucleus ? kIStHadronInTheNucleus : kIStStableFinalState ;
  

  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);


  GHepParticle * neutrino  = event->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  GHepParticle * p = 0;
  TIter particle_iter(particle_list);
  while ((p = (GHepParticle *) particle_iter.Next())) 
  {

    int pdgc = p -> Pdg() ;
        
    // set the proper status according to a number of things:
    // interaction on a nucleaus or nucleon, particle type
    GHepStatus_t ist = ( p -> Status() ==1 ) ? istfin : kIStDISPreFragmHadronicState;
    // not that this status here ^ is the pythia state

    // handle gammas, and leptons that might come from internal pythia decays
    // mark them as final state particles
    bool not_hadron = ( pdgc == kPdgGamma ||
			pdg::IsNeutralLepton(pdgc) ||
			pdg::IsChargedLepton(pdgc) ) ;

    if( not_hadron )  { ist = kIStStableFinalState; }
    p -> SetStatus( ist ) ;
 
    p->SetFirstMother(mom + p->FirstMother() );
    p->SetLastMother( -1 );
    
    // In Pythia6 having no daughters means daughter == 0 hnce the following check
    int ifd = (p->FirstDaughter() <= 0 ) ? -1 : mom  + p->FirstDaughter();
    int ild = (p->LastDaughter()  <= 0 ) ? -1 : mom  + p->LastDaughter();
    p->SetFirstDaughter(ifd);
    p->SetLastDaughter (ild);
    
    // the Pythia particle position is overridden
    p -> SetPosition( vtx ) ; 

    event->AddParticle(*p);
  }

  particle_list -> Delete() ; 
  delete particle_list ;
}
//____________________________________________________________________________
TClonesArray * 
  Pythia6Hadronization::Hadronize(
         const Interaction * interaction) const
{
  LOG("Pythia6Had", pNOTICE) << "Running PYTHIA hadronizer";

  if(!this->AssertValidity(interaction)) {
     LOG("Pythia6Had", pERROR) << "Returning a null particle list!";
     return 0;
  }

  // get kinematics / init-state / process-info

  const Kinematics &   kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  assert(target.HitQrkIsSet());

  double W = kinematics.W();

  int  probe       = init_state.ProbePdg();
  int  hit_nucleon = target.HitNucPdg();
  int  hit_quark   = target.HitQrkPdg();
  bool from_sea    = target.HitSeaQrk();

  LOG("Pythia6Had", pNOTICE)
          << "Hit nucleon pdgc = " << hit_nucleon << ", W = " << W;
  LOG("Pythia6Had", pNOTICE)
            << "Selected hit quark pdgc = " << hit_quark
                           << ((from_sea) ? "[sea]" : "[valence]");

  // check hit-nucleon assignment, input neutrino & interaction type
  bool isp  = pdg::IsProton           (hit_nucleon);
  bool isn  = pdg::IsNeutron          (hit_nucleon);
  bool isv  = pdg::IsNeutrino         (probe);
  bool isvb = pdg::IsAntiNeutrino     (probe);
//bool isl  = pdg::IsNegChargedLepton (probe);
//bool islb = pdg::IsPosChargedLepton (probe);
  bool iscc = proc_info.IsWeakCC      ();
  bool isnc = proc_info.IsWeakNC      ();
  bool isdm = proc_info.IsDarkMatter  ();
  bool isem = proc_info.IsEM          ();
  bool isu  = pdg::IsUQuark           (hit_quark);
  bool isd  = pdg::IsDQuark           (hit_quark);
  bool iss  = pdg::IsSQuark           (hit_quark);
  bool isub = pdg::IsAntiUQuark       (hit_quark);
  bool isdb = pdg::IsAntiDQuark       (hit_quark);
  bool issb = pdg::IsAntiSQuark       (hit_quark);

  //
  // Generate the quark system (q + qq) initiating the hadronization
  //

  int  final_quark = 0; // leading quark (hit quark after the interaction)
  int  diquark     = 0; // remnant diquark (xF<0 at hadronic CMS)

  // Figure out the what happens to the hit quark after the interaction
  if (isnc || isem || isdm) {
    // NC, EM
    final_quark = hit_quark;
  } else {
    // CC
    if      (isv  && isd ) final_quark = kPdgUQuark;
    else if (isv  && iss ) final_quark = kPdgUQuark;
    else if (isv  && isub) final_quark = kPdgAntiDQuark;
    else if (isvb && isu ) final_quark = kPdgDQuark;
    else if (isvb && isdb) final_quark = kPdgAntiUQuark;
    else if (isvb && issb) final_quark = kPdgAntiUQuark;
    else {
      LOG("Pythia6Had", pERROR)
        << "Not allowed mode. Refused to make a final quark assignment!";
      return 0;
    }
  }//CC

  // Figure out what the remnant diquark is.
  // Note from Hugh, following a conversation with his local HEP theorist
  // (Gary Goldstein): "I am told that the probability of finding the diquark
  // in the singlet vs. triplet states is 50-50."

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

    // The following section needs revisiting.

    // The lepton is scattered off a sea antiquark, materializing its quark
    // partner and leaving me with a 5q system ( <qbar + q> + qqq(valence) )
    // I will force few qbar+q annhilations below to get my quark/diquark system
    // Probably it is best to leave the qqq system in the final state and then
    // just do the fragmentation of the qbar q system? But how do I figure out
    // how to split the available energy?

    /* bar{u} (-> bar{d}) + u uud => u + uu */
    if(isp && isub && iscc)         {final_quark = kPdgUQuark; diquark = kPdgUUDiquarkS1;}
    /* bar{u} (-> bar{u}) + u uud => u + ud */
    if(isp && isub && (isnc||isem||isdm)) {final_quark = kPdgUQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{d} (-> bar{u}) + d uud => d + ud */
    if(isp && isdb && iscc)         {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{d} (-> bar{d}) + d uud => d + uu */
    if(isp && isdb && (isnc||isem||isdm)) {final_quark = kPdgDQuark; diquark = kPdgUUDiquarkS1;}
    /* bar{u} (-> bar{d}) + u udd => u + ud */
    if(isn && isub && iscc)         {final_quark = kPdgUQuark; diquark = kPdgUDDiquarkS1;}
    /* bar{u} (-> bar{u}) + u udd => u + dd */
    if(isn && isub && (isnc||isem||isdm)) {final_quark = kPdgUQuark; diquark = kPdgDDDiquarkS1;}
    /* bar{d} (-> bar{u}) + d udd => d + dd */
    if(isn && isdb && iscc)         {final_quark = kPdgDQuark; diquark = kPdgDDDiquarkS1;}
    /* bar{d} (-> bar{d}) + d udd => d + ud */
    if(isn && isdb && (isnc||isem||isdm)) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}

    // The neutrino is scatterred off s or sbar sea quarks
    // For the time being I will handle s like d and sbar like dbar (copy & paste
    // from above) so that I conserve charge.

    if(iss || issb) {
       LOG("Pythia6Had", pNOTICE) 
                 << "Can not really handle a hit s or sbar quark / Faking it";

       if(isp && iss) { diquark = kPdgUUDiquarkS1; }
       if(isn && iss) { diquark = kPdgUDDiquarkS1; }

       if(isp && issb && iscc)         {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
       if(isp && issb && (isnc||isem||isdm)) {final_quark = kPdgDQuark; diquark = kPdgUUDiquarkS1;}
       if(isn && issb && iscc)         {final_quark = kPdgDQuark; diquark = kPdgDDDiquarkS1;}
       if(isn && issb && (isnc||isem||isdm)) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
    }

    // if the diquark is a ud, switch it to the singlet state with 50% probability
    if(diquark == kPdgUDDiquarkS1) {
      RandomGen * rnd = RandomGen::Instance();
      double Rqq = rnd->RndHadro().Rndm();
      if(Rqq<0.5) diquark = kPdgUDDiquarkS0;
    }
  }
  assert(diquark!=0);

  //
  // PYTHIA -> HADRONIZATION
  //

  LOG("Pythia6Had", pNOTICE)
        << "Fragmentation / Init System: "
        << "q = " << final_quark << ", qq = " << diquark;
  int ip = 0;

  // Determine how jetset treats un-stable particles appearing in hadronization

  int pi0_decflag = fPythia->GetMDCY(fPythia->Pycomp(kPdgPi0),              1);
  int K0_decflag  = fPythia->GetMDCY(fPythia->Pycomp(kPdgK0),               1);
  int K0b_decflag = fPythia->GetMDCY(fPythia->Pycomp(kPdgAntiK0),           1);
  int L0_decflag  = fPythia->GetMDCY(fPythia->Pycomp(kPdgLambda),           1);
  int L0b_decflag = fPythia->GetMDCY(fPythia->Pycomp(kPdgAntiLambda),       1);
  int Dm_decflag  = fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),  1);
  int D0_decflag  = fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),  1);
  int Dp_decflag  = fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),  1);
  int Dpp_decflag = fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP), 1);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for pi0           =  " << pi0_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for K0            =  " << K0_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for \bar{K0}      =  " << K0b_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for Lambda        =  " << L0_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for \bar{Lambda0} =  " << L0b_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for D-            =  " << Dm_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for D0            =  " << D0_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for D+            =  " << Dp_decflag;
  LOG("Pythia6Had", pDEBUG) << "Original decay flag for D++           =  " << Dpp_decflag;
#endif

  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),               1,0); // don't decay pi0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),                1,0); // don't decay K0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),            1,0); // don't decay \bar{K0}
  fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),            1,0); // don't decay Lambda0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda),        1,0); // don't decay \bar{Lambda0}
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),   1,1); // decay Delta-
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),   1,1); // decay Delta0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),   1,1); // decay Delta+
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP),  1,1); // decay Delta++

  // -- hadronize --
  py2ent_(&ip, &final_quark, &diquark, &W); // hadronizer

  // restore pythia decay settings so as not to interfere with decayer
  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),              1, pi0_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),               1, K0_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),           1, K0b_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),           1, L0_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda),       1, L0b_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),  1, Dm_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),  1, D0_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),  1, Dp_decflag);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP), 1, Dpp_decflag);

  // get LUJETS record
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles =
       (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
  particle_list->SetOwner(true);

  unsigned int i = 0;
  TMCParticle * p = 0;
  TIter particle_iter(pythia_particles);

  // Hadronic 4vec
  TLorentzVector p4Had = kinematics.HadSystP4();

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  while( (p = (TMCParticle *) particle_iter.Next()) ) {
     // The fragmentation products are generated in the hadronic CM frame
     // where the z>0 axis is the \vec{phad} direction. For each particle 
     // returned by the hadronizer:
     // - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
     // - rotate its 3-momentum from LAB' to LAB
     TLorentzVector p4o(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());
     p4o.Boost(beta); 
     TVector3 p3 = p4o.Vect();
     p3.RotateUz(unitvq); 
     TLorentzVector p4(p3,p4o.Energy());

     // Convert from TMCParticle to GHepParticle
     GHepParticle particle = GHepParticle(
         p->GetKF(),                // pdg
         GHepStatus_t(p->GetKS()),  // status
         p->GetParent(),            // first parent
         -1,                        // second parent
         p->GetFirstChild(),        // first daughter
         p->GetLastChild(),         // second daughter
         p4.Px(),                   // px
         p4.Py(),                   // py
         p4.Pz(),                   // pz
         p4.Energy(),               // e
         p->GetVx(),                // x
         p->GetVy(),                // y
         p->GetVz(),                // z
         p->GetTime()               // t
     );

     LOG("Pythia6Had", pDEBUG)
          << "Adding final state particle pdgc = " << particle.Pdg() 
          << " with status = " << particle.Status();

     if(particle.Status() == 1) {
       if( pdg::IsQuark  (particle.Pdg()) || 
	   pdg::IsDiQuark(particle.Pdg()) ) {
	 LOG("Pythia6Had", pERROR)
	   << "Hadronization failed! Bare quark/di-quarks appear in final state!";
	 
	 particle_list->Delete();
	 delete particle_list;
	 return 0;
       }
     }
     
     // insert the particle in the list
     new ( (*particle_list)[i++] ) GHepParticle(particle);
  }

  utils::fragmrec::Print(particle_list);

  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * 
   Pythia6Hadronization::SelectParticles(
            const Interaction * interaction) const
{
// Works the opposite way (compared with the KNO hadronization model)
// Rather than having this method as one of the hadronization model components,
// we extract the list of particles from the fragmentation record after the
// hadronization has been completed.

  TClonesArray * particle_list = this->Hadronize(interaction);

  if(!particle_list) return 0;

  bool allowdup=true;
  PDGCodeList * pdgcv = new PDGCodeList(allowdup);
  pdgcv->reserve(particle_list->GetEntries());

  GHepParticle * particle = 0;
  TIter particle_iter(particle_list);

  while ((particle = (GHepParticle *) particle_iter.Next())) 
  {
    if (particle->Status()==kIStStableFinalState) pdgcv->push_back(particle->Pdg());
  }
  particle_list->Delete();
  delete particle_list;

  return pdgcv;
}
//____________________________________________________________________________
TH1D * Pythia6Hadronization::MultiplicityProb( const Interaction * interaction ) const
{
// Similar comments apply as in SelectParticles()

  if(!this->AssertValidity(interaction)) {
     LOG("Pythia6Had", pWARN) 
                << "Returning a null multipicity probability distribution!";
     return 0;
  }
  double maxmult   = this->MaxMult(interaction);
  TH1D * mult_prob = this->CreateMultProbHist(maxmult);

  const int nev=500;
  GHepParticle * particle = 0;

  for(int iev=0; iev<nev; iev++) {

     TClonesArray * particle_list = this->Hadronize(interaction);
     double         weight        = this->Weight();

     if(!particle_list) { iev--; continue; }

     int n = 0;
     TIter particle_iter(particle_list);
     while ((particle = (GHepParticle *) particle_iter.Next())) 
     {
       if (particle->Status()==kIStStableFinalState) n++;
     }   
     particle_list->Delete();
     delete particle_list;
     mult_prob->Fill( (double)n, weight);
  }

  double integral = mult_prob->Integral("width");
  if(integral>0) {
    // Normalize the probability distribution
    mult_prob->Scale(1.0/integral);
  } else {
    SLOG("Pythia6Had", pWARN) << "probability distribution integral = 0";
    return mult_prob;
  }

  return mult_prob;
}
//____________________________________________________________________________
double Pythia6Hadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
bool Pythia6Hadronization::AssertValidity(const Interaction * interaction) const {

  // check that there is no charm production
  // (GENIE uses a special model for these cases)
  if(interaction->ExclTag().IsCharmEvent()) {
    LOG("Hadronization", pWARN) << "Can't hadronize charm events";
    return false;
  }
  // check the available mass
  double W = utils::kinematics::W(interaction);
  if(W < this->Wmin()) {
    LOG("Hadronization", pWARN) << "Low invariant mass, W = "
				<< W << " GeV!!";
    return false;
  }

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  if( ! target.HitQrkIsSet() ) {
    LOG("Hadronization", pWARN) << "Hit quark was not set!";
    return false;
  }

  int  probe       = init_state.ProbePdg();
  int  hit_nucleon = target.HitNucPdg();
  int  hit_quark   = target.HitQrkPdg();
  //bool from_sea    = target.HitSeaQrk();

  // check hit-nucleon assignment, input neutrino & weak current
  bool isp  = pdg::IsProton           (hit_nucleon);
  bool isn  = pdg::IsNeutron          (hit_nucleon);
  bool isv  = pdg::IsNeutrino         (probe);
  bool isvb = pdg::IsAntiNeutrino     (probe);
  bool isdm = pdg::IsDarkMatter         (probe);
  bool isl  = pdg::IsNegChargedLepton (probe);
  bool islb = pdg::IsPosChargedLepton (probe);
  bool iscc = proc_info.IsWeakCC      ();
  bool isnc = proc_info.IsWeakNC      ();
  bool isdmi = proc_info.IsDarkMatter  ();
  bool isem = proc_info.IsEM          ();
  if( !(iscc||isnc||isem||isdmi) ) {
    LOG("Hadronization", pWARN)
      << "Can only handle electro-weak interactions";
    return false;
  }
  if( !(isp||isn) || !(isv||isvb||isl||islb||isdm) ) {
    LOG("Hadronization", pWARN)
      << "Invalid initial state: probe = "
      << probe << ", hit_nucleon = " << hit_nucleon;
    return false;
  }

  // assert that the interaction mode is allowed
  bool isu  = pdg::IsUQuark     (hit_quark);
  bool isd  = pdg::IsDQuark     (hit_quark);
  bool iss  = pdg::IsSQuark     (hit_quark);
  bool isub = pdg::IsAntiUQuark (hit_quark);
  bool isdb = pdg::IsAntiDQuark (hit_quark);
  bool issb = pdg::IsAntiSQuark (hit_quark);

  bool allowed = (iscc && isv  && (isd||isub||iss))  ||
    (iscc && isvb && (isu||isdb||issb)) ||
    (isnc && (isv||isvb) && (isu||isd||isub||isdb||iss||issb)) ||
    (isdmi && isdm && (isu||isd||isub||isdb||iss||issb)) ||
    (isem && (isl||islb) && (isu||isd||isub||isdb||iss||issb));
  if(!allowed) {
    LOG("Hadronization", pWARN)
      << "Impossible interaction type / probe / hit quark combination!";
    return false;
  }

  return true;
}
//____________________________________________________________________________
double Pythia6Hadronization::MaxMult(const Interaction * interaction) const
{
  double W = interaction->Kine().W();

  double maxmult = TMath::Floor(1 + (W-kNeutronMass)/kPionMass);
  return maxmult;
}
//____________________________________________________________________________
TH1D * Pythia6Hadronization::CreateMultProbHist(double maxmult) const
{
  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  TH1D * mult_prob = new TH1D("mult_prob",
			      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
  mult_prob->SetDirectory(0);

  return mult_prob;
}
//____________________________________________________________________________
double Pythia6Hadronization::Wmin(void) const {

  return (kNucleonMass+kPionMass);
}
//____________________________________________________________________________
void Pythia6Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadronization::LoadConfig(void)
{
  // the configurable PYTHIA parameters used here are the ones used by NUX 
  // (see A.Rubbia's talk @ NuINT-01)
  // The defaults are the values used by PYTHIA
  // Use the NUX config set to set the tuned values as used in NUX.

  GetParam( "PYTHIA-SSBarSuppression", fSSBarSuppression ) ;
  GetParam( "PYTHIA-GaussianPt2",      fGaussianPt2      ) ;
  GetParam( "PYTHIA-NonGaussianPt2Tail", fNonGaussianPt2Tail  ) ;
  GetParam( "PYTHIA-RemainingEnergyCutoff", fRemainingECutoff ) ;

  GetParam( "PYTHIA-DiQuarkSuppression", fDiQuarkSuppression ) ;
  GetParam( "PYTHIA-LightVMesonSuppression", fLightVMesonSuppression ) ;
  GetParam( "PYTHIA-SVMesonSuppression", fSVMesonSuppression ) ;
  GetParam( "PYTHIA-Lunda", fLunda ) ;
  GetParam( "PYTHIA-Lundb", fLundb ) ;
  GetParam( "PYTHIA-LundaDiq", fLundaDiq ) ;

  fPythia->SetPARJ(1,  fDiQuarkSuppression ) ;
  fPythia->SetPARJ(11, fLightVMesonSuppression ) ;
  fPythia->SetPARJ(12, fSVMesonSuppression ) ;
  fPythia->SetPARJ(41, fLunda ) ;
  fPythia->SetPARJ(42, fLundb ) ;
  fPythia->SetPARJ(45, fLundaDiq ) ;
  
  fPythia->SetPARJ(2,  fSSBarSuppression);
  fPythia->SetPARJ(21, fGaussianPt2);
  fPythia->SetPARJ(23, fNonGaussianPt2Tail);
  fPythia->SetPARJ(33, fRemainingECutoff);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  double Wcut, Wmin ;
  GetParam( "Wcut",            Wcut ) ;
  GetParam( "KNO2PYTHIA-Wmin", Wmin ) ;

  if ( Wcut > Wmin ) {
    LOG("Pythia6Hadronization", pFATAL) << "Wcut value too high and in conflict with the KNO2PYTHIA-Wmin" ;
    LOG("Pythia6Hadronization", pFATAL) << "Wcut = " << Wcut ;
    LOG("Pythia6Hadronization", pFATAL) << "KNO2PYTHIA-Wmin = " << Wmin ;
  }

  LOG("Pythia6Hadronization", pDEBUG) << GetConfig() ;
}

