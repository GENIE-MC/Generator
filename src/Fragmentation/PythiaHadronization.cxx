//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 04, 2007 - CA
   Handle very rare failure mode where a bare quark or di-quark appears in
   the final state
*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TClonesArray.h>
#include <TMath.h>
#include <TH1D.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Decay/DecayModelI.h"
#include "Fragmentation/PythiaHadronization.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/FragmRecUtils.h"

using namespace genie;
using namespace genie::constants;

//-- the actual PYTHIA call

extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
PythiaHadronization::PythiaHadronization() :
HadronizationModelBase("genie::PythiaHadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaHadronization::PythiaHadronization(string config) :
HadronizationModelBase("genie::PythiaHadronization", config)
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

  if(!this->AssertValidity(interaction)) {
     LOG("PythiaHad", pWARN) << "Returning a null particle list!";
     return 0;
  }
  this->SyncSeeds();

  //-- get kinematics / init-state / process-info

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
  bool iss  = pdg::IsSQuark     (hit_quark);
  bool isub = pdg::IsAntiUQuark (hit_quark);
  bool isdb = pdg::IsAntiDQuark (hit_quark);
  bool issb = pdg::IsAntiSQuark (hit_quark);

  bool allowed = (iscc && isv  && (isd||isub||iss))  ||
                 (iscc && isvb && (isu||isdb||issb)) ||
                 (isnc && (isv||isvb) && (isu||isd||isub||isdb||iss||issb));
  assert(allowed);

  //-- generate the quark system (q + qq) initiating the hadronization

  int  final_quark = 0; // leading quark (hit quark after the interaction)
  int  diquark     = 0; // remnant diquark (xF<0 at hadronic CMS)

  //-- figure out the what happens to the hit quark after the interaction

  if (proc_info.IsWeakNC()) final_quark = hit_quark;
  else {
    if      (isv  && isd ) final_quark = kPdgUQuark;
    else if (isv  && iss ) final_quark = kPdgUQuark;
    else if (isv  && isub) final_quark = kPdgAntiDQuark;
    else if (isvb && isu ) final_quark = kPdgDQuark;
    else if (isvb && isdb) final_quark = kPdgAntiUQuark;
    else if (isvb && issb) final_quark = kPdgAntiUQuark;
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

    // The following section needs revisiting.

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

    // The neutrino is scatterred off s or sbar sea quarks 
    // For the time being I will handle s like d and sbar like dbar (copy & paste
    // from above) so that I conserve charge. 

    if(iss || issb) {
       LOG("PythiaHad", pNOTICE) 
                 << "Can not really handle a hit s or sbar quark / Faking it";

       if(isp && iss) { diquark = kPdgUUDiquarkS1; }
       if(isn && iss) { diquark = kPdgUDDiquarkS1; }
       if(isp && issb && iscc) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
       if(isp && issb && isnc) {final_quark = kPdgDQuark; diquark = kPdgUUDiquarkS1;}
       if(isn && issb && iscc) {final_quark = kPdgDQuark; diquark = kPdgDDDiquarkS1;}
       if(isn && issb && isnc) {final_quark = kPdgDQuark; diquark = kPdgUDDiquarkS1;}
    }
 
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

  // Determine how jetset treats unsstable particles appearning in hadronization
  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),        1,0); // don't decay pi0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),         1,0); // don't decay K0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),     1,0); // don't decay \bar{K0}
  fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),     1,0); // don't decay Lambda0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda), 1,0); // don't decay \bar{Lambda0}
  fPythia->SetMDCY(fPythia->Pycomp(2214),           1,1); // decay Delta+
  fPythia->SetMDCY(fPythia->Pycomp(2224),           1,1); // decay Delta++

  // -- hadronize --
  py2ent_(&ip, &final_quark, &diquark, &W); // hadronizer

  // restore pythia settings
  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),        1,1);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),         1,1);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),     1,1);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),     1,1); 
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda), 1,1); 

  // get LUJETS record
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles =
                      (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("TMCParticle", np);
  particle_list->SetOwner(true);

  register unsigned int i = 0;
  TMCParticle * particle = 0;
  TIter particle_iter(pythia_particles);

  while( (particle = (TMCParticle *) particle_iter.Next()) ) {
     LOG("PythiaHad", pDEBUG)
          << "Adding final state particle pdgc = " << particle->GetKF() 
          << " with status = " << particle->GetKS();

     if(particle->GetKS() == 1) {
	if( pdg::IsQuark  (particle->GetKF()) || 
            pdg::IsDiQuark(particle->GetKF()) ) {
                LOG("PythiaHad", pERROR)
                  << "Hadronization failed! Bare quark/di-quarks appear in final state!";
            particle_list->Delete();
            delete particle_list;
            return 0;            
        }
     }

     // fix numbering scheme used for mother/daughter assignments
     particle->SetParent     (particle->GetParent()     - 1);
     particle->SetFirstChild (particle->GetFirstChild() - 1);
     particle->SetLastChild  (particle->GetLastChild()  - 1);

     // insert the particle in the list
     new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  utils::fragmrec::Print(particle_list);
  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * PythiaHadronization::SelectParticles(
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

  TMCParticle * particle = 0;
  TIter particle_iter(particle_list);

  while ((particle = (TMCParticle *) particle_iter.Next())) 
  {
    if (particle->GetKS()==1) pdgcv->push_back(particle->GetKF());
  }
  particle_list->Delete();
  delete particle_list;

  return pdgcv;
}
//____________________________________________________________________________
TH1D * PythiaHadronization::MultiplicityProb(
                       const Interaction * interaction, Option_t * opt) const
{
// Similar comments apply as in SelectParticles()

  if(!this->AssertValidity(interaction)) {
     LOG("PythiaHad", pWARN) 
                << "Returning a null multipicity probability distribution!";
     return 0;
  }
  double maxmult   = this->MaxMult(interaction);
  TH1D * mult_prob = this->CreateMultProbHist(maxmult);

  const int nev=500;
  TMCParticle * particle = 0;

  for(int iev=0; iev<nev; iev++) {

     TClonesArray * particle_list = this->Hadronize(interaction);
     double         weight        = this->Weight();

     if(!particle_list) { iev--; continue; }

     int n = 0;
     TIter particle_iter(particle_list);
     while ((particle = (TMCParticle *) particle_iter.Next())) 
     {
       if (particle->GetKS()==1) n++;
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
    SLOG("PythiaHad", pWARN) << "probability distribution integral = 0";
    return mult_prob;
  }

  string option(opt);

  bool apply_neugen_Rijk = option.find("+LowMultSuppr") != string::npos;
  bool renormalize       = option.find("+Renormalize")  != string::npos;

  // Apply the NeuGEN probability scaling factors -if requested-
  if(apply_neugen_Rijk) {
    SLOG("KNOHad", pINFO) << "Applying NeuGEN scaling factors";
     // Only do so for W<Wcut
     const Kinematics & kinematics = interaction->Kine();
     double W = kinematics.W();
     if(W<fWcut) {
       this->ApplyRijk(interaction, renormalize, mult_prob);
     } else {
        SLOG("PythiaHad", pDEBUG)
              << "W = " << W << " < Wcut = " << fWcut
                                << " - Will not apply scaling factors";
     }//<wcut?
  }//apply?

  return mult_prob;
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
              "SSBarSuppression", gc->GetDouble("PYTHIA-SSBarSuppression"));
  fGaussianPt2 = fConfig->GetDoubleDef(
                        "GaussianPt2", gc->GetDouble("PYTHIA-GaussianPt2"));
  fNonGaussianPt2Tail = fConfig->GetDoubleDef(
          "NonGaussianPt2Tail", gc->GetDouble("PYTHIA-NonGaussianPt2Tail"));
  fRemainingECutoff = fConfig->GetDoubleDef(
    "RemainingEnergyCutoff", gc->GetDouble("PYTHIA-RemainingEnergyCutoff"));

  fPythia->SetPARJ(2,  fSSBarSuppression);
  fPythia->SetPARJ(21, fGaussianPt2);
  fPythia->SetPARJ(23, fNonGaussianPt2Tail);
  fPythia->SetPARJ(33, fRemainingECutoff);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));

  // decayer
  fDecayer = 0;
  if(fConfig->Exists("Decayer")) {
     fDecayer = dynamic_cast<const DecayModelI *> (this->SubAlg("Decayer"));
     assert(fDecayer);
  }

  // Load NEUGEN multiplicity probability scaling parameters Rijk
  fRvpCCm2  = fConfig->GetDoubleDef(
                      "R-vp-CC-m2", gc->GetDouble("DIS-HMultWgt-vp-CC-m2"));
  fRvpCCm3  = fConfig->GetDoubleDef(
                      "R-vp-CC-m3", gc->GetDouble("DIS-HMultWgt-vp-CC-m3"));
  fRvpNCm2  = fConfig->GetDoubleDef(
                      "R-vp-NC-m2", gc->GetDouble("DIS-HMultWgt-vp-NC-m2"));
  fRvpNCm3  = fConfig->GetDoubleDef(
                      "R-vp-NC-m3", gc->GetDouble("DIS-HMultWgt-vp-NC-m3"));
  fRvnCCm2  = fConfig->GetDoubleDef(
                      "R-vn-CC-m2", gc->GetDouble("DIS-HMultWgt-vn-CC-m2"));
  fRvnCCm3  = fConfig->GetDoubleDef(
                      "R-vn-CC-m3", gc->GetDouble("DIS-HMultWgt-vn-CC-m3"));
  fRvnNCm2  = fConfig->GetDoubleDef(
                      "R-vn-NC-m2", gc->GetDouble("DIS-HMultWgt-vn-NC-m2"));
  fRvnNCm3  = fConfig->GetDoubleDef(
                      "R-vn-NC-m3", gc->GetDouble("DIS-HMultWgt-vn-NC-m3"));
  fRvbpCCm2 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m2",gc->GetDouble("DIS-HMultWgt-vbp-CC-m2"));
  fRvbpCCm3 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m3",gc->GetDouble("DIS-HMultWgt-vbp-CC-m3"));
  fRvbpNCm2 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m2",gc->GetDouble("DIS-HMultWgt-vbp-NC-m2"));
  fRvbpNCm3 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m3",gc->GetDouble("DIS-HMultWgt-vbp-NC-m3"));
  fRvbnCCm2 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m2",gc->GetDouble("DIS-HMultWgt-vbn-CC-m2"));
  fRvbnCCm3 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m3",gc->GetDouble("DIS-HMultWgt-vbn-CC-m3"));
  fRvbnNCm2 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m2",gc->GetDouble("DIS-HMultWgt-vbn-NC-m2"));
  fRvbnNCm3 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m3",gc->GetDouble("DIS-HMultWgt-vbn-NC-m3"));

  LOG("PythiaHad", pDEBUG) << *fConfig;
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
bool PythiaHadronization::AssertValidity(const Interaction * interaction) const
{
  if(interaction->ExclTag().IsCharmEvent()) {
     LOG("PythiaHad", pWARN) << "Can't hadronize charm events";
     return false;
  }
  double W = utils::kinematics::W(interaction);
  if(W < this->Wmin()) {
     LOG("PythiaHad", pWARN) << "Low invariant mass, W = " << W << " GeV!!";
     return false;
  }
  return true;
}
//____________________________________________________________________________
void PythiaHadronization::SwitchDecays(int pdgc, bool on_off) const
{
  LOG("PythiaHad", pNOTICE)
     << "Switching " << ((on_off) ? "ON" : "OFF")
                     << " all PYTHIA decay channels for particle = " << pdgc;

  int flag     = (on_off) ? 1 : 0;
  int kc       = fPythia->Pycomp(pdgc);
  int first_ch = fPythia->GetMDCY(kc,2);
  int last_ch  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  for(int ich = first_ch; ich < last_ch; ich++) fPythia->SetMDME(ich,1,flag);
}
//____________________________________________________________________________
void PythiaHadronization::HandleDecays(TClonesArray * plist) const
{
// Handle decays of unstable particles if requested through the XML config.
// The default is not to decay the particles at this stage (during event
// generation, the UnstableParticleDecayer event record visitor decays what
// is needed to be decayed later on). But, when comparing various models
// (eg PYTHIA vs KNO) independently and not within the full MC simulation
// framework it might be necessary to force the decays at this point.

  if(!fDecayer) {
    LOG("PythiaHad", pWARN) << "No decayer was specified!";
    return;
  }

  this->SwitchDecays(kPdgLambda,     true); // decay Lambda
  this->SwitchDecays(kPdgAntiLambda, true); // decay \bar{Lambda}
  this->SwitchDecays(kPdgSigmaP,     true); // decay Sigma+
  this->SwitchDecays(kPdgSigma0,     true); // decay Sigma0
  this->SwitchDecays(kPdgSigmaM,     true); // decay Sigma-
  this->SwitchDecays(kPdgAntiSigmaP, true); // decay Sigma+
  this->SwitchDecays(kPdgAntiSigma0, true); // decay Sigma0
  this->SwitchDecays(kPdgAntiSigmaM, true); // decay Sigma-
  this->SwitchDecays(kPdgXi0,        true); // decay Xi0
  this->SwitchDecays(kPdgXiM,        true); // decay Xi-
  this->SwitchDecays(kPdgAntiXi0,    true); // decay \bar{Xi0}
  this->SwitchDecays(kPdgAntiXiP,    true); // decay \bar{Xi+}
  this->SwitchDecays(kPdgOmegaM,     true); // decay Omega-
  this->SwitchDecays(kPdgAntiOmegaP, true); // decay \bar{Omega+}

  int mstj21 = fPythia->GetMSTJ(21);
  fPythia->SetMSTJ(21,1); 
  fPythia->SetMSTJ(22,2);                  
  fPythia->SetPARJ(71,100);                  

  //-- loop through the fragmentation event record & decay unstables
  int idecaying   = -1; // position of decaying particle
  TMCParticle * p =  0; // current particle

  TIter piter(plist);
  while ( (p = (TMCParticle *) piter.Next()) ) {
     idecaying++;
     int status = p->GetKS();
     int pdg    = p->GetKF();

     bool decay_it = (status<10) && 
                     ( pdg == kPdgLambda ||
                       pdg == kPdgAntiLambda ||
                       pdg == kPdgSigmaP ||
                       pdg == kPdgSigma0 ||
                       pdg == kPdgSigmaM ||
                       pdg == kPdgAntiSigmaP ||
                       pdg == kPdgAntiSigma0 ||
                       pdg == kPdgAntiSigmaM ||
                       pdg == kPdgXi0 ||
                       pdg == kPdgXiM ||
                       pdg == kPdgAntiXi0 ||
                       pdg == kPdgAntiXiP ||
                       pdg == kPdgOmegaM  ||
                       pdg == kPdgAntiOmegaP );

     // bother for final state particle only
     if(decay_it) {

          LOG("PythiaHad", pINFO)
                     << "Decaying particle with pdgc = " << p->GetKF();

          DecayerInputs_t dinp;

          TLorentzVector p4;
          p4.SetPxPyPzE(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());

          dinp.PdgCode = p->GetKF();
          dinp.P4      = &p4;

          TClonesArray * decay_products = fDecayer->Decay(dinp);
          if(decay_products) {
                  //--  mark the parent particle as decayed & set daughters
                  p->SetKS(11);

                  int nfp = plist->GetEntries();          // n. fragm. products
                  int ndp = decay_products->GetEntries(); // n. decay products

                  p->SetFirstChild ( nfp );          // decay products added at
                  p->SetLastChild  ( nfp + ndp -1 ); // the end of the fragm.rec.

                  //--  add decay products to the fragmentation record
                  TMCParticle * dp = 0;
                  TIter dpiter(decay_products);

                  while ( (dp = (TMCParticle *) dpiter.Next()) ) {
  	  	     if(dp->GetKS()>10) continue;
                     dp->SetParent(idecaying);
                     new ( (*plist)[plist->GetEntries()] ) TMCParticle(*dp);
                  }

                  //-- clean up decay products
                  decay_products->Delete();
                  delete decay_products;
           }

     } // KS < 10 : final state particle (as in PYTHIA LUJETS record)
  } // particles in fragmentation record

  fPythia->SetMSTJ(21,mstj21); // restore mstj(21)
}
//____________________________________________________________________________

