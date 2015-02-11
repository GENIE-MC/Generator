//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 21, 2009 - CA
   Remove SyncSeeds() function. Now the GENIE/PYTHIA6 random number generator
   seed is synchonized at the genie::RandomGen() initialization.
 @ Oct 02, 2009 - CA
   Re-organize code and implement the `UnInhibitDecay(int,TDecayChannel*)
   const' and `InhibitDecay(int,TDecayChannel*) const' methods.
   Test/fix the code to match a ROOT TDecayChannel to a PYTHIA6 decay channel.
   Decay() returns null if decay is inhibited or if the sum{branching ratios}
   for all enabled decay channels is non-positive. In case of inhibited decay
   channels, a weight is calculated as w = 1./sum{BR for enabled channels}.
 @ Feb 04, 2010 - CA
   Comment out (unused) code using the fForceDecay flag

*/
//____________________________________________________________________________

#include <vector>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TDecayChannel.h>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Units.h"
#include "Conventions/Constants.h"
#include "Decay/PythiaDecayer.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"

using std::vector;

using namespace genie;

//the actual PYTHIA calls

extern "C" void py1ent_(int *,  int *, double *, double *, double *);
extern "C" void pydecy_(int *);

//____________________________________________________________________________
PythiaDecayer::PythiaDecayer() :
DecayModelI("genie::PythiaDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::PythiaDecayer(string config) :
DecayModelI("genie::PythiaDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::~PythiaDecayer() 
{ 

}
//____________________________________________________________________________
bool PythiaDecayer::IsHandled(int code) const
{
// does not handle requests to decay baryon resonances
  
  if( utils::res::IsBaryonResonance(code) ) {
     LOG("PythiaDec", pINFO)
       << "This algorithm can not decay particles with PDG code = " << code;
     return false;
  } else return true;
}
//____________________________________________________________________________
void PythiaDecayer::Initialize(void) const
{
  fPythia = TPythia6::Instance();
  fWeight = 1.;

  // sync GENIE/PYTHIA6 seeds
  RandomGen::Instance();
}
//____________________________________________________________________________
TClonesArray * PythiaDecayer::Decay(const DecayerInputs_t & inp) const
{
  fWeight = 1.; // reset weight

  int pdgc = inp.PdgCode;

  if ( ! this->IsHandled(pdgc) ) return 0;
  
  int kc   = fPythia->Pycomp(pdgc);
  int mdcy = fPythia->GetMDCY(kc, 1);
  if(mdcy == 0) {
    LOG("PythiaDec", pNOTICE)
       << (PDGLibrary::Instance())->Find(pdgc)->GetName() 
       << " decays are inhibited!";
    return 0;
  }

  double sumbr = this->SumBR(kc);
  if(sumbr <= 0) {
    LOG("PythiaDec", pNOTICE)
       << "The sum of enabled "
       << (PDGLibrary::Instance())->Find(pdgc)->GetName() 
       << " decay channel branching rations is non-positive!";
    return 0;
  }

  fWeight = 1./sumbr; // update weight to account for inhibited channels

  int    ip    = 0;
  double E     = inp.P4->Energy();
  double Theta = inp.P4->Theta();
  double Phi   = inp.P4->Phi();

  fPythia->SetMSTJ(22,1);

  py1ent_(&ip, &pdgc, &E, &Theta, &Phi);

/*
  //-- check whether we are asked to force the decay 
  if(fForceDecay) {
    int F = 1;
    pydecy_(&F); // FORCE DECAY
  }
*/
  
  //-- get decay products
  fPythia->GetPrimaries();
  TClonesArray * impl = (TClonesArray *) fPythia->ImportParticles("All");

  if(!impl) return 0;
  
  //-- copy PYTHIA container to a new TClonesArray so as to transfer ownership
  //   of the container and of its elements to the calling method
  TClonesArray * pl = new TClonesArray("TMCParticle", impl->GetEntries());

  register unsigned int i = 0;
  TMCParticle * p = 0;
  TIter particle_iter(impl);

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

    string type = (p->GetKS()==11) ? "mother: " : "+ daughter: ";
    SLOG("PythiaDec", pINFO)
       << type << p->GetName() << " (pdg-code = "
          << p->GetKF() << ", m = " << p->GetMass() 
             << ", E = " << p->GetEnergy() << ")";

    p->SetLifetime (p->GetLifetime() * units::mm/constants::kLightSpeed);
    p->SetTime     (p->GetTime()     * units::mm/constants::kLightSpeed);
    p->SetVx       (p->GetVx()       * units::mm);
    p->SetVy       (p->GetVy()       * units::mm);
    p->SetVz       (p->GetVz()       * units::mm);

    new ( (*pl)[i++] ) TMCParticle(*p);
  }

  //-- transfer ownership and return
  pl->SetOwner(true);    
  return pl;
}
//____________________________________________________________________________
double PythiaDecayer::Weight(void) const 
{
  return fWeight; 
}
//____________________________________________________________________________
void PythiaDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return; 

  int kc = fPythia->Pycomp(pdgc);

  if(!dc) {
    LOG("PythiaDec", pINFO)
       << "Switching OFF ALL decay channels for particle = " << pdgc;
    fPythia->SetMDCY(kc, 1,0); 
    return;
  }

  LOG("PythiaDec", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdgc;

  int ichannel = this->FindPythiaDecayChannel(kc, dc);
  if(ichannel != -1) {
    fPythia->SetMDME(ichannel,1,0); // switch-off
  }
}
//____________________________________________________________________________
void PythiaDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return; 

  int kc = fPythia->Pycomp(pdgc);

  if(!dc) {
    LOG("PythiaDec", pINFO)
      << "Switching ON all PYTHIA decay channels for particle = " << pdgc;

    fPythia->SetMDCY(kc, 1,1); 

    int first_channel = fPythia->GetMDCY(kc,2);
    int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

    for(int ichannel = first_channel; 
            ichannel < last_channel; ichannel++) {
         fPythia->SetMDME(ichannel,1,1); // switch-on
    }
    return;
  }//!dc

  LOG("PythiaDec", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdgc;

  int ichannel = this->FindPythiaDecayChannel(kc, dc);
  if(ichannel != -1) {
    fPythia->SetMDME(ichannel,1,1); // switch-on
  }
}
//____________________________________________________________________________
double PythiaDecayer::SumBR(int kc) const
{
// Sum of branching ratios for enabled channels
//
  double sumbr=0.;

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  bool has_inhibited_channels=false;

  // loop over pythia decay channels
  for(int ichannel = first_channel; 
          ichannel < last_channel; ichannel++) {

     bool enabled = (fPythia->GetMDME(ichannel,1) == 1);
     if (!enabled) { 
       has_inhibited_channels = true; 
     } else {
       sumbr += fPythia->GetBRAT(ichannel);
     }

/*
     LOG("PythiaDec", pDEBUG) 
        << "ich: " << ichannel << ", " << ((enabled)? "ON" : "OFF")
        << ", br = " << fPythia->GetBRAT(ichannel)
        << ", curr_sum{br}[enabled] = " << sumbr;
*/
  }

  if(!has_inhibited_channels) return 1.;

  LOG("PythiaDec", pINFO) << "Sum{BR} = " << sumbr;

  return sumbr;
}
//____________________________________________________________________________
int PythiaDecayer::FindPythiaDecayChannel(int kc, TDecayChannel* dc) const
{	
  if(!dc) return -1;

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  bool found_match = false;

  // loop over pythia decay channels
  for(int ichannel = first_channel; 
          ichannel < last_channel; ichannel++) {

     // does the  current pythia channel matches the input TDecayChannel?
     LOG("PythiaDec", pINFO)
         << "\nComparing PYTHIA's channel = " << ichannel
         << " with TDecayChannel = " << dc->Number();
 
     found_match = this->MatchDecayChannels(ichannel,dc);
     if(found_match) {
         LOG("PythiaDec", pNOTICE)
            << " ** TDecayChannel id = " << dc->Number()
            << " corresponds to PYTHIA6 channel id = " << ichannel;
         return ichannel;
     }//match?
  }//loop pythia decay ch.

  LOG("PythiaDec", pWARN)
     << " ** No PYTHIA6 decay channel match found for "
     << "TDecayChannel id = " << dc->Number();

  return -1;
}
//____________________________________________________________________________
bool PythiaDecayer::MatchDecayChannels(int ichannel, TDecayChannel* dc) const
{
  // num. of daughters in the input TDecayChannel & the input PYTHIA ichannel
  int nd = dc->NDaughters(); 

  int py_nd = 0;
  for (int i = 1; i <= 5; i++) {
     if(fPythia->GetKFDP(ichannel,i) != 0) py_nd++;
  }

  LOG("PythiaDec", pDEBUG)
    << "NDaughters: PYTHIA = " << py_nd << ", ROOT's TDecayChannel = " << nd;

  if(nd != py_nd) return false;
  
  //
  // if the two channels have the same num. of daughters, then compare them
  //
 
  // store decay daughters for the input TDecayChannel
  vector<int> dc_daughter(nd); 
  int k=0;
  for( ; k < nd; k++) {
     dc_daughter[k] = dc->DaughterPdgCode(k);
  }
  // store decay daughters for the input PYTHIA's ichannel
  vector<int> py_daughter(nd); 
  k=0;
  for(int i = 1; i <= 5; i++) {
     if(fPythia->GetKFDP(ichannel,i) == 0) continue;
     py_daughter[k] = fPythia->GetKFDP(ichannel,i);
     k++;
  }

  // sort both daughter lists 
  sort( dc_daughter.begin(), dc_daughter.end() ); 
  sort( py_daughter.begin(), py_daughter.end() );

  // compare  
  for(int i = 0; i < nd; i++) {
    if(dc_daughter[i] != py_daughter[i]) return false;
  }
  return true;
}
//____________________________________________________________________________
void PythiaDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaDecayer::LoadConfig(void)
{
// Read configuration options or set defaults
/*
  // check whether we are asked to force the decay / default = false
  fForceDecay = fConfig->GetBoolDef("ForceDecay", false);
*/
}
//____________________________________________________________________________
