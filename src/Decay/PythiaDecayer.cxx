//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - June 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <vector>

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
  //delete fPythia;
}
//____________________________________________________________________________
bool PythiaDecayer::IsHandled(int code) const
{
// does not handle requests to decay baryon resonances
  
  if( utils::res::IsBaryonResonance(code) ) {
     LOG("Decay", pINFO)
       << "This algorithm can not decay particles with PDG code = " << code;
     return false;
  } else return true;
}
//____________________________________________________________________________
void PythiaDecayer::Initialize(void) const
{
  fPythia = TPythia6::Instance();
}
//____________________________________________________________________________
TClonesArray * PythiaDecayer::Decay(const DecayerInputs_t & inp) const
{
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;
  
  this->SyncSeeds();

  //-- check whether we should inhibit some channels
  if(inp.InhibitedChannels)
       this->SwitchOffInhibitedChannels(inp.PdgCode, inp.InhibitedChannels);

  int    ip    = 0;
  int    pdgc  = inp.PdgCode;
  double E     = inp.P4->Energy();
  double Theta = inp.P4->Theta();
  double Phi   = inp.P4->Phi();

  fPythia->SetMSTJ(22,1);

  py1ent_(&ip, &pdgc, &E, &Theta, &Phi);

  //-- check whether we are asked to force the decay 
  if(fForceDecay) {
    int F = 1;
    pydecy_(&F); // FORCE DECAY
  }
  
  //-- get decay products
  fPythia->GetPrimaries();
  TClonesArray * impl = (TClonesArray *) fPythia->ImportParticles("All");

  if(!impl) return 0;

  //-- if we switched of some channels, now we should restore TPythia's state
  if(inp.InhibitedChannels) this->SwitchOnAllChannels(inp.PdgCode);
  
  //-- copy PYTHIA container to a new TClonesArray so as to transfer ownership
  //   of the container and of its elements to the calling method
  TClonesArray * pl = new TClonesArray("TMCParticle", impl->GetEntries());

  register unsigned int i = 0;
  TMCParticle * p = 0;
  TIter particle_iter(impl);

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

    string type = (p->GetKS()==11) ? "mother: " : "+ daughter: ";
    SLOG("Decay", pINFO)
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
  return 1; // does not generate weighted decays
}
//____________________________________________________________________________
void PythiaDecayer::SwitchOnAllChannels(int pdgc) const
{
  LOG("Decay", pINFO)
         << "Switching ON all PYTHIA decay channels for particle = " << pdgc;

  int kc = fPythia->Pycomp( pdgc );

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  for(int ichannel = first_channel; ichannel < last_channel; ichannel++) {

     fPythia->SetMDME(ichannel,1,1); // switch-on
  }
}
//____________________________________________________________________________
void PythiaDecayer::SwitchOffInhibitedChannels(
                               int pdgc, const TClonesArray * inhibited) const
{
  LOG("Decay", pINFO)
          << "Switching OFF inhibited decay channels for particle = " << pdgc;

  this->SwitchOnAllChannels(pdgc);

  int kc = fPythia->Pycomp( pdgc );

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  for(int ichannel = first_channel; ichannel < last_channel; ichannel++) {

     TDecayChannel * dc = 0;  

     TIter dciter(inhibited);

     //-- loop over the input inhibited decay channels and check whether
     //   one can be matched with the current decay channel

     //   The channels are matched if all decay products are matched
     
     while( (dc = (TDecayChannel *) dciter.Next()) ) {

        LOG("Decay", pINFO)
              << "\nComparing PYTHIA's channel = " << ichannel
                                  << " with TDecayChannel = " << dc->Number();
 
        if( this->MatchDecayChannel(ichannel,*dc) )
                                fPythia->SetMDME(ichannel,1,0); // switch-off
     }
  }
}
//____________________________________________________________________________
bool PythiaDecayer::MatchDecayChannel(int ichannel, TDecayChannel & dc) const
{
  // num. of daughters in the input TDecayChannel & the input PYTHIA ichannel

  int nd = dc.NDaughters(); 

  int py_nd = 0;
  for (int i = 0; i < 5; i++) if( fPythia->GetKFDP(ichannel,i) ) py_nd++;

  LOG("Decay", pINFO)
    << "NDaughters: PYTHIA = " << py_nd << ", ROOT's TDecayChannel = " << nd;

  if(nd != py_nd) return false;
  
  // if the two channels have the same num. of daughters, then compare them

  vector<int> dc_daughter(nd); // daughters for the input TDecayChannel
  vector<int> py_daughter(nd); // daughters for the input PYTHIA's ichannel

  for(int i = 0; i < nd; i++) dc_daughter[i] = dc.DaughterPdgCode(i);
  for(int i = 0; i < nd; i++) py_daughter[i] = fPythia->GetKFDP(ichannel,i);

  // sort both daughter lists 
  sort( dc_daughter.begin(), dc_daughter.end() ); 
  sort( py_daughter.begin(), py_daughter.end() );

  bool channel_matched = true;
  
  for(int i = 0; i < nd; i++)
      channel_matched = channel_matched && (dc_daughter[i] == py_daughter[i]);

  if( channel_matched ) { LOG("Decay", pINFO) << " *** channels matched"; }
  
  return channel_matched;
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

  //-- check whether we are asked to force the decay / default = false
  fForceDecay = fConfig->GetBoolDef("ForceDecay", false);
}
//____________________________________________________________________________
void PythiaDecayer::SyncSeeds(void) const
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

