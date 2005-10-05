//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayers.

          The PythiaDecayer is a concrete implementation of the DecayModelI
          interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 20, 2004
 
*/
//____________________________________________________________________________

#include <vector>

#include <TMCParticle6.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Decay/PythiaDecayer.h"
#include "Messenger/Messenger.h"

using std::vector;

using namespace genie;

//the actual PYTHIA calls

extern "C" void py1ent_(int *,  int *, double *, double *, double *);
extern "C" void pydecy_(int *);

//____________________________________________________________________________
PythiaDecayer::PythiaDecayer() :
DecayModelI()
{
  fName     = "genie::PythiaDecayer";
  fParamSet = "Default";

  FindConfig();

  fPythia = new TPythia6();
}
//____________________________________________________________________________
PythiaDecayer::PythiaDecayer(const char * param_set) :
DecayModelI(param_set)
{
  fName = "genie::PythiaDecayer";

  FindConfig();

  fPythia = new TPythia6();
}
//____________________________________________________________________________
PythiaDecayer::~PythiaDecayer() 
{ 
  delete fPythia;
}
//____________________________________________________________________________
bool PythiaDecayer::IsHandled(int pdg_code) const
{
// does not handle requests to decay baryon resonances
  
  if( utils::res::IsBaryonResonance(pdg_code) ) {

    LOG("Decay", pINFO)
         << "\n *** The particle with PDG-Code = "
                           << pdg_code << " is not decayed by this algorithm";
    return false;

  } else return true;
}
//____________________________________________________________________________
void PythiaDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * PythiaDecayer::Decay(const DecayerInputs_t & inp) const
{
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;

  
  // check whether we should inhibit some channels

  if(inp.InhibitedChannels)
       this->SwitchOffInhibitedChannels(inp.PdgCode, inp.InhibitedChannels);

  int ip = 0;

  int    pdgc  = inp.PdgCode;
  double E     = inp.P4->Energy();
  double Theta = inp.P4->Theta();
  double Phi   = inp.P4->Phi();

  fPythia->SetMSTJ(22,1);

  py1ent_(&ip, &pdgc, &E, &Theta, &Phi);

  // check whether we are asked to force the decay / default = false
 
  bool force_decay = (fConfig->Exists("force-decay")) ? 
                                    fConfig->GetBool("force-decay") : false;
  if(force_decay) {
    int F = 1;
    pydecy_(&F); // FORCE DECAY
  }

  
  // get decay products

  fPythia->GetPrimaries();
  
  TClonesArray * pythia_particles =
                           (TClonesArray *) fPythia->ImportParticles("All");

                           
  //-- if we switched of some channels, now we should restore TPythia's state

  if(inp.InhibitedChannels) this->SwitchOnAllChannels(inp.PdgCode);
  

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  TClonesArray * particle_list = new TClonesArray(
                             "TMCParticle", pythia_particles->GetEntries() );

  register unsigned int i = 0;
  TMCParticle * particle = 0;

  TIter particle_iter(pythia_particles);

  while( (particle = (TMCParticle *) particle_iter.Next()) ) {

       LOG("Decay", pINFO)
               << "Adding decay product with PDGC = " << particle->GetKF();

       new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  particle_list->SetOwner(true);
    
  return particle_list;
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
