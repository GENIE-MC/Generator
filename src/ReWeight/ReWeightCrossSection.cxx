//____________________________________________________________________________
/*!

\class   genie::ReWeightCrossSection

\brief   A cross section model reweighting engine.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#include "Algorithm/AlgCmp.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "ReWeight/ReWeightCrossSection.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
ReWeightCrossSection::ReWeightCrossSection(void)
{
  this->Initialize();
}
//___________________________________________________________________________
ReWeightCrossSection::~ReWeightCrossSection(void)
{

}
//___________________________________________________________________________
void ReWeightCrossSection::Initialize(void)
{
  this->DiffCrossSecType( kScQuasiElastic,  kPSQ2fE  );
  this->DiffCrossSecType( kScDeepInelastic, kPSxyfE  );
  this->DiffCrossSecType( kScResonant,      kPSWQ2fE );
  this->DiffCrossSecType( kScCoherent,      kPSxyfE  );
}
//___________________________________________________________________________
void ReWeightCrossSection::HandleInitState(const InitialState & init_state)
{
  GEVGDriver * evg_driver = fGPool.FindDriver(init_state);
  if(!evg_driver) {
    evg_driver = new GEVGDriver;
    evg_driver->Configure(init_state);
    fGPool.insert( GEVGPool::value_type(init_state.AsString(), evg_driver) );
  }
}
//___________________________________________________________________________
void ReWeightCrossSection::DontReweight(const Interaction & interaction)
{
  Interaction * in = new Interaction(interaction);
  fNoRewProc.push_back(in);
}
//___________________________________________________________________________
void ReWeightCrossSection::DiffCrossSecType(
                     ScatteringType_t sct, KinePhaseSpace_t kps)
{
  fCrossSecModelPhSp.insert(
               map<ScatteringType_t,KinePhaseSpace_t>::value_type(sct,kps));
}
//___________________________________________________________________________
double ReWeightCrossSection::NewWeight(const EventRecord & event)
{
  // Get event summary (Interaction) from the input event
  Interaction & interaction = * event.Summary();
/*
  if(!interaction) {
      LOG("ReWeight", pWARN) << "Null interaction!!";
      return 0;
  }
*/

  // Reweight?
  InteractionList::const_iterator iiter = fNoRewProc.begin();
  for( ; iiter != fNoRewProc.end(); ++iiter) {
	const Interaction & norewint = **iiter;
	if(interaction == norewint) {
            LOG("ReWeight", pDEBUG) << "Skipping reweighting: " << interaction;
            return 1;  
        }
  }

  LOG("ReWeight", pDEBUG) << "Computing new weight for: \n" << interaction;

  // Find the event generation driver that handles the given initial state
  const InitialState & init_state = interaction.InitState();
  GEVGDriver * evg_driver = fGPool.FindDriver(init_state);
  if(!evg_driver) {
    LOG("ReWeight", pERROR)
       << "No event generator driver for init state: " << init_state.AsString();
    return 0;
  }

  // Find the event generation thread that handles the given interaction
  const EventGeneratorI * evg_thread = evg_driver->FindGenerator(&interaction);
  if(!evg_thread) {
    LOG("ReWeight", pERROR)
      << "No event generator thread for interaction: " << interaction;
    return 0;
  }

  // Get the cross section model associated with that thread
  const XSecAlgorithmI * xsec_model = evg_thread->CrossSectionAlg();
  if(!xsec_model) {
    LOG("ReWeight", pERROR)
      << "No cross section model for interaction: " << interaction;
    return 0;
  }

  // Get the kinematical phase space used for computing the differential
  // cross sections stored in the event
  ScatteringType_t sct = interaction.ProcInfo().ScatteringTypeId();
  map<ScatteringType_t, KinePhaseSpace_t>::const_iterator 
                               kpsi = fCrossSecModelPhSp.find(sct); 

  if(kpsi == fCrossSecModelPhSp.end()) return 1;
  KinePhaseSpace_t kps = kpsi->second;
  if(kps==kPSNull) return 1;

  // Copy the 'selected' kinematics into the 'running' kinematics
  interaction.KinePtr()->UseSelectedKinematics();

  // tweak
  interaction.SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = xsec_model->XSec(&interaction,kps);
  double new_weight = old_weight * (new_xsec/old_xsec);

  // Clear the 'running' kinematics buffer
  interaction.KinePtr()->ClearRunningValues();

  // un-tweak
  interaction.ResetBit(kIAssumeFreeNucleon);

  LOG("ReWeight", pINFO)
     << "Event d{xsec}/dK : " << old_xsec   << " --> " << new_xsec;
  LOG("ReWeight", pINFO)
     << "Event weight     : " << old_weight << " ---> " << new_weight;

  return new_weight;
}
//___________________________________________________________________________

