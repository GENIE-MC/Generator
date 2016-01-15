//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2007 - CA
   This file was added in 2.0.1
 @ Sep 08, 2009 - CA
   Renamed from ReWeightCrossSection to GReWeightNuXSecHelper and included in
   the genie::rew namespace. Integrated with new event reweighting framework.
 @ Apr 27, 2010 - CA
   Added option to reweight differential cross sections normalizing to a const
   integral (shape only effect of tweaked physics parameter)
*/
//____________________________________________________________________________

#include "Algorithm/AlgCmp.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "ReWeight/GReWeightNuXSecHelper.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::rew;

//___________________________________________________________________________
GReWeightNuXSecHelper::GReWeightNuXSecHelper(void)
{
  this->Initialize();
}
//___________________________________________________________________________
GReWeightNuXSecHelper::~GReWeightNuXSecHelper(void)
{

}
//___________________________________________________________________________
void GReWeightNuXSecHelper::Initialize(void)
{
  this->DiffCrossSecType( kScQuasiElastic,    kPSQ2fE  );
  this->DiffCrossSecType( kScDeepInelastic,   kPSxyfE  );
  this->DiffCrossSecType( kScResonant,        kPSWQ2fE );
  this->DiffCrossSecType( kScCoherent,        kPSxyfE  );
}
//___________________________________________________________________________
void GReWeightNuXSecHelper::HandleInitState(const InitialState & is)
{
  // form initial state filtering out any unwanted info
  InitialState init_state(is.TgtPdg(), is.ProbePdg()); 

  // check for an event generation configured for that initial state
  GEVGDriver * evg_driver = fGPool.FindDriver(init_state);

  // if none was found  then create/configure/store one now
  if(!evg_driver) {
    LOG("ReW", pNOTICE) 
        << "Adding event generation driver for initial state = " 
        << init_state.AsString();
    evg_driver = new GEVGDriver;
    evg_driver->Configure(init_state);
    fGPool.insert( GEVGPool::value_type(init_state.AsString(), evg_driver) );
  }
}
//___________________________________________________________________________
void GReWeightNuXSecHelper::DiffCrossSecType(
        ScatteringType_t sct, KinePhaseSpace_t kps)
{
  fCrossSecModelPhSp.insert(
     map<ScatteringType_t,KinePhaseSpace_t>::value_type(sct,kps));
}
//___________________________________________________________________________
double GReWeightNuXSecHelper::NewWeight(
  const EventRecord & event, bool shape_only)
{
  // Get event summary (Interaction) from the input event
  assert(event.Summary());
  Interaction & interaction = * event.Summary();

  //LOG("ReW", pDEBUG) << "Computing new weight for: \n" << interaction;

  // Find the event generation driver that handles the given initial state
  const InitialState & init_state = interaction.InitState();
  GEVGDriver * evg_driver = fGPool.FindDriver(init_state);
  if(!evg_driver) {
    LOG("ReW", pINFO)
      << "Adding generator driver for init state: " << init_state.AsString();
    evg_driver = new GEVGDriver;
    evg_driver->Configure(init_state);
    fGPool.insert( GEVGPool::value_type(init_state.AsString(), evg_driver) );
  }
  assert(evg_driver);

  // Find the event generation thread that handles the given interaction
  const EventGeneratorI * evg_thread = evg_driver->FindGenerator(&interaction);
  if(!evg_thread) {
    LOG("ReW", pERROR)
      << "No event generator thread for interaction: " << interaction;
    return 0;
  }

  // Get the cross section model associated with that thread
  const XSecAlgorithmI * xsec_model = evg_thread->CrossSectionAlg();
  if(!xsec_model) {
    LOG("ReW", pERROR)
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

  // hack to match what was stored during event generation
  // -- is currently revisited -- 
  if(interaction.ProcInfo().IsQuasiElastic()) 
		interaction.SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = xsec_model->XSec(&interaction,kps);
  double new_weight = old_weight * (new_xsec/old_xsec);

  if(shape_only) {
    double old_integrated_xsec = event.XSec();
    double new_integrated_xsec = xsec_model->Integral(&interaction);
    assert(new_integrated_xsec > 0);
    new_weight *= (old_integrated_xsec/new_integrated_xsec);
  }

  // hack - closing parenthesis
  if(interaction.ProcInfo().IsQuasiElastic()) 
  		interaction.ResetBit(kIAssumeFreeNucleon);

  // Clear the 'running' kinematics buffer
  interaction.KinePtr()->ClearRunningValues();

  LOG("ReW", pINFO)
     << "Event d{xsec}/dK : " << old_xsec   << " --> " << new_xsec;
  LOG("ReW", pINFO)
     << "Event weight     : " << old_weight << " ---> " << new_weight;

  return new_weight;
}
//___________________________________________________________________________

