//____________________________________________________________________________
/*
 Copyright (c) 2023-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp() :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp")
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp(string config) :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp", config)
{

}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::~MAIDHelicityAmplModelEMp()
{

}
//____________________________________________________________________________
const MAIDHelicityAmpl & MAIDHelicityAmplModelEMp::Compute( const Interaction * interaction ) const {

  Resonance_t res = interaction->ExclTag().Resonance();

  switch(res) {

  }//switch

  return fAmpl;
}
//____________________________________________________________________________
