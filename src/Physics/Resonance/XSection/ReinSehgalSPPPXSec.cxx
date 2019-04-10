//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/ReinSehgalSPPPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSehgalSPPPXSec::ReinSehgalSPPPXSec() :
XSecAlgorithmI("genie::ReinSehgalSPPPXSec")
{

}
//____________________________________________________________________________
ReinSehgalSPPPXSec::ReinSehgalSPPPXSec(string config) :
XSecAlgorithmI("genie::ReinSehgalSPPPXSec", config)
{

}
//____________________________________________________________________________
ReinSehgalSPPPXSec::~ReinSehgalSPPPXSec()
{

}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalSpp", pDEBUG)
                 << "Computing a cross section for " << *interaction;
#endif
  //-- Check whether a resonance has been specified
  //   If yes, compute only the contribution of this resonance at the 
  //   specified exclusive state  

  Resonance_t inpres = interaction->ExclTag().Resonance();
  if(inpres != kNoResonance) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
	LOG("ReinSehgalSpp", pDEBUG) 
               << "Computing only the contribution from: " << utils::res::AsString(inpres);
#endif        
    if(!fResList.Find(inpres)) {
       LOG("ReinSehgalSpp", pWARN) 
           << "Resonance: " << utils::res::AsString(inpres) << " was not found in my list";
       return 0;
    }
    //-- Compute the contribution of this resonance
    //-- Get the Breit-Wigner weighted xsec for exciting the resonance
    
    return fSingleResXSecModel->XSec(interaction,kps);
  }

  //-- Loop over the specified list of baryon resonances and compute
  //   the cross section for the input exclusive channel
  
  return this->XSecNRES(interaction,kps);
}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::XSecNRES(
                const Interaction * interaction, KinePhaseSpace_t kps) const
{
// computes the 1pi cros section taking into account the contribution of all
// specified baryon resonances

  unsigned int nres = fResList.NResonances();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__ 
  LOG("ReinSehgalSpp", pDEBUG)
    << "Computing SPP cross section using " << nres << " resonances";
#endif
  
  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalSpp", pDEBUG)
              << "SPP channel " << SppChannel::AsString(spp_channel);
#endif
              
  double xsec = 0;
  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Get next resonance from the resonance list
     Resonance_t res = fResList.ResonanceId(ires);

     //-- Set current resonance to interaction object
     interaction->ExclTagPtr()->SetResonance(res);

     //-- Get the BR for the (resonance) -> (exclusive final state)
	 double br = SppChannel::BranchingRatio(spp_channel, res);

	 //-- Get the Isospin Clebsch-Gordon coefficient for the given resonance
	 //   and exclusive final state
	 double igg = SppChannel::IsospinWeight(spp_channel, res);

	 //-- Compute the weighted xsec
	 //  (total weight = Breit-Wigner * BR * isospin Clebsch-Gordon)
	 double res_xsec_contrib = fSingleResXSecModel->XSec(interaction,kps)*br*igg;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
	 LOG("ReinSehgalSpp", pDEBUG)
     << "Contrib. from [" << utils::res::AsString(res) << "] = "
     << "<Clebsch-Gordon = " << igg
     << "> * <BR(->1pi) = " << br
     << "> * <Breit-Wigner * d^nxsec/dK^n = " << rxsec
     << "> = " << res_xsec_contrib;
#endif

     //-- Add contribution of this resonance to the cross section
     xsec += res_xsec_contrib;
  }

  //-- delete the resonance from the input interaction
  interaction->ExclTagPtr()->SetResonance(kNoResonance);

  return xsec;
}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::Integral(const Interaction * interaction) const
{
  return fXSecIntegrator->Integrate(this,interaction);
}
//____________________________________________________________________________
bool ReinSehgalSPPPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    LOG("ReinSehgalSpp", pERROR)
            << "\n *** Insufficient SPP exclusive final state information!";
    return false;
  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalSpp", pDEBUG)
                       << "Reaction: " << SppChannel::AsString(spp_channel);
#endif
  return true;
}
//____________________________________________________________________________
void ReinSehgalSPPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalSPPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalSPPPXSec::LoadConfig(void)
{
// load the single resonance cross section algorithm specified in the config.

  fSingleResXSecModel =
   dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("SingleRESDiffXSecAlg"));
  assert(fSingleResXSecModel);

  //-- Create a BaryonResList by decoding the resonance list from
  //   the XML input
  //   The list of resonances can be specified as a string with
  //   comma separated resonance names (eg "P33(1233),S11(1535),D13(1520)")
  //   The BaryonResList can also decode lists of pdg-codes or
  //   resonance-ids (Resonance_t enumerations).
  //   Support for this will be added here as well.

  fResList.Clear();


  string resonances ;
  GetParam( "ResonanceNameList", resonances ) ;
  fResList.DecodeFromNameList(resonances);

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
