//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Base/XSecIntegratorI.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "ReinSehgal/ReinSehgalSPPPXSec.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

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

  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  LOG("ReinSehgalSpp", pINFO)
                 << "Computing a cross section for " << *interaction;
  LOG("ReinSehgalSpp", pNOTICE)
              << "SPP channel " << SppChannel::AsString(spp_channel);

  //-- Check whether a resonance has been specified
  //   If yes, compute only the contribution of this resonance at the 
  //   specified exclusive state  

  Resonance_t inpres = interaction->ExclTag().Resonance();
  if(inpres != kNoResonance) {

    string rname = utils::res::AsString(inpres);

    LOG("ReinSehgalSpp", pNOTICE) 
               << "Computing only the contribution from: " << rname;
    if(!fResList.Find(inpres)) {
       LOG("ReinSehgalSpp", pWARN) 
           << "Resonance: " << rname << " was not found in my list";
       return 0;
    }
    //-- Compute the contribution of this resonance
    double xsec1 = this->XSec1RES(interaction,kps);
    LOG("ReinSehgalSpp", pNOTICE) << "d^nxsec/ dK^n = " << xsec1;
    return xsec1;
  }

  //-- Loop over the specified list of baryon resonances and compute
  //   the cross section for the input exclusive channel

  double xsecN = this->XSecNRES(interaction,kps);
  LOG("ReinSehgalSpp", pNOTICE) << "d^nxsec/ dK^n = " << xsecN;
  return xsecN;
}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::XSecNRES(
                const Interaction * interaction, KinePhaseSpace_t kps) const
{
// computes the 1pi cros section taking into account the contribution of all
// specified baryon resonances

  unsigned int nres = fResList.NResonances();
  LOG("ReinSehgalSpp", pNOTICE)
    << "Computing SPP cross section using " << nres << " resonances";

  double xsec = 0;
  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Get next resonance from the resonance list
     Resonance_t res = fResList.ResonanceId(ires);

     //-- Set current resonance to interaction object
     interaction->ExclTagPtr()->SetResonance(res);

     //-- Compute the contribution of this resonance
     double res_xsec_contrib = this->XSec1RES(interaction,kps);

     //-- Add contribution of this resonance to the cross section
     xsec += res_xsec_contrib;
  }

  //-- delete the resonance from the input interaction
  interaction->ExclTagPtr()->SetResonance(kNoResonance);

  return xsec;
}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::XSec1RES(
               const Interaction * interaction, KinePhaseSpace_t kps) const
{
// computes the contribution of a resonance to a 1pi exlusive reaction

  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  Resonance_t res = interaction->ExclTag().Resonance();

  //-- Get the Breit-Wigner weighted xsec for exciting the resonance
  double rxsec = fSingleResXSecModel->XSec(interaction,kps);

  //-- Get the BR for the (resonance) -> (exclusive final state)
  double br = SppChannel::BranchingRatio(spp_channel, res);

  //-- Get the Isospin Glebsch-Gordon coefficient for the given resonance
  //   and exclusive final state
  double igg = SppChannel::IsospinWeight(spp_channel, res);

  //-- Compute the weighted xsec
  //  (total weight = Breit-Wigner * BR * isospin Glebsch-Gordon)
  double res_xsec_contrib = rxsec*br*igg;

  SLOG("ReinSehgalSpp", pINFO)
     << "Contrib. from [" << utils::res::AsString(res) << "] = "
     << "<Glebsch-Gordon = " << igg
     << "> * <BR(->1pi) = " << br
     << "> * <Breit-Wigner * d^nxsec/dK^n = " << rxsec
     << "> = " << res_xsec_contrib;

  return res_xsec_contrib;
}
//____________________________________________________________________________
double ReinSehgalSPPPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
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
  LOG("ReinSehgalSpp", pINFO)
                       << "Reaction: " << SppChannel::AsString(spp_channel);
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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  string resonances = fConfig->GetStringDef(
                     "ResonanceNameList", gc->GetString("ResonanceNameList"));
  fResList.DecodeFromNameList(resonances);

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
