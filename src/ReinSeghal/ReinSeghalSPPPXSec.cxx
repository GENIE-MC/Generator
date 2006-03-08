//____________________________________________________________________________
/*!

\class    genie::ReinSeghalSPPPXSec

\brief    Computes the differential resonance cross section  for an exclusive
          1pi reaction from resonance neutrinoproduction.

          The computed xsec is the double differential d^2 xsec/ dQ^2 dW \n
          where \n
          \li Q^2 : momentum transfer ^ 2
          \li W   : invariant mass of the final state hadronic system

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Seghal single resonance cross sections
          weighted:

          \li With the value of their Breit-Wigner distributions at the given
              W,Q^2 (The code for BW weighting is included in the single
              resonance cross section algorithm. The user needs to make sure
              that he does not run the single resonance cross section code with
              a configuration that inhibits weighting).

          \li With the isospin Glebsch-Gordon coefficient determining the
              contribution of each resonance to the exclusive final state.

          \li With the BR for the produced resonance to decay into the given
              exclusive final state.

          In this algorithm we follow the non-coherent approach: we sum
          the weighted resonance production cross sections rather than the
          resonance production amplitudes.

          Is a concrete implementation of the XSecAlgorithmI initerface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "ReinSeghal/ReinSeghalSPPPXSec.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalSPPPXSec::ReinSeghalSPPPXSec() :
XSecAlgorithmI("genie::ReinSeghalSPPPXSec")
{

}
//____________________________________________________________________________
ReinSeghalSPPPXSec::ReinSeghalSPPPXSec(string config) :
XSecAlgorithmI("genie::ReinSeghalSPPPXSec", config)
{

}
//____________________________________________________________________________
ReinSeghalSPPPXSec::~ReinSeghalSPPPXSec()
{

}
//____________________________________________________________________________
double ReinSeghalSPPPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalSpp", pINFO) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  LOG("ReinSeghalSpp", pINFO)
                 << "Computing a cross section for " << *interaction;
  LOG("ReinSeghalSpp", pNOTICE)
              << "SPP channel " << SppChannel::AsString(spp_channel);

  //-- Check whether a resonance has been specified
  //   If yes, compute only the contribution of this resonance at the 
  //   specified exclusive state  

  Resonance_t inpres = interaction->GetExclusiveTag().Resonance();
  if(inpres != kNoResonance) {
    LOG("ReinSeghalSpp", pNOTICE) 
      << "Computing only the contribution from: " 
                                     << utils::res::AsString(inpres);
    if(!fResList.Find(inpres)) {
       LOG("ReinSeghalSpp", pWARN) 
            << "The resonance: " << utils::res::AsString(inpres)
                            << " was not found in my resonance list";
       return 0;
    }
    //-- Compute the contribution of this resonance
    double res_xsec_contrib = this->XSec1RES(interaction);
    LOG("ReinSeghalSpp", pNOTICE) 
                       << "d^2 xsec/ dQ^2 dW = " << res_xsec_contrib;
    return res_xsec_contrib;
  }

  //-- Loop over the specified list of baryon resonances and compute
  //   the cross section for the input exlusive channel

  unsigned int nres = fResList.NResonances();
  LOG("ReinSeghalSpp", pNOTICE)
    << "Computing SPP cross section using " << nres << " resonances";

  double xsec = 0;
  for(unsigned int ires = 0; ires < nres; ires++) {

     //-- Get next resonance from the resonance list
     Resonance_t res = fResList.ResonanceId(ires);

     //-- Set current resonance to interaction object
     interaction->GetExclusiveTagPtr()->SetResonance(res);

     //-- Compute the contribution of this resonance
     double res_xsec_contrib = this->XSec1RES(interaction);

     //-- Add contribution of this resonance to the cross section
     xsec += res_xsec_contrib;
  }

  //-- delete the resonance from the input interaction
  interaction->GetExclusiveTagPtr()->SetResonance(kNoResonance);

  LOG("ReinSeghalSpp", pNOTICE) << "d^2 xsec/ dQ^2 dW = " << xsec;
  return xsec;
}
//____________________________________________________________________________
double ReinSeghalSPPPXSec::XSec1RES(const Interaction * interaction) const
{
// computes the contribution of a resonance to a 1pi exlusive reaction

  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  Resonance_t res = interaction->GetExclusiveTag().Resonance();

  //-- Get the Breit-Wigner weighted xsec for exciting the resonance
  double rxsec = fSingleResXSecModel->XSec(interaction);

  //-- Get the BR for the (resonance) -> (exclusive final state)
  double br = SppChannel::BranchingRatio(spp_channel, res);

  //-- Get the Isospin Glebsch-Gordon coefficient for the given resonance
  //   and exclusive final state
  double igg = SppChannel::IsospinWeight(spp_channel, res);

  //-- Compute the weighted xsec
  //  (total weight = Breit-Wigner * BR * isospin Glebsch-Gordon)
  double res_xsec_contrib = rxsec*br*igg;

  SLOG("ReinSeghalSpp", pINFO)
     << "Contrib. from [" << utils::res::AsString(res) << "] = "
     << "<Glebsch-Gordon = " << igg
     << "> * <BR(->1pi) = " << br
     << "> * <Breit-Wigner * d^2xsec/dQ^2dW = " << rxsec
     << "> = " << res_xsec_contrib;

  return res_xsec_contrib;
}
//____________________________________________________________________________
bool ReinSeghalSPPPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    LOG("ReinSeghalSpp", pERROR)
            << "\n *** Insufficient SPP exclusive final state information!";
    return false;
  }
  LOG("ReinSeghalSpp", pINFO)
                       << "Reaction: " << SppChannel::AsString(spp_channel);
  return true;
}
//____________________________________________________________________________
bool ReinSeghalSPPPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double W    = kinematics.W();
  double q2   = kinematics.q2();

  //-- Check energy threshold & kinematical limits in q2, W
  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(E <= EvThr) {
    LOG("ReinSeghalSpp", pINFO) << "E  = " << E << " < Ethr = " << EvThr;
    return false;
  }

  //-- Check against physical range in W and Q2
  Range1D_t rW  = utils::kinematics::WRange(interaction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction);

  bool in_physical_range = utils::math::IsWithinLimits(W, rW) &&
                           utils::math::IsWithinLimits(-q2, rQ2);
  if(!in_physical_range) return false;

  return true;
}
//____________________________________________________________________________
void ReinSeghalSPPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->GetResonanceList();
}
//____________________________________________________________________________
void ReinSeghalSPPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->GetResonanceList();
}
//____________________________________________________________________________
void ReinSeghalSPPPXSec::LoadSubAlg(void)
{
// load the single resonance cross section algorithm specified in the config.

  fSingleResXSecModel = 0;

  fSingleResXSecModel =
       dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                    "single-res-xsec-alg-name", "single-res-xsec-param-set"));
  assert(fSingleResXSecModel);
}
//____________________________________________________________________________
void ReinSeghalSPPPXSec::GetResonanceList(void)
{
// create the baryon resonance list specified in the config.

  fResList.Clear();

  assert( fConfig->Exists("resonance-name-list") );
  string resonanes = fConfig->GetString("resonance-name-list");

  //-- Create a BaryonResList by decoding the resonance list from
  //   the XML input
  //   The list of resonances can be specified as a string with
  //   comma separated resonance names (eg "P33(1233),S11(1535),D13(1520)")
  //   The BaryonResList can also decode lists of pdg-codes or
  //   resonance-ids (Resonance_t enumerations).
  //   Support for this will be added here as well.

  fResList.DecodeFromNameList(resonanes);
}
//____________________________________________________________________________

